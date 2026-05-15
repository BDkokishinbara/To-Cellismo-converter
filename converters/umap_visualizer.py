"""
h5mu (MuData) ファイルから UMAP を生成するモジュール

指定されたモダリティに対して scanpy の標準前処理パイプラインを実行し、
UMAP の散布図を PNG として保存する。ユーザーがボタン 1 つでデータの
構造を素早く確認できるようにすることが目的。
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
# Flask のリクエストスレッドからも安全に動作させるため、非対話バックエンドを強制する。
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import mudata as md
import scanpy as sc


def umap_from_h5mu(
    h5mu_path: str,
    output_png_path: str,
    modality: str | None = None,
    min_cells: int = 3,
    min_genes: int = 200,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    color_by: str | None = None,
    do_clustering: bool = True,
    random_state: int = 0,
    progress_cb=None,
) -> dict:
    """h5mu ファイルを読み込み、標準前処理 → UMAP → プロットの一連を実行し、
    結果の散布図を ``output_png_path`` に保存する。

    ``progress_cb`` は各処理ステージ毎に短い説明文字列を受け取るコールバック。
    UI 側で進捗をリアルタイム表示する用途に使う（省略可）。

    Returns
    -------
    dict
        生成結果のサマリ。
    """
    def _step(msg: str):
        """ログとコールバックの両方に進捗を流す内部ヘルパー。"""
        print(f"[UMAP] {msg}")
        if progress_cb is not None:
            try:
                progress_cb(msg)
            except Exception:
                # UI コールバックの失敗で計算側を壊さないよう例外を握りつぶす。
                pass

    h5mu_path = str(h5mu_path)
    output_png_path = str(output_png_path)

    _step('h5mu を読み込み中...')
    mdata = md.read_h5mu(h5mu_path)

    # 使用するモダリティを決める。'rna' があれば優先、それ以外は最初のもの。
    if modality is None:
        if 'rna' in mdata.mod:
            modality = 'rna'
        else:
            modality = next(iter(mdata.mod.keys()))
    if modality not in mdata.mod:
        raise ValueError(f"Modality '{modality}' not found in file. Available: {list(mdata.mod.keys())}")

    adata = mdata[modality].copy()

    # データが極端に小さい場合は早い段階で分かりやすいメッセージで失敗させる。
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError(
            f"データが小さすぎて UMAP を作成できません: cells={adata.n_obs}, genes={adata.n_vars}"
        )

    # scanpy が期待する一意な文字列名にそろえる
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # ── 高速パス ────────────────────────────────────────────────────────
    # 既に obsm['X_umap'] が格納されていれば、再計算せずそれをそのまま使う。
    # CXG キュレーション済みの h5ad は X_umap を事前計算した状態で配布される
    # ことが多い。
    if 'X_umap' in getattr(adata, 'obsm', {}) and adata.obsm['X_umap'] is not None:
        _step('既存の UMAP (obsm[X_umap]) を再利用します — PCA/近傍探索/UMAP 計算をスキップ')
        chosen_color = _pick_color_column(adata, color_by, do_clustering, random_state, _step)
        return _save_umap_plot(adata, modality, chosen_color, output_png_path, _step)

    # ── 標準的な scanpy 前処理パイプライン ─────────────────────────────
    _step(f'細胞・遺伝子をフィルタ中... ({adata.n_obs} cells × {adata.n_vars} genes)')
    sc.pp.filter_cells(adata, min_genes=min(min_genes, adata.n_vars - 1))
    sc.pp.filter_genes(adata, min_cells=min(min_cells, adata.n_obs - 1))

    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError(
            "フィルタリング後に十分な細胞/遺伝子が残りませんでした。データを確認してください。"
        )

    _step('正規化中... (normalize_total + log1p)')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 高変動遺伝子 (HVG) を選択して以降の解析対象を絞る
    _step('高変動遺伝子 (HVG) を選択中...')
    try:
        n_top = min(n_top_genes, adata.n_vars)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat')
        if 'highly_variable' in adata.var.columns and adata.var['highly_variable'].any():
            adata = adata[:, adata.var['highly_variable']].copy()
    except Exception:
        # 極端にスパース・小規模なデータで HVG 計算が失敗することがあるので、
        # その場合はスキップして先に進む。
        pass

    # PCA — n_comps は (n_obs, n_vars) の小さい方より小さくする必要がある
    n_comps = max(2, min(n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    _step(f'PCA 計算中... (n_comps={n_comps})')
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)

    # 近傍グラフの構築（UMAP/Leiden の前提）
    n_nb = max(2, min(n_neighbors, adata.n_obs - 1))
    _step(f'近傍グラフ構築中... (n_neighbors={n_nb})')
    sc.pp.neighbors(adata, n_neighbors=n_nb, n_pcs=n_comps, random_state=random_state)

    # UMAP 本体
    _step('UMAP 計算中... (数十秒〜数分かかる場合があります)')
    sc.tl.umap(adata, random_state=random_state)

    chosen_color = _pick_color_column(adata, color_by, do_clustering, random_state, _step)
    return _save_umap_plot(adata, modality, chosen_color, output_png_path, _step)


# UMAP プロットの着色に使う候補列。優先順は「細胞型らしいもの」を上位に。
# 一般的に最も情報量が多いため。
_PREFERRED_COLOR_COLUMNS = (
    'celltype', 'cell_type', 'cell_types', 'CellType', 'Cell_Type',
    'celltype_l1', 'celltype_l2', 'celltype_l3',
    'sorted_celltype', 'annotation', 'cluster', 'clusters',
    'leiden', 'louvain',
    'n_genes_by_counts', 'n_counts', 'total_counts', 'n_genes',
)


def _pick_color_column(adata, color_by, do_clustering, random_state, step_fn):
    """UMAP の着色に使う obs 列を決める。"""
    if color_by is not None:
        return color_by

    # 既存の細胞型アノテーションがあればそれを優先
    for col in _PREFERRED_COLOR_COLUMNS:
        if col in adata.obs.columns:
            return col

    # 生物学的アノテーションが無い → クラスタリングを試みる
    if do_clustering:
        try:
            step_fn('Leiden クラスタリング中...')
            sc.tl.leiden(adata, random_state=random_state)
            return 'leiden'
        except Exception:
            pass

    return None


def _save_umap_plot(adata, modality, chosen_color, output_png_path, step_fn):
    """共通の描画＆保存処理。高速パス・通常パスの両方から呼ばれる。"""
    Path(output_png_path).parent.mkdir(parents=True, exist_ok=True)
    step_fn('プロットを描画中...')
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    try:
        sc.pl.umap(
            adata,
            color=chosen_color,
            ax=ax,
            show=False,
            frameon=False,
            legend_fontsize=7,
            size=12,
        )
        title = f"UMAP — modality: {modality}"
        if chosen_color:
            title += f"  (color: {chosen_color})"
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(output_png_path, dpi=140, bbox_inches='tight')
    finally:
        plt.close(fig)

    step_fn('完了')
    return {
        'modality': modality,
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'color_by': chosen_color,
        'output': output_png_path,
    }
