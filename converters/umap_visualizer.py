"""
UMAP visualization from an h5mu (MuData) file.

Runs a standard scanpy preprocessing pipeline on the chosen modality and
saves a UMAP scatter plot as a PNG so the user can quickly inspect the
structure of their data.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
# Force non-interactive backend so this works in a Flask request thread.
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
    """
    Read an h5mu file, run a standard preprocessing pipeline, compute UMAP
    and save the plot to ``output_png_path``.

    ``progress_cb`` is an optional callable that receives short stage
    descriptions so a caller can surface live progress to a UI.

    Returns a small dict summarising what was produced.
    """
    def _step(msg: str):
        print(f"[UMAP] {msg}")
        if progress_cb is not None:
            try:
                progress_cb(msg)
            except Exception:
                # never let UI-callback errors break the pipeline
                pass

    h5mu_path = str(h5mu_path)
    output_png_path = str(output_png_path)

    _step('h5mu を読み込み中...')
    mdata = md.read_h5mu(h5mu_path)

    # Pick a modality. Default to 'rna' if present, otherwise the first one.
    if modality is None:
        if 'rna' in mdata.mod:
            modality = 'rna'
        else:
            modality = next(iter(mdata.mod.keys()))
    if modality not in mdata.mod:
        raise ValueError(f"Modality '{modality}' not found in file. Available: {list(mdata.mod.keys())}")

    adata = mdata[modality].copy()

    # Guard against degenerate inputs early with a clear message.
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError(
            f"データが小さすぎて UMAP を作成できません: cells={adata.n_obs}, genes={adata.n_vars}"
        )

    # Coerce var/obs names to unique strings for scanpy.
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # ── Fast path: if a UMAP embedding is already stored on the object, reuse
    # it instead of re-computing from scratch. CXG-curated h5ad files
    # typically ship with obsm['X_umap'] pre-computed.
    if 'X_umap' in getattr(adata, 'obsm', {}) and adata.obsm['X_umap'] is not None:
        _step('既存の UMAP (obsm[X_umap]) を再利用します — PCA/近傍探索/UMAP 計算をスキップ')
        chosen_color = _pick_color_column(adata, color_by, do_clustering, random_state, _step)
        return _save_umap_plot(adata, modality, chosen_color, output_png_path, _step)

    # ── Standard scanpy preprocessing ──────────────────────────────────────
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

    # HVG → keep
    _step('高変動遺伝子 (HVG) を選択中...')
    try:
        n_top = min(n_top_genes, adata.n_vars)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat')
        if 'highly_variable' in adata.var.columns and adata.var['highly_variable'].any():
            adata = adata[:, adata.var['highly_variable']].copy()
    except Exception:
        # If HVG step fails (e.g. extremely sparse / small data) skip it.
        pass

    # PCA: cap n_comps to min(n_obs, n_vars) - 1
    n_comps = max(2, min(n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    _step(f'PCA 計算中... (n_comps={n_comps})')
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)

    # Neighbours
    n_nb = max(2, min(n_neighbors, adata.n_obs - 1))
    _step(f'近傍グラフ構築中... (n_neighbors={n_nb})')
    sc.pp.neighbors(adata, n_neighbors=n_nb, n_pcs=n_comps, random_state=random_state)

    # UMAP
    _step('UMAP 計算中... (数十秒〜数分かかる場合があります)')
    sc.tl.umap(adata, random_state=random_state)

    chosen_color = _pick_color_column(adata, color_by, do_clustering, random_state, _step)
    return _save_umap_plot(adata, modality, chosen_color, output_png_path, _step)


# Candidate obs columns that we prefer (in order) when picking a color for UMAP.
# Cell-type-like columns first because they are usually the most informative.
_PREFERRED_COLOR_COLUMNS = (
    'celltype', 'cell_type', 'cell_types', 'CellType', 'Cell_Type',
    'celltype_l1', 'celltype_l2', 'celltype_l3',
    'sorted_celltype', 'annotation', 'cluster', 'clusters',
    'leiden', 'louvain',
    'n_genes_by_counts', 'n_counts', 'total_counts', 'n_genes',
)


def _pick_color_column(adata, color_by, do_clustering, random_state, step_fn):
    """Pick an obs column to color the UMAP scatter with."""
    if color_by is not None:
        return color_by

    # Prefer cell-type-like annotations already present.
    for col in _PREFERRED_COLOR_COLUMNS:
        if col in adata.obs.columns:
            return col

    # No biological annotation available — try clustering.
    if do_clustering:
        try:
            step_fn('Leiden クラスタリング中...')
            sc.tl.leiden(adata, random_state=random_state)
            return 'leiden'
        except Exception:
            pass

    return None


def _save_umap_plot(adata, modality, chosen_color, output_png_path, step_fn):
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
