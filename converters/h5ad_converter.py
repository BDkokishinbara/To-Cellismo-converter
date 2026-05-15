"""
h5ad → h5mu 変換モジュール

AnnData (.h5ad) を 1 つのモダリティ ``rna`` を持つ MuData にラップして
h5mu 形式で書き出す。CXG 形式の h5ad では ``var.index`` が Ensembl ID
（``ENSG...``）になっていることが多いため、遺伝子記号列が見つかれば
``var.index`` に差し替える処理も行う。
"""
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd


# 遺伝子記号（読みやすい名前）が格納されていそうな var の列名候補。
# 上から順に検索し、最初に見つかった列を採用する。
_GENE_SYMBOL_COLUMNS = ('feature_name', 'gene_symbol', 'gene_symbols', 'gene_name', 'gene_names', 'Symbol')


def _remap_one_var_frame(var_df):
    """1 つの var データフレームに対して、Ensembl ID → 遺伝子記号への
    差し替え用 Index を計算して返す。

    Parameters
    ----------
    var_df : pandas.DataFrame
        ``adata.var`` または ``adata.raw.var`` のどちらか。

    Returns
    -------
    (new_index, source_col) : (pandas.Index, str) または (None, None)
        差し替え可能なら新しい Index と参照した列名を返す。
        既に記号らしい場合や、記号列が見つからない場合は (None, None)。

    Notes
    -----
    引数の DataFrame は変更しない（純関数）。差し替えは呼び出し元で行う。
    """
    # 先頭 50 件をサンプリングして、半数以上が "ENS" で始まるなら Ensembl ID と判定。
    sample = var_df.index[:50].astype(str).tolist()
    looks_like_ensembl = sum(1 for s in sample if s.startswith('ENS')) >= len(sample) * 0.5
    if not looks_like_ensembl:
        return None, None

    # 候補列を順に試して、最初に見つかったものを採用。
    for col in _GENE_SYMBOL_COLUMNS:
        if col in var_df.columns:
            symbols = var_df[col].astype(str)
            ensembl_strs = var_df.index.astype(str)
            # 空文字 / NaN の記号は元の Ensembl ID にフォールバック。
            blank = symbols.isin(['', 'nan', 'None', 'NaN'])
            symbols = symbols.where(~blank, pd.Series(ensembl_strs, index=var_df.index))
            return pd.Index(symbols, name='gene'), col
    return None, None


def _remap_var_index_to_symbols(adata):
    """``var.index`` と（存在すれば）``raw.var.index`` を Ensembl ID から
    遺伝子記号に差し替える。下流のツール（Cellismo など）で読みやすい
    名前が表示されるようにするための処理。

    元の Ensembl ID は ``var`` / ``raw.var`` の ``ensembl_id`` 列に保存される。
    既に記号らしい場合は何もしない。

    Returns
    -------
    (remapped, source_column) : (bool, str | None)
        差し替えを行ったか、参照した列名は何かを返す。
    """
    # メインの var を処理
    new_index, src = _remap_one_var_frame(adata.var)
    remapped = False
    if new_index is not None:
        if 'ensembl_id' not in adata.var.columns:
            adata.var['ensembl_id'] = adata.var.index.astype(str)
        adata.var.index = new_index
        adata.var_names_make_unique()
        remapped = True

    # raw が存在する場合はそちらも処理。
    # CXG 形式の h5ad は ``adata.raw`` に正規化前のカウントが格納されており、
    # ``raw.var.index`` も Ensembl ID のまま。Cellismo はこの raw の方を
    # 参照することがあるため、両方差し替える必要がある。
    if getattr(adata, 'raw', None) is not None:
        raw_new_index, raw_src = _remap_one_var_frame(adata.raw.var)
        if raw_new_index is not None:
            # ``adata.raw`` は読み取り専用なので、いったん AnnData に戻して
            # 編集し、再代入する形で差し替える。
            raw_adata = adata.raw.to_adata()
            if 'ensembl_id' not in raw_adata.var.columns:
                raw_adata.var['ensembl_id'] = raw_adata.var.index.astype(str)
            raw_adata.var.index = raw_new_index
            raw_adata.var_names_make_unique()
            adata.raw = raw_adata
            remapped = True
            src = src or raw_src

    return remapped, src


def h5ad_to_h5mu(h5ad_path, output_path, modality='rna'):
    """h5ad (AnnData) ファイルを h5mu (MuData) 形式に変換する。

    Parameters
    ----------
    h5ad_path : str
        入力 .h5ad ファイルのパス。
    output_path : str
        出力 .h5mu ファイルのパス。
    modality : str
        MuData 内でこの AnnData を格納するモダリティ名。
        本プロジェクトの他のコンバータと統一するため、既定は ``'rna'``。

    Returns
    -------
    str
        作成された h5mu ファイルのパス。
    """
    try:
        print(f"Reading h5ad file: {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path)
        print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # var.index を読みやすい遺伝子記号に差し替える（可能なら）
        remapped, source_col = _remap_var_index_to_symbols(adata)
        if remapped:
            print(f"Remapped var_names to gene symbols using column '{source_col}' (Ensembl IDs preserved in 'ensembl_id')")

        # 基本的なメタデータを追記（既に存在する場合は上書きしない）
        if 'n_counts' not in adata.obs.columns:
            adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        if 'n_cells' not in adata.var.columns:
            adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        # 1 モダリティの MuData として包む
        mdata = md.MuData({modality: adata})

        print(f"Saving to: {output_path}")
        # gzip 圧縮を使って入力 h5ad に近いサイズに保つ。
        # 圧縮なしで書き出すとスパース行列が数倍に膨らんでしまうので
        # この対応は重要。
        try:
            mdata.write(output_path, compression='gzip')
        except TypeError:
            # 古い mudata では compression 引数を受け付けない場合があるので
            # 圧縮なしで書き出す。
            mdata.write(output_path)
        print("Conversion completed successfully!")
        return output_path

    except Exception as e:
        raise Exception(f"Error converting h5ad to h5mu: {str(e)}")


def validate_h5ad(h5ad_path):
    """h5ad ファイルの概要情報を返す軽量バリデータ。

    細胞数・遺伝子数・obs/var の主な列名だけを返し、ファイル全体は
    メモリに読み込まない（backed='r' でディスク参照）。
    """
    try:
        adata = ad.read_h5ad(h5ad_path, backed='r')
        info = {
            'n_obs': int(adata.n_obs),
            'n_vars': int(adata.n_vars),
            'obs_columns': list(adata.obs.columns)[:20],
            'var_columns': list(adata.var.columns)[:20],
        }
        try:
            adata.file.close()
        except Exception:
            pass
        return info
    except Exception as e:
        raise Exception(f"Error validating h5ad: {str(e)}")
