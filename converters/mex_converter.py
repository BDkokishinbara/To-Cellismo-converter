"""
MEX (10x Genomics) → h5mu 変換モジュール

10x Genomics の MEX 形式（matrix.mtx + features.tsv + barcodes.tsv の
3 ファイルからなるスパース行列）を読み込んで h5mu (MuData) 形式に変換する。
複数の feature type が混在する場合（例: Gene Expression + Antibody
Capture）は、それぞれ別のモダリティとして MuData に格納する。
"""
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from scipy.io import mmread
from pathlib import Path
import gzip


def mex_to_h5mu(mex_dir, output_path, genome=None):
    """MEX 形式（10x Genomics）を h5mu 形式に変換する。

    MEX 形式は以下の 3 ファイルで構成される:
      - matrix.mtx (または matrix.mtx.gz) : Matrix Market 形式のスパース行列
      - features.tsv (または genes.tsv とその .gz 版) : 遺伝子・フィーチャー情報
      - barcodes.tsv (または barcodes.tsv.gz) : 細胞バーコード情報

    Parameters
    ----------
    mex_dir : str
        MEX 形式のファイル群を含むディレクトリのパス。
    output_path : str
        出力 h5mu ファイルのパス。
    genome : str, optional
        複数ゲノムが含まれる場合のゲノム名（例: 'GRCh38'）。
        ※ 現状未使用。将来の拡張用。

    Returns
    -------
    str
        作成された h5mu ファイルのパス。
    """
    try:
        print(f"Reading MEX files from: {mex_dir}")
        mex_path = Path(mex_dir)

        # 行列ファイルを探す
        matrix_file = find_file(mex_path, ['matrix.mtx.gz', 'matrix.mtx'])
        if not matrix_file:
            raise FileNotFoundError("Matrix file (matrix.mtx or matrix.mtx.gz) not found")

        # フィーチャー / 遺伝子ファイルを探す（10x v3 は features、v2 は genes）
        features_file = find_file(mex_path, ['features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv'])
        if not features_file:
            raise FileNotFoundError("Features/genes file not found")

        # バーコードファイルを探す
        barcodes_file = find_file(mex_path, ['barcodes.tsv.gz', 'barcodes.tsv'])
        if not barcodes_file:
            raise FileNotFoundError("Barcodes file not found")

        print(f"Found files:")
        print(f"  Matrix: {matrix_file.name}")
        print(f"  Features: {features_file.name}")
        print(f"  Barcodes: {barcodes_file.name}")

        # 行列の読み込み
        print("Reading matrix...")
        matrix = mmread(matrix_file)
        # 行方向の操作を高速化するため CSR 形式に変換
        matrix = matrix.tocsr()
        print(f"Matrix shape: {matrix.shape}")

        # フィーチャー（遺伝子）情報の読み込み
        print("Reading features...")
        features_df = read_tsv(features_file)

        # 列数に応じて 10x のバージョンを判定
        if features_df.shape[1] >= 3:
            # 10x Genomics v3 形式: gene_id, gene_name, feature_type
            features_df.columns = ['gene_id', 'gene_name', 'feature_type'][:features_df.shape[1]]
        elif features_df.shape[1] == 2:
            # 10x Genomics v2 形式: gene_id, gene_name
            features_df.columns = ['gene_id', 'gene_name']
        else:
            # 1 列のみ: gene_id だけ
            features_df.columns = ['gene_id']

        # var_names として gene_name を優先（無ければ gene_id）
        if 'gene_name' in features_df.columns:
            var_names = features_df['gene_name'].astype(str).values
        else:
            var_names = features_df['gene_id'].astype(str).values

        # バーコードの読み込み
        print("Reading barcodes...")
        barcodes_df = read_tsv(barcodes_file, header=None)
        barcodes = barcodes_df[0].astype(str).values

        # AnnData を作成
        print("Creating AnnData object...")
        # MEX の行列は (遺伝子 × 細胞) なので、AnnData が期待する (細胞 × 遺伝子) に転置する
        adata = ad.AnnData(X=matrix.T)

        # obs（細胞）の名前を設定
        adata.obs_names = barcodes

        # var（遺伝子）の名前を設定
        adata.var_names = var_names

        # フィーチャーのメタデータを var に追記
        for col in features_df.columns:
            adata.var[col] = features_df[col].values

        # 基本的な QC 指標を追加
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
        adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # 複数の feature_type が含まれている場合はモダリティ別に分割する
        if 'feature_type' in adata.var.columns:
            feature_types = adata.var['feature_type'].unique()
            print(f"Feature types found: {feature_types}")

            modalities = {}

            for feat_type in feature_types:
                mask = adata.var['feature_type'] == feat_type
                adata_subset = adata[:, mask].copy()

                # feature_type を MuData のモダリティ名にマッピング
                if feat_type == 'Gene Expression':
                    mod_name = 'rna'
                elif feat_type == 'Antibody Capture':
                    mod_name = 'adt'
                elif feat_type == 'CRISPR Guide Capture':
                    mod_name = 'crispr'
                else:
                    mod_name = feat_type.lower().replace(' ', '_')

                modalities[mod_name] = adata_subset
                print(f"  {feat_type}: {adata_subset.shape[1]} features -> '{mod_name}' modality")

            # 複数モダリティを持つ MuData を作成
            mdata = md.MuData(modalities)
        else:
            # 単一モダリティ（RNA）として保存
            mdata = md.MuData({'rna': adata})

        # h5mu として書き出す（ファイルサイズ抑制のため gzip 圧縮）
        print(f"Saving to: {output_path}")
        try:
            mdata.write(output_path, compression='gzip')
        except TypeError:
            mdata.write(output_path)

        print("Conversion completed successfully!")
        return output_path

    except Exception as e:
        raise Exception(f"Error converting MEX to h5mu: {str(e)}")


def find_file(directory, filenames):
    """候補ファイル名のリストから、最初に見つかったファイルのパスを返す。

    Parameters
    ----------
    directory : Path
        検索対象のディレクトリ。
    filenames : list
        候補ファイル名のリスト（優先順）。

    Returns
    -------
    Path or None
        見つかったファイルのパス。どれも無ければ None。
    """
    for filename in filenames:
        filepath = directory / filename
        if filepath.exists():
            return filepath
    return None


def read_tsv(filepath, header=None):
    """TSV ファイルを読み込む。プレーン / gzip 圧縮の両方に対応。

    Parameters
    ----------
    filepath : Path
        TSV ファイルのパス。
    header : int or None
        ヘッダーとして扱う行番号。なしなら None。

    Returns
    -------
    pandas.DataFrame
        読み込んだデータフレーム。
    """
    if filepath.suffix == '.gz':
        with gzip.open(filepath, 'rt') as f:
            return pd.read_csv(f, sep='\t', header=header)
    else:
        return pd.read_csv(filepath, sep='\t', header=header)


def validate_mex_directory(mex_dir):
    """MEX ディレクトリの構造をチェックする。

    Parameters
    ----------
    mex_dir : str
        MEX ファイル群を含むディレクトリのパス。

    Returns
    -------
    dict
        各ファイルの有無と検出ファイル名を含む情報辞書。
    """
    try:
        mex_path = Path(mex_dir)

        if not mex_path.exists():
            raise FileNotFoundError(f"Directory not found: {mex_dir}")

        if not mex_path.is_dir():
            raise NotADirectoryError(f"Not a directory: {mex_dir}")

        # 必要な 3 種のファイルがそろっているかチェック
        matrix_file = find_file(mex_path, ['matrix.mtx.gz', 'matrix.mtx'])
        features_file = find_file(mex_path, ['features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv'])
        barcodes_file = find_file(mex_path, ['barcodes.tsv.gz', 'barcodes.tsv'])

        info = {
            'has_matrix': matrix_file is not None,
            'has_features': features_file is not None,
            'has_barcodes': barcodes_file is not None,
            'matrix_file': matrix_file.name if matrix_file else None,
            'features_file': features_file.name if features_file else None,
            'barcodes_file': barcodes_file.name if barcodes_file else None,
        }

        return info

    except Exception as e:
        raise Exception(f"Error validating MEX directory: {str(e)}")
