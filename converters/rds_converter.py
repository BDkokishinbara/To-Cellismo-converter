"""
RDS → h5mu 変換モジュール

R で保存された RDS ファイル（Seurat オブジェクト または
SingleCellExperiment オブジェクト）を h5mu (MuData) 形式に変換する。
本モジュールは ``rpy2`` 経由で R を呼び出すため、R 本体と rpy2 のインストールが必要。
"""
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from pathlib import Path

# rpy2 がインストールされていない環境でも import エラーにならないように
# try/except で囲み、利用可否フラグを立てておく。
try:
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()
    numpy2ri.activate()
    RPY2_AVAILABLE = True
except ImportError:
    RPY2_AVAILABLE = False


def rds_to_h5mu(rds_path, output_path, object_type='auto'):
    """RDS ファイル（Seurat または SingleCellExperiment）を h5mu 形式に変換する。

    Parameters
    ----------
    rds_path : str
        入力 RDS ファイルのパス。
    output_path : str
        出力 h5mu ファイルのパス。
    object_type : str
        R オブジェクトの型: ``'seurat'`` / ``'sce'`` / ``'auto'`` （既定）。
        ``'auto'`` の場合は class() で自動判定する。

    Returns
    -------
    str
        作成された h5mu ファイルのパス。
    """
    if not RPY2_AVAILABLE:
        raise ImportError(
            "rpy2 is not installed. Install it with: pip install rpy2\n"
            "Also make sure R is installed on your system."
        )

    try:
        print(f"Reading RDS file: {rds_path}")

        # R オブジェクトとしてロード
        readRDS = robjects.r['readRDS']
        r_object = readRDS(rds_path)

        # 型を自動判定
        if object_type == 'auto':
            object_class = robjects.r['class'](r_object)[0]
            print(f"Detected R object class: {object_class}")

            if 'Seurat' in object_class:
                object_type = 'seurat'
            elif 'SingleCellExperiment' in object_class or 'SCE' in object_class:
                object_type = 'sce'
            else:
                raise ValueError(
                    f"Unknown R object class: {object_class}. "
                    "Please specify object_type as 'seurat' or 'sce'"
                )

        # 型に応じて AnnData へ変換
        if object_type == 'seurat':
            adata = convert_seurat_to_anndata(r_object)
        elif object_type == 'sce':
            adata = convert_sce_to_anndata(r_object)
        else:
            raise ValueError(f"Unknown object_type: {object_type}")

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # 1 モダリティの MuData として包む
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
        raise Exception(f"Error converting RDS to h5mu: {str(e)}")


def convert_seurat_to_anndata(seurat_obj):
    """Seurat オブジェクトを AnnData に変換する。

    Parameters
    ----------
    seurat_obj : rpy2 R オブジェクト
        R 側の Seurat オブジェクト。

    Returns
    -------
    AnnData
        変換後の AnnData オブジェクト。
    """
    try:
        # Seurat 変換関数の利用準備
        base = importr('base')

        # カウント行列の取得
        # Seurat のバージョンによって API が異なるので、上から順に試す。
        try:
            # Seurat v5
            robjects.r('suppressMessages(library(Seurat))')
            counts = robjects.r('GetAssayData')(seurat_obj, slot='counts')
        except:
            try:
                # Seurat v3 / v4
                counts = robjects.r('GetAssayData')(seurat_obj, assay='RNA', slot='counts')
            except:
                # フォールバック: スロットを直接参照
                counts = seurat_obj.slots['assays'].slots['RNA'].slots['counts']

        # スパース行列（dgCMatrix）なら密行列に変換
        if robjects.r('inherits')(counts, 'dgCMatrix')[0]:
            counts = robjects.r('as.matrix')(counts)

        # numpy 配列に変換し、AnnData が期待する (細胞 × 遺伝子) に転置
        counts_np = np.array(counts).T

        # 細胞メタデータ（meta.data）を取得
        try:
            metadata = robjects.r('as.data.frame')(seurat_obj.slots['meta.data'])
            obs = pandas2ri.rpy2py(metadata)
        except:
            obs = pd.DataFrame(index=range(counts_np.shape[0]))

        # 遺伝子名を取得
        try:
            var_names = list(robjects.r('rownames')(counts))
            var = pd.DataFrame(index=var_names)
        except:
            var = pd.DataFrame(index=range(counts_np.shape[1]))

        # AnnData として組み立てる
        adata = ad.AnnData(X=counts_np, obs=obs, var=var)

        return adata

    except Exception as e:
        raise Exception(f"Error converting Seurat object: {str(e)}")


def convert_sce_to_anndata(sce_obj):
    """SingleCellExperiment オブジェクトを AnnData に変換する。

    Parameters
    ----------
    sce_obj : rpy2 R オブジェクト
        R 側の SingleCellExperiment オブジェクト。

    Returns
    -------
    AnnData
        変換後の AnnData オブジェクト。
    """
    try:
        # SingleCellExperiment のインポート
        robjects.r('suppressMessages(library(SingleCellExperiment))')

        # カウント行列を取得
        counts = robjects.r('counts')(sce_obj)

        # スパース（dgCMatrix）なら密行列に変換
        if robjects.r('inherits')(counts, 'dgCMatrix')[0]:
            counts = robjects.r('as.matrix')(counts)

        # numpy 配列に変換し、(細胞 × 遺伝子) に転置
        counts_np = np.array(counts).T

        # 細胞メタデータ（colData）を取得
        try:
            col_data = robjects.r('colData')(sce_obj)
            col_data = robjects.r('as.data.frame')(col_data)
            obs = pandas2ri.rpy2py(col_data)
        except:
            obs = pd.DataFrame(index=range(counts_np.shape[0]))

        # 遺伝子メタデータ（rowData）を取得
        try:
            row_data = robjects.r('rowData')(sce_obj)
            row_data = robjects.r('as.data.frame')(row_data)
            var = pandas2ri.rpy2py(row_data)
        except:
            # rowData が取れなければ、せめて遺伝子名（行名）だけは取得
            try:
                var_names = list(robjects.r('rownames')(sce_obj))
                var = pd.DataFrame(index=var_names)
            except:
                var = pd.DataFrame(index=range(counts_np.shape[1]))

        # AnnData として組み立てる
        adata = ad.AnnData(X=counts_np, obs=obs, var=var)

        return adata

    except Exception as e:
        raise Exception(f"Error converting SingleCellExperiment object: {str(e)}")


def check_rpy2_available():
    """rpy2 が利用可能かどうかを返す。"""
    return RPY2_AVAILABLE
