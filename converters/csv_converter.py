"""
CSV → h5mu 変換モジュール

シングルセル解析データのカウントマトリックス CSV を読み込んで、
h5mu (MuData) 形式に変換する。SeqGeq 形式（メタデータセクション付き）
の CSV と、行/列の自動転置判定にも対応している。
"""
import re
import pandas as pd
import anndata as ad
import mudata as md
import numpy as np
from pathlib import Path


# ── 行/列のどちらに遺伝子名があるかを推定するためのヒューリスティック ──

# 遺伝子の記号 / ID らしく見えるパターン
_GENE_PATTERNS = [
    # ヒトの大文字表記の遺伝子記号: ACTB, CD3D, HLA-DRA, MT-ATP6, RPS4Y1, IFITM3 など
    re.compile(r'^[A-Z][A-Z0-9._-]{0,19}$'),
    # 各種生物の Ensembl ID: ENSG..., ENSMUSG..., ENSDARG... など
    re.compile(r'^ENS[A-Z]{0,3}[GTPR]\d{6,}'),
    # マウス由来の混合大小文字表記: Actb, Cd3d, Ms4a1, Hbb-b1 など
    re.compile(r'^[A-Z][a-z][a-zA-Z0-9.-]{1,18}$'),
    # MGI（Mouse Genome Informatics）のアクセション ID
    re.compile(r'^MGI:\d+$'),
]

# 細胞バーコードらしく見えるパターン（＝該当軸は細胞であって遺伝子ではない）
_BARCODE_PATTERNS = [
    # 10x Genomics 形式のバーコード。任意で -1 等のサンプルサフィックス付き。
    re.compile(r'^[ACGTN]{12,24}(-\d+)?$'),
    # サンプル名プレフィックス付きのバーコード: Sample1_AAACCTGAG..., Patient3-AAACCT...
    re.compile(r'^[A-Za-z0-9]+[_-][ACGTN]{12,24}(-\d+)?$'),
]


def _fraction_matching(labels, patterns):
    """``labels`` のうち、``patterns`` のいずれかにマッチするものの割合を返す。"""
    if not labels:
        return 0.0
    n = sum(1 for s in (str(x) for x in labels) if any(p.match(s) for p in patterns))
    return n / len(labels)


def detect_transpose_needed(df, sample_size=500):
    """``df`` を AnnData / scanpy / MuData が期待する形（細胞 × 遺伝子）に
    そろえるために転置が必要かどうかを判定する。

    Returns
    -------
    (transpose, info) : (bool, dict)
        ``transpose`` が True のときは、入力データが行=遺伝子の並びで、
        列方向にひっくり返す必要がある（＝結果として行=細胞、列=遺伝子になる）。
        ``info`` には判定根拠となったスコアが入っている。UI 側で
        「なぜそう判断したか」を表示するために利用する。
    """
    rows = df.index[:sample_size].astype(str).tolist()
    cols = df.columns[:sample_size].astype(str).tolist()

    row_gene = _fraction_matching(rows, _GENE_PATTERNS)
    col_gene = _fraction_matching(cols, _GENE_PATTERNS)
    row_bc = _fraction_matching(rows, _BARCODE_PATTERNS)
    col_bc = _fraction_matching(cols, _BARCODE_PATTERNS)

    info = {
        'shape': tuple(df.shape),
        'row_gene_score': round(row_gene, 3),
        'col_gene_score': round(col_gene, 3),
        'row_barcode_score': round(row_bc, 3),
        'col_barcode_score': round(col_bc, 3),
    }

    # 行に遺伝子名が明確に多い → 転置して行=細胞にする
    if row_gene > 0.30 and row_gene > col_gene + 0.10:
        info['reason'] = '行に遺伝子名が多く検出されたため転置（行=遺伝子 → 行=細胞）'
        return True, info

    # 列に遺伝子名が明確に多い → そのまま使う
    if col_gene > 0.30 and col_gene > row_gene + 0.10:
        info['reason'] = '列に遺伝子名が多く検出されたため転置せず'
        return False, info

    # バーコードからも推定できる場合
    if col_bc > 0.50 and row_bc < 0.30:
        info['reason'] = '列に細胞バーコードが検出されたため転置（列=細胞 → 行=細胞）'
        return True, info
    if row_bc > 0.50 and col_bc < 0.30:
        info['reason'] = '行に細胞バーコードが検出されたため転置せず'
        return False, info

    # 明確な手掛かりが得られない場合は安全側（既に細胞×遺伝子の並びと仮定）にフォールバック。
    # 多くの前処理済み出力はこの並びになっているため、デフォルトとして妥当。
    info['reason'] = '明確な手掛かりが得られなかったため既定（cells × genes）として扱い、転置せず'
    return False, info


def detect_seqgeq_format(csv_path):
    """SeqGeq 形式の CSV か（[Metadata] と [Data] セクションを持つか）を判定する。

    Parameters
    ----------
    csv_path : str
        判定対象 CSV のパス。

    Returns
    -------
    int or None
        SeqGeq 形式の場合は、本体データの直前までスキップすべき行数を返す。
        SeqGeq 形式でない場合は None。
    """
    try:
        with open(csv_path, 'r') as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.strip())
                # 最初の 20 行だけ見れば判定には十分。
                if i >= 20:
                    break

            if len(lines) > 0 and lines[0] == '[Metadata]':
                for j, line in enumerate(lines):
                    if line == '[Data]':
                        # [Data] 行を含めてスキップする行数を返す
                        return j + 1
        return None
    except Exception:
        return None


def csv_to_h5mu(csv_path, output_path, transpose='auto', has_header=True, has_index=True):
    """CSV ファイルを h5mu 形式に変換する。

    Parameters
    ----------
    csv_path : str
        入力 CSV ファイルのパス。
    output_path : str
        出力 h5mu ファイルのパス。
    transpose : 'auto' | bool
        ``'auto'`` （既定）の場合、行・列のラベルから遺伝子軸を自動判定する。
        ``True`` を渡すと強制的に転置、``False`` を渡すと転置しない。
    has_header : bool
        True なら 1 行目をヘッダー（細胞名）として扱う。
    has_index : bool
        True なら 1 列目をインデックス（遺伝子名）として扱う。

    Returns
    -------
    dict
        以下のキーを持つ辞書:
          - ``output_path``: 書き出した h5mu のパス
          - ``transpose_applied``: 実際に転置を適用したか
          - ``auto_detected``: 自動判定モードだったか
          - ``detect_info``: 自動判定時の根拠スコアと理由
    """
    try:
        # SeqGeq 形式ならメタデータ行をスキップする
        skiprows = detect_seqgeq_format(csv_path)
        if skiprows:
            print(f"Detected SeqGeq format, skipping {skiprows} metadata rows")

        print(f"Reading CSV file: {csv_path}")
        if has_index:
            df = pd.read_csv(csv_path, index_col=0, skiprows=skiprows)
        else:
            df = pd.read_csv(csv_path, skiprows=skiprows)

        print(f"CSV shape: {df.shape}")

        # 転置の要否を決定する
        auto_detected = False
        detect_info = {}
        if transpose == 'auto' or transpose is None:
            auto_detected = True
            transpose_needed, detect_info = detect_transpose_needed(df)
            print(
                "自動転置判定: "
                f"transpose={transpose_needed} ("
                f"row_gene={detect_info['row_gene_score']}, "
                f"col_gene={detect_info['col_gene_score']}, "
                f"row_bc={detect_info['row_barcode_score']}, "
                f"col_bc={detect_info['col_barcode_score']})"
            )
            print(f"  理由: {detect_info['reason']}")
        else:
            transpose_needed = bool(transpose)

        if transpose_needed:
            df = df.T
            print(f"Transposed to: {df.shape}")

        # AnnData オブジェクトを作成（細胞 × 遺伝子の並び）
        adata = ad.AnnData(X=df.values.astype(np.float32))
        adata.obs_names = df.index.astype(str)
        adata.var_names = df.columns.astype(str)

        # 基本的なメタデータを追加
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # 1 モダリティの MuData として包む
        mdata = md.MuData({'rna': adata})

        print(f"Saving to: {output_path}")
        # gzip 圧縮でファイルサイズを抑える。古い mudata は compression
        # 引数を取らないので、その場合は無圧縮で書き出す。
        try:
            mdata.write(output_path, compression='gzip')
        except TypeError:
            mdata.write(output_path)

        print("Conversion completed successfully!")

        return {
            'output_path': output_path,
            'transpose_applied': bool(transpose_needed),
            'auto_detected': auto_detected,
            'detect_info': detect_info,
        }

    except Exception as e:
        raise Exception(f"Error converting CSV to h5mu: {str(e)}")


def validate_csv(csv_path):
    """CSV ファイルの形式を簡易チェックし、概要情報を返す。"""
    try:
        skiprows = detect_seqgeq_format(csv_path)
        df_preview = pd.read_csv(csv_path, nrows=5, skiprows=skiprows)
        info = {
            'rows': len(df_preview),
            'columns': len(df_preview.columns),
            'preview_columns': list(df_preview.columns[:10]),
            'has_numeric': df_preview.select_dtypes(include=[np.number]).shape[1] > 0,
        }
        return info
    except Exception as e:
        raise Exception(f"Error validating CSV: {str(e)}")
