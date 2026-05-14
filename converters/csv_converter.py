"""
CSV to h5mu converter for single-cell analysis data
"""
import re
import pandas as pd
import anndata as ad
import mudata as md
import numpy as np
from pathlib import Path


# ── Heuristics for auto-detecting which axis holds gene names ─────────────

# Patterns that look like gene symbols / IDs.
_GENE_PATTERNS = [
    # Uppercase human-style symbols: ACTB, CD3D, HLA-DRA, MT-ATP6, RPS4Y1, IFITM3
    re.compile(r'^[A-Z][A-Z0-9._-]{0,19}$'),
    # Ensembl IDs across species: ENSG..., ENSMUSG..., ENSDARG..., etc.
    re.compile(r'^ENS[A-Z]{0,3}[GTPR]\d{6,}'),
    # Mouse-style mixed case: Actb, Cd3d, Ms4a1, Hbb-b1
    re.compile(r'^[A-Z][a-z][a-zA-Z0-9.-]{1,18}$'),
    # MGI accession IDs
    re.compile(r'^MGI:\d+$'),
]

# Patterns that look like single-cell barcodes (i.e. the axis is cells, not genes).
_BARCODE_PATTERNS = [
    # 10x Genomics style barcodes, optionally with -1 sample suffix
    re.compile(r'^[ACGTN]{12,24}(-\d+)?$'),
    # Sample-prefixed barcodes: Sample1_AAACCTGAG..., Patient3-AAACCT...
    re.compile(r'^[A-Za-z0-9]+[_-][ACGTN]{12,24}(-\d+)?$'),
]


def _fraction_matching(labels, patterns):
    if not labels:
        return 0.0
    n = sum(1 for s in (str(x) for x in labels) if any(p.match(s) for p in patterns))
    return n / len(labels)


def detect_transpose_needed(df, sample_size=500):
    """
    Decide whether ``df`` should be transposed to reach the canonical
    (cells x genes) orientation expected by AnnData / scanpy / MuData.

    Returns
    -------
    (transpose : bool, info : dict)
        ``transpose`` is True when the input has genes laid out as rows and
        must be flipped so genes end up as columns (= ``var`` axis in AnnData).
        ``info`` carries the diagnostic scores so the UI can show the user
        why a decision was made.
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

    # Clear "genes in rows" signal: transpose to put cells in rows.
    if row_gene > 0.30 and row_gene > col_gene + 0.10:
        info['reason'] = '行に遺伝子名が多く検出されたため転置（行=遺伝子 → 行=細胞）'
        return True, info

    # Clear "genes in columns" signal: keep as-is.
    if col_gene > 0.30 and col_gene > row_gene + 0.10:
        info['reason'] = '列に遺伝子名が多く検出されたため転置せず'
        return False, info

    # Barcode signal can also disambiguate.
    if col_bc > 0.50 and row_bc < 0.30:
        info['reason'] = '列に細胞バーコードが検出されたため転置（列=細胞 → 行=細胞）'
        return True, info
    if row_bc > 0.50 and col_bc < 0.30:
        info['reason'] = '行に細胞バーコードが検出されたため転置せず'
        return False, info

    # Fallback: assume the file is already in (cells x genes) form. This
    # matches the most common pre-processed output and is the safer default.
    info['reason'] = '明確な手掛かりが得られなかったため既定（cells × genes）として扱い、転置せず'
    return False, info


def detect_seqgeq_format(csv_path):
    """
    Detect if CSV is in SeqGeq format (with [Metadata] and [Data] sections)

    Parameters
    ----------
    csv_path : str
        Path to CSV file

    Returns
    -------
    int or None
        Number of rows to skip if SeqGeq format, None otherwise
    """
    try:
        with open(csv_path, 'r') as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.strip())
                if i >= 20:
                    break

            if len(lines) > 0 and lines[0] == '[Metadata]':
                for j, line in enumerate(lines):
                    if line == '[Data]':
                        return j + 1
        return None
    except Exception:
        return None


def csv_to_h5mu(csv_path, output_path, transpose='auto', has_header=True, has_index=True):
    """
    Convert CSV file to h5mu format.

    Parameters
    ----------
    csv_path : str
        Path to input CSV file.
    output_path : str
        Path to output h5mu file.
    transpose : 'auto' | bool
        ``'auto'`` (default): auto-detect orientation from row/column labels.
        ``True`` forces a transpose, ``False`` forces no transpose.
    has_header : bool
        Treat first row as header (cell names) if True.
    has_index : bool
        Treat first column as index (gene names) if True.

    Returns
    -------
    dict
        Keys:
          - ``output_path``: written file path
          - ``transpose_applied``: bool of whether a transpose was applied
          - ``auto_detected``: bool, True when 'auto' was used
          - ``detect_info``: diagnostic info from detection (when auto)
    """
    try:
        skiprows = detect_seqgeq_format(csv_path)
        if skiprows:
            print(f"Detected SeqGeq format, skipping {skiprows} metadata rows")

        print(f"Reading CSV file: {csv_path}")
        if has_index:
            df = pd.read_csv(csv_path, index_col=0, skiprows=skiprows)
        else:
            df = pd.read_csv(csv_path, skiprows=skiprows)

        print(f"CSV shape: {df.shape}")

        # Decide transpose.
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

        # Create AnnData object (cells x genes)
        adata = ad.AnnData(X=df.values.astype(np.float32))
        adata.obs_names = df.index.astype(str)
        adata.var_names = df.columns.astype(str)

        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        mdata = md.MuData({'rna': adata})

        print(f"Saving to: {output_path}")
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
    """Validate CSV file and return information about it."""
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
