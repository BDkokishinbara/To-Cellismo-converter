"""
CSV to h5mu converter for single-cell analysis data
"""
import pandas as pd
import anndata as ad
import mudata as md
import numpy as np
from pathlib import Path


def detect_seqgeq_format(csv_path):
    """
    Detect if CSV is in SeqGeq format (with [Metadata] and [Data] sections)

    Parameters:
    -----------
    csv_path : str
        Path to CSV file

    Returns:
    --------
    int or None : Number of rows to skip if SeqGeq format, None otherwise
    """
    try:
        with open(csv_path, 'r') as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.strip())
                # Only check first 20 lines to avoid reading entire file
                if i >= 20:
                    break

            # Check for SeqGeq format
            if len(lines) > 0 and lines[0] == '[Metadata]':
                # Look for [Data] section
                for j, line in enumerate(lines):
                    if line == '[Data]':
                        # Skip up to and including [Data] line
                        return j + 1
        return None
    except Exception:
        return None


def csv_to_h5mu(csv_path, output_path, transpose=False, has_header=True, has_index=True):
    """
    Convert CSV file to h5mu format

    Parameters:
    -----------
    csv_path : str
        Path to input CSV file (genes x cells or cells x genes matrix)
    output_path : str
        Path to output h5mu file
    transpose : bool
        If True, transpose the matrix (useful if CSV is genes x cells)
    has_header : bool
        If True, first row is header (cell names)
    has_index : bool
        If True, first column is index (gene names)

    Returns:
    --------
    str : Path to the created h5mu file
    """
    try:
        # Check for SeqGeq format
        skiprows = detect_seqgeq_format(csv_path)
        if skiprows:
            print(f"Detected SeqGeq format, skipping {skiprows} metadata rows")

        # Read CSV file
        print(f"Reading CSV file: {csv_path}")

        if has_index:
            df = pd.read_csv(csv_path, index_col=0, skiprows=skiprows)
        else:
            df = pd.read_csv(csv_path, skiprows=skiprows)

        print(f"CSV shape: {df.shape}")

        # Transpose if needed (typically we want cells x genes)
        if transpose:
            df = df.T
            print(f"Transposed to: {df.shape}")

        # Create AnnData object
        # AnnData expects cells x genes format
        adata = ad.AnnData(X=df.values.astype(np.float32))

        # Set observations (cells) and variables (genes)
        adata.obs_names = df.index.astype(str)
        adata.var_names = df.columns.astype(str)

        # Add basic metadata
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # Create MuData object (can contain multiple modalities)
        mdata = md.MuData({'rna': adata})

        # Save as h5mu
        print(f"Saving to: {output_path}")
        mdata.write(output_path)

        print("Conversion completed successfully!")
        return output_path

    except Exception as e:
        raise Exception(f"Error converting CSV to h5mu: {str(e)}")


def validate_csv(csv_path):
    """
    Validate CSV file and return information about it

    Returns:
    --------
    dict : Information about the CSV file
    """
    try:
        # Check for SeqGeq format
        skiprows = detect_seqgeq_format(csv_path)

        # Read first few rows to check format
        df_preview = pd.read_csv(csv_path, nrows=5, skiprows=skiprows)

        info = {
            'rows': len(df_preview),
            'columns': len(df_preview.columns),
            'preview_columns': list(df_preview.columns[:10]),
            'has_numeric': df_preview.select_dtypes(include=[np.number]).shape[1] > 0
        }

        return info

    except Exception as e:
        raise Exception(f"Error validating CSV: {str(e)}")
