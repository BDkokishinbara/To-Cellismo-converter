"""
h5ad to h5mu converter for single-cell analysis data
Wraps an AnnData file as a single 'rna' modality MuData and writes h5mu.
"""
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd


# CXG / standard h5ad columns that typically hold the human-readable gene
# symbol when var.index is an Ensembl ID.
_GENE_SYMBOL_COLUMNS = ('feature_name', 'gene_symbol', 'gene_symbols', 'gene_name', 'gene_names', 'Symbol')


def _remap_one_var_frame(var_df):
    """Return (new_index, source_col) if a remap is possible, else (None, None).

    ``var_df`` is either ``adata.var`` or ``adata.raw.var``. The function does
    NOT mutate the frame; the caller installs the new index.
    """
    sample = var_df.index[:50].astype(str).tolist()
    looks_like_ensembl = sum(1 for s in sample if s.startswith('ENS')) >= len(sample) * 0.5
    if not looks_like_ensembl:
        return None, None

    for col in _GENE_SYMBOL_COLUMNS:
        if col in var_df.columns:
            symbols = var_df[col].astype(str)
            ensembl_strs = var_df.index.astype(str)
            blank = symbols.isin(['', 'nan', 'None', 'NaN'])
            symbols = symbols.where(~blank, pd.Series(ensembl_strs, index=var_df.index))
            return pd.Index(symbols, name='gene'), col
    return None, None


def _remap_var_index_to_symbols(adata):
    """Swap ``var.index`` (and ``raw.var.index`` if any) from Ensembl IDs to
    gene symbols so downstream tools like Cellismo display readable names.

    Original Ensembl IDs are preserved as the ``ensembl_id`` column on both
    ``var`` and ``raw.var``. No-op when the index already looks like symbols.
    Returns (remapped: bool, source_column: str | None).
    """
    # Main var
    new_index, src = _remap_one_var_frame(adata.var)
    remapped = False
    if new_index is not None:
        if 'ensembl_id' not in adata.var.columns:
            adata.var['ensembl_id'] = adata.var.index.astype(str)
        adata.var.index = new_index
        adata.var_names_make_unique()
        remapped = True

    # Raw var (if present). CXG h5ad files keep the un-normalized counts in
    # adata.raw and its var.index is also Ensembl IDs — Cellismo may display
    # these instead of adata.var_names.
    if getattr(adata, 'raw', None) is not None:
        raw_new_index, raw_src = _remap_one_var_frame(adata.raw.var)
        if raw_new_index is not None:
            # adata.raw is read-only; rebuild it by hand-promoting to AnnData.
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
    """
    Convert h5ad (AnnData) file to h5mu (MuData) format.

    Parameters
    ----------
    h5ad_path : str
        Path to input .h5ad file.
    output_path : str
        Path to output .h5mu file.
    modality : str
        Modality name under which the AnnData is stored in the MuData
        container. Defaults to 'rna' for consistency with the other
        converters in this project.

    Returns
    -------
    str
        Path to the created h5mu file.
    """
    try:
        print(f"Reading h5ad file: {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path)
        print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # Make sure var.index is human-readable gene symbols when possible.
        remapped, source_col = _remap_var_index_to_symbols(adata)
        if remapped:
            print(f"Remapped var_names to gene symbols using column '{source_col}' (Ensembl IDs preserved in 'ensembl_id')")

        if 'n_counts' not in adata.obs.columns:
            adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        if 'n_cells' not in adata.var.columns:
            adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        mdata = md.MuData({modality: adata})

        print(f"Saving to: {output_path}")
        # Use gzip compression to keep file size close to the original h5ad.
        # Without this, sparse h5ad inputs can balloon several-fold when
        # written back out as h5mu.
        try:
            mdata.write(output_path, compression='gzip')
        except TypeError:
            # Older mudata versions may not accept compression kwarg.
            mdata.write(output_path)
        print("Conversion completed successfully!")
        return output_path

    except Exception as e:
        raise Exception(f"Error converting h5ad to h5mu: {str(e)}")


def validate_h5ad(h5ad_path):
    """Return a small summary dict describing the h5ad file."""
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
