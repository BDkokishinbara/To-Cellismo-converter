"""
MEX (10x Genomics) to h5mu converter for single-cell analysis data
"""
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from scipy.io import mmread
from pathlib import Path
import gzip


def mex_to_h5mu(mex_dir, output_path, genome=None):
    """
    Convert MEX format (10x Genomics) to h5mu format

    MEX format consists of:
    - matrix.mtx (or matrix.mtx.gz): sparse matrix in Matrix Market format
    - genes.tsv (or features.tsv.gz): gene/feature information
    - barcodes.tsv (or barcodes.tsv.gz): cell barcode information

    Parameters:
    -----------
    mex_dir : str
        Path to directory containing MEX files
    output_path : str
        Path to output h5mu file
    genome : str, optional
        Genome name (e.g., 'GRCh38') if multiple genomes present

    Returns:
    --------
    str : Path to the created h5mu file
    """
    try:
        print(f"Reading MEX files from: {mex_dir}")
        mex_path = Path(mex_dir)

        # Find matrix file
        matrix_file = find_file(mex_path, ['matrix.mtx.gz', 'matrix.mtx'])
        if not matrix_file:
            raise FileNotFoundError("Matrix file (matrix.mtx or matrix.mtx.gz) not found")

        # Find features/genes file
        features_file = find_file(mex_path, ['features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv'])
        if not features_file:
            raise FileNotFoundError("Features/genes file not found")

        # Find barcodes file
        barcodes_file = find_file(mex_path, ['barcodes.tsv.gz', 'barcodes.tsv'])
        if not barcodes_file:
            raise FileNotFoundError("Barcodes file not found")

        print(f"Found files:")
        print(f"  Matrix: {matrix_file.name}")
        print(f"  Features: {features_file.name}")
        print(f"  Barcodes: {barcodes_file.name}")

        # Read matrix
        print("Reading matrix...")
        matrix = mmread(matrix_file)
        # Convert to CSR format for efficient row operations
        matrix = matrix.tocsr()
        print(f"Matrix shape: {matrix.shape}")

        # Read features/genes
        print("Reading features...")
        features_df = read_tsv(features_file)

        # Handle different feature file formats
        if features_df.shape[1] >= 3:
            # 10x Genomics v3 format: gene_id, gene_name, feature_type
            features_df.columns = ['gene_id', 'gene_name', 'feature_type'][:features_df.shape[1]]
        elif features_df.shape[1] == 2:
            # 10x Genomics v2 format: gene_id, gene_name
            features_df.columns = ['gene_id', 'gene_name']
        else:
            # Single column: gene_id only
            features_df.columns = ['gene_id']

        # Use gene_name as index if available, otherwise gene_id
        if 'gene_name' in features_df.columns:
            var_names = features_df['gene_name'].astype(str).values
        else:
            var_names = features_df['gene_id'].astype(str).values

        # Read barcodes
        print("Reading barcodes...")
        barcodes_df = read_tsv(barcodes_file, header=None)
        barcodes = barcodes_df[0].astype(str).values

        # Create AnnData object
        print("Creating AnnData object...")
        # Matrix is genes x cells, need to transpose to cells x genes
        adata = ad.AnnData(X=matrix.T)

        # Set observation (cell) names
        adata.obs_names = barcodes

        # Set variable (gene) names
        adata.var_names = var_names

        # Add feature metadata
        for col in features_df.columns:
            adata.var[col] = features_df[col].values

        # Add basic QC metrics
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
        adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # Filter by feature type if multiple types present
        if 'feature_type' in adata.var.columns:
            feature_types = adata.var['feature_type'].unique()
            print(f"Feature types found: {feature_types}")

            # Create separate modalities for different feature types
            modalities = {}

            for feat_type in feature_types:
                mask = adata.var['feature_type'] == feat_type
                adata_subset = adata[:, mask].copy()

                # Map feature type to modality name
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

            # Create MuData with multiple modalities
            mdata = md.MuData(modalities)
        else:
            # Single modality (RNA)
            mdata = md.MuData({'rna': adata})

        # Save as h5mu
        print(f"Saving to: {output_path}")
        mdata.write(output_path)

        print("Conversion completed successfully!")
        return output_path

    except Exception as e:
        raise Exception(f"Error converting MEX to h5mu: {str(e)}")


def find_file(directory, filenames):
    """
    Find first matching file from a list of possible filenames

    Parameters:
    -----------
    directory : Path
        Directory to search in
    filenames : list
        List of possible filenames

    Returns:
    --------
    Path or None : Path to found file or None
    """
    for filename in filenames:
        filepath = directory / filename
        if filepath.exists():
            return filepath
    return None


def read_tsv(filepath, header=None):
    """
    Read TSV file (handles both plain and gzipped files)

    Parameters:
    -----------
    filepath : Path
        Path to TSV file
    header : int or None
        Row number to use as header

    Returns:
    --------
    DataFrame : Loaded dataframe
    """
    if filepath.suffix == '.gz':
        with gzip.open(filepath, 'rt') as f:
            return pd.read_csv(f, sep='\t', header=header)
    else:
        return pd.read_csv(filepath, sep='\t', header=header)


def validate_mex_directory(mex_dir):
    """
    Validate MEX directory structure

    Parameters:
    -----------
    mex_dir : str
        Path to directory containing MEX files

    Returns:
    --------
    dict : Information about the MEX files
    """
    try:
        mex_path = Path(mex_dir)

        if not mex_path.exists():
            raise FileNotFoundError(f"Directory not found: {mex_dir}")

        if not mex_path.is_dir():
            raise NotADirectoryError(f"Not a directory: {mex_dir}")

        # Check for required files
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
