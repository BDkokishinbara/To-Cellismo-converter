"""
RDS to h5mu converter for single-cell analysis data
Supports Seurat and SingleCellExperiment objects
"""
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from pathlib import Path

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
    """
    Convert RDS file (Seurat or SingleCellExperiment) to h5mu format

    Parameters:
    -----------
    rds_path : str
        Path to input RDS file
    output_path : str
        Path to output h5mu file
    object_type : str
        Type of R object: 'seurat', 'sce', or 'auto' (default)

    Returns:
    --------
    str : Path to the created h5mu file
    """
    if not RPY2_AVAILABLE:
        raise ImportError(
            "rpy2 is not installed. Install it with: pip install rpy2\n"
            "Also make sure R is installed on your system."
        )

    try:
        print(f"Reading RDS file: {rds_path}")

        # Load R object
        readRDS = robjects.r['readRDS']
        r_object = readRDS(rds_path)

        # Detect object type
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

        # Convert based on object type
        if object_type == 'seurat':
            adata = convert_seurat_to_anndata(r_object)
        elif object_type == 'sce':
            adata = convert_sce_to_anndata(r_object)
        else:
            raise ValueError(f"Unknown object_type: {object_type}")

        print(f"Created AnnData object: {adata.shape[0]} cells x {adata.shape[1]} genes")

        # Create MuData object
        mdata = md.MuData({'rna': adata})

        # Save as h5mu
        print(f"Saving to: {output_path}")
        mdata.write(output_path)

        print("Conversion completed successfully!")
        return output_path

    except Exception as e:
        raise Exception(f"Error converting RDS to h5mu: {str(e)}")


def convert_seurat_to_anndata(seurat_obj):
    """
    Convert Seurat object to AnnData

    Parameters:
    -----------
    seurat_obj : rpy2 R object
        Seurat object from R

    Returns:
    --------
    AnnData : Converted AnnData object
    """
    try:
        # Try to use Seurat conversion functions
        base = importr('base')

        # Get count matrix
        # Try different Seurat versions
        try:
            # Seurat v5
            robjects.r('suppressMessages(library(Seurat))')
            counts = robjects.r('GetAssayData')(seurat_obj, slot='counts')
        except:
            try:
                # Seurat v3/v4
                counts = robjects.r('GetAssayData')(seurat_obj, assay='RNA', slot='counts')
            except:
                # Fallback
                counts = seurat_obj.slots['assays'].slots['RNA'].slots['counts']

        # Convert sparse matrix to dense (if needed)
        if robjects.r('inherits')(counts, 'dgCMatrix')[0]:
            counts = robjects.r('as.matrix')(counts)

        # Convert to numpy array
        counts_np = np.array(counts).T  # Transpose to cells x genes

        # Get cell metadata
        try:
            metadata = robjects.r('as.data.frame')(seurat_obj.slots['meta.data'])
            obs = pandas2ri.rpy2py(metadata)
        except:
            obs = pd.DataFrame(index=range(counts_np.shape[0]))

        # Get gene names
        try:
            var_names = list(robjects.r('rownames')(counts))
            var = pd.DataFrame(index=var_names)
        except:
            var = pd.DataFrame(index=range(counts_np.shape[1]))

        # Create AnnData object
        adata = ad.AnnData(X=counts_np, obs=obs, var=var)

        return adata

    except Exception as e:
        raise Exception(f"Error converting Seurat object: {str(e)}")


def convert_sce_to_anndata(sce_obj):
    """
    Convert SingleCellExperiment object to AnnData

    Parameters:
    -----------
    sce_obj : rpy2 R object
        SingleCellExperiment object from R

    Returns:
    --------
    AnnData : Converted AnnData object
    """
    try:
        # Import SingleCellExperiment functions
        robjects.r('suppressMessages(library(SingleCellExperiment))')

        # Get count matrix
        counts = robjects.r('counts')(sce_obj)

        # Convert to dense matrix if sparse
        if robjects.r('inherits')(counts, 'dgCMatrix')[0]:
            counts = robjects.r('as.matrix')(counts)

        # Convert to numpy array
        counts_np = np.array(counts).T  # Transpose to cells x genes

        # Get cell metadata (colData)
        try:
            col_data = robjects.r('colData')(sce_obj)
            col_data = robjects.r('as.data.frame')(col_data)
            obs = pandas2ri.rpy2py(col_data)
        except:
            obs = pd.DataFrame(index=range(counts_np.shape[0]))

        # Get gene metadata (rowData)
        try:
            row_data = robjects.r('rowData')(sce_obj)
            row_data = robjects.r('as.data.frame')(row_data)
            var = pandas2ri.rpy2py(row_data)
        except:
            # At least get gene names
            try:
                var_names = list(robjects.r('rownames')(sce_obj))
                var = pd.DataFrame(index=var_names)
            except:
                var = pd.DataFrame(index=range(counts_np.shape[1]))

        # Create AnnData object
        adata = ad.AnnData(X=counts_np, obs=obs, var=var)

        return adata

    except Exception as e:
        raise Exception(f"Error converting SingleCellExperiment object: {str(e)}")


def check_rpy2_available():
    """Check if rpy2 is available"""
    return RPY2_AVAILABLE
