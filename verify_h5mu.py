"""
Verify the converted h5mu file
"""
import mudata as md

h5mu_path = "/home/bdkoki/デスクトップ/csv_h5mu/SeqGeq_demo.h5mu"

print("Loading h5mu file...")
mdata = md.read_h5mu(h5mu_path)

print("\n=== MuData Summary ===")
print(f"Modalities: {list(mdata.mod.keys())}")
print(f"Number of observations (cells): {mdata.n_obs}")
print(f"Number of variables (genes): {mdata.n_vars}")

print("\n=== RNA Modality ===")
adata = mdata['rna']
print(f"Shape: {adata.shape}")
print(f"Observations (cells): {adata.n_obs}")
print(f"Variables (genes): {adata.n_vars}")

print("\n=== First 5 cell IDs ===")
print(adata.obs_names[:5].tolist())

print("\n=== First 5 gene names ===")
print(adata.var_names[:5].tolist())

print("\n=== Data matrix info ===")
print(f"Data type: {type(adata.X)}")
print(f"Data dtype: {adata.X.dtype}")
print(f"Data shape: {adata.X.shape}")

print("\n=== Observation metadata ===")
print(adata.obs.columns.tolist())
print(f"Mean counts per cell (first 5): {adata.obs['n_counts'].head().tolist()}")

print("\n=== Variable metadata ===")
print(adata.var.columns.tolist())
print(f"Number of cells expressing each gene (first 5): {adata.var['n_cells'].head().tolist()}")

print("\n✓ h5mu file is valid and can be read successfully!")
