"""
Verify the converted h5mu file from tar.gz MEX
"""
import mudata as md

h5mu_path = "/home/bdkoki/デスクトップ/csv_h5mu/outputs/convert_GSM4630028_ccRCC1.h5mu"

print("Verifying converted h5mu file...")
print(f"File: {h5mu_path}\n")

# Load h5mu file
mdata = md.read_h5mu(h5mu_path)

print("=== MuData Summary ===")
print(f"Modalities: {list(mdata.mod.keys())}")
print(f"Number of observations (cells): {mdata.n_obs}")
print(f"Number of variables (genes): {mdata.n_vars}")

print("\n=== RNA Modality ===")
adata = mdata['rna']
print(f"Shape: {adata.shape}")
print(f"Observations (cells): {adata.n_obs}")
print(f"Variables (genes): {adata.n_vars}")

print("\n=== First 5 cell barcodes ===")
print(adata.obs_names[:5].tolist())

print("\n=== First 5 gene names ===")
print(adata.var_names[:5].tolist())

print("\n=== Observation metadata columns ===")
print(adata.obs.columns.tolist())

print("\n=== Variable metadata columns ===")
print(adata.var.columns.tolist())

print("\n=== Data matrix info ===")
print(f"Data type: {type(adata.X)}")
print(f"Data shape: {adata.X.shape}")

print("\n✓ Verification complete! The h5mu file is valid.")
