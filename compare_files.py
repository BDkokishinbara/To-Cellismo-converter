"""
Compare input CSV and output h5mu to verify cell and gene counts match
"""
import pandas as pd
import mudata as md
import numpy as np

csv_path = "/home/bdkoki/デスクトップ/csv_h5mu/SeqGeq_demo.csv"
h5mu_path = "/home/bdkoki/デスクトップ/csv_h5mu/SeqGeq_demo.h5mu"

print("="*80)
print("Reading input CSV file (with skiprows=6 for SeqGeq format)...")
print("="*80)

# Read CSV with SeqGeq metadata rows skipped
df = pd.read_csv(csv_path, skiprows=6, index_col=0)

print(f"CSV shape: {df.shape}")
print(f"CSV rows (cells): {df.shape[0]}")
print(f"CSV columns (genes): {df.shape[1]}")
print(f"CSV index (cell IDs) - first 5: {df.index[:5].tolist()}")
print(f"CSV columns (gene names) - first 5: {df.columns[:5].tolist()}")

print("\n" + "="*80)
print("Reading output h5mu file...")
print("="*80)

mdata = md.read_h5mu(h5mu_path)
adata = mdata['rna']

print(f"h5mu shape: {adata.shape}")
print(f"h5mu rows (cells): {adata.n_obs}")
print(f"h5mu columns (genes): {adata.n_vars}")
print(f"h5mu obs_names (cell IDs) - first 5: {adata.obs_names[:5].tolist()}")
print(f"h5mu var_names (gene names) - first 5: {adata.var_names[:5].tolist()}")

print("\n" + "="*80)
print("COMPARISON RESULTS")
print("="*80)

# Compare dimensions
cells_match = df.shape[0] == adata.n_obs
genes_match = df.shape[1] == adata.n_vars

print(f"\n✓ Number of cells match: {cells_match}")
print(f"  - CSV: {df.shape[0]}")
print(f"  - h5mu: {adata.n_obs}")

print(f"\n✓ Number of genes match: {genes_match}")
print(f"  - CSV: {df.shape[1]}")
print(f"  - h5mu: {adata.n_vars}")

# Compare cell IDs
cell_ids_match = (df.index == adata.obs_names).all()
print(f"\n✓ Cell IDs match: {cell_ids_match}")
if not cell_ids_match:
    print(f"  First mismatched cell: CSV={df.index[0]}, h5mu={adata.obs_names[0]}")

# Compare gene names
gene_names_match = (df.columns == adata.var_names).all()
print(f"\n✓ Gene names match: {gene_names_match}")
if not gene_names_match:
    print(f"  First mismatched gene: CSV={df.columns[0]}, h5mu={adata.var_names[0]}")

# Compare data values (sample check - first 10x10)
print(f"\n✓ Checking data values (first 10x10 submatrix)...")
csv_sample = df.iloc[:10, :10].values.astype(np.float32)
h5mu_sample = adata.X[:10, :10]
data_match = np.allclose(csv_sample, h5mu_sample, rtol=1e-5, atol=1e-8)
print(f"  Data values match: {data_match}")

if data_match:
    print(f"  Sample CSV value [0,0]: {csv_sample[0,0]}")
    print(f"  Sample h5mu value [0,0]: {h5mu_sample[0,0]}")

print("\n" + "="*80)
if cells_match and genes_match and cell_ids_match and gene_names_match and data_match:
    print("✓✓✓ ALL CHECKS PASSED! Input and output are consistent. ✓✓✓")
else:
    print("⚠ WARNING: Some checks failed. Please review the differences above.")
print("="*80)
