import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import random

# --- Configuration ---
random.seed(42)
np.random.seed(42)

n_donors = 30
cells_per_donor = 20
n_genes = 20000
n_cells = n_donors * cells_per_donor  # 10,000 cells
n_cell_types = 10
k_neighbors = 15  # Standard default for k-NN

# --- 1. Generate Gene Expression Data (Matrix X) ---
print("Generating expression matrix...")
X = np.random.normal(loc=0.0, scale=1.0, size=(n_cells, n_genes))

# --- 2. Generate Metadata (Donor & Pheno) ---
print("Generating metadata...")
donor_ids = [f"Donor_{i}" for i in range(n_donors)]
cell_donor_ids = np.repeat(donor_ids, cells_per_donor)

# Phenotypes (Normal dist, mean 0, var 1) distributed across donors
donor_pheno_values = np.random.normal(loc=0.0, scale=1.0, size=n_donors)
donor_to_pheno_map = dict(zip(donor_ids, donor_pheno_values))
cell_phenotypes = [donor_to_pheno_map[d] for d in cell_donor_ids]

# Define 10 distinct cell type labels
cell_type_labels = [f"CellType_{i}" for i in range(n_cell_types)]
cell_assigned_types = np.random.choice(cell_type_labels, size=n_cells)

# --- 3. Construct AnnData ---
obs = pd.DataFrame({'Donor': cell_donor_ids, 'Pheno': cell_phenotypes, 'CellType': cell_assigned_types})
obs.index = [f"Cell_{i}" for i in range(n_cells)]
var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])

adata = ad.AnnData(X=X, obs=obs, var=var)

# --- 4. Generate k-NN Graph ---
print("Computing PCA...")
# Reduce dimensions first for efficiency (standard workflow)
# We use n_comps=50 as a standard default
sc.pp.pca(adata, n_comps=50)

print("Computing Neighbors...")
# Compute the k-NN graph
# n_neighbors sets the 'k' in k-NN
sc.pp.neighbors(adata, n_neighbors=k_neighbors, n_pcs=50)

# --- 5. Map to Requested Location ---
# Modern Scanpy stores connectivities in adata.obsp['connectivities'].
# We manually map it to adata.uns['neighbors']['connectivities'] to match your specific request.
adata.uns["neighbors"]["connectivities"] = adata.obsp["connectivities"]

# --- 6. Verify ---
print("\nData Summary:")
print(adata)
print("\nCheck connectivity matrix location:")
print(f"Graph shape in .uns: {adata.uns['neighbors']['connectivities'].shape}")

# --- 7. Save single cell data ---
adata.write_h5ad("./input/test_scdata.h5ad")

# --- 8. Generate test MAGMA data ---
df_magma = pd.DataFrame({'GENE': adata.var.index, 'ZSTAT': np.random.normal(size=n_genes)})
df_magma['ZSTAT'].values[0:1000] = np.random.normal(loc=10, size=1000)
df_magma.to_csv("./input/test_magma.txt", index=False, sep="\t")