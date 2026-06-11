import scanpy as sc
import scanpy.external as sce
import numpy as np

# Load single-cell data in h5ad format
in_fnm = '<input single-cell data file name>'
adata = sc.read_h5ad(in_fnm)

# Remove batches with few cells
batch_nm = '<batch variable name>'
batch_cell_count = adata.obs[batch_nm].value_counts()
keep_batch = batch_cell_count[batch_cell_count>=3].index.values.tolist()
adata = adata[adata.obs[batch_nm].isin(keep_batch)]

# Perform QC and log-normalization
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Calculate PCA
sc.pp.highly_variable_genes(adata, n_top_genes=2000, inplace=True)
sc.pp.pca(adata, use_highly_variable=True) 

# Integrate cells in the PCA space
sce.pp.harmony_integrate(adata, batch_nm, basis='X_pca', max_iter_harmony=50)
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# Save single-cell data
out_fnm = '<output single-cell data file name>'        
adata.write_h5ad(out_fnm)

# Regress out batch
sc.pp.regress_out(adata, keys=[batch_nm])
adata.X = adata.X.astype(np.float32)

# Save regressed out single-cell data
out_fnm_regressout = '<output regressed out single-cell data file name>'        
adata.write_h5ad(out_fnm_regressout)

# Save single-cell data without X
del adata.X
out_fnm_noX = '<output single-cell data (without adata.X) file name>'        
adata.write_h5ad(out_fnm_noX)
