import sys
sys.path.append('/path/to/sceps_tool')

from src.sceps_core import *

# Create args with required arguments and default settings to run scEPS
args = create_sceps_default_args(adata='/path/to/the/single/cell/data',
                                gene_list='/path/to/the/gene/level/MAGMA/association/statistics')

# Specify additional settings as needed
args.cell_id_col = '<cell ID column in adata.obs>'
args.gene_id_col = '<gene ID column in adata.var>'
args.donor_id_col = '<donor ID column in adata.obs>'
args.pheno_file = '/path/to/an/external/phenotype/file'
args.pheno = '<phenotype column in adata.obs>'
args.covar = '<a list of covariates to adjust, separated by commas>'
args.out = '<output file name>'

# Set seed for random number generator
np.random.seed(args.seed)
random.seed(args.seed)

# Load the single-cell data
adata = load_single_cell_data(args)

# Load the list of genes with MAGMA Z-stat
gene_list = load_gene_stats(args)

# Intersect the genes in the adata with genes with MAGMA Z-stat
adata = subset_adata_genes(args, adata, gene_list)

# Determine the list of GWAS genes
gwas_genes = select_disease_gwas_genes(args, gene_list)

# Run scEPS for the first 10 cells in the single-cell data
start_idx, stop_idx = 0, 10
sceps_out = get_sceps_stats(args, adata, gwas_genes, start_idx=start_idx, stop_idx=stop_idx)

# Save results
out_suffix = '{}-{}'.format(start_idx, stop_idx)
save_sceps_stats(args, sceps_out, out_suffix)