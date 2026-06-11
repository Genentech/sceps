import sys, os, pickle, gc, random, time, warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import scipy
import scipy.stats

from src.sceps_core import *
from tqdm import tqdm

def main():
   
    # Get command line input
    args = get_command_line()

    # Convert args to a DataFrame (saved after command finishes)
    df_args = pd.DataFrame.from_dict(args.__dict__, orient='index', columns=['VALUE']).reset_index()
    df_args.rename(columns={'index':'ARGUMENT'}, inplace=True)
    df_args = df_args[['ARGUMENT', 'VALUE']]

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

    # Load the set of focal cell neighborhoods to analyze (this is empty, if not specified in args)
    focal_cells = load_focal_cells(args)

    # Calculate the the start and stop indices of the cells to analyze
    start_idx, stop_idx = get_start_stop_index(args, adata.shape[0])

    # Create the output prefix based on start_idx and stop_idx
    out_suffix = '{}-{}'.format(start_idx, stop_idx)

    # Check if scEPS results already exist
    if (check_output(args, out_suffix) == True) and (args.force_rerun == False):
        logging.info('Output already exists')
        sys.exit()

    # Estimate scEPS statistics
    sceps_out = get_sceps_stats(args, adata, gwas_genes, start_idx, stop_idx, focal_cells=focal_cells)

    # Save scEPS output
    save_sceps_stats(args, sceps_out, out_suffix)

    # Save the arguments as a DataFrame in a log file
    df_args.to_csv(args.out+'.{}.log'.format(out_suffix), sep='\t', na_rep='NA', index=False)


if __name__ == "__main__":
    main()
