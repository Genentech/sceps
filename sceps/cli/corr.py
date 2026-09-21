import sys, os, pickle, gc, random, time, argparse, logging

import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import scipy.stats

from src.utils import get_start_stop_index
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

EPS = 1e-16

estimand_col = ['OMEGA_GWAS', 'OMEGA_CONTROL', 'OMEGA_REST', 'OMEGA_OVERALL', 'OMEGA_DIFF']

def main():
   
    # Get command line input
    args = get_command_line()

    # Set seed for the random number generator
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Step if results already exists
    if os.path.exists(args.out):
        sys.exit()

    # Load scEPS results
    logging.info("Loading scEPS results from {}".format(args.sceps_result))
    df_sceps = pd.read_table(args.sceps_result)        

    # Load single-cell data
    adata = load_single_cell_data(args)

    # Intersect with scEPS results to include only shared cells
    shared_cell = np.intersect1d(adata.obs[args.cell_id_col], df_sceps['CELL'])
    adata = adata[adata.obs[args.cell_id_col].isin(shared_cell)]
    df_sceps = df_sceps[df_sceps['CELL'].isin(shared_cell)]

    # Define neighborhood clusters if number of bootstrap samples is non-zero
    if args.neighborhood_clusters != '':
        logging.info("Using pre-computed clusters of cell neighborhoods from {}".format(args.neighborhood_clusters))
        df_neighborhood_cluster = pd.read_table(args.neighborhood_clusters)
        cell2cluster = dict(zip(df_neighborhood_cluster[args.cell_id_col], df_neighborhood_cluster[args.block_bootstrap]))
        df_sceps['sceps.neighborhood_cluster'] = df_sceps['CELL'].map(cell2cluster)
    
    # Calculate the the start and stop indices of the cells to analyze
    start_idx, stop_idx = get_start_stop_index(args, adata.shape[1])

    # Calculate and save the correlation
    df_sceps = df_sceps.set_index('CELL')
    df_sceps = df_sceps.loc[adata.obs[args.cell_id_col]]
    df_out = calc_corr_all_gene(args, df_sceps, adata, start_idx=start_idx, stop_idx=stop_idx,
        min_num_nonzero=args.min_num_nonzero)
    df_out.to_csv('{}.{}-{}.txt.gz'.format(args.out, start_idx, stop_idx),
        sep='\t', index=False, float_format='%.5g')


# load single cell data
def load_single_cell_data(args):

    # load h5ad file
    adata = sc.read_h5ad(args.adata)
        
    # If args.cell_id_col not in adata.obs.columns, we use adata.obs.index as replacement
    obs = adata.obs.copy()
    if args.cell_id_col not in obs.columns:
        if args.cell_id_col == '':
            args.cell_id_col = 'sceps.cell_index'
        obs[args.cell_id_col] = obs.index.copy()
        adata.obs = obs
    
    # Check if gene ID column is already in adata.var.columns
    if args.gene_id_col not in adata.var.columns:
        if args.gene_id_col == '':
            args.gene_id_col = 'sceps.gene_index'
        adata.var[args.gene_id_col] = adata.var.index.copy()

    # select the focal cell type, if specified
    if args.cell_type_col != '':
        adata = adata[adata.obs[args.cell_type_col]==args.focal_cell_type]

    # select cells in the focal_cells file
    if args.focal_cells != '':
        df_focal_cells = pd.read_table(args.focal_cells)
        focal_cells = set(df_focal_cells[args.cell_id_col].values.tolist())
        adata = adata[adata.obs[args.cell_id_col].isin(focal_cells)]

    # Apply log normalization if specified
    if args.lognorm_count == True:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    return adata


# calculating correlation between sceps statistics and the expression of all gene
def calc_corr_all_gene(args, df_sceps, adata, start_idx=None, stop_idx=None,
    nonzero_other=True, min_num_nonzero=100):

    # Set the block bootstrap column
    block = None
    if (args.block_bootstrap != '') and (args.block_bootstrap in df_sceps.columns):
        block = df_sceps[args.block_bootstrap].values

    # Set start_idx and stop_idx
    if start_idx is None or stop_idx is None:
        start_idx, stop_idx = 0, adata.shape[1]

    # Extract the genes
    X = adata.X[:, start_idx:stop_idx]
    gene_list = adata.var[args.gene_id_col].values[start_idx:stop_idx]
    num_genes = gene_list.shape[0]
    
    # Iterate through the genes
    df_out = []
    for i in tqdm(range(num_genes)):
        vec_expr = X[:,i].toarray().flatten()
        if (nonzero_other == True) and (np.sum(vec_expr != 0) < min_num_nonzero):
            continue
        df_row = [gene_list[i]]
        for col in estimand_col:
            vec_sceps = df_sceps[col].values
            corr, se_corr, z_corr, p_corr = calc_corr(vec_sceps, vec_expr, block, nbs=args.num_bootstrap)
            df_row = df_row + [corr, se_corr, z_corr, p_corr]
        df_out.append(pd.DataFrame(df_row).transpose())

    # Create the output data frame
    df_out = pd.concat(df_out, ignore_index=True)

    # Rename the columns
    out_cols = ['GENE']
    for col in estimand_col:
        out_cols.append('R_{}'.format(col))
        out_cols.append('SE_R_{}'.format(col))
        out_cols.append('Z_R_{}'.format(col))
        out_cols.append('P_R_{}'.format(col))
    df_out.columns = out_cols

    # Remove any nan entries
    df_out = df_out.dropna(axis=1, how='all').reset_index(drop=True)

    return df_out


# calculating correlation between 2 vectors using bootstrap for testing
def calc_corr(vec_sceps_, vec_other_, block_, corr_type='spearman', nonzero_other=True, nbs=1000):

    nentry = vec_sceps_.shape[0]
    use_idx = np.array(range(nentry))
    if nonzero_other == True:
        use_idx = np.where(vec_other_ != 0)[0]
    vec_sceps = vec_sceps_[use_idx]
    vec_other = vec_other_[use_idx]

    if corr_type == 'pearson':
        corr = np.corrcoef(vec_sceps, vec_other)[0,1] 
    else:
        corr = scipy.stats.spearmanr(vec_sceps, vec_other).statistic
    all_corr_bs = []

    # block bootstrap based on the block
    if block_ is not None:
        block = block_[use_idx]
        all_block = np.array(np.unique(block))
        nblock = all_block.shape[0]
        for _ in range(nbs):
            use_block = np.random.choice(all_block, size=nblock, replace=True)
            use_idx = np.concatenate([np.where(block == blk)[0] for blk in use_block])
            if corr_type == 'pearson':
                corr_bs = np.corrcoef(vec_sceps[use_idx], vec_other[use_idx])[0,1] 
            else:
                corr_bs = scipy.stats.spearmanr(vec_sceps[use_idx], vec_other[use_idx]).statistic
            all_corr_bs.append(corr_bs)
    
    # ordinary bootstrap
    else:
        ncell = vec_sceps.shape[0]
        all_index = np.array(range(ncell))
        for _ in range(nbs):
            use_idx = np.random.choice(all_index, size=ncell, replace=True)
            if corr_type == 'pearson':
                corr_bs = np.corrcoef(vec_sceps[use_idx], vec_other[use_idx])[0,1] 
            else:
                corr_bs = scipy.stats.spearmanr(vec_sceps[use_idx], vec_other[use_idx]).statistic
            all_corr_bs.append(corr_bs)

    # get test statistics
    se_corr, z_corr, p_corr = np.nan, np.nan, np.nan
    if len(all_corr_bs) > 0:
        se_corr = np.std(all_corr_bs)
        z_corr = corr / (se_corr + EPS)
        p_corr = (1-scipy.stats.norm.cdf(np.fabs(z_corr)))*2.0

    return corr, se_corr, z_corr, p_corr


def get_command_line():
 
    # Create the parser
    parser = argparse.ArgumentParser(description="""This tool calculates the correlation between """ \
    """scEPS statistics and gene expression across cells.""")


    # input-related command line argument
    parser.add_argument('--sceps-result', type=str, required=False,
        help="""Used to specify a text file containing the data frame for estimates of scEPS """ \
        """statistics across all cell neighborhoods.""")

    parser.add_argument('--adata', type=str, required=False, default='',
        help="""Used to specify the path to the single-cell data containing log-noramlized gene expression matrix.""")

    parser.add_argument('--cell-id-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents cell IDs in the adata.obs data frame """ \
        """of the single-cell data. If left empty, scEPS will use what's in adata.obs.index as cell IDs.""")
    
    parser.add_argument('--gene-id-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents gene symbols/IDs in the adata.var data """ \
        """frame of the single-cell data. If left empty, scEPS will use what's in adata.var.index as gene symbols/IDs.""")
    
    parser.add_argument('--lognorm-count', default=False, required=False, action='store_true',
        help="""If specified, the tool will apply log-normalization on the single-cell data. This flag should """ \
        """only be specified if the data contains raw read count.""")

    parser.add_argument('--cell-type-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents cell types in the adata.obs data frame of """ \
        """the single-cell data. This is an emtpy string by default.""")

    parser.add_argument('--focal-cell-type', type=str, required=False, default='',
        help="""Used to specify the focal cell type to analyze. If specified, scEPS will only analyze cell """ \
        """neighborhoods in the specified cell type. This is an emtpy string by default.""")
    
    parser.add_argument('--focal-cells', type=str, required=False, default='',
        help="""Used to specify a text file containing a list of focal cells. """ \
        """If specified, the tool will calculate correlations across cells in this files.""")


    # paralleilization related argument
    parser.add_argument('--total-num-job', type=int, required=False, default=None,
        help="""Used to specify the total number of parallel jobs used for the scEPS analysis.""")
    
    parser.add_argument('--job-idx', type=int, required=False, default=None,
        help="""Used to specify the index of the parallel job for analyzing a subset of the single-cell data.""")
    
    parser.add_argument('--start-idx', type=int, required=False, default=None,
        help="""Used to specify the starting index (inclusive) of the cell neighborhood to analyze.""")
    
    parser.add_argument('--stop-idx', type=int, required=False, default=None,
        help="""Used to specify the stopping index (exclusive) of the cell neighborhood to analyze.""")


    # testing related command line argument
    parser.add_argument('--corr-type', type=str, required=False, default='spearman',
        help="""Used to specify the correlation type. The user can choose from "spearman" and "pearson".""")

    parser.add_argument('--min-num-nonzero', type=int, required=False, default=100,
        help="""Used to specify the minimum number of cells the gene needs to be expressed in (100 by default).""")

    parser.add_argument('--neighborhood-clusters', type=str, required=False, default='',
        help="""Used to specify the text file from step 2, representing a pre-computed mapping of cell """ \
        """neighborhoods to approximately independent blocks of cell neighborhoods.""")
    
    parser.add_argument('--num-bootstrap', type=int, required=False, default=0,
        help="""Used to specify the number of bootstrap samples (0 by default, not calculating test statistics).""")

    parser.add_argument('--block-bootstrap', type=str, required=False, default='sceps.neighborhood_cluster',
        help="""Used to specify the column that represents approximately independent cell neighborhood blocks """ \
        """(sceps.neighborhood_cluster by default). If left empty, the tool will fall back to regular bootstrap """ \
        """across individual (instead of blocks of) cell neighborhoods.""")
    
    parser.add_argument('--seed', type=int, required=False, default=0,
        help="""Used to specify the seed for the random number generator (0 by default).""")


    # output related command line argument
    parser.add_argument('--out', type=str, required=False,
        help="""Used to specify the prefix of the output file name.""")


    # Execute the parse_args() method
    args = parser.parse_args()
    
    return args


if __name__ == "__main__":
    main()
