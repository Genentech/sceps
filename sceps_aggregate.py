import argparse, sys, glob, os
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
from tqdm import tqdm
import logging, random
import scanpy as sc

from src.utils import EPS
from src.stats_test import get_weighted_mean

from src.scdata import *
from sklearn.cluster import MiniBatchKMeans
import statsmodels.stats.multitest as smm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

estimand_col = ['OMEGA_GWAS', 'OMEGA_CONTROL', 'OMEGA_REST', 'OMEGA_OVERALL', 'OMEGA_DIFF']
info_col = ['NEIGHBORHOOD_SIZE', 'NUM_DONOR', \
            'MEAN_MEAN_EXPR_GWAS', 'MEAN_MEAN_EXPR_CONTROL', 'MEAN_MEAN_EXPR_REST', 'MEAN_MEAN_EXPR_ALL', \
            'MEAN_VAR_EXPR_GWAS', 'MEAN_VAR_EXPR_CONTROL', 'MEAN_VAR_EXPR_REST', 'MEAN_VAR_EXPR_ALL', \
            'VAR_OUTCOME', 'NUM_GENE_GWAS', 'NUM_GENE_CONTROL', 'NUM_GENE_REST']

def main():
    
    # get command line
    args = get_command_line()

    # set seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    # append prefix to the estimand columns if using weighted estimands
    if args.use_sigma:
        for i in range(len(estimand_col)):
            estimand_col[i] = estimand_col[i].replace('OMEGA', 'SIGMA')

    # load the single-cell data
    logging.info("Loading single-cell data from {}".format(args.adata))
    adata = sc.read_h5ad(args.adata)

    # If args.cell_id_col not in adata.obs.columns, we use adata.obs.index as replacement
    obs = adata.obs.copy()
    if args.cell_id_col not in obs.columns:
        if (args.cell_id_col is None) or (args.cell_id_col == ''):
            args.cell_id_col = 'sceps.cell_index'
        obs[args.cell_id_col] = obs.index.copy()
        adata.obs = obs

    # define neighborhood clusters for bootstrap
    if args.neighborhood_clusters != '':
        logging.info("Using pre-computed clusters of cell neighborhoods from {}".format(args.neighborhood_clusters))
        df_neighborhood_cluster = pd.read_table(args.neighborhood_clusters)
        adata = adata[adata.obs[args.cell_id_col].isin(df_neighborhood_cluster[args.cell_id_col])].copy()
        cell2cluster = dict(zip(df_neighborhood_cluster[args.cell_id_col], df_neighborhood_cluster[args.block_bootstrap]))
        adata.obs[args.block_bootstrap] = adata.obs[args.cell_id_col].map(cell2cluster)

    # aggregate results across batches of sceps runs
    logging.info("Aggregating scEPS results across batches of runs")
    df_sceps = aggregate_sceps_batches(args)
    adata = adata[adata.obs[args.cell_id_col].isin(df_sceps['CELL'])].copy()
    obs = adata.obs.copy()

    # perform testing at cell type level
    if args.cell_type_col != '':
        all_cell_type_col = args.cell_type_col.split(',')
        for cell_type_col in all_cell_type_col:

            logging.info("Testing at cell type level {}".format(cell_type_col))

            # create output file name
            if args.use_sigma == True:
                out_fnm = '{}.{}.sceps.sigma.celltype.txt'.format(args.out, cell_type_col)
            else:
                out_fnm = '{}.{}.sceps.omega.celltype.txt'.format(args.out, cell_type_col)
            
            # test all cell types in the cell_type_col column
            df_test_ct = test_all_cell_types(args, cell_type_col, obs, df_sceps)

            # save results
            df_test_ct.to_csv(out_fnm, sep='\t', index=False, float_format='%.5g')
    else:
        # no cell type column specified, aggregate across all cells
        if args.use_sigma == True:
            out_fnm = '{}.sceps.sigma.celltype.txt'.format(args.out)
        else:
            out_fnm = '{}.sceps.omega.celltype.txt'.format(args.out)

        # create a temporary cell type column
        obs['scEPS_tmp_cell_type_col'] = 'All'
        df_test_ct = test_all_cell_types(args, 'scEPS_tmp_cell_type_col', obs, df_sceps)

        # save results
        df_test_ct.to_csv(out_fnm, sep='\t', index=False, float_format='%.5g')


def aggregate_sceps_batches(args):
    """
    Aggregate scEPS results across batches of cells
    """

    # load the list of neighborhoods to use
    df_usenb = None
    if args.use_neighborhoods != '':
        df_usenb = pd.read_table(args.use_neighborhoods, header=None)

    # load sceps output files
    if args.prefix.endswith('.txt.gz') == False:
        if args.use_sigma == False:
            all_out_ = glob.glob('{}*sceps.omega.txt.gz'.format(args.prefix))
        else:
            all_out_ = glob.glob('{}*sceps.sigma.txt.gz'.format(args.prefix))
    else:
        all_out_ = glob.glob(args.prefix)
    
    # get output file name
    if args.use_sigma == False:
        out_fnm = '{}.sceps.omega.txt.gz'.format(args.out)
    else:
        out_fnm = '{}.sceps.sigma.txt.gz'.format(args.out)

    # exclude files that contains specific string
    all_out = []
    if args.exclude_files_with_str is not None:
        for fnm in all_out_:
            if fnm.find(args.exclude_files_with_str) == -1:
                all_out.append(fnm)
    else:
        all_out = all_out_

    # aggregate the results
    logging.info('Found {} scEPS score files'.format(len(all_out)))
    df_sceps = []
    for fnm in tqdm(all_out):
        df_sceps_batch = pd.read_table(fnm)
        if df_usenb is not None:
            df_sceps_batch = df_sceps_batch[df_sceps_batch['CELL'].isin(df_usenb[0])].reset_index(drop=True)
        df_sceps.append(df_sceps_batch)
    df_sceps = pd.concat(df_sceps, ignore_index=True)

    # apply testing at cell neighborhood level
    logging.info("Testing at cell neighborhood level")
    test_cell(df_sceps)

    # save the aggregated sceps score file
    df_sceps.to_csv(out_fnm, sep='\t', index=False, float_format='%.5g')

    return df_sceps

def test_cell(df_sceps):

    ncell = df_sceps.shape[0]

    # iterate through estimand
    for estimand in estimand_col:

        pval = df_sceps['P_Z_{}'.format(estimand)]
        est_val = df_sceps[estimand].values

        signif_fdr5 = np.zeros(ncell, dtype=bool)
        signif_fdr10 = np.zeros(ncell, dtype=bool)
        signif_fdr20 = np.zeros(ncell, dtype=bool)

        signif_fdr5 = (est_val > 0) & (smm.multipletests(pval, alpha=0.05, method='fdr_bh')[0])
        signif_fdr10 = (est_val > 0) & (smm.multipletests(pval, alpha=0.10, method='fdr_bh')[0])
        signif_fdr20 = (est_val > 0) & (smm.multipletests(pval, alpha=0.20, method='fdr_bh')[0])

        if estimand.find('DIFF') > 0:
            estimand_gwas = estimand.replace('DIFF', 'GWAS')
            signif_fdr5 = signif_fdr5 & (df_sceps[estimand_gwas].values > 0)
            signif_fdr10 = signif_fdr10 & (df_sceps[estimand_gwas].values > 0)
            signif_fdr20 = signif_fdr20 & (df_sceps[estimand_gwas].values > 0)

        df_sceps['SIGNIF_FDR5_{}'.format(estimand)] = signif_fdr5
        df_sceps['SIGNIF_FDR10_{}'.format(estimand)] = signif_fdr10
        df_sceps['SIGNIF_FDR20_{}'.format(estimand)] = signif_fdr20

    return df_sceps


def test_all_cell_types(args, cell_type_col, obs, df_sceps):
    """
    Test all the cell types
    """

    # get all cell type
    all_ct = pd.unique(obs[cell_type_col])

    # perform testing for each cell type
    df_test_out = []
    for ct in tqdm(all_ct):

        # get test stats
        df_test_out_ct = test_cell_type(args, cell_type_col, obs, ct, df_sceps)
        if df_test_out_ct is None:
            continue

        # add additional information
        df_test_out_ct['CELL_TYPE'] = ct
        if args.add_column is not None:
            if (args.add_column[0] in df_test_out_ct.columns) == False:
                df_test_out_ct[args.add_column[0]] = args.add_column[1]

        # append to list
        df_test_out.append(df_test_out_ct)

    # aggregate results
    if len(df_test_out) > 0:
        df_test_out = pd.concat(df_test_out)
        df_test_out = df_test_out.dropna(axis=1, how='all')
        ct_col = df_test_out.pop('CELL_TYPE')
        df_test_out.insert(0, 'CELL_TYPE', ct_col) 

        return df_test_out

    # return none if nothing to return
    return None


def test_cell_type(args, cell_type_col, obs, ct, df_sceps):
    """
    Test a particular cell type
    """

    # extract cells from the cell type
    ct_cells = obs[obs[cell_type_col]==ct][args.cell_id_col]
    df_sceps_ct = df_sceps[df_sceps['CELL'].isin(ct_cells)].copy()

    # add block information if using block bootstrap
    if args.block_bootstrap != '':
        cell2block = dict(zip(obs[args.cell_id_col], obs[args.block_bootstrap]))
        df_sceps_ct['BLOCK'] = df_sceps_ct['CELL'].map(cell2block)

    # check if dataframe empty
    num_cells = df_sceps_ct.shape[0]
    if num_cells == 0:
        return None

    # get summary info
    df_all_test_out = []
    for col in info_col:
        if col in df_sceps_ct.columns:
            df_all_test_out.append(pd.DataFrame({'MEAN_'+col: [np.mean(df_sceps_ct[col])]}))

    # get bootstrap test statistics
    use_inv_var_wgt = False
    if args.use_inverse_variance_weights == True:
        use_inv_var_wgt = True
    for col in estimand_col:
        test_out = get_bootstrap_test_stats(df_sceps_ct, col,
            nbs=args.num_bootstrap, use_inv_var_wgt=use_inv_var_wgt)
        df_all_test_out.append(pd.DataFrame(test_out))
    df_all_test_out = pd.concat(df_all_test_out, axis=1)
    
    # add number of neighborhoods
    df_all_test_out['NUM_CELL'] = num_cells

    return df_all_test_out


def get_bootstrap_test_stats(df, col, nbs=1000, use_inv_var_wgt=False, count_signif=True):

    # prepare bootstrap
    ncell = df.shape[0]
    stats_vec = df[col].values
    se_vec = df['SE_'+col].values
    var_vec = np.square(se_vec)
    inv_var_vec = 1.0 / (var_vec + EPS)
    weight_vec = np.ones(ncell) / ncell
    if use_inv_var_wgt == True:
        weight_vec = inv_var_vec / np.sum(inv_var_vec)
    mean_stats = np.sum(stats_vec * weight_vec)
    all_mean_stats_bs = []

    # get number of significant cell neighborhood
    if count_signif == True:
        num_signif_fdr5 = np.sum(df['SIGNIF_FDR5_{}'.format(col)])
        num_signif_fdr10 = np.sum(df['SIGNIF_FDR10_{}'.format(col)])
        num_signif_fdr20 = np.sum(df['SIGNIF_FDR20_{}'.format(col)])

    # standard bootstrap
    if 'BLOCK' not in df.columns:
        all_idx = np.array(range(ncell))
        weight_vec_bs = weight_vec.copy()
        for _ in nbs:
            # bootstrap the cells
            use_idx = np.random.choice(all_idx, size=ncell, replace=True)
            if use_inv_var_wgt == True:
                weight_vec_bs = inv_var_vec[use_idx] / np.sum(inv_var_vec[use_idx])
            all_mean_stats_bs.append(np.sum(stats_vec[use_idx] * weight_vec_bs[use_idx]))
    # block bootstrap
    else:
        block_val = df['BLOCK'].values
        all_block = np.array(pd.unique(df['BLOCK']))
        nblock = all_block.shape[0]
        weight_vec_bs = weight_vec.copy()
        # bootstrap the blocks
        for _ in range(nbs):
            use_block = np.random.choice(all_block, size=nblock, replace=True)
            use_idx = np.concatenate([np.where(block_val == blk)[0] for blk in use_block])
            if use_inv_var_wgt == True:
                weight_vec_bs = inv_var_vec[use_idx] / np.sum(inv_var_vec[use_idx])
            all_mean_stats_bs.append(np.sum(stats_vec[use_idx] * weight_vec_bs[use_idx]))
    
    # get test statistics
    se_mean_stats = np.std(all_mean_stats_bs)
    z_mean_stats = mean_stats / (se_mean_stats + EPS)
    p_mean_stats = (1-scipy.stats.norm.cdf(np.fabs(z_mean_stats)))*2.0

    # create out dict
    if count_signif == True:
        out = {'NUM_SIGNIF_FDR5_{}'.format(col): [num_signif_fdr5],
            'NUM_SIGNIF_FDR10_{}'.format(col): [num_signif_fdr10],
            'NUM_SIGNIF_FDR20_{}'.format(col): [num_signif_fdr20],
            'MEAN_{}'.format(col): [mean_stats],
            'SE_MEAN_{}'.format(col): [se_mean_stats],
            'Z_MEAN_{}'.format(col): [z_mean_stats],
            'P_Z_MEAN_{}'.format(col): [p_mean_stats]}
    else:
        out = {'MEAN_{}'.format(col): [mean_stats],
               'SE_MEAN_{}'.format(col): [se_mean_stats],
               'Z_MEAN_{}'.format(col): [z_mean_stats],
               'P_Z_MEAN_{}'.format(col): [p_mean_stats]}

    return out

def get_command_line():
 
    # Create the parser
    parser = argparse.ArgumentParser(description="Summarize the results")

    # Input regulated command line arguments
    parser.add_argument('--prefix', type=str, required=False,
        help="""Used to specify a regular expression for the file names of the output """ \
        """for individual cell neighborhoods from step 1.""")

    parser.add_argument('--adata', type=str, required=False,
        help="""Used to specify the same single-cell data used for obtaining scEPS statistics """ \
        """at individual cell neighborhood level, with cell ID column specified by the --cell-id-col flag.""")

    parser.add_argument('--cell-id-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents cell IDs in the adata.obs data """ \
        """frame of the single-cell data. If left empty, scEPS will use what's in adata.obs.index as cell IDs.""")

    parser.add_argument('--exclude-files-with-str', type=str, default=None, required=False,
        help="""Used to filter out files with specific strings in their file names. This is """ \
        """primarily used for debugging purposes.""")

    parser.add_argument('--neighborhood-clusters', type=str, required=False, default='',
        help="""Used to specify the text file from step 2, representing a pre-computed mapping of cell """ \
        """neighborhoods to approximately independent blocks of cell neighborhoods""")
    
    # Testing related command line argument
    parser.add_argument('--cell-type-col', type=str, required=False, default='',
        help="""Used to specify a list of column names (e.g., at different resolutions) in adata.obs """ \
        """representing cell types. The list of column names need to be separated by commas. The tool will """ \
        """calculate the average scEPS statistics for each cell type under each cell type column, in separate files. """ \
        """If this flag is not specified, scEPS will aggregate results across all cell neighborhoods.""")

    parser.add_argument('--use-neighborhoods', type=str, required=False, default='',
        help="""Used to specify a text file containing a list of cell IDs representing cell neighborhoods. """ \
        """By default, this is an empty string, and the tool aggregates results from all cell neighborhoods. """ \
        """If this is non-empty, scEPS will calculate aggregated statistics using only the cell neighborhoods """ \
        """listed in the text file.""")

    parser.add_argument('--use-inverse-variance-weights', default=False, required=False, action='store_true',
        help="""If specified (not recommended), the tool will use inverse variance weighted average as """ \
        """the aggregated statistics for each cell type.""")

    parser.add_argument('--use-sigma', default=False, required=False, action='store_true',
        help="""If specified, the tool will aggregate the scEPS SIGMA statistics instead of the OMEGA """ \
        """statistics. If this flag is specified, the --prefix flag should specify file names with """ \
        """.sceps.sigma.txt.gz as suffix.""")

    parser.add_argument('--block-bootstrap', type=str, required=False, default='sceps.neighborhood_cluster',
        help="""Used to specify the column in that represents approximately independent cell """  \
        """neighborhood blocks (sceps.neighborhood_cluster by default). If left empty, """ \
        """the tool will fall back to regular bootstrap across individual (instead of blocks of) """ \
        """cell neighborhoods.""")

    parser.add_argument('--num-bootstrap', type=int, required=False, default=1000,
        help="""Used to specify the number of bootstrap samples (1,000 by default).""")

    parser.add_argument('--seed', type=int, required=False, default=0,
        help="""Used to specify the seed for the random number generator (0 by default).""")

    # Output related command line argument
    parser.add_argument('--add-column', type=str, required=False, nargs=2, default=None,
        help="""Used to specify the column name and a value for the additional column to add.""")

    parser.add_argument('--out', type=str, required=False,
        help="""Used to specify the prefix of the output file name.""")

    # Execute the parse_args() method
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    main()
