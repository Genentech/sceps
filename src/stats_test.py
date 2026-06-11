import argparse, sys, glob, os
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
from tqdm import tqdm
import logging
import statsmodels as sm
import statsmodels.stats.meta_analysis

from scipy.stats import beta
from src.utils import EPS

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def test_cell(df_sceps, df_sceps_testing, prefix='OMEGA_'):
    """
    Test the results at cell level
    """

    # get number of donors and parameters
    num_param = df_sceps.groupby('CELL', sort=False)['NUM_PARAM'].mean().values[0]
    num_donor = df_sceps.groupby('CELL', sort=False)['NUM_DONOR'].mean().values[0]

    # get mean of mean of expression
    mean_mean_gwas = df_sceps.groupby('CELL', sort=False)['MEAN_MEAN_EXPR_GWAS'].mean().values
    mean_mean_ctrl = df_sceps.groupby('CELL', sort=False)['MEAN_MEAN_EXPR_CONTROL'].mean().values
    mean_mean_rest = None
    if 'MEAN_MEAN_EXPR_REST' in df_sceps.columns:
        mean_mean_rest = df_sceps.groupby('CELL', sort=False)['MEAN_MEAN_EXPR_REST'].mean().values
    mean_mean_all = df_sceps.groupby('CELL', sort=False)['MEAN_MEAN_EXPR_ALL'].mean().values

    # get mean of variance of expression
    mean_var_gwas = df_sceps.groupby('CELL', sort=False)['MEAN_VAR_EXPR_GWAS'].mean().values
    mean_var_ctrl = df_sceps.groupby('CELL', sort=False)['MEAN_VAR_EXPR_CONTROL'].mean().values
    mean_var_rest = None
    if 'MEAN_VAR_EXPR_REST' in df_sceps.columns:
        mean_var_rest = df_sceps.groupby('CELL', sort=False)['MEAN_VAR_EXPR_REST'].mean().values
    mean_var_all = df_sceps.groupby('CELL', sort=False)['MEAN_VAR_EXPR_ALL'].mean().values

    # get p-values for the statistics
    mean_sigma_gwas, p_sigma_gwas_nominal, z_sigma_gwas, se_sigma_gwas = get_test_stats(df_sceps, df_sceps_testing, 'GWAS', prefix=prefix)
    mean_sigma_ctrl, p_sigma_ctrl_nominal, z_sigma_ctrl, se_sigma_ctrl = get_test_stats(df_sceps, df_sceps_testing, 'CONTROL', prefix=prefix)
    mean_sigma_rest, p_sigma_rest_nominal, z_sigma_rest, se_sigma_rest = get_test_stats(df_sceps, df_sceps_testing, 'REST', prefix=prefix)
    mean_sigma_overall, p_sigma_overall_nominal, z_sigma_overall, se_sigma_overall = get_test_stats(df_sceps, \
        df_sceps_testing, 'OVERALL', prefix=prefix)
    mean_diff, p_diff_nominal, z_diff, se_diff = get_test_stats(df_sceps, df_sceps_testing, 'DIFF', prefix=prefix)

    # prepare the output data frame
    df_sceps_cell = df_sceps.drop_duplicates(subset=['CELL']).copy()

    df_sceps_cell.loc[:,prefix+'GWAS'] = mean_sigma_gwas
    df_sceps_cell.loc[:,'SE_{}GWAS'.format(prefix)] = se_sigma_gwas
    df_sceps_cell.loc[:,'P_{}GWAS_PERMUTE'.format(prefix)] = p_sigma_gwas_nominal
    df_sceps_cell.loc[:,'Z_{}GWAS'.format(prefix)] = z_sigma_gwas
    df_sceps_cell.loc[:,'P_Z_{}GWAS'.format(prefix)] = (1-scipy.stats.t.cdf(np.fabs(z_sigma_gwas), num_donor-num_param))*2.0
    
    df_sceps_cell.loc[:,prefix+'CONTROL'] = mean_sigma_ctrl
    df_sceps_cell.loc[:,'SE_{}CONTROL'.format(prefix)] = se_sigma_ctrl
    df_sceps_cell.loc[:,'P_{}CONTROL_PERMUTE'.format(prefix)] = p_sigma_ctrl_nominal
    df_sceps_cell.loc[:,'Z_{}CONTROL'.format(prefix)] = z_sigma_ctrl
    df_sceps_cell.loc[:,'P_Z_{}CONTROL'.format(prefix)] = (1-scipy.stats.t.cdf(np.fabs(z_sigma_ctrl), num_donor-num_param))*2.0
    
    df_sceps_cell.loc[:,prefix+'REST'] = mean_sigma_rest
    df_sceps_cell.loc[:,'SE_{}REST'.format(prefix)] = se_sigma_rest
    df_sceps_cell.loc[:,'P_{}REST_PERMUTE'.format(prefix)] = p_sigma_rest_nominal
    df_sceps_cell.loc[:,'Z_{}REST'.format(prefix)] = z_sigma_rest
    df_sceps_cell.loc[:,'P_Z_{}REST'.format(prefix)] = (1-scipy.stats.t.cdf(np.fabs(z_sigma_rest), num_donor-num_param))*2.0

    df_sceps_cell.loc[:,prefix+'OVERALL'] = mean_sigma_overall
    df_sceps_cell.loc[:,'SE_{}OVERALL'.format(prefix)] = se_sigma_overall
    df_sceps_cell.loc[:,'P_{}OVERALL_PERMUTE'.format(prefix)] = p_sigma_overall_nominal
    df_sceps_cell.loc[:,'Z_{}OVERALL'.format(prefix)] = z_sigma_overall
    df_sceps_cell.loc[:,'P_Z_{}OVERALL'.format(prefix)] = (1-scipy.stats.t.cdf(np.fabs(z_sigma_overall), num_donor-num_param-2))*2.0

    df_sceps_cell.loc[:,prefix+'DIFF'] = mean_diff
    df_sceps_cell.loc[:,'SE_{}DIFF'.format(prefix)] = se_diff
    df_sceps_cell.loc[:,'P_{}DIFF_PERMUTE'.format(prefix)] = p_diff_nominal
    df_sceps_cell.loc[:,'Z_{}DIFF'.format(prefix)] = z_diff
    df_sceps_cell.loc[:,'P_Z_{}DIFF'.format(prefix)] = (1-scipy.stats.t.cdf(np.fabs(z_diff), num_donor-num_param-1))*2.0
   
    df_sceps_cell.loc[:,'MEAN_MEAN_EXPR_GWAS'] = mean_mean_gwas
    df_sceps_cell.loc[:,'MEAN_MEAN_EXPR_CONTROL'] = mean_mean_ctrl
    df_sceps_cell.loc[:,'MEAN_MEAN_EXPR_ALL'] = mean_mean_all

    df_sceps_cell.loc[:,'MEAN_VAR_EXPR_GWAS'] = mean_var_gwas
    df_sceps_cell.loc[:,'MEAN_VAR_EXPR_CONTROL'] = mean_var_ctrl
    df_sceps_cell.loc[:,'MEAN_VAR_EXPR_ALL'] = mean_var_all

    # keep relevant columns
    if (mean_var_rest is not None) and (mean_mean_rest is not None):
        out_cols = ['CELL', 'DONOR_ID', 'STEP_SIZE', 'TRANSIT_PROP_THRES', 'NEIGHBORHOOD_SIZE', 'NUM_DONOR', 'ELAPSED_TIME', \
                    'MEAN_MEAN_EXPR_GWAS', 'MEAN_MEAN_EXPR_CONTROL', 'MEAN_MEAN_EXPR_REST', 'MEAN_MEAN_EXPR_ALL', \
                    'MEAN_VAR_EXPR_GWAS', 'MEAN_VAR_EXPR_CONTROL', 'MEAN_VAR_EXPR_REST', 'MEAN_VAR_EXPR_ALL', \
                    'VAR_OUTCOME', 'NUM_GENE_GWAS', 'NUM_GENE_CONTROL', 'NUM_GENE_REST']
    else:
        out_cols = ['CELL', 'DONOR_ID', 'STEP_SIZE', 'TRANSIT_PROP_THRES', 'NEIGHBORHOOD_SIZE', 'NUM_DONOR', 'ELAPSED_TIME', \
                    'MEAN_MEAN_EXPR_GWAS', 'MEAN_MEAN_EXPR_CONTROL', 'MEAN_MEAN_EXPR_ALL', \
                    'MEAN_VAR_EXPR_GWAS', 'MEAN_VAR_EXPR_CONTROL', 'MEAN_VAR_EXPR_ALL', \
                    'VAR_OUTCOME', 'NUM_GENE_GWAS', 'NUM_GENE_CONTROL', 'NUM_GENE_REST']
    for stats in ['GWAS', 'CONTROL', 'REST', 'OVERALL', 'DIFF']:
        out_cols = out_cols + [prefix+stats, 'SE_{}{}'.format(prefix, stats), 'Z_{}{}'.format(prefix, stats)]
        out_cols = out_cols + ['P_{}{}_PERMUTE'.format(prefix, stats), 'P_Z_{}{}'.format(prefix, stats)]
    
    # intersect with what's in the data frame
    out_cols_use = []
    for c in out_cols:
        if c in df_sceps_cell.columns:
            out_cols_use.append(c)
    df_sceps_cell = df_sceps_cell[out_cols_use]

    # drop columns with all nan values
    df_sceps_cell = df_sceps_cell.dropna(axis=1, how='all')

    return df_sceps_cell


def get_weighted_mean(df_sceps, col, wgt):
    """
    Get weighted mean
    """

    df_sceps.loc[:,col+'_TMP'] = df_sceps[col] * df_sceps[wgt]
    numerator = df_sceps.groupby('CELL', sort=False)[col+'_TMP'].sum().values
    denom = df_sceps.groupby('CELL', sort=False)[wgt].sum().values
    mean = numerator / denom

    return mean


def get_test_stats(df_sceps, df_sceps_testing, col, prefix=''):

    # prepare for testing
    all_cells = pd.unique(df_sceps['CELL'])
    ncell = all_cells.shape[0]
    cells_val = df_sceps['CELL'].values
    if df_sceps_testing is not None:
        cells_testing_val = df_sceps_testing['CELL'].values

    stats_val = df_sceps[prefix+col].values
    se_stats_val = df_sceps['SE_{}{}'.format(prefix, col)].values
    if df_sceps_testing is not None:
        stats_testing_val = df_sceps_testing[prefix+col].values

    # iterate through cells
    all_p_nominal = []
    all_stats = []
    all_se_stats = []
    for i in range(ncell):

        # get data for cell i
        cell_i = all_cells[i]
        idx = np.where(cells_val == cell_i)[0]
        if df_sceps_testing is not None:
            idx_testing = np.where(cells_testing_val == cell_i)[0]

        # prepare for testing
        stats_i = stats_val[idx]
        se_stats_i = se_stats_val[idx]
        if df_sceps_testing is not None:
            stats_testing_i = stats_testing_val[idx_testing]

        # perform testing
        p_nominal = np.nan
        if df_sceps_testing is not None:
            p_nominal_raw = (np.fabs(np.mean(stats_i)) < np.fabs(stats_testing_i))
            p_nominal = (np.sum(p_nominal_raw) + 1) / (p_nominal_raw.shape[0] + 1)
        
        all_p_nominal.append(p_nominal)
        idx_sorted = np.argsort(stats_i)
        idx_pct50 = idx_sorted[int(np.floor(float(stats_i.shape[0])/2.0))]
        all_stats.append(stats_i[idx_pct50])
        all_se_stats.append(se_stats_i[idx_pct50])

    # convert to numpy array
    all_stats = np.array(all_stats)
    all_se_stats = np.array(all_se_stats)
    all_p_nominal = np.array(all_p_nominal)
    all_z_stats = all_stats / (all_se_stats + EPS)

    return all_stats, all_p_nominal, all_z_stats, all_se_stats
