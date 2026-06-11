import sys, gc, time

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

import scipy
import scipy.linalg
import scipy.stats as st
import scipy.sparse as sp
from scipy.sparse import csr_matrix

import statsmodels.api as sm

from src.utils import *
from src.scdata import *

import matplotlib.pyplot as plt
import seaborn as sns


def permute_donor_id(args, donor_id_vec):

    permuted_donor_id = []

    all_donor_id = np.unique(donor_id_vec)
    for _ in range(args.num_control_gene_set_testing):
        np.random.shuffle(all_donor_id)
        permuted_donor_id.append(all_donor_id.copy())

    return permuted_donor_id


def prep_expression(args, X, mean_X, std_X, gene_idx):
    """
    prepare the expression data
    """

    X = X[:, gene_idx]
    std_X_gene = std_X[gene_idx]
    mean_X_gene = std_X[gene_idx]

    return X, np.mean(mean_X_gene), np.mean(np.square(std_X_gene))


def prep_regression(args, adata, genes, permuted_donor_id):
    """
    prepare expression, covariates, and phenotype data for regression
    """

    # remove missing rows
    obs = adata.obs.copy()
    nan_rows = obs.isna().any(axis=1)
    obs = obs[~nan_rows].reset_index(drop=True)
    X_all_gene = adata.X.copy()
    std_X_all_gene = adata.var['std'].values
    mean_X_all_gene = adata.var['mean'].values
    X_all_gene = X_all_gene[~nan_rows, :]
    if np.sum(nan_rows) == obs.shape[0]:
        logging.warning('Covariate data is missing for all donors')
        return None

    # prepare expression
    gene_idx = adata.var.loc[genes]['GENE_INDEX'].values
    X, mean_mean_X, mean_var_X = prep_expression(args, X_all_gene, mean_X_all_gene, std_X_all_gene, gene_idx)

    # prepare phenotype
    df_pheno = obs[[args.donor_id_col]+[args.pheno]]
    df_pheno = df_pheno.set_index(args.donor_id_col)
    if pd.unique(df_pheno[df_pheno.columns[0]]).shape[0] == 1:
        logging.warning('No variation in phenotype')
    if args.check_pheno_in_neighborhood == True:
        if pd.unique(df_pheno[df_pheno.columns[0]]).shape[0] == 1:
            return None
    if pd.unique(df_pheno[df_pheno.columns[0]]).shape[0] == 1:
        df_pheno = pd.get_dummies(df_pheno).astype(float)
    else:
        df_pheno = pd.get_dummies(df_pheno, drop_first=True).astype(float)
    df_pheno.loc[:, 'DONOR_IDX'] = range(df_pheno.shape[0])
    y = df_pheno[df_pheno.columns[0]].values
    y = y - np.mean(y)

    # prepare index information
    if permuted_donor_id is not None:
        permuted_donor_idx = []
        donor_id_set = set(df_pheno.index.tolist())
        for i in range(len(permuted_donor_id)):
            donor_id_shared = [d for d in permuted_donor_id[i] if d in donor_id_set]
            permuted_donor_idx.append(df_pheno.loc[donor_id_shared,'DONOR_IDX'].values)

    # prepare batch
    batch_val = None
    if args.batch_key != "":
        df_batch = obs[[args.donor_id_col]+[args.batch_key]]
        df_batch = df_batch.set_index(args.donor_id_col)
        batch_val = df_batch[df_batch.columns[0]].values

    # prepare covar
    covar = []
    if args.covar != '':
        for c in args.covar.split(','):
            if pd.unique(obs[c]).shape[0] > 1:
                covar.append(c)

    # prepare covariates
    covar_mat = np.array([])
    cond_covar_mat = 1.0
    if len(covar) > 0:
        df_covar = obs[[args.donor_id_col]+covar]
        df_covar = df_covar.set_index(args.donor_id_col)
        df_covar = pd.get_dummies(df_covar, drop_first=True).astype(float)
        covar_mat = df_covar.values.reshape((df_covar.shape[0], df_covar.shape[1]))
        covar_mat = sm.add_constant(covar_mat)
        cond_covar_mat = np.linalg.cond(np.dot(covar_mat.T, covar_mat))
    else:
        covar_mat = np.ones((df_pheno.shape[0],1))

    # prepare output dictionary
    out_dict = dict()
    out_dict['X'] = X
    out_dict['y'] = y
    out_dict['covar_mat'] = covar_mat
    out_dict['cond_covar_mat'] = cond_covar_mat
    out_dict['mean_var_X'] = mean_var_X
    out_dict['mean_mean_X'] = mean_mean_X
    out_dict['X_all_gene'] = X_all_gene
    out_dict['std_X_all_gene'] = std_X_all_gene
    out_dict['mean_X_all_gene'] = mean_X_all_gene
    out_dict['gwas_gene_idx'] = gene_idx
    if batch_val is not None:
        out_dict['batch'] = batch_val
    if permuted_donor_id is not None:
        out_dict['permuted_donor_idx'] = permuted_donor_idx

    # record prop cell info
    covar_mat_with_prop_cell = np.array([])
    covar_with_prop_cell = covar.copy()
    if 'sceps.prop_cells_local' not in covar_with_prop_cell:
        covar_with_prop_cell.insert(0, 'sceps.prop_cells_local')
    else:
        covar_with_prop_cell.remove('sceps.prop_cells_local')
        covar_with_prop_cell.insert(0, 'sceps.prop_cells_local')
    
    df_covar = obs[[args.donor_id_col]+covar_with_prop_cell]
    df_covar = df_covar.set_index(args.donor_id_col)
    df_covar = pd.get_dummies(df_covar, drop_first=True).astype(float)
    covar_mat_with_prop_cell = df_covar.values.reshape((df_covar.shape[0], df_covar.shape[1]))
    covar_mat_with_prop_cell = sm.add_constant(covar_mat_with_prop_cell, prepend=False)
    out_dict['covar_mat_with_prop_cell'] = covar_mat_with_prop_cell

    return out_dict


def covmat_to_vec(cov_mat, offset=1):
    """
    Convert a covariance matrix to a vector
    """

    nsample = cov_mat.shape[0]
    idx = np.triu_indices(nsample, k=offset)
    vec = cov_mat[idx]
    
    return vec


def get_cross_prod(vec, offset=1):
    """
    Obtain cross product of a vector
    """

    nsample = vec.shape[0]
    cp_mat = np.outer(vec, vec)

    idx = np.triu_indices(nsample, k=offset)
    cp = cp_mat[idx]

    return cp


def get_reg_wgt(rel_mat_all, y, mean_var_X_all, offset):
    """
    estimate the regression weights for wls
    """

    # get initial estimate of the parameter
    cp = covmat_to_vec(np.outer(y,y))
    M = covmat_to_vec(rel_mat_all).flatten()
    var_y = np.var(y)

    param = M.dot(cp) / M.dot(M)
    param_low, param_high = 0.0, var_y/mean_var_X_all
    if param < param_low:
        param = param_low
    if param > param_high:
        param = param_high
    
    # estimate the variance of the cross product
    pred_cov = param * rel_mat_all
    pred_var = 2.0 * np.square(pred_cov) + np.square(var_y) - np.square(pred_cov)
    np.fill_diagonal(pred_var, 2.0*np.square(var_y))

    # get reg wgt
    reg_wgt = 1 / (covmat_to_vec(pred_var, offset=offset) + EPS)

    return reg_wgt


def get_bootstrap_disattenuation_factor(M, all_M_bs, num_bins=1):
    """
    get bootstrap disattenuation factor using bootstrapped samples
    """

    # create an array for the original M and bootstrapped M
    all_M = np.zeros((all_M_bs.shape[0]+1, all_M_bs.shape[1]))
    all_M[0,:] = M
    all_M[1:,:] = all_M_bs

    # find the indices of the M vector used for estimating disattenuation factor
    M_sq = np.square(M)
    idx_sort = np.argsort(M_sq)
    idx_sort_bins = np.array_split(idx_sort, num_bins)

    # estimate the disattenuation factor
    all_disatt_factor = []
    all_weight = []
    for idx in idx_sort_bins:
        mean_overall = np.mean(all_M[:,idx])
        mean_M = all_M[:,idx].mean(axis=0)
        sb, sw = np.var(mean_M - mean_overall), np.var(all_M[:,idx] - mean_M)
        disatt_factor_bin = 1.0 + sw / sb
        all_disatt_factor.append(disatt_factor_bin)
        all_weight.append(np.sum(M_sq[idx]))
    all_disatt_factor = np.array(all_disatt_factor)
    all_weight = np.array(all_weight)

    # final disatt factor is weighted average of each bin
    disatt_factor = np.sum(all_disatt_factor * all_weight) / np.sum(all_weight)

    return disatt_factor


def estimate_disattenuation_factor(args, X, M, all_M_bs, offset):
    """
    calculate the disattenuation factor for regression dilution
    """

    if args.no_disattenuation == True:
        return 1.0

    # disattenuation factor without the diag
    if offset == 1:
        disatt_factor = get_bootstrap_disattenuation_factor(M, all_M_bs)
    else:
        # get diagonal and off diagonal indices
        ndonor = X.shape[0]
        rows, cols = np.triu_indices(ndonor)
        idx_offdiag = np.where(rows!=cols)[0]
        idx_diag = np.where(rows==cols)[0]

        # get disatt factor for the diagonals and off diagonals separately
        disatt_factor_offdiag = get_bootstrap_disattenuation_factor(M[idx_offdiag], all_M_bs[:, idx_offdiag])
        disatt_factor_diag = get_bootstrap_disattenuation_factor(M[idx_diag], all_M_bs[:, idx_diag])

        # get average disatt factor
        disatt_factor = disatt_factor_offdiag * idx_offdiag.shape[0] + disatt_factor_diag * idx_diag[1]
        disatt_factor = disatt_factor / M.shape[0]

    return disatt_factor


def get_rel_mat(X, offset, nbs=1000):

    ngene = X.shape[1]
    rel_mat = np.dot(X, X.T) / (ngene-1.0)
    rel_mat = covmat_to_vec(rel_mat, offset)
    all_rel_mat_bs = []

    all_idx = np.array(range(ngene))
    for _ in range(nbs):
        use_idx = np.random.choice(all_idx, size=ngene, replace=True)
        X_bs = X[:,use_idx]
        rel_mat_bs = np.dot(X_bs, X_bs.T) / (ngene-1.0)
        rel_mat_bs = covmat_to_vec(rel_mat_bs, offset)
        all_rel_mat_bs.append(rel_mat_bs)
    all_rel_mat_bs = np.array(all_rel_mat_bs)

    rel_mat = 2.0 * rel_mat - np.mean(all_rel_mat_bs, axis=0)

    return rel_mat, all_rel_mat_bs


def wls_bootstrap(cp, M, all_M_bs, wgt, nbs=1000):

    # get estimates using all data
    wls_model = sm.WLS(cp, M, weights=wgt)
    results = wls_model.fit()
    params = results.params

    # get bootstrapped params over samples
    nrow = cp.shape[0]
    all_idx = np.array(range(nrow))
    params_bs_sample = []
    for _ in range(nbs):
        idx = np.random.choice(all_idx, size=nrow, replace=True)
        wls_model_bs = sm.WLS(cp[idx], M[idx], weights=wgt[idx])
        results_bs = wls_model_bs.fit()
        params_bs_sample.append(results_bs.params)
    params_bs_sample = np.array(params_bs_sample)

    # get bootstrapped parameters
    cov_params_sample = np.cov(params_bs_sample.T)
    cov_params = cov_params_sample

    # correct for biases
    params = 2.0 * params - np.mean(params_bs_sample, axis=0)

    return params, cov_params


def regression(args, adata, genes1, genes2_generator, shuffle_pheno=False,
    permuted_donor_id=None):
    """
    use method of moment to fit the variance component model
    """

    # prepare input
    prep_out1 = prep_regression(args, adata, genes1, permuted_donor_id)
    if prep_out1 is None:
        return None

    # parse out results
    X1 = prep_out1['X']
    y = prep_out1['y']
    covar_mat = prep_out1['covar_mat']
    covar_mat_with_prop_cell = prep_out1['covar_mat_with_prop_cell']
    mean_mean1 = prep_out1['mean_mean_X']
    mean_var1 = prep_out1['mean_var_X']
    X_all_gene = prep_out1['X_all_gene']
    std_X_all_gene = prep_out1['std_X_all_gene']
    mean_X_all_gene = prep_out1['mean_X_all_gene']
    mean_var_all = np.mean(np.square(std_X_all_gene))
    mean_mean_all = np.mean(mean_X_all_gene)
    gwas_gene_idx = prep_out1['gwas_gene_idx']

    # get permutation information
    if permuted_donor_id is not None:
        permuted_donor_idx = prep_out1['permuted_donor_idx']

    # get batch information
    batch_val = None
    uniq_batch = None
    if 'batch' in prep_out1:
        batch_val = prep_out1['batch']
        uniq_batch = np.unique(batch_val)

    # construct relatedness matrix -- double check sample size constraints unless for simulation
    nsample, ngene1 = X1.shape
    if (nsample < args.min_num_donor) and (args.sample_frac_cell_in_neighborhood is None):
        return None

    # prepare phenotypes
    resid = y.copy()
    if covar_mat.size > 0:
        fit = sm.OLS(y, covar_mat).fit()
        resid = fit.resid
    resid = resid - np.mean(resid)
    if args.scale_pheno_neighborhood == True:
        resid = resid / (np.std(resid)+EPS)
    var_pheno = np.var(resid)

    # check if diagonal component is included
    offset = 0
    if args.exclude_diag == True:
        offset = 1

    # prepare regression design matrix
    dim = int(nsample*(nsample-1) / 2) + (1-offset)*nsample
    M = np.ones((dim, args.num_var_comp + (1-offset) + int(not args.exclude_intercept)))
    if args.exclude_diag == False:
        M[:,args.num_var_comp] = covmat_to_vec(np.eye(nsample), offset)
    M1, M1_bs = get_rel_mat(X1, offset, nbs=args.num_bootstrap_disattenuation)
    M[:,0] = M1

    # initialize bootstrap for M
    if (args.use_analytical_stderr == True) and (args.no_disattenuation == True):
        args.num_bootstrap_disattenuation = 0
    all_M_bs = []
    for bs_idx in range(args.num_bootstrap_disattenuation):
        M_bs = M.copy()
        M_bs[:,0] = M1_bs[bs_idx,:]
        all_M_bs.append(M_bs)
    
    # prepare regression response vector
    cp = get_cross_prod(resid, offset)
    
    # estimate disattenuation factor
    disatt_vec = np.ones(M.shape[1])
    disatt_vec[0] = estimate_disattenuation_factor(args, X1, M1, M1_bs, offset)

    # get rel_mat for all genes
    ngene_all = X_all_gene.shape[1]
    rel_mat_all = np.dot(X_all_gene, X_all_gene.T)
    rel_mat_all = rel_mat_all / (ngene_all-1.0)

    # get regression weights
    if args.use_ols == True:
        reg_wgt = np.ones(M.shape[0])
    else:
        reg_wgt = get_reg_wgt(rel_mat_all, resid, mean_var_all, offset)

    # iterate through control genes
    all_out = []
    control_idx = 0
    all_gene_idx = np.array(range(ngene_all))
    for shuf_idx, genes2_idx in enumerate(genes2_generator):

        # shuffle pheno
        if (shuffle_pheno == True) and (permuted_donor_idx is not None):
            resid_shuf = resid[permuted_donor_idx[shuf_idx]]
            cp = get_cross_prod(resid_shuf, offset)

        # construct related mat for control genes
        X2, mean_mean2, mean_var2 = prep_expression(args, X_all_gene, mean_X_all_gene, std_X_all_gene, genes2_idx)
        _, ngene2 = X2.shape
        M2, M2_bs = get_rel_mat(X2, offset, nbs=args.num_bootstrap_disattenuation)
        M[:,1] = M2
        disatt_vec[1] = estimate_disattenuation_factor(args, X2, M2, M2_bs, offset)
        for bs_idx in range(args.num_bootstrap_disattenuation):
            M_bs = all_M_bs[bs_idx]
            M_bs[:,1] = M2_bs[bs_idx,:]
            all_M_bs[bs_idx] = M_bs

        # constract related mat for remaining genes
        if args.num_var_comp == 3:

            rest_idx = np.setdiff1d(all_gene_idx, genes2_idx)
            rest_idx = np.setdiff1d(rest_idx, gwas_gene_idx)
            X3 = X_all_gene[:, rest_idx]
            M3, M3_bs = get_rel_mat(X3, offset, nbs=args.num_bootstrap_disattenuation)
            M[:,2] = M3
            for bs_idx in range(args.num_bootstrap_disattenuation):
                M_bs = all_M_bs[bs_idx]
                M_bs[:,2] = M3_bs[bs_idx,:]
                all_M_bs[bs_idx] = M_bs

            disatt_vec[2] = estimate_disattenuation_factor(args, X2, M3, M3_bs, offset)
            mean_var_rest = np.sum(np.square(std_X_all_gene)) - ngene1*mean_var1 - ngene2*mean_var2
            mean_var_rest = mean_var_rest / (ngene_all-ngene2-ngene1-1.0)
            mean_mean_rest = np.sum(mean_X_all_gene) - ngene1*mean_mean1 - ngene2*mean_mean2
            mean_mean_rest = mean_mean_rest / (ngene_all-ngene2-ngene1-1.0)

        # get regression results
        if args.use_analytical_stderr == True:
            wls_model = sm.WLS(cp, M, weights=reg_wgt)
            results = wls_model.fit()
            fit_params = results.params
            cov_fit_params = results.normalized_cov_params
        else:
            fit_params, cov_fit_params = wls_bootstrap(cp, M, all_M_bs, reg_wgt, nbs=args.num_bootstrap_regression)

        # obtain the params
        params = fit_params * disatt_vec
        mean_var = np.array([mean_var1, mean_var2])
        est = params[0:2] / np.array([ngene1, ngene2])
        est_weighted = est * mean_var
        diff = est[0] - est[1]
        diff_weighted = est_weighted[0] - est_weighted[1]

        var_params = (np.diag(disatt_vec).dot(cov_fit_params)).dot(np.diag(disatt_vec))
        var_params_sub = var_params[0:2,:][:,0:2]

        gwas_vec = np.array([1/ngene1, 0])
        gwas_vec_weighted = gwas_vec * mean_var
        se_gwas = np.sqrt(var_params_sub.dot(gwas_vec).dot(gwas_vec))
        se_gwas_weighted = np.sqrt(var_params_sub.dot(gwas_vec_weighted).dot(gwas_vec_weighted))

        ctrl_vec = np.array([0, 1/ngene2])
        ctrl_vec_weighted = ctrl_vec * mean_var
        se_ctrl = np.sqrt(var_params_sub.dot(ctrl_vec).dot(ctrl_vec))
        se_ctrl_weighted = np.sqrt(var_params_sub.dot(ctrl_vec_weighted).dot(ctrl_vec_weighted))

        diff_vec = np.array([1/ngene1, -1/ngene2])
        diff_vec_weighted = diff_vec * mean_var
        se_diff = np.sqrt(var_params_sub.dot(diff_vec).dot(diff_vec))
        se_diff_weighted = np.sqrt(var_params_sub.dot(diff_vec_weighted).dot(diff_vec_weighted))

        # estimate the sum and rest
        if args.num_var_comp == 3:
            
            overall_vec = np.array([1/ngene_all, 1/ngene_all, 1/ngene_all])
            sigma_overall = params[0:3].dot(overall_vec)
            sigma_overall_weighted = sigma_overall * mean_var_all

            var_params_sub = var_params[0:3,:][:,0:3]
            se_sigma_overall = np.sqrt(var_params_sub.dot(overall_vec).dot(overall_vec))
            se_sigma_overall_weighted = se_sigma_overall * mean_var_all

            ngene_rest = ngene_all - ngene1 - ngene2
            sigma_rest = params[2] / ngene_rest
            se_sigma_rest = np.sqrt(var_params_sub[2,2]) / ngene_rest
            
            sigma_rest_weighted = sigma_rest * mean_var_rest
            se_sigma_rest_weighted = se_sigma_rest * mean_var_rest

        # create data frame
        out = {'CONTROL_IDX': [control_idx],
               'SIGMA_GWAS': [est[0]],
               'OMEGA_GWAS': [est_weighted[0]],
               'SIGMA_CONTROL': [est[1]],
               'OMEGA_CONTROL': [est_weighted[1]],
               'SIGMA_DIFF': [diff],
               'OMEGA_DIFF': [diff_weighted],
               'SE_SIGMA_GWAS': [se_gwas],
               'SE_SIGMA_CONTROL': [se_ctrl],
               'SE_SIGMA_DIFF': [se_diff],
               'SE_OMEGA_GWAS': [se_gwas_weighted],
               'SE_OMEGA_CONTROL': [se_ctrl_weighted],
               'SE_OMEGA_DIFF': [se_diff_weighted],
               'MEAN_MEAN_EXPR_GWAS': [mean_mean1],
               'MEAN_MEAN_EXPR_CONTROL': [mean_mean2],
               'MEAN_MEAN_EXPR_ALL': [mean_mean_all],
               'MEAN_VAR_EXPR_GWAS': [mean_var1],
               'MEAN_VAR_EXPR_CONTROL': [mean_var2],
               'MEAN_VAR_EXPR_ALL': [mean_var_all],
               'NUM_GENE_GWAS': [ngene1],
               'NUM_GENE_CONTROL': [ngene2],
               'VAR_PHENO': [var_pheno],
               'NUM_PARAM': [M.shape[1]]}

        if args.num_var_comp == 3:
            out['NUM_GENE_REST'] = ngene_rest
            out['MEAN_MEAN_EXPR_REST'] = [mean_mean_rest]
            out['MEAN_VAR_EXPR_REST'] = [mean_var_rest]

            out['SIGMA_REST'] = sigma_rest
            out['SE_SIGMA_REST'] = se_sigma_rest
            out['OMEGA_REST'] = sigma_rest_weighted
            out['SE_OMEGA_REST'] = se_sigma_rest_weighted

            out['SIGMA_OVERALL'] = sigma_overall
            out['SE_SIGMA_OVERALL'] = se_sigma_overall
            out['OMEGA_OVERALL'] = sigma_overall_weighted
            out['SE_OMEGA_OVERALL'] = se_sigma_overall_weighted

        all_out.append(pd.DataFrame(out))

        control_idx += 1

    # return the result
    if len(all_out) > 0:
        all_out = pd.concat(all_out, ignore_index=True)
        return all_out

    return None
