import sys, gc, collections
from packaging import version

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import statsmodels.api as sm

import scipy
import scipy.stats as st
import scipy.sparse as sp

from sklearn.utils.sparsefuncs import inplace_row_scale
from src.utils import *


def get_pseudobulk(args, adata):
    """
    Get pseudobulk from single-cell data
    """

    # get donor id's
    obs = adata.obs.copy()
    obs = obs.drop_duplicates(subset=[args.donor_id_col], keep='first').reset_index(allow_duplicates=True)
    all_donor_id = pd.unique(obs[args.donor_id_col])

    # create pseudobulk
    ndonors = all_donor_id.shape[0]
    ngenes = adata.var.shape[0]
    pseudobulk_X = np.zeros((ndonors, ngenes))
    num_cells = np.zeros(ndonors)

    # iterate through donor id
    X = adata.X
    donor_id_val = adata.obs[args.donor_id_col].values
    for i, donor_id in enumerate(all_donor_id):

        # extract data for donor i
        cell_idx = np.where(donor_id_val==donor_id)[0]
        adata_donor_X = X[cell_idx, :]

        # save pseudobulk data
        pseudobulk_X_donor = np.mean(adata_donor_X, axis=0)
        pseudobulk_X[i,:] = pseudobulk_X_donor

        # save number of cells
        num_cells[i] = cell_idx.shape[0]

    # create adata
    adata_pseudobulk = ad.AnnData(pseudobulk_X, dtype=np.float32)
    adata_pseudobulk.obs = obs.copy()
    adata_pseudobulk.var = adata.var.copy()
    adata_pseudobulk.obs['sceps.num_cells'] = num_cells
    adata_pseudobulk.obs['sceps.prop_cells_local'] = num_cells / np.sum(num_cells)

    # filter out samples with small number of cells
    adata_pseudobulk = adata_pseudobulk[adata_pseudobulk.obs['sceps.num_cells']>=args.min_num_cell]

    # filter out genes with std less than threshold
    adata_pseudobulk.var['mean'] = adata_pseudobulk.X.mean(axis=0)
    adata_pseudobulk.var['std'] = adata_pseudobulk.X.std(axis=0)
    expr_std = adata_pseudobulk.var['std'].values
    adata_pseudobulk = adata_pseudobulk[:, (np.square(expr_std) > args.min_expr_var_thres) &\
                                           (np.square(expr_std) < args.max_expr_var_thres)]

    return adata_pseudobulk


def next_control_genes(args, adata, target, nbins=10, nset=500):
    """
    Select mean matched control genes
    """

    # create bins and count how many top genes are in each bin
    adata_var = adata.var.copy()

    # if control gene list not empty, use predetermined control genes
    if args.control_gene_list != '':
        controls = pd.read_table(args.control_gene_list, header=None)[0]
        controls_index = adata_var[adata_var[args.gene_id_col].isin(controls)]['GENE_INDEX'].values
        for _ in range(nset):
            yield controls_index
    
    # use randomly selected mean expression matched control genes
    else:
        out_mean = pd.qcut(adata_var['mean'].values, nbins, duplicates='drop')
        adata_var['bin_mean'] = out_mean.astype(str)
        adata_var['bin_std'] = ''

        if args.match_expr_std:
            out_std = pd.qcut(adata_var['std'].values, nbins, duplicates='drop')
            adata_var['bin_std'] = out_std.astype(str)

        adata_var['bin'] = adata_var['bin_mean'] + ' | ' + adata_var['bin_std']

        # count occurance of GWAS genes in each bin
        adata_var_gwas = adata_var[adata_var[args.gene_id_col].isin(target)].copy()
        adata_var['weight'] = 0.0
        gwas_cnt = adata_var_gwas['bin'].value_counts()
        for idx, cnt in enumerate(gwas_cnt):

            gene_idx = adata_var['bin']==gwas_cnt.index[idx]
            adata_var.loc[gene_idx,'weight'] = cnt

            # enrich for low gene-level association statistics
            if args.enrich_low_score_control_genes == True:
                gene_score_bin = adata_var.loc[gene_idx,'sceps.SCORE'].values
                score_rank = (np.argsort(np.argsort(-1*gene_score_bin)) + 1.0) / gene_score_bin.shape[0]
                adata_var.loc[gene_idx,'weight'] = adata_var.loc[gene_idx,'weight'] + score_rank

        # remove target gene from adata_var
        if args.include_gwas_genes_in_control_genes == False:
            adata_var = adata_var[~adata_var[args.gene_id_col].isin(target)]
        adata_var = adata_var.reset_index(drop=True)

        # sample control genes according to the weights
        for _ in range(nset):
            controls = adata_var.sample(n=target.shape[0], weights=adata_var['weight'])
            controls_index = controls['GENE_INDEX'].values
            yield controls_index


def prep_adata_obs(args, adata):

    # retain only relevant data
    obs = adata.obs.copy()

    # If args.cell_id_col not in adata.obs.columns, we use adata.obs.index as replacement
    if args.cell_id_col not in obs.columns:
        if args.cell_id_col == '':
            args.cell_id_col = 'sceps.cell_index'
        obs[args.cell_id_col] = obs.index.copy()
        adata.obs = obs
    
    # calculate global cell proportions
    obs_donor_cell = obs[[args.donor_id_col, args.cell_id_col]].copy()
    prop_cells_global = obs_donor_cell.groupby(args.donor_id_col).count()
    prop_cells_global = prop_cells_global / obs.shape[0]
    prop_cells_global = prop_cells_global.reset_index()
    prop_cells_global.columns = [args.donor_id_col, 'sceps.prop_cells_global']
    donor_id_prop_cell_global = prop_cells_global.set_index(args.donor_id_col)['sceps.prop_cells_global'].to_dict()
    obs['sceps.prop_cells_global'] = obs[args.donor_id_col].map(donor_id_prop_cell_global)

    # check if cell id cols is in obs again
    if args.cell_id_col not in obs.columns.tolist():
        logging.info('{} not found in adata.obs'.format(args.cell_id_col))
        sys.exit()

    # prepare the covariate to adjust list
    covar_list = []
    if args.covar != '':
        covar_list = args.covar.split(',')
        if args.batch_key != '':
            covar_list = covar_list + [args.batch_key]
        if '' in covar_list:
            covar_list.remove('')
        if 'sceps.prop_cells_local' in covar_list:
            covar_list.remove('sceps.prop_cells_local') # sceps.prop_cells_local will be added later
        covar_list = np.unique(np.array(covar_list)).tolist()
        obs = obs[[args.cell_id_col, args.donor_id_col, args.pheno]+covar_list]
    else:
        if args.batch_key == '':
            obs = obs[[args.cell_id_col, args.donor_id_col, args.pheno]]
        else:
            obs = obs[[args.cell_id_col, args.donor_id_col, args.pheno, args.batch_key]]

    # remove samples with missing entries
    nan_rows = obs.isna().any(axis=1)
    obs = obs[~nan_rows].reset_index(drop=True)
    if obs.shape[0] == 0:
        logging.info('Data is empty')
        sys.exit()

    # prepare covariates to adjust
    covar_to_adj = []
    for c in covar_list:
        if pd.unique(obs[c]).shape[0] > 1:
            covar_to_adj.append(c)

    # get unique phenotypes and covariates information
    df_pheno_covar = obs[[args.donor_id_col]+[args.pheno]+covar_to_adj].copy()
    df_pheno_covar = df_pheno_covar.drop_duplicates()

    # prepare covariates adjustment
    covar_mat = covar_mat = np.ones((df_pheno_covar.shape[0],1))

    # construct covariates matrix to adjust
    df_covar = df_pheno_covar[[args.donor_id_col]+covar_to_adj]
    df_covar = df_covar.set_index(args.donor_id_col)
    df_pheno_covar = df_pheno_covar.set_index(args.donor_id_col)
    if df_covar.empty == False:
        df_covar = pd.get_dummies(df_covar, drop_first=True).astype(float)
        covar_mat = df_covar.values.reshape((df_covar.shape[0], df_covar.shape[1]))
        covar_mat = sm.add_constant(covar_mat)

    # perform the adjustment
    y = df_pheno_covar[args.pheno].values 
    if pd.unique(df_pheno_covar[args.pheno]).shape[0] <= 2:
        y = pd.get_dummies(df_pheno_covar[args.pheno], drop_first=True).astype(float)
    fit = sm.OLS(y, covar_mat).fit()
    resid = fit.resid
    if args.scale_pheno == True:
        resid = resid / (np.std(resid)+EPS)

    # update the phenotype column
    df_resid = pd.DataFrame({args.donor_id_col: df_pheno_covar.index, \
                                args.pheno: resid})
    
    # update the obs
    obs = obs[obs[args.donor_id_col].isin(df_resid[args.donor_id_col])]
    donor_pheno_dict = dict(zip(df_resid[args.donor_id_col], df_resid[args.pheno]))
    obs[args.pheno] = obs[args.donor_id_col].map(donor_pheno_dict)

    # retain num_cells only if specified
    if 'sceps.prop_cells_local' in args.covar.split(','):
        args.covar = 'sceps.prop_cells_local'
    else:
        args.covar = ''

    obs.index = obs.index.astype(str)

    return obs


def load_single_cell_data(args):

    # Load the phenotype data
    df_pheno = None
    if args.pheno_file != '':
        df_pheno = pd.read_table(args.pheno_file)

    # load the single-cell data
    adata = sc.read_h5ad(args.adata)

    # select the focal cell type, if specified
    if args.cell_type_col != '':
        adata = adata[adata.obs[args.cell_type_col]==args.focal_cell_type]
    logging.info('Loaded {} cells for {} samples'.format(adata.shape[0],
        pd.unique(adata.obs[args.donor_id_col]).shape[0]))

    # If args.gene_id_col not in adata.var.columns, we use adata.var.index as replacement
    if args.gene_id_col not in adata.var.columns:
        if args.gene_id_col == '':
            args.gene_id_col = 'sceps.gene_index'
        adata.var[args.gene_id_col] = adata.var.index.copy()

    # update adata.obs with phenotype
    if df_pheno is not None:
        adata_donor_id = pd.unique(adata.obs[args.donor_id_col])
        df_pheno = df_pheno[df_pheno[args.donor_id_col].isin(adata_donor_id)].reset_index(drop=True)
        adata = adata[adata.obs[args.donor_id_col].isin(df_pheno[args.donor_id_col])]
        adata_donor_id = pd.unique(adata.obs[args.donor_id_col])
        assert adata_donor_id.shape[0] == df_pheno.shape[0]
        adata.obs = adata.obs.reset_index()
        adata.obs = adata.obs.merge(df_pheno, on=[args.donor_id_col], how="left")
        adata.obs = adata.obs.set_index(args.cell_id_col)
        
    obs_updt = prep_adata_obs(args, adata)
    adata = adata[adata.obs[args.cell_id_col].isin(obs_updt[args.cell_id_col])]
    adata.obs = obs_updt

    return adata


def subset_adata_genes(args, adata, gene_list):

    adata = adata[:,adata.var[args.gene_id_col].isin(gene_list[args.gene_id_col])]
    adata.var = adata.var.merge(gene_list, on=[args.gene_id_col])
    adata.var.index = adata.var.index.astype(str)
    gc.collect()
    logging.info('{} genes after intersecting gene list'.format(adata.shape[1]))

    return adata


def extract_donor_pheno(args, adata):

    # Get the vector of each individual cell's corresponding donor ID
    donor_id_vec = adata.obs[args.donor_id_col].values

    # Get the donor phenotype
    df_donor_pheno = adata.obs[[args.donor_id_col, args.pheno]].copy().drop_duplicates()
    if pd.unique(df_donor_pheno[args.pheno]).shape[0] <= 2:
        df_donor_pheno[args.pheno] = pd.get_dummies(df_donor_pheno[args.pheno],\
                                                drop_first=True).astype(float)
    ndonor = df_donor_pheno[args.donor_id_col].unique().shape[0]
    df_donor_pheno = df_donor_pheno.set_index(args.donor_id_col, drop=False)
    df_donor_pheno['INDEX'] = np.array(range(df_donor_pheno.shape[0]))

    # Check to make sure that the donor ID column of df_donor_pheno does not have duplicate
    if ndonor > np.unique(donor_id_vec).shape[0]:
        logging.info('Donor not unique in df_donor_pheno')
        sys.exit()

    return donor_id_vec, df_donor_pheno
