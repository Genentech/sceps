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


def get_connectivity(data):
    """
    Get the adjacency matrix from the anndata
    """

    av = ad.__version__
    if type(av) == str:
        av = version.parse(av)
    if av < version.parse("0.7.2"):
        return data.uns["neighbors"]["connectivities"]
    else:
        return data.obsp["connectivities"]


def get_transition_matrix(adj_mat):
    """
    Calculate transition matrix from adjacency matrix
    """

    ncells = adj_mat.shape[0]
    trans_mat = adj_mat + sp.identity(ncells, format=adj_mat.getformat())
    colsums = np.array(trans_mat.sum(axis=0)).flatten()
    inplace_row_scale(trans_mat, 1.0/colsums)

    return trans_mat


def diffuse_stepwise(trans_mat, s, maxnsteps=15):
    """
    Calculate the nam matrix in a stepwise manner
    """

    # do diffusion
    for _ in range(maxnsteps):

        # 2nd term is to keep s as a dataframe
        s = trans_mat.dot(s) + 0.0*s
        yield s


def choose_random_walk_nsteps(args, trans_mat, obs, nsteps=None, maxnsteps=15):
    """
    Determine the optimal number of random walk steps
    """

    def R(A, B):
        return ((A - A.mean(axis=0))*(B - B.mean(axis=0))).mean(axis=0) \
            / A.std(axis=0) / B.std(axis=0)

    S = pd.get_dummies(obs[args.donor_id_col])
    C = S.sum(axis=0)

    prevmedkurt = np.inf
    old_s = np.zeros(S.shape)
    for i, s in enumerate(diffuse_stepwise(trans_mat, S, maxnsteps=maxnsteps)):
        medkurt = np.median(st.kurtosis(s/C, axis=1))
        #R2 = R(s, old_s)**2
        old_s = s
        if nsteps is None:
            if prevmedkurt - medkurt < 3 and i+1 >= 3:
                break
            prevmedkurt = medkurt
        elif i+1 == nsteps:
            break
        gc.collect()

    init_nsteps = i+1

    return init_nsteps, s


def get_target_nam(args, adata):

    adj_mat = get_connectivity(adata)

    obs = adata.obs.copy()
    trans_mat = get_transition_matrix(adj_mat)

    init_nsteps, nam_matrix = choose_random_walk_nsteps(args, trans_mat, obs)

    if args.verbose == True:
        logging.info('Determined random walk step size is {}'.format(init_nsteps))
    if args.save_target_nam == True:
        nam_matrix.to_csv(args.out + '.nam.txt.gz', sep='\t')
    
    return trans_mat, init_nsteps, nam_matrix


def prob_i_reach_all_j(trans_mat, idx_i, nsteps):
    """
    Compute the probability that node i can reach all other node j
    """

    ncells = trans_mat.shape[0]

    vec = np.zeros(ncells)
    vec[idx_i] = 1

    for _ in range(nsteps):
        vec = trans_mat.transpose().dot(vec)

    vec = vec.flatten()

    return vec


def get_neighborhood_cells_thres(prob_vec, prob_thres=EPS):
    """
    Return the indices of the cells in the neighborhood
    """

    idx = np.where(prob_vec >= prob_thres)[0]
    
    return idx


def get_neighborhood_cells_nam(args, prob_vec, donor_id_vec, df_donor_pheno,
    nam_vec, ntimes=1000, check_pheno_in_neighborhood=False):
    """
    Return the indices of the cells in the neighborhood
    """

    # set up the data frame to record the result
    nam_vec_tmp = pd.DataFrame(nam_vec)
    nam_vec_tmp.columns = ['expected_nam']
    nam_vec_tmp.loc[:, 'num_cells'] = 0.0
    nam_vec_tmp.loc[:, 'realized_nam'] = 0.0

    # sort all values based on transition probability
    df_sc = pd.DataFrame({'PROB': prob_vec, 'DONOR': donor_id_vec})
    df_sc = df_sc.sort_values(by=['PROB'], ascending=False)
    prob_vec_sorted = df_sc.loc[:,'PROB'].values
    donor_id_vec_sorted = df_sc.loc[:,'DONOR'].values
    num_donor = pd.unique(df_sc['DONOR']).shape[0]

    # initialize the grid search
    has_multiple_pheno = False
    thres_diff_norm, pheno_ok_donor = [], []
    idx_sorted, last_idx_sorted, nbhood_size, num_ok_donor = 0, 0, 0, 0
    all_thres = np.flip(np.linspace(0.0, 0.01, num=ntimes))
    for thres in all_thres:

        # find the index of the prob just passing the thres
        while idx_sorted < prob_vec.shape[0]:
            if prob_vec_sorted[idx_sorted] <= thres:
                break
            idx_sorted += 1
        
        # stop if trying to include everything
        if prob_vec_sorted[min(idx_sorted, prob_vec_sorted.shape[0]-1)] == 0.0:
            break

        # calculate the number of donors to add
        donor_to_add = donor_id_vec_sorted[last_idx_sorted:idx_sorted]
        ncell_nbhood, last_idx_sorted = idx_sorted, idx_sorted

        # increment neighborhood size
        counter = collections.Counter(donor_to_add)
        for donor in counter:
            nam_vec_tmp.loc[donor, 'num_cells'] += counter[donor]
            nbhood_size += counter[donor]

        # record no. donors that pass minimum requirement
        num_ok_donor = np.sum(nam_vec_tmp['num_cells'] > args.min_num_cell)

        # check if there's variation in phenotype
        if check_pheno_in_neighborhood == True:
            if has_multiple_pheno == False:
                df_ok_donor = nam_vec_tmp[nam_vec_tmp['num_cells'] > args.min_num_cell]
                pheno_ok_donor = df_donor_pheno.loc[df_ok_donor.index, args.pheno]
                has_multiple_pheno = (np.unique(pheno_ok_donor).shape[0] > 1)
            if has_multiple_pheno == False:
                continue

        # calculate the difference between expected and realized nam
        nam_vec_tmp.loc[:, 'realized_nam'] = nam_vec_tmp.loc[:,'num_cells'] / (ncell_nbhood + EPS)
        diff_norm = np.linalg.norm(nam_vec_tmp['expected_nam'] - nam_vec_tmp['realized_nam'])
        
        # record the result
        thres_diff_norm.append([thres, diff_norm, num_ok_donor, nbhood_size])

    # check if any result
    if len(thres_diff_norm) == 0:
        return None

    # select the optimal threshold
    df_thres_diff_norm = pd.DataFrame(thres_diff_norm)
    df_thres_diff_norm.columns = ['THRES', 'DIFF_NORM', 'NUM_OK_DONOR', 'NEIGHBORHOOD_SIZE']
    min_num_donor = max(args.min_num_donor, int(np.ceil(args.min_frac_donor * float(num_donor))))
    df_thres_diff_norm = df_thres_diff_norm[df_thres_diff_norm['NUM_OK_DONOR']>min_num_donor]

    # check if any result
    if df_thres_diff_norm.shape[0] == 0:
        return None

    # select optimal threshold
    df_thres_diff_norm['NEG_NUM_OK_DONOR'] = -df_thres_diff_norm['NUM_OK_DONOR']
    df_thres_diff_norm = df_thres_diff_norm.sort_values(by=['DIFF_NORM', 'NEG_NUM_OK_DONOR', 'NEIGHBORHOOD_SIZE'])
    opt_thres = df_thres_diff_norm.iloc[0]['THRES']

    # get the cell neighborhood using the optimal threshold
    cell_idx_nbhood = get_neighborhood_cells_thres(prob_vec, prob_thres=opt_thres)

    return cell_idx_nbhood, opt_thres


def get_neighborhood_cells(args, idx, trans_mat, nstep, target_nam, donor_id_vec, df_donor_pheno):

    # Initialize neighborhood cell index
    neighborhood_cell_idx, prob_thres = None, None

    # Calculate the probability that cell i can reach all other cell j in nstep based on the trans_mat
    prob_i = prob_i_reach_all_j(trans_mat, idx, nstep)

    # Define cell neighborhood using the "nam" method
    if args.neighborhood_definition_method == 'nam':
        target_nam_row = target_nam.iloc[idx][:]
        nam_neighborhood = get_neighborhood_cells_nam(args, prob_i, donor_id_vec,
            df_donor_pheno, target_nam_row, check_pheno_in_neighborhood=args.check_pheno_in_neighborhood)
        
        # Defining cell neighborhoods failed using the "nam" method
        if nam_neighborhood is None:
            if args.verbose == True:
                logging.info('Neighborhood definition failed, incrementing nstep')
            return None, None
        
        neighborhood_cell_idx = nam_neighborhood[0]
        prob_thres = nam_neighborhood[1]
    
    # Define cell neighborhood using the "soft threshold" method
    elif args.neighborhood_definition_method == 'soft threshold':
        prob_thres = nstep / adata.shape[0]
        neighborhood_cell_idx = get_neighborhood_cells_thres(prob_i, prob_thres=prob_thres)

    # Define cell neighborhood using the "hard threshold" method
    elif args.neighborhood_definition_method == 'hard threshold':
        prob_thres = args.prob_thres
        neighborhood_cell_idx = get_neighborhood_cells_thres(prob_i, prob_thres=prob_thres)
    
    # The neighborhood definition method doesn't exist
    else:
        logging.warning('Neighborhood definition method does not exit')
        sys.exit()
    
    return neighborhood_cell_idx, prob_thres