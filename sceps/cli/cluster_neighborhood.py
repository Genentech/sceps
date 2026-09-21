import argparse, sys, glob, os
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
from tqdm import tqdm
import logging, random
import scanpy as sc
import anndata as ad
from sklearn.cluster import MiniBatchKMeans

from src.utils import EPS
from src.neighborhood import *

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    
    # get command line
    args = get_command_line()

    # set seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    # define clusters for bootstrap
    logging.info("Loading single-cell data from {}".format(args.adata))
    adata = sc.read_h5ad(args.adata)
    
    # If args.cell_id_col not in adata.obs.columns, we use adata.obs.index as replacement
    obs = adata.obs.copy()
    if args.cell_id_col not in obs.columns:
        if args.cell_id_col == '':
            args.cell_id_col = 'sceps.cell_index'
        obs[args.cell_id_col] = obs.index.copy()
        adata.obs = obs

    # define clusters for bootstrap
    logging.info("Identifying independent groups of cell neighborhoods")
    adata.obs['sceps.neighborhood_cluster'] = get_nbhood_clusters(args, adata)
    
    # Create and save a data frame mapping cells to their assigned cluster
    df_out = adata.obs[[args.cell_id_col, 'sceps.neighborhood_cluster']]
    df_out.to_csv(args.out+'.txt.gz', sep='\t', index=False)


def get_nbhood_clusters(args, adata):

    # check if neighborhood is already calculated
    has_knn = False
    av = ad.__version__
    if type(av) == str:
        av = version.parse(av)
    if av < version.parse("0.7.2"):
        if "neighbors" in adata.uns:
            if "connectivities" in adata.uns["neighbors"]:
                has_knn = True
    else:
        if "connectivities" in adata.obsp:
            has_knn = True

    # calculate knn only when necessary
    if has_knn == False:
        sc.pp.neighbors(adata, use_rep=args.neighbors_use_rep)
    
    # calculate nam matrix
    adj_mat = get_connectivity(adata)
    trans_mat = get_transition_matrix(adj_mat)
    _, nam_matrix = choose_random_walk_nsteps(args, trans_mat, adata.obs)
    nam_matrix = nam_matrix.values
    nam_matrix = (nam_matrix - nam_matrix.mean(axis=0)) / (nam_matrix.std(axis=0)+EPS)

    # run kmeans
    kmeans = MiniBatchKMeans(n_clusters=args.num_kmeans_cluster,
        batch_size=int(0.05*adata.shape[0]), random_state=args.seed, n_init='auto')
    
    return kmeans.fit(nam_matrix).labels_


def get_command_line():
 
    # Create the parser
    parser = argparse.ArgumentParser(description="This tool aggregates scEPS statistics across groups of cell neighborhoods.")

    # Add the arguments
    parser.add_argument('--adata', type=str, required=False,
        help="""Used to specify the input single-cell RNA-seq data in h5ad format.""" \
        """This should be same single-cell data as analyzed by scEPS. However, the user """ \
        """may remove adata.X to reduce memory usage, as the clustering tool only requires k-NN """ \
        """graph for the cells.""")

    parser.add_argument('--cell-id-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents cell IDs in the adata.obs """ \
        """data frame of the single-cell data. If left empty, scEPS will use what's in """ \
        """adata.obs.index as cell IDs.""")

    parser.add_argument('--donor-id-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents donor IDs in the adata.obs """ \
        """data frame of the single-cell data.""")

    parser.add_argument('--neighbors-use-rep', type=str, required=False, default='X_pca_harmony',
        help="""Used to specify the cell embedding (e.g., PCA, scVI embeddings, etc.) used to """ \
        """construct the k-NN graph.""")

    parser.add_argument('--num-kmeans-cluster', type=int, required=False, default=50,
        help="""Used to specify the desired number of clusters (i.e., approximately independent """ \
        """blocks of cell neighborhoods).""")

    parser.add_argument('--seed', type=int, required=False, default=0,
        help="""Used to specify the seed for the random number generator. (This is set to 0, by default)""")

    parser.add_argument('--out', type=str, required=False,
        help="""Used to specify the output file name.""")

    # Execute the parse_args() method
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    main()
