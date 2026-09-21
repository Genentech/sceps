import sys, os, pickle, gc, random, time, warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import scipy
import scipy.stats

from .utils import *
from .scdata import *
from .neighborhood import *
from .estimation import *
from .stats_test import *

from tqdm import tqdm


def select_disease_gwas_genes(args, gene_list):
    """
    Selects disease GWAS genes for scEPS analysis. Apply automatic gene selection based
    on MAGMA FDR, unless the user specifies manual selection and provides the number
    of GWAS genes to select.

    Args:
        args (argparse.Namespace): An object containing user specifications for selecting
        disease GWAS genes.

        gene_list (pd.DataFrame): A DataFrame containing gene symbols and MAGMA Z-stat.

    Returns:
        pd.Series: A pandas Series containing the selected disease GWAS genes.
    """

    # Sort the gene list based on MAGMA Z-stat in descending order 
    gene_list = gene_list.sort_values(by=['ZSTAT'], ascending=False)
    gene_list = gene_list.reset_index(drop=True)

    # If the user chooses to manually select the genes, return the GWAS genes specified by the user
    if args.auto_gene_selection == False:
        return gene_list.iloc[0:args.num_gwas_genes][args.gene_id_col]

    # Calculate FDR based on MAGMA Z-stat
    gene_list['P'] = 1.0 - scipy.stats.norm.cdf(gene_list['ZSTAT'])
    gene_list['FDR'] = scipy.stats.false_discovery_control(gene_list['P'].values)

    # Find the indices of genes passing the MAGMA FDR threshold
    use_idx = np.where(gene_list['FDR']<args.magma_fdr_thres)[0]
    
    # Cap the maximum number of GWAS genes
    if use_idx.shape[0] > args.max_num_gwas_genes:
        use_idx = range(0, args.max_num_gwas_genes)
    
    # Cap the minimum number of GWAS genes
    elif use_idx.shape[0] < args.min_num_gwas_genes:
        use_idx = range(0, args.min_num_gwas_genes)

    # Return the gene names
    selected_genes = gene_list.iloc[use_idx][args.gene_id_col]
    logging.info('Selected {} genes for analysis'.format(selected_genes.shape[0]))

    # Check if the list of selected genes is empty
    if selected_genes.shape[0] == 0:
        logging.info('Empty list of GWAS genes returned')
        sys.exit()

    return selected_genes


def get_sceps_stats(args, adata, gwas_genes, start_idx=None, stop_idx=None, focal_cells=set()):

    # Set start_idx and stop_idx
    if start_idx is None or stop_idx is None:
        start_idx, stop_idx = 0, adata.shape[0]

    # Get the target neighborhood abundance matrix
    trans_mat, init_nsteps, target_nam = get_target_nam(args, adata)

    # Extract the list of donors and donor phenotypes from adata
    donor_id_vec, df_donor_pheno = extract_donor_pheno(args, adata)

    # Generate a random permutation of samples, if the user chooses to use permutation for testing
    permuted_donor_id = permute_donor_id(args, donor_id_vec)

    # Initialize the results
    all_reg_out = []                # Estimated scEPS statistics
    all_reg_out_testing = []        # Permutation based test statistics, if specified by the user
    all_neighborhood_cells = []     # Record the cells in each neighborhood, for diagnosis purposes
    all_neighborhood_donors = []    # Record the donors in each neighborhood, for diagnosis purposes
    nam_realized = []               # Record the realized NAM, for diagnosis purposes

    # iterate through cell neighborhood
    for idx in tqdm(range(start_idx, stop_idx)):

        # get cell id
        cell_id = adata.obs[args.cell_id_col].values[idx]
        cell_donor_id = adata.obs[args.donor_id_col].values[idx]

        # If focal cell not specified, sample cells based on sample_frac
        if len(focal_cells) == 0:
            rand_num = random.random()
            if rand_num > args.sample_frac:
                continue
        # Skip cells not in the set of focal cells
        else:
            if cell_id not in focal_cells:
                continue

        # Log the current cell neighborhood
        if args.verbose == True:
            logging.info('Processing cell {}'.format(cell_id))

        # Get scEPS statistics for a single cell neighborhood
        sceps_stats_nbhood = get_sceps_stats_neighborhood(args, idx, adata, gwas_genes,
            donor_id_vec, df_donor_pheno, trans_mat, init_nsteps, target_nam)
        reg_out = sceps_stats_nbhood['reg_out']
        reg_out_testing = sceps_stats_nbhood['reg_out_testing']

        # Append output
        if reg_out is not None:
            reg_out['CELL_IDX'] = idx
            reg_out['CELL'] = cell_id
            reg_out['DONOR_ID'] = cell_donor_id
            neighborhood_cell_idx = sceps_stats_nbhood['neighborhood_cell_idx']
            reg_out['NEIGHBORHOOD_SIZE'] = neighborhood_cell_idx.shape[0]
            reg_out['NUM_DONOR'] = sceps_stats_nbhood['num_donor']
            reg_out['STEP_SIZE'] = sceps_stats_nbhood['nstep']
            reg_out['ELAPSED_TIME'] = sceps_stats_nbhood['elapsed_time']
            if sceps_stats_nbhood['prob_thres'] is not None:
                reg_out['TRANSIT_PROB_THRES'] = sceps_stats_nbhood['prob_thres']
            all_reg_out.append(pd.DataFrame(reg_out))

            if args.verbose == True:
                logging.info('Elapsed time: {}s'.format(sceps_stats_nbhood['elapsed_time']))
                logging.info('Current OMEGA_DIFF: {}'.format(np.mean(reg_out['OMEGA_DIFF'])))

            # Save scEPS statistics obtained based on permutation
            if reg_out_testing is not None:
                reg_out_testing['CELL_IDX'] = idx
                reg_out_testing['CELL'] = cell_id
                reg_out_testing['DONOR_ID'] = cell_donor_id
                reg_out_testing['NEIGHBORHOOD_SIZE'] = neighborhood_cell_idx.shape[0]
                reg_out_testing['NUM_DONOR'] = sceps_stats_nbhood['num_donor']
                reg_out_testing['STEP_SIZE'] = sceps_stats_nbhood['nstep']
                all_reg_out_testing.append(pd.DataFrame(reg_out_testing))

            # Save realized nam
            nam_realized_row = np.zeros(df_donor_pheno.shape[0])
            donor_idx = df_donor_pheno.loc[sceps_stats_nbhood['neighborhood_donors']]['INDEX']
            nam_realized_row[donor_idx] = sceps_stats_nbhood['neighborhood_abundance']
            nam_realized.append([cell_id] + nam_realized_row.tolist())

            # Save cells and donors in the neighborhood
            neighborhood_idx_str = ','.join([str(nci) for nci in neighborhood_cell_idx.tolist()])
            neighborhood_cells = sceps_stats_nbhood['neighborhood_cells']
            neighborhood_cell_str = ','.join([str(nc) for nc in neighborhood_cells.tolist()])
            all_neighborhood_cells.append([idx, adata.obs.index[idx], \
                                            neighborhood_idx_str, neighborhood_cell_str])
            neighborhood_donors = sceps_stats_nbhood['neighborhood_donors']
            neighborhood_donor_str = ','.join(str(nd) for nd in neighborhood_donors.tolist())
            all_neighborhood_donors.append([idx, adata.obs.index[idx], neighborhood_donor_str])
        else:
            if args.verbose == True:
                logging.info('scEPS failed for cell {}'.format(idx))

    # Convert the results to DataFrame
    if len(all_reg_out) > 0:

        # Convert the scEPS results to DataFrame
        all_reg_out = pd.concat(all_reg_out, ignore_index=True)
        if len(all_reg_out_testing) > 0:
            all_reg_out_testing = pd.concat(all_reg_out_testing, ignore_index=True)
        else:
            all_reg_out_testing = None
        
        # Perform statistics testing at cell neighborhood level
        df_sceps_sigma = test_cell(all_reg_out, all_reg_out_testing, prefix='SIGMA_')
        df_sceps_omega = test_cell(all_reg_out, all_reg_out_testing, prefix='OMEGA_')

        # Convert the realized NAM to DataFrame
        nam_realized = pd.DataFrame(nam_realized)
        nam_realized.columns = ['CELL'] + df_donor_pheno[args.donor_id_col].tolist()
        nam_realized = nam_realized.set_index('CELL')

        # Convert the neighborhood cells to DataFrame
        all_neighborhood_cells = pd.DataFrame(all_neighborhood_cells)
        all_neighborhood_cells.columns = ['CELL_IDX', 'CELL', 'NEIGHBORHOOD_CELL_IDX', 'NEIGHBORHOOD_CELL']

        # Convert the neighborhood donors to DataFrame
        all_neighborhood_donors = pd.DataFrame(all_neighborhood_donors)
        all_neighborhood_donors.columns = ['CELL_IDX', 'CELL', 'DONOR_ID']

        # Make output dictionary
        out_dict = dict()
        out_dict['df_sceps_omega'] = df_sceps_omega                     # DataFrame storing the scEPS omega stats
        out_dict['df_sceps_sigma'] = df_sceps_sigma                     # DataFrame storing the scEPS sigma stats (for diagnosis)
        out_dict['all_reg_out'] = all_reg_out                           # DataFrame for raw scEPS stats (for diagnosis)
        out_dict['all_reg_out_testing'] = all_reg_out_testing           # DataFrame for raw scEPS stats based on permutation (for diagnosis)
        out_dict['nam_realized'] = nam_realized                         # DataFrame storing the realized NAM (for diagnosis)
        out_dict['all_neighborhood_cells'] = all_neighborhood_cells     # DataFrame storing the cells in neighborhoods (for diagnosis)
        out_dict['all_neighborhood_donors'] = all_neighborhood_donors   # DataFrame storing the donors in neighborhoods (for diagnosis)

        return out_dict
    
    # Return None if the result is empty
    return None


def get_sceps_stats_neighborhood(args, idx, adata, gwas_genes, donor_id_vec,
    df_donor_pheno, trans_mat, init_nsteps, target_nam):
    
    # Initialize results for the current cell neighborhood        
    reg_out, reg_out_testing, neighborhood_cell_idx = None, None, None
    start_time = time.time()

    # Increment nstep until max step size is reached
    for nstep in range(init_nsteps, args.max_step_size+1):

        # Try defining neighborhood using random walk step size, nstep
        if args.verbose == True:
            logging.info('Using step size: {}'.format(nstep))
        neighborhood_cell_idx, prob_thres = get_neighborhood_cells(args, idx,
            trans_mat, nstep, target_nam, donor_id_vec, df_donor_pheno)
        
        # Retry neighborhood definition if nstep fails
        if (neighborhood_cell_idx is None) or (prob_thres is None):
            continue
        
        # get data for the neighborhood
        neighborhood_cells = adata.obs.index[neighborhood_cell_idx]
        if neighborhood_cell_idx.shape[0] == 0:
            if args.verbose == True:
                logging.info('Empty neighborhood at cell {}, incrementing nstep'.format(idx))
            continue
        adata_neighbor = adata[neighborhood_cell_idx, :]
        if args.verbose == True:
            logging.info('Neighborhood size {}'.format(neighborhood_cell_idx.shape[0]))

        # get pseudobulk
        adata_pseudobulk = get_pseudobulk(args, adata_neighbor)
        neighborhood_donors = adata_pseudobulk.obs[args.donor_id_col]
        if args.verbose == True:
            logging.info('Number of donors {} for {} genes'.format(adata_pseudobulk.shape[0], \
                            adata_pseudobulk.shape[1]))
        if (adata_pseudobulk.shape[0]) < args.min_num_donor or (adata_pseudobulk.shape[1]==0):
            if args.verbose == True:
                logging.info('Number of donors too small at cell {}'.format(idx))
            continue

        # subsample cells if specified -- used for simulation purpose only
        if args.sample_frac_cell_in_neighborhood is not None:
            sc.pp.subsample(adata_neighbor, fraction=args.sample_frac_cell_in_neighborhood)
            adata_pseudobulk = get_pseudobulk(args, adata_neighbor)
            neighborhood_donors = adata_pseudobulk.obs[args.donor_id_col]
            if args.verbose == True:
                logging.info('*Simulation* {} cells from {} donors for {} genes, after subsampling cells'.format(\
                    adata_neighbor.shape[0], adata_pseudobulk.shape[0], adata_pseudobulk.shape[1]))

        # prepare pseudobulk adata
        adata_pseudobulk.X = adata_pseudobulk.X - np.mean(adata_pseudobulk.X, axis=0)
        adata_pseudobulk.var['GENE_INDEX'] = np.array(range(adata_pseudobulk.var.shape[0])).astype(int)
        adata_pseudobulk.var = adata_pseudobulk.var.set_index(args.gene_id_col, drop=False)

        # create control gene set generators
        gwas_genes_nb = gwas_genes[gwas_genes.isin(adata_pseudobulk.var[args.gene_id_col])]
        control_genes_generator = next_control_genes(args, adata_pseudobulk, \
                                    gwas_genes_nb, nbins=args.num_expr_bins, nset=args.num_control_gene_set)

        # create control gene set for permutation based testing
        if args.num_control_gene_set_testing > 0:
            control_genes_generator_testing = next_control_genes(args, adata_pseudobulk, \
                                    gwas_genes_nb, nbins=args.num_expr_bins, nset=args.num_control_gene_set_testing)

        # using method of moment to estimate parameters
        if args.skip_regression == False:
            reg_out = regression(args, adata_pseudobulk, gwas_genes_nb, control_genes_generator)

            # get sceps statistics for permutation based testing
            if args.num_control_gene_set_testing > 0:
                reg_out_testing = regression(args, adata_pseudobulk, gwas_genes_nb,
                    control_genes_generator_testing, shuffle_pheno=True, permuted_donor_id=permuted_donor_id)

        # regression succeeds     
        if reg_out is not None:
            break
    
    # Record the time used to analyze one cell neighborhood
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Make output dictionary
    out_dict = dict()
    out_dict['reg_out'] = reg_out
    out_dict['reg_out_testing'] = reg_out_testing
    out_dict['neighborhood_cell_idx'] = neighborhood_cell_idx
    
    # Analysis failed, all other variables not created
    if reg_out is None:
        return out_dict

    # Analysis succeeded, all other variables created
    out_dict['neighborhood_cells'] = neighborhood_cells
    out_dict['num_donor'] = adata_pseudobulk.shape[0]
    out_dict['elapsed_time'] = elapsed_time
    out_dict['nstep'] = nstep
    out_dict['prob_thres'] = prob_thres
    out_dict['neighborhood_donors'] = adata_pseudobulk.obs[args.donor_id_col]
    out_dict['neighborhood_abundance'] = adata_pseudobulk.obs['sceps.num_cells'].values

    return out_dict


def save_sceps_stats(args, sceps_out, out_suffix):

    # Directly return nothing if sceps_out is None
    if sceps_out is None:
        return

    # Save the main scEPS omega stats
    sceps_out['df_sceps_omega'].to_csv(args.out+'.{}.sceps.omega.txt.gz'.format(out_suffix),
        sep='\t', index=False, na_rep='NA', float_format='%.5g')

    # Save the scEPS sigma stats
    if args.save_sceps_sigma:
        sceps_out['df_sceps_sigma'].to_csv(args.out+'.{}.sceps.sigma.txt.gz'.format(out_suffix),
            sep='\t', index=False, na_rep='NA', float_format='%.5g')

    # Save the raw scEPS stats and scEPS stats for testing based on permutation
    if args.save_sceps_raw == True:
        sceps_out['all_reg_out'].to_csv(args.out+'.{}.sceps.raw.txt.gz'.format(out_suffix),
                sep='\t', index=False, na_rep='NA', float_format='%.5g')
        if sceps_out['all_reg_out_testing'] is not None:
            sceps_out['all_reg_out_testing'].to_csv(args.out+'.{}.sceps.raw.testing.txt.gz'.format(out_suffix),
                sep='\t', index=False, na_rep='NA', float_format='%.5g')

    # Save the realized NAM
    if args.save_realized_nam == True:
        sceps_out['nam_realized'].to_csv(args.out + '.{}.realized_nam.txt.gz'.format(out_suffix), sep='\t')

    # Save the cells in the neighborhoods
    if args.save_neighborhood_cells == True:
        sceps_out['all_neighborhood_cells'].to_csv(args.out + '.{}.neighborhood_cells.txt.gz'.format(out_suffix),
            sep='\t', index=False, na_rep='NA')

    # Save the donors in the neighborhoods
    if args.save_neighborhood_donors == True:
        sceps_out['all_neighborhood_donors'].to_csv(args.out + '.{}.neighborhood_donors.txt.gz'.format(out_suffix),
            sep='\t', index=False, na_rep='NA')


def create_sceps_default_args(adata, gene_list):

    # Calling get_command_line() to assign default values to args
    args = get_command_line()
    
    # Assign user input values for required arguments to args
    args.adata = adata
    args.gene_list = gene_list

    # Also set a few recommended flag to be true
    args.scale_pheno = True
    args.scale_pheno_neighborhood = True
    args.auto_gene_selection = True

    return args