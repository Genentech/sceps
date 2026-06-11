import argparse, logging, gzip, os, sys
import scanpy as sc
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

EPS = 1e-16
NUM_DISATT_FACTOR_BIN = 5

def get_command_line():
 
    # Create the parser
    parser = argparse.ArgumentParser(description="""scEPS integrates GWAS and single-cell disease """ \
    """cell atlas data to identify disease-associated cell neighborhoods.""")


    # Input related options
    parser.add_argument('--adata', type=str, required=False,
        help="""Used to specify the input single-cell RNA-seq data in h5ad format. """ \
        """The single-cell data should consist of multiple donors, showing variations in their disease phenotypes.""")

    parser.add_argument('--cell-type-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents cell types in """ \
        """the adata.obs data frame of the single-cell data. This is an emtpy string by default.""")

    parser.add_argument('--focal-cell-type', type=str, required=False, default='',
        help="""Used to specify the focal cell type to analyze. If specified, scEPS will """ \
        """only analyze cell neighborhoods in the specified cell type. This is an emtpy string by default.""")

    parser.add_argument('--gene-id-col', type=str, required=False, default='',
        help="""Used to specify the name of the column that represents gene symbols/IDs in the adata.var data """ \
        """frame of the single-cell data. If left empty, scEPS will use what's in adata.var.index as gene symbols/IDs.""")
    
    parser.add_argument('--cell-id-col', type=str, required=False, default='', 
        help="""Used to specify the name of the column that represents cell ID in the adata.obs data """ \
        """frame of the single-cell data. If left empty, scEPS will use what's in adata.obs.index as cell IDs.""")
    
    parser.add_argument('--donor-id-col', type=str, required=False,
        help="""Used to specify the name of the column that represents donor IDs in the adata.obs """ \
        """data frame of the single-cell data.""")

    parser.add_argument('--pheno-file', type=str, required=False, default='',
        help="""Used to specify an external text file containing donor phenotypes. The text file """ \
        """should contain a donor ID column, same as the one used in the single-cell data, and one or """ \
        """multiple columns representing donor phenotypes.""")

    parser.add_argument('--batch-key', type=str, required=False, default='',
        help="""Used to specify batch information for the donors. If specified, scEPS will permute the """ \
        """donors within each batch, if using permutation-based testing procedure.""")

    parser.add_argument('--gene-list', type=str, required=False,
        help="""Used to specify a text file containing the gene-level MAGMA associaiton statistics. """ \
        """The text file should contain a GENE column and a ZSTAT column.""")

    parser.add_argument('--control-gene-list', type=str, required=False, default='',
        help='Used to specify a text file containing a list of pre-determined control gene list. """ \
        """By default, scEPS will automatically randomly select a list of control genes.')

    parser.add_argument('--focal-cells', type=str, required=False, default='',
        help="""Used to specify a text file containing a list of focal cell neighborhoods to analyze. """ \
        """If specified, scEPS will only analyze cells in this file.""")
    

    # Parallelization related options
    parser.add_argument('--total-num-job', type=int, required=False, default=None,
        help="""Used to specify the total number of parallel jobs used for the scEPS analysis.""")
    
    parser.add_argument('--job-idx', type=int, required=False, default=None,
        help="""Used to specify the index of the parallel job for analyzing a subset of the single-cell data.""")

    parser.add_argument('--start-idx', type=int, required=False, default=None,
        help="""Used to specify the starting index (inclusive) of the cell neighborhood to analyze.""")
    
    parser.add_argument('--stop-idx', type=int, required=False, default=None,
        help='Used to specify the stopping index (exclusive) of the cell neighborhood to analyze.')


    # QC related options
    parser.add_argument('--min-expr-var-thres', type=float, required=False, default=1e-4,
        help='Used to specify the minimum threshold on the variance of gene expression across donors. """ \
        """By default, genes with variance across donors less than 1e-4 are filtered.')

    parser.add_argument('--max-expr-var-thres', type=float, required=False, default=100.0,
        help='Used to specify the maximum threshold on the variance of gene expression across donors. """ \
        """By default, genes with variance across donors greater than 100.0 are filtered.')


    # Neighborhood definition related options
    parser.add_argument('--min-num-donor', type=int, required=False, default=8,
        help="""Used to specify the minimum number of donors required in each cell neighborhood. """ \
        """By default, we set this number to 8.""")

    parser.add_argument('--min-frac-donor', type=float, required=False, default=0.333333,
        help="""Used the specify the fraction of total number of donors in the single-cell data """ \
        """required to be represented in each cell neighborhood. By default, we set this to 0.333333 """ \
        """(i.e., approximately 1/3).""")

    parser.add_argument('--min-num-cell', type=int, required=False, default=5,
        help="""Used to specify the minimum number of cells required in each cell neighborhood. """ \
        """By default, we set this to 5.""")

    parser.add_argument('--neighborhood-definition-method', type=str, required=False, default='nam',
        help="""Used to specify the neighborhood definition method. The user can choose from "nam",""" \
        """ "soft threshold", and "hard threshold". By default scEPS uses the "nam" approach, defining """ \
        """cell neighborhoods based on the neighborhood abundance matrix.""")

    parser.add_argument('--check-pheno-in-neighborhood', required=False, default=False, action='store_true',
        help="""Is specified, scEPS will ensure that the variance of the phenotype across donors in the """ \
        """cell neighborhood is non-zero.""")

    parser.add_argument('--max-step-size', type=int, required=False, default=15,
        help="""Used to specify the maximum number of random walk step size to test. By default, """ \
        """this is set to 15.""")

    parser.add_argument('--prob-thres', type=float, required=False, default=0.0001,
        help="""Used to specify the threshold for the transition probability, if --neighborhood-definition-method """ \
        """is set to "hard threshold".""")


    # GWAS gene selction related options
    parser.add_argument('--auto-gene-selection', required=False, default=False, action='store_true',
        help="""Is specified, scEPS will automatically select GWAS genes based on the gene-level MAGMA association """ \
        """statistics. The user can also bypass auto gene selection, using the --num-gwas-genes flag to specify the number """
        """of top GWAS genes.""")

    parser.add_argument('--num-gwas-genes', type=int, required=False, default=500,
        help="""Used to specify the number of top GWAS genes, if the user chooses to manually select the top GWAS """ \
        """genes to analyze.""")
    
    parser.add_argument('--magma-fdr-thres', type=float, required=False, default=0.05,
        help='Used to specify the FDR threshold to select GWAS genes, if the user chooses to use the automatically """ \
        """approach for selecting GWAS genes.')

    parser.add_argument('--min-num-gwas-genes', type=int, required=False, default=500,
        help='Used to specify the minimum number of GWAS genes included in the scEPS analysis. By default, """ \
        """this is set to 500.')

    parser.add_argument('--max-num-gwas-genes', type=int, required=False, default=2000,
        help="""Used to specify the maximum number of GWAS genes included in the scEPS analysis. """ \
        """By default, this is set to 2,000.""")


    # Control gene selection related options
    parser.add_argument('--enrich-low-score-control-genes', required=False, default=False, action='store_true',
        help="""If specified, scEPS will preferentially select control genes with the lowest gene-level """ \
        """MAGMA association statistics.""")

    parser.add_argument('--include-gwas-genes-in-control-genes', required=False, default=False, action='store_true',
        help="""If specified, the set of GWAS genes will also be included in the sampling of control genes.""")

    parser.add_argument('--num-expr-bins', type=int, required=False, default=10,
        help="""Used to specify the number of bins for binning all the genes based on their mean expression. """\
        """The default is 10.""")

    parser.add_argument('--num-control-gene-set', type=int, required=False, default=1,
        help="""Used to specify the number of control gene sets. This should almost always be set to 1.""")
    
    parser.add_argument('--num-control-gene-set-testing', type=int, required=False, default=0,
        help="""Used to specify the number of control gene sets use for permutation-based testing. """ \
        """By default, this number is set to 0.""")

    parser.add_argument('--match-expr-std', required=False, default=False, action='store_true',
        help="""If specified, scEPS will also match the standard deviation of the gene expression of """ \
        """the control genes with that of GWAS genes.""")

    
    # Estimation related options
    parser.add_argument('--pheno', type=str, required=False,
        help="""Used to specify the column that represents phenotype in the adata.obs data frame of """ \
        """the single-cell data or in the text file external to the single-cell data.""")

    parser.add_argument('--covar', type=str, required=False, default='', 
        help="""Used to specify a list of columns in adata.obs that represent covariates to adjust. """\
        """Multiple covariates names should be separated by commas. We recommend to include sceps.prop_cells_local """ \
        """in this flag, to adjust for the impact of neighborhood abundance on the phenotype. An example specification """ \
        """for --covar could be --covar "age,gender,sceps.prop_cells_local".""")

    parser.add_argument('--num-var-comp', type=int, required=False, default=3,
        help="""Used to specify the number of variance components. By default, this is set to 3. And scEPS will """ \
        """decompose the phenotypic variance into GWAS genes, control genes, and remaining (non-GWAS/control) genes. """ \
        """If this is set to 2 (not recommended), scEPS will decompose the phenotypic variance into GWAS genes and """ \
        """control genes.""")

    parser.add_argument('--num-bootstrap-disattenuation', type=int, required=False, default=100,
        help="""Used to specify the number of bootstap samples used to estimate the disattenuation factor. """ \
        """By default, this is set to 100.""")

    parser.add_argument('--num-bootstrap-regression', type=int, required=False, default=1000,
        help="""Used to specify the number of bootstrap samples used to estimate the standard errors for the scEPS """ \
        """statistics. By default, this is set to 1,000.""")

    parser.add_argument('--use-ols', required=False, default=False, action='store_true',
        help="""If specified, scEPS will use ordinary least square regression to estimate the variance components. """ \
        """By default, scEPS uses weighted least square regression.""")

    parser.add_argument('--use-analytical-stderr', required=False, default=False, action='store_true',
        help="""If specified, scEPS will use analytical approach to calculate standard errors. This is a lot """ \
        """faster than using the bootstrap approach, and can be used if computational resource is limited.""")

    parser.add_argument('--scale-pheno', required=False, default=False, action='store_true',
        help="""If specified (recommended), scEPS will mean-center and scale the phenotype to have unit variance, """ \
        """at the global level (i.e., before variance decomposition at individual cell neighborhood level). """ \
        """Otherwise, scEPS will only mean-center the phenotype.""")

    parser.add_argument('--scale-pheno-neighborhood', required=False, default=False, action='store_true',
        help="""If specified (recommended), scEPS will mean-center and scale the phenotype to have unit variance, """ \
        """at each cell neighborhood, on top of the global level standardization. Otherwise, scEPS will only """ \
        """mean-center the phenotype.""")

    parser.add_argument('--exclude-diag', required=False, default=False, action='store_true',
        help="""If specified (not recommended), scEPS will only use cross products of the phenotypes between """ \
        """different individuals to estimate the variance components. By default, scEPS includes squared phenotypes """ \
        """of the same individuals.""")

    parser.add_argument('--exclude-intercept', required=False, default=False, action='store_true',
        help="""If specified (not recommended), scEPS will not include an intercept term in the regression. """ \
        """By default, an intercept term is included in the regression.""")

    parser.add_argument('--no-disattenuation', required=False, default=False, action='store_true',
        help="""If specified (not recommended), scEPS will not apply disattenuation factor on the estimates.""")

    parser.add_argument('--skip-regression', required=False, default=False, action='store_true',
        help="""If specified, scEPS will skip the estimation procedure, but will record summary statistics """ \
        """(e.g., number of cells/donors) for each neighborhood. This is typically used for diagnosis purposes only.""")

    parser.add_argument('--permute-within-batch', required=False, default=False, action='store_true',
        help="""If specified, scEPS will run permutation within each batch of donors. This is only used if the """ \
        """user chooses to use permutation to obtain p-values.""")


    # Output related options
    parser.add_argument('--save-neighborhood-cells', required=False, default=False, action='store_true',
        help="""If specified, scEPS will save the cells in each neighborhood as a text file.""")
    
    parser.add_argument('--save-neighborhood-donors', required=False, default=False, action='store_true',
        help="""If specified, scEPS will save the donors in each neighborhood as a text file.""")

    parser.add_argument('--save-target-nam', required=False, default=False, action='store_true',
        help="""If specififed, scEPS will save the target neighborhood abundance matrix as a text file.""")

    parser.add_argument('--save-realized-nam', required=False, default=False, action='store_true',
        help="""If specified, scEPS will save the realized neighborhood abundance matrix as a text file.""")

    parser.add_argument('--save-sceps-raw', required=False, default=False, action='store_true',
        help="""If specified, scEPS will save the intermediate raw results as a text file.""")

    parser.add_argument('--save-sceps-sigma', required=False, default=False, action='store_true',
        help="""If specified, scEPS will save the sigma statistics as a text file.""")

    parser.add_argument('--out', type=str, required=False,
        help="""Used to specify the prefix for the output file.""")


    # Miscellaneous options
    parser.add_argument('--seed', type=int, required=False, default=0,
        help="""Used to specify the seed for the random number generator. By default, this is set to 0.""")

    parser.add_argument('--sample-frac', type=float, required=False, default=1.0,
        help="""Used to down-sample the number of cell neighborhoods to analyze. By default, this is set to 1.0.""")

    parser.add_argument('--force-rerun', required=False, default=False, action='store_true',
        help="""Used to force rerun the analysis. By default, if scEPS detects the output exists, it will not """ \
        """perform the analysis.""")
    
    parser.add_argument('--verbose', required=False, default=False, action='store_true',
        help="""If specified, scEPS will print intermediate output to the screen.""")


    # Simulation related options -- should not be used in analysis of real data
    parser.add_argument('--sample-frac-cell-in-neighborhood', type=float, required=False, default=None,
        help="""Used to down-sample cells within a neighborhood. By default, this is set to 1.0. This flag """ \
        """is used to assess the robustness of scEPS in simulations, and should not be used in analysis of real traits.""")


    # execute the parse_args() method
    args = parser.parse_args()
    
    return args


def check_gzip(fnm):

    chunksize=10000000 # 10 Mbytes

    if os.path.exists(fnm) == False:
        return False

    ok = True
    with gzip.open(fnm, 'rb') as f:
        try:
            while f.read(chunksize) != b'':
                pass
        # the file is not a gzip file.
        except gzip.BadGzipFile:
            ok = False
        # EOFError: Compressed file ended before the end-of-stream marker was reached
        # a truncated gzip file.
        except EOFError:
            ok = False

    return ok


def check_output(args, out_suffix):
    """
    Return true if output already exists
    """
    ok = check_gzip(args.out+'{}.sceps.omega.txt.gz'.format(out_suffix))

    return ok


def get_start_stop_index(args, ntasks):

    # By default, perform all the tasks
    start_idx, stop_idx = 0, ntasks

    # Calculate the starting and stopping indices based on the number of jobs provided
    if (args.total_num_job is not None) and ((args.job_idx is not None)):
        job_batch_size = ntasks // args.total_num_job
        remainder = ntasks % args.total_num_job
        if args.job_idx < remainder:
            start_idx = args.job_idx * (job_batch_size + 1)
            stop_idx = start_idx + job_batch_size + 1
        else:
            start_idx = args.job_idx * job_batch_size + remainder
            stop_idx = start_idx + job_batch_size
    
    # Use user specified starting and stopping indices
    if (args.start_idx is not None) and (args.stop_idx is not None):
        start_idx = max(args.start_idx, 0)
        stop_idx = min(args.stop_idx, ntasks)
    
    # Check to make sure start_idx is less than stop_idx
    if start_idx >= stop_idx:
        sys.exit()
    
    return start_idx, stop_idx


def load_gene_stats(args):

    gene_list = pd.read_table(args.gene_list, delim_whitespace=True)

    if ('GENE' not in gene_list.columns) or ('ZSTAT' not in gene_list.columns):
        logging.info('Missing GENE or ZSTAT column in the gene list file')
        sys.exit()

    gene_list = gene_list[['GENE', 'ZSTAT']]
    gene_list.columns = [args.gene_id_col, 'ZSTAT']

    return gene_list


def load_focal_cells(args):

    focal_cells = set()
    if args.focal_cells != '':
        df_focal_cells = pd.read_table(args.focal_cells)
        focal_cells = set(df_focal_cells[args.cell_id_col].values.tolist())
    
    return focal_cells