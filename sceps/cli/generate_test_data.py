import argparse, os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import random


def main():

    # Get command line input
    args = get_command_line()

    # --- Configuration ---
    random.seed(args.seed)
    np.random.seed(args.seed)

    n_donors = args.num_donors
    cells_per_donor = args.cells_per_donor
    n_genes = args.num_genes
    n_cells = n_donors * cells_per_donor
    n_cell_types = args.num_cell_types
    k_neighbors = args.num_neighbors
    n_comps = args.num_pcs
    n_gwas_genes = args.num_gwas_genes

    # Create the output directory if it does not already exist
    os.makedirs(args.out_dir, exist_ok=True)
    out_adata = os.path.join(args.out_dir, args.adata_file_name)
    out_magma = os.path.join(args.out_dir, args.magma_file_name)

    # --- 1. Generate Gene Expression Data (Matrix X) ---
    print("Generating expression matrix...")
    X = np.random.normal(loc=0.0, scale=1.0, size=(n_cells, n_genes))

    # --- 2. Generate Metadata (Donor & Pheno) ---
    print("Generating metadata...")
    donor_ids = [f"Donor_{i}" for i in range(n_donors)]
    cell_donor_ids = np.repeat(donor_ids, cells_per_donor)

    # Phenotypes (Normal dist, mean 0, var 1) distributed across donors
    donor_pheno_values = np.random.normal(loc=0.0, scale=1.0, size=n_donors)
    donor_to_pheno_map = dict(zip(donor_ids, donor_pheno_values))
    cell_phenotypes = [donor_to_pheno_map[d] for d in cell_donor_ids]

    # Define the distinct cell type labels
    cell_type_labels = [f"CellType_{i}" for i in range(n_cell_types)]
    cell_assigned_types = np.random.choice(cell_type_labels, size=n_cells)

    # --- 3. Construct AnnData ---
    obs = pd.DataFrame({'Donor': cell_donor_ids, 'Pheno': cell_phenotypes, 'CellType': cell_assigned_types})
    obs.index = [f"Cell_{i}" for i in range(n_cells)]
    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # --- 4. Generate k-NN Graph ---
    print("Computing PCA...")
    # Reduce dimensions first for efficiency (standard workflow)
    sc.pp.pca(adata, n_comps=n_comps)

    print("Computing Neighbors...")
    # Compute the k-NN graph; n_neighbors sets the 'k' in k-NN
    sc.pp.neighbors(adata, n_neighbors=k_neighbors, n_pcs=n_comps)

    # --- 5. Map to Requested Location ---
    # Modern Scanpy stores connectivities in adata.obsp['connectivities'].
    # We manually map it to adata.uns['neighbors']['connectivities'], which is
    # the layout scEPS falls back to for older AnnData versions.
    adata.uns["neighbors"]["connectivities"] = adata.obsp["connectivities"]

    # --- 6. Verify ---
    print("\nData Summary:")
    print(adata)
    print("\nCheck connectivity matrix location:")
    print(f"Graph shape in .uns: {adata.uns['neighbors']['connectivities'].shape}")

    # --- 7. Save single cell data ---
    adata.write_h5ad(out_adata)

    # --- 8. Generate test MAGMA data ---
    df_magma = pd.DataFrame({'GENE': adata.var.index, 'ZSTAT': np.random.normal(size=n_genes)})
    df_magma['ZSTAT'].values[0:n_gwas_genes] = np.random.normal(loc=10, size=n_gwas_genes)
    df_magma.to_csv(out_magma, index=False, sep="\t")

    # --- 9. Report what was written ---
    print("\nWrote:")
    print("  {} -- simulated single-cell data with {} cells of {} cell types "
          "across {} genes for {} donors".format(out_adata, n_cells, n_cell_types, n_genes, n_donors))
    print("  {} -- simulated MAGMA gene-level association statistics".format(out_magma))


def get_command_line():

    # Create the parser
    parser = argparse.ArgumentParser(description="""This tool generates a small simulated single-cell """ \
        """data set and a matching simulated MAGMA gene-level association statistics file, for testing """ \
        """the scEPS workflow. The default settings reproduce the test data described in the scEPS """ \
        """repository. Note that the expression values and phenotypes are drawn independently, so no """ \
        """cell neighborhood is expected to show a genuine disease association.""")

    # Output related command line arguments
    parser.add_argument('--out-dir', type=str, required=False, default='./input',
        help="""Used to specify the directory to write the test data to (./input by default). """ \
        """The directory is created if it does not already exist.""")
    parser.add_argument('--adata-file-name', type=str, required=False, default='test_scdata.h5ad',
        help="""Used to specify the file name of the simulated single-cell data """ \
        """(test_scdata.h5ad by default).""")
    parser.add_argument('--magma-file-name', type=str, required=False, default='test_magma.txt',
        help="""Used to specify the file name of the simulated MAGMA association statistics """ \
        """(test_magma.txt by default).""")

    # Simulation related command line arguments
    parser.add_argument('--num-donors', type=int, required=False, default=30,
        help="""Used to specify the number of donors to simulate (30 by default). scEPS expects at """ \
        """least 20 donors showing variations in their disease phenotype.""")
    parser.add_argument('--cells-per-donor', type=int, required=False, default=20,
        help="""Used to specify the number of cells to simulate per donor (20 by default). The total """ \
        """number of cells is --num-donors multiplied by --cells-per-donor.""")
    parser.add_argument('--num-genes', type=int, required=False, default=20000,
        help="""Used to specify the number of genes to simulate (20,000 by default).""")
    parser.add_argument('--num-cell-types', type=int, required=False, default=10,
        help="""Used to specify the number of distinct cell type labels to assign (10 by default).""")
    parser.add_argument('--num-gwas-genes', type=int, required=False, default=1000,
        help="""Used to specify the number of genes given inflated MAGMA association statistics """ \
        """(1,000 by default).""")
    parser.add_argument('--num-neighbors', type=int, required=False, default=15,
        help="""Used to specify the number of nearest neighbors used to build the k-NN graph """ \
        """(15 by default).""")
    parser.add_argument('--num-pcs', type=int, required=False, default=50,
        help="""Used to specify the number of principal components used to build the k-NN graph """ \
        """(50 by default).""")
    parser.add_argument('--seed', type=int, required=False, default=42,
        help="""Used to specify the seed for the random number generator (42 by default).""")

    # Execute the parse_args() method
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
