# Path to the MAGMA tool
MAGMA=<path to the MAGMA package>/magma_v1.10_static/magma


# Mapping SNP to genes using MAGMA
WINDOW=<upstream window size>,<downstream window size>	# By default, WINDOW=10,10
SNP_LOC=<prefix to a reference panel plink file across all SNPs>
GENE_LOC=<a gene annotation file>			# Provide Entrez ID, chromosome number, start base pair, stop base pair, strand, gene symbol
ANNOT_OUT=<MAGMA annotation output file name>

$MAGMA \
    --annotate window=$WINDOW \
    --snp-loc $SNP_LOC.bim \
    --gene-loc $GENE_LOC \
    --out $ANNOT_OUT


# Get gene-level association statistics
PVAL=<GWAS summary statistics data file>		# Path to the SNP-level GWAS summary statistics file
USE=<SNP ID column name>,<p-value column name>		# Provide the SNP ID and p-value column in the GWAS summary statsitics file
NCOL=<sample size colname name>				# Provide the sample size column in the GWAS summary statistics file
GS_OUT=<output file name>

$MAGMA \
    --bfile $SNP_LOC \
    --pval $PVAL use=$USE ncol=$NCOL \
    --gene-annot $ANNOT_OUT.genes.annot \
    --out $GS_OUT
