python ../sceps.py \
    --adata ./input/test_scdata.h5ad \
    --donor-id-col Donor \
    --gene-list ./input/test_magma.txt \
    --auto-gene-selection \
    --pheno Pheno \
    --scale-pheno \
    --scale-pheno-neighborhood \
    --out ./output/step1

