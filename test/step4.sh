python ../sceps_corr.py \
    --adata ./input/test_scdata.h5ad \
    --sceps-result ./output/step3.sceps.omega.txt.gz \
    --cell-type-col CellType \
    --focal-cell-type CellType_1 \
    --min-num-nonzero 3 \
    --out ./output/step4.CellType_1

python ../sceps_corr.py \
    --adata ./input/test_scdata.h5ad \
    --sceps-result ./output/step3.sceps.omega.txt.gz \
    --min-num-nonzero 3 \
    --out ./output/step4
