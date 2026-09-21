sceps-aggregate \
    --prefix "./output/step1.*.txt.gz" \
    --adata ./input/test_scdata.h5ad \
    --neighborhood-clusters ./output/step2.txt.gz \
    --cell-type-col CellType \
    --out ./output/step3

sceps-aggregate \
    --prefix "./output/step1.*.txt.gz" \
    --adata ./input/test_scdata.h5ad \
    --neighborhood-clusters ./output/step2.txt.gz \
    --out ./output/step3
