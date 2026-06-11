python ../sceps_aggregate.py \
    --prefix "./output/step1.*.txt.gz" \
    --adata ./input/test_scdata.h5ad \
    --neighborhood-clusters ./output/step2.txt.gz \
    --cell-type-col CellType \
    --out ./output/step3

python ../sceps_aggregate.py \
    --prefix "./output/step1.*.txt.gz" \
    --adata ./input/test_scdata.h5ad \
    --neighborhood-clusters ./output/step2.txt.gz \
    --out ./output/step3
