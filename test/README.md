To run the test case, please first run
```
python generate_test_data.py
```

This will create 2 files under ./input:
1. test_scdata.h5ad -- a simulated single-cell data with 600 cells of 10 different cell types across 20000 genes for 30 donors
2. test_magma.txt -- a text file containing a simulated MAGMA result

Please then run:
1. step1.sh -- for estimating scEPS statistics for individual cell neighborhoods
2. step2.sh -- for defining approximately independently cell neighborhood blocks
3. step3.sh -- for aggregating scEPS statistics across cell types and across all cells
4. step4.sh -- for calculating the correlation between scEPS statistics and gene expression across CellType_1 and across all cells

In ./output, please find:
1. Output from running step1.sh: step1.0-600.sceps.omega.txt.gz, step1.0-600.log
2. Output from running step2.sh: step2.txt.gz
3. Output from running step3.sh: step3.CellType.sceps.omega.celltype.txt, step3.sceps.omega.celltype.txt, step3.sceps.omega.txt.gz
4. Output from running step4.sh: step4.CellType_1.0-20000.txt.gz, step4.0-20000.txt.gz
