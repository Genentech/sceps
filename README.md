# scEPS
**scEPS** (single-cell Expression exPlainability Statistics)

This repo contains the code of the method, **scEPS**, for integrating GWAS and single-cell disease cell atlas data to identify disease-associated cell neighborhoods. scEPS calculates a $d$ statistic at each cell neighborhood, representing the difference between the variance in disease explained by variations in the expression of each GWAS vs. each mean-expression matched control gene. An illustration of the scEPS method is shown below:

![scEPS illustration](https://raw.githubusercontent.com/Genentech/sceps/master/img/scEPS_overview.png  "Overview of the scEPS method")

# Reference

The current draft of the manuscript is available [here](https://www.medrxiv.org/content/10.64898/2026.06.26.26356714v1). The code we used to create the figures in the manuscript is available [here](https://github.com/Genentech/sceps_manuscript). We also implemented [CNA*](https://github.com/Genentech/cna_star), a simple extension of [CNA](https://github.com/immunogenomics/cna), for estimating the variance in disease attributable to variations in cell abundance at each cell neighborhood.

We also provide a web UI for visualizing the results in the scEPS manuscript [here](https://scepsresultexplorer.streamlit.app/).

# Manual
We provide a detailed manual of scEPS in the [Wiki page](https://github.com/Genentech/sceps/wiki).

# Installation

## Option 1: using pip

The easiest way to install scEPS is from [PyPI](https://pypi.org/project/sceps/):
```shell
pip install sceps
```

This installs the `sceps` Python package along with the four command-line tools described under [Usage](#usage).

## Option 2: using Anaconda or Miniforge
scEPS may also be installed into a dedicated environment through [Anaconda](https://www.anaconda.com/download) or [Miniforge](https://github.com/conda-forge/miniforge). To do this, please first install Anaconda or Miniforge on your machine. You may then install scEPS using the following commands:
``` shell
git clone git@github.com:Genentech/sceps.git
cd sceps
conda env create -f sceps.yml
conda activate sceps
pip install .
```

The `sceps.yml` file installs the dependencies through conda; the final `pip install .` installs scEPS itself and its command-line tools. Use `pip install -e .` instead if you intend to modify the scEPS source.

## Option 3: manually install required packages

The user may also manually install the required packages to run scEPS. scEPS requires Python 3.9 or newer and the following packages:

| Package | Minimum | Version pinned in `sceps.yml` |
| --- | --- | --- |
| [numpy](https://numpy.org/) | 1.23 | 1.26.2 |
| [pandas](https://pandas.pydata.org/) | 1.5 | 1.5.3 |
| [scipy](https://scipy.org/) | 1.9 | 1.13.1 |
| [anndata](https://anndata.readthedocs.io/) | 0.10 | 0.10.7 |
| [scanpy](https://scanpy.readthedocs.io/) | 1.10 | 1.10.3 |
| [scikit-learn](https://scikit-learn.org/) | 1.1 | 1.3.2 |
| [statsmodels](https://www.statsmodels.org/) | 0.13 | 0.14.5 |
| [tqdm](https://tqdm.github.io/) | 4.60 | 4.67.1 |
| [packaging](https://packaging.pypa.io/) | 20 | 25.0 |
| [matplotlib](https://matplotlib.org/) | 3.6 | 3.9.4 |
| [seaborn](https://seaborn.pydata.org/) | 0.12 | 0.13.2 |

These can be installed with a single command:
```shell
conda install -c conda-forge python=3.9 numpy=1.26.2 pandas=1.5.3 scipy=1.13.1 \
    anndata=0.10.7 scanpy=1.10.3 scikit-learn=1.3.2 statsmodels=0.14.5 \
    tqdm=4.67.1 packaging=25.0 matplotlib-base=3.9.4 seaborn=0.13.2
```

The pinned versions are those used for the analyses in the manuscript, and `sceps.yml` reproduces that environment exactly. The minimums are the floors declared in `pyproject.toml`; scEPS has also been verified to reproduce identical output on numpy 2.x, pandas 2.x, anndata 0.12 and scanpy 1.11.

The optional preprocessing helper script `misc/preprocess_scdata.py` additionally requires [harmonypy](https://github.com/slowkow/harmonypy) for batch integration. This is also available as an extra:
```shell
pip install "sceps[preprocess]"
```

Once the required packages to run scEPS are installed, the user may then install scEPS using:
```shell
git clone git@github.com:Genentech/sceps.git
cd sceps
pip install --no-deps .
```

# Usage

Installing scEPS provides the following command-line tools. The first four correspond to the four steps of the scEPS workflow:

| Command | Purpose |
| --- | --- |
| `sceps` | Estimate scEPS statistics for individual cell neighborhoods |
| `sceps-cluster-neighborhood` | Define approximately independent cell neighborhood blocks |
| `sceps-aggregate` | Aggregate scEPS statistics across cell types and across all cells |
| `sceps-corr` | Correlate scEPS statistics with gene expression |
| `sceps-generate-test-data` | Simulate a small test data set for trying out the workflow |

Pass `--help` to any of them for the full list of options, e.g. `sceps --help`. A detailed description of each step is available in the [Wiki page](https://github.com/Genentech/sceps/wiki).

scEPS can also be driven from Python rather than the command line:
```python
from sceps.sceps_core import *
```
See [misc/run_sceps_from_python.py](https://github.com/Genentech/sceps/blob/master/misc/run_sceps_from_python.py) for a worked example.

# Testing scEPS

The `sceps-generate-test-data` command simulates a small data set so that the whole workflow can be exercised from a plain `pip install`, without cloning this repository:

```shell
sceps-generate-test-data
```

This writes `test_scdata.h5ad` (600 cells of 10 cell types across 20,000 genes for 30 donors) and `test_magma.txt` into `./input`. The four workflow steps can then be run against it:

```shell
sceps --adata ./input/test_scdata.h5ad --donor-id-col Donor \
    --gene-list ./input/test_magma.txt --auto-gene-selection \
    --pheno Pheno --scale-pheno --scale-pheno-neighborhood --out ./output/step1
sceps-cluster-neighborhood --adata ./input/test_scdata.h5ad --donor-id-col Donor \
    --neighbors-use-rep X_pca --out ./output/step2
sceps-aggregate --prefix "./output/step1.*.txt.gz" --adata ./input/test_scdata.h5ad \
    --neighborhood-clusters ./output/step2.txt.gz --cell-type-col CellType --out ./output/step3
sceps-corr --adata ./input/test_scdata.h5ad --sceps-result ./output/step3.sceps.omega.txt.gz \
    --min-num-nonzero 3 --out ./output/step4
```

Create the `./output` directory first. The first step analyzes every cell neighborhood in turn and takes a few minutes on the simulated data; add `--start-idx 0 --stop-idx 40` to run a subset instead.

Because the simulated expression values and phenotypes are drawn independently, no cell neighborhood is expected to show a genuine disease association. The purpose of the test data is to confirm that the workflow runs end to end and to illustrate the format of each output file.

The equivalent shell scripts, together with reference output files to compare against, are also provided [here](https://github.com/Genentech/sceps/tree/master/test).

# Contact

Please create a GitHub issue if you experience any issue with running scEPS.
