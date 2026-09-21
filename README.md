# scEPS
**scEPS** (single-cell Expression exPlainability Statistics)

This repo contains the code of the method, **scEPS**, for integrating GWAS and single-cell disease cell atlas data to identify disease-associated cell neighborhoods. scEPS calculates a $d$ statistic at each cell neighborhood, representing the difference between the variance in disease explained by variations in the expression of each GWAS vs. each mean-expression matched control gene. An illustration of the scEPS method is shown below:

![scEPS illustration](https://github.com/Genentech/sceps/blob/master/img/scEPS_overview.png  "Overview of the scEPS method")

# Reference

The current draft of the manuscript is available [here](https://www.medrxiv.org/content/10.64898/2026.06.26.26356714v1). The code we used to create the figures in the manuscript is available [here](https://github.com/Genentech/sceps_manuscript). We also implemented [CNA*](https://github.com/Genentech/cna_star), a simple extension of [CNA](https://github.com/immunogenomics/cna), for estimating the variance in disease attributable to variations in cell abundance at each cell neighborhood.

We also provide a web UI for visualizing the results in the scEPS manuscript [here](https://scepsresultexplorer.streamlit.app/).

# Manual
We provide a detailed manual of scEPS in the [Wiki page](https://github.com/Genentech/sceps/wiki).

# Installation

## Option 1: using Anaconda or Miniforge
The easiest way to install scEPS is by creating a dedicated environment through [Anaconda](https://www.anaconda.com/download) or [Miniforge](https://github.com/conda-forge/miniforge). To do this, please first install Anaconda or Miniforge on your machine. You may then install scEPS using the following commands:
``` shell
git clone git@github.com:Genentech/sceps.git
cd sceps
conda env create -f sceps.yml
conda activate sceps
```

## Option 2: manually install required packages

The user may also manually install the required packages to run scEPS. scEPS requires Python 3.9 and the following packages:

| Package | Version tested |
| --- | --- |
| [numpy](https://numpy.org/) | 1.26.2 |
| [pandas](https://pandas.pydata.org/) | 1.5.3 |
| [scipy](https://scipy.org/) | 1.13.1 |
| [anndata](https://anndata.readthedocs.io/) | 0.10.7 |
| [scanpy](https://scanpy.readthedocs.io/) | 1.10.3 |
| [scikit-learn](https://scikit-learn.org/) | 1.3.2 |
| [statsmodels](https://www.statsmodels.org/) | 0.14.5 |
| [tqdm](https://tqdm.github.io/) | 4.67.1 |
| [packaging](https://packaging.pypa.io/) | 25.0 |
| [matplotlib](https://matplotlib.org/) | 3.9.4 |
| [seaborn](https://seaborn.pydata.org/) | 0.13.2 |

These can be installed with a single command:
```shell
conda install -c conda-forge python=3.9 numpy=1.26.2 pandas=1.5.3 scipy=1.13.1 \
    anndata=0.10.7 scanpy=1.10.3 scikit-learn=1.3.2 statsmodels=0.14.5 \
    tqdm=4.67.1 packaging=25.0 matplotlib-base=3.9.4 seaborn=0.13.2
```

The versions above are the ones scEPS has been tested with; other recent versions are likely to work as well.

The optional preprocessing helper script `misc/preprocess_scdata.py` additionally requires [harmonypy](https://github.com/slowkow/harmonypy) for batch integration:
```shell
conda install -c conda-forge harmonypy
```

Once the required packages to run scEPS are installed, the user may then install scEPS using:
```shell
git clone git@github.com:Genentech/sceps.git
```

# Testing scEPS

We provide examples script to test the scEPS workflow [here](https://github.com/Genentech/sceps/tree/master/test).

# Contact

Please create a GitHub issue if you experience any issue with running scEPS.
