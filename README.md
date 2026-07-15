# scEPS
**scEPS** (single-cell Expression exPlainability Statistics)

This repo contains the code of the method, **scEPS**, for integrating GWAS and single-cell disease cell atlas data to identify disease-associated cell neighborhoods. scEPS calculates a $d$ statistic at each cell neighborhood, representing the difference between the variance in disease explained by variations in the expression of each GWAS vs. each mean-expression matched control gene. An illustration of the scEPS method is shown below:

![scEPS illustration](https://github.com/Genentech/sceps/blob/master/img/scEPS_overview.png  "Overview of the scEPS method")

# Reference

The current draft of the manuscript is available [here](https://www.medrxiv.org/content/10.64898/2026.06.26.26356714v1). The code we used to create the figures in the manuscript is available [here](https://github.com/Genentech/sceps_manuscript).

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

The user may also manually install the required packages to run scEPS using the following commands:
```shell
conda install pandas=1.5.3
conda install numpy=1.26.2
conda install scipy=1.13.1
conda install scanpy=1.10.3
conda install anndata=0.10.7
conda install scikit-learn=1.3.2
```

Once the required packages to run scEPS are installed, the user may then install scEPS using:
```shell
git clone git@github.com:Genentech/sceps.git
```

# Testing scEPS

We provide examples script to test the scEPS workflow [here](https://github.com/Genentech/sceps/tree/master/test).

# Contact

Please create a GitHub issue if you experience any issue with running scEPS.
