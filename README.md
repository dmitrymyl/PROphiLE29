# PROphiLE29 by OMG

![OMG](./omg.jpg)

A pipeline for prediction of &phi;29 polymerase-like activity in a set of proteins.
Use a linux-based system for installation and usage. `conda` package managment system
is required for dependency installation.

## Overview

The pipeline consists of 3 steps.

### Prefiltering

Proteins are prefiltered with `hmmer` to detect putative DNA polymerases type B
(PF00136).

### Feature extraction

Aminoacid residues are extracted from filtered proteins based on their positions
in reference alignment.

### Prediction

Prediction is made with pretrained model for CatBoost. Details of training the model
are available in the dedicated Jupyter notebook in `scripts` directory.

## Installation

Just clone this repo and install `conda` environment from file.

```bash
$ git clone https://github.com/dmitrymyl/PROphiLE29.git
$ cd PROphiLE29
$ conda env create -f environment.yml
```

## Usage

Activate conda environment `prophile29` and run python script `predict.py` with
arguments listed below.

```bash
$ conda activate prophile29
$ python path/to/repo/predict.py -h
usage: predict.py [-h] --query [<>.fasta] [--threshold [x]] --output [<>.csv]

optional arguments:
  -h, --help          show this help message and exit
  --query [<>.fasta]  query protein sequences in fasta format (default: None)
  --threshold [x]     HMM e-value threshold for prefiltering (default: 5)
  --output [<>.csv]   output csv file with predictions ranked by proba (default: None)
```

Output is a whitespace-delimited .csv file with two columns: protein IDs from input .fasta file
and probabilities of exhibiting &phi;29 polymerase-like activity.