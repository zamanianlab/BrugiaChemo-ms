# BrugiaChemo-ms
This repository contains all the Bash, Python, and R scripts necessary to reproduce the figures published in Wheeler et al. 2020.

## Description of scripts in `bash/`:

All Bash scripts use a convention for organizing data in local directories outside of GitHub version control. In our local `.bash_profile`, we add the following two lines to define local variables:

```bash
export GIT_PATH="$HOME/GitHub"
export GIT_DATA="$HOME/Box/GHdata"
```

In the header of each Bash script, we define `proj` as the repo name and reference to the GitHub repo as well as a local directory into which the GitHub repo will write out all resulting files, for example:

```bash
proj="BrugiaChemo-ms"

gh_dir="${GIT_PATH}"/"${proj}"
local_dir="${GIT_DATA}"/"${proj}"
```

For this repo to be cloned/forked and reproduced properly, the local variables `GIT_DATA` and `GIT_PATH` will need to be defined according to the user's organization.

**NOTE:** All bash scripts use GNU Coreutils. If using a Mac, you will need to install these via Homebrew. If a command isn't working try preceding the command with a 'g' (e.g. `ggrep` instead of `grep`).

**CNG_analysis.sh**: Bash script of all command line programs and parameters used for cyclic-nucleotide gated channel gene identification and annotation.

**ChemoR_analysis.sh**: Bash script of all command line programs and parameters used for chemoreceptor (GPCR) gene identification and annotation.

**TRP_analysis.sh**: Bash script of all command line programs and parameters used for transient receptor potential channel gene identification and annotation.

**WBP_download.sh**: Bash script used to populate a local directory with all genome files from the WormBase ParaSite FTP.

## Description of scripts in `R/`:

This directory contains a subdirectory for each gene family (chemoreceptors, CNGs, and TRPs), which contains R scripts for plotting trees and other figures where appropriate, and an `experiments` subdirectory. Everything in `experiments` relates to benchwork, while everything in the gene family subdirectories relates to computational work. Data used to create figures are held in `data` and output plots are written to `plots`.

## Description of files in `auxillary/`:

**CNG, ChemoR, and TRP**: includes all seed sequences and alignments used for gene identification and annotation.

**pfam_HMMs**: includes hidden Markov models curated by Pfam.

**species_info.csv**: a spreadsheet of all the genomes analyzed in the comparative pipeline, their BioProject accessions, and their nematode clade designation.

**species.txt**: a list of species analyzed, used in all Bash scripts for iterating through species genome directories
