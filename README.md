# BrugiaChemo-ms
This repository contains all the Bash, Python, and R scripts necessary to reproduce the data published in Wheeler et al. 2019.

## Description of scripts in `bash/`:

All Bash scripts use a convention

**CNG_analysis.sh**: Bash script of all command line programs and parameters used for cyclic-nucleotide gated channel gene identification and annotation.

**ChemoR_analysis.sh**: Bash script of all command line programs and parameters used for chemoreceptor (GPCR) gene identification and annotation.

**TRP_analysis.sh**: Bash script of all command line programs and parameters used for transient receptor potential channel gene identification and annotation.

**WBP_download.sh**: Bash script used to populate a local directory with all genome files from the WormBase ParaSite FTP.

## Description of scripts in `R/`:



## Description of files in `auxillary/`:

**CNG, ChemoR, and TRP**: includes all seed sequences and alignment used for gene identification and annotation.

**pfam_HMMs**: includes hidden Markov models curated by Pfam.

**species_info.csv**: a spreadsheet of all the genomes analyzed in the comparative pipeline, their BioProject accessions, and their nematode clade designation.

**species.txt**: a list of species analyzed, used in all Bash scripts for iterating through species genome directories
