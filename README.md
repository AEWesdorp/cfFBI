# *cf*FBI


We developed a **F**oal *c*ell-*f*ree DNA sequencing for **B**acterial **I**dentification (*cf***FBI**) workflow, integrating wetlab and computational protocols to detect increased bacterial cfDNA abundance in blood. Specifically, our workflow first subtracts host-derived reads by mapping to the foal reference genome(s), then applies stringent bacterial read classification, and finally performs extensive *in silico* decontamination to remove residual contaminants. This repository is associated with: "Bacterial cell-free DNA profiling reveals co-elevation of multiple bacteria in newborn foals with suspected sepsis."

For details on the Snakemake pipeline, used for the processing of paired-end Illumina sequencing data see: [cfFBI-pipeline](https://github.com/AEWesdorp/cfFBI/tree/main/pipeline).

## Table of Contents
1. [Introduction](#introduction)
2. [Repository Content](#repository-content)
4. [Usage](#usage)
5. [License](#license)
6. [Contact](#contact)

## Introduction

cfFBI leverages paired-end Illumina sequencing of plasma cell-free DNA to identify bacteria, providing an end-to-end pipeline from raw data to final analysis. This repository was created by Li-Ting Chen & Emmy Wesdorp from [De Ridder lab](https://www.deridderlab.nl/) at the Center of Molecular Medicine, University Medical Center Utrecht, the Netherlands.

You can find analysis pipelines and scripts for generating figures related to the paper: "Bacterial cell-free DNA profiling reveals co-elevation of multiple bacteria in newborn foals with suspected sepsis."

## Repository Content
The various parts of the analyses are organized into different folders within the main directory. Each folder contains a `README.md`  file with specific details relevant to the analyses conducted/information provided within that folder.

#### cfFBI-pipeline
Details of the pipeline, which processes paired-end Illumina sequencing data to identify pathogenic species, with optimization for detecting bacterial genera and species in foals: [cfFBI](https://github.com/AEWesdorp/cfFBI/tree/main/pipeline).


#### microbial
... see: [microbial](https://github.com/AEWesdorp/cfFBI/tree/main/microbial)

#### fragmentomics
For details on the fragmentomic characteristics of the horse-derived cfDNA molecules analyzed in this study, see: [fragmentomics](https://github.com/AEWesdorp/cfFBI/tree/main/fragmentomics).

#### figures
Details on data processing and figure generation: [figures](https://github.com/AEWesdorp/cfFBI/tree/main/figures).

## Usage
A `README.md` file could be found in each folder concerning relevant analyses.

## License
This project is licensed under the GNU GENERAL PUBLIC LICENSE. See the LICENSE file for more details.

## Contact
Please contact the authors and create an issue on github to get help.
