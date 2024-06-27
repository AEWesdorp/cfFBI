# Fragmentomics
In this project, we investigated the fragmentomics characteristics of horse-derived cfDNA molecules. Specifically, after mapping to the host genome, we analyzed either the cfDNA template length or the cfDNA end-motif frequency. After **processing** the BAM files with `Snakemake`, we performed further **post-processing** in `R`, resulting in the following files:
`library_prep_comparison_EquCabAll_nonMT_IS_meta.csv` & `../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif.csv`


## `Snakemake` for processing BAM files 
For instructions on how to get started with Snakemake, we refer to `../pipeline/README.md`. 
#### Comparing Library Preparations in a Foal's Liquid Biopsy Sample: SRSLY retention protocol (moderate/small) & bead-based size-selection (yes/no)
```bash
snakemake --configfile ./config/fragmentomics_library_prep_comparison_EquCab3.yaml \
--snakefile  ./processing/Snakemake_fragmentomics_host --profile ./profile/slurm --conda-frontend conda --use-conda
```

#### Computational Analysis Cohort Samples: Mapping to All Equus caballus Reference Genomes in NCBI RefSeq
```bash
snakemake --configfile ./config/fragmentomics_foal_cohort_EquCabAll.yaml \
--snakefile ./processing/Snakemake_fragmentomics_host --profile ./profile/slurm --conda-frontend conda --use-conda
```

## Post-processing in `R`
After processing the BAM files, we performed further post-processing. 

#### Prerequisites
Post-processing was conducted using R version 4.2.0 along with the following packages:
- Hmisc: 5.1.0
- dplyr: 1.1.2
- R.utils: 2.12.3
- Biostrings: 2.66.0
- data.table: 1.14.8
Ensure all these packages are installed before proceeding.

#### Post-processing Workflow
Navigate to the `post_processing/` directory:
```bash
cd post_processing/
```
#### Post-processing comparison library preparations
```bash
Rscript pp_library_prep_comparison_EquCabAll.R
```
The resulting `../output/04_fragmentomics/library_prep_comparison_EquCabAll_nonMT_IS_meta.csv` file is intended for visualization purposes.

#### Post-processingCohort Samples
```bash
Rscript pp_foal_cohort_EquCabAll_host_motif.R
```
The resulting `../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif.csv` file is intended for visualization purposes. 
