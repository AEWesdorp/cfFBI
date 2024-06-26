# Fragmentomics
In this project, we investigated the fragmentomics characteristics of horse-derived cfDNA molecules. Specifically, after mapping to the host genome, we analyzed either the cfDNA template length or the cfDNA end-motif frequency. After **processing** the BAM files with `Snakemake`, we performed further **post-processing** in `R`, resulting in the following files:
- !!!
- !!!
- !!!

## `Snakemake` for processing BAM files 
#### Comparing Library Preparations in a Foal's Liquid Biopsy Sample: SRSLY retention protocol (moderate/small) & bead-based size-selection (yes/no)
```bash
snakemake --configfile ./config/fragmentomics_library_prep_comparison_EquCabAll.yaml \
--snakefile  ./processing/Snakemake_fragmentomics_host --profile ./profile/slurm --conda-frontend conda --use-conda

snakemake --configfile ./config/fragmentomics_library_prep_comparison_EquCab3.yaml \
--snakefile  ./processing/Snakemake_fragmentomics_host --profile ./profile/slurm --conda-frontend conda --use-conda
```

#### Computational Analysis Cohort Samples: Mapping to All Equus caballus Reference Genomes in NCBI RefSeq
```bash
!!! running for Equ3, foal cohort!!!!! 
```

## Post-processing in `R`
After processing the BAM files, we performed further post-processing using `R`:
```bash
!!! 
```
The resulting .csv files where used for visualization purposes. 
