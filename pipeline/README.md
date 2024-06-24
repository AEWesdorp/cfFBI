
#### Comparing Library Preparations in a Foal's Liquid Biopsy Sample: Moderate/Small Retention; Bead-Based Selection Yes/No

```bash
snakemake --configfile ./config/config_library_prep_comparison_EquCabAll.yaml \
--snakefile ./workflow/Snakefile_FOALS --profile ./profile/slurm --conda-frontend conda --use-conda
```

#### Computational Analysis Cohort Samples: Mapping to All Equus caballus Reference Genomes in NCBI RefSeq

```bash
snakemake --configfile ./config/config_foal_cohort_EquCabAll.yaml \
--snakefile ./workflow/Snakefile_FOALS --profile ./profile/slurm --conda-frontend conda --use-conda
```

#### Computational Analysis Cohort Samples: Mapping to the NCBI RefSeq Equus caballus Reference Genome (EquCab3.0)
```bash
snakemake --configfile ./config/config_foal_cohort_EquCab3.yaml \
--snakefile ./workflow/Snakefile_FOALS --profile ./profile/slurm --conda-frontend conda --use-conda --until h_kraken2
```

