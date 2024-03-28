

```bash
snakemake \
--configfile ./config/config_foals_sepsis_EquCabAll.yaml  \
--snakefile workflow/Snakefile_FOALS_NEW  \
--profile ./profile/slurm --conda-frontend conda --use-conda --dry-run
```
