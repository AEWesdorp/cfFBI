# FBI-pipeline
This repository contains the computational pipeline for cell-**f**ree DNA sequencing for **b**acteria **i**dentification, also known as **FBI**, an open-source data analysis Snakemake workflow designed for bacterial profiling of foals' liquid biopsy plasma samples. The FBI-pipeline processes paired-end Illumina sequencing data.

#### FBI-Pipeline Overview:
1. **Removal of Synthetic Spike-in Reads:** Reads from synthetic spike-in sequences (50, 100, and 150 bp) are removed using *bbduk.sh* from *BBMap*.
2. **Duplicate Removal:** Duplicates are removed using *Nubeam*.
3. **Quality Control and Filtering:** High-quality sequencing data is generated using *fastp*, which removes low-quality reads and applies a low complexity filter. Adapter sequences and short reads (<35 bp) are removed using *AdapterRemoval*.
4. **Host Sequence Subtraction:** Host sequences are subtracted by mapping to the human reference genome using *Bowtie2*.
5. **Taxonomic Classification:** Remaining paired-end reads are taxonomically classified using *Kraken2*.
6. **Taxonomic Abundance Estimation:** For a select set of positive control samples, taxonomic abundance is estimated using *Bracken*.

## Getting started
To use the FBI-pipeline, follow these steps: 

### Prerequisites
1. Ensure you have either `conda` installed on your system.

### Installation
1. Clone the GitHub Repository:
    ```bash
    git clone https://github.com/AEWesdorp/FBI.git
    ```
2. Install Snakemake:
   follow the installation instructions on the [Snakemake website](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
3. Create and Activate a Conda Environment:
    ```bash
    # Create a new empty environment called "FBI_env"
    conda create -c conda-forge -c bioconda -n FBI_env snakemake
    # Activate the environment "FBI_env"
    conda activate FBI_env
    ```

### Create a samplesheet and configfile 
Create a samplesheet and a configfile in folder `configs/` accordingly. 
An example samplesheet can be found at **`config/units_foal_cohort.txt`** and an example configfile at **`config/config_foal_cohort_EquCabAll.yaml`**. 

In **`config/units_foal_cohort.txt`**, ensure to specify the following: 
- *sample_name*: Specify name of your sample.
- *library_prep*: Specify the library preparation you have used (options: `SRSLY` or `KAPA`)
- *adapter_type*: Specify which adapter were used during Illumina library contstruction (options: `SRSLY_dual_index`, `KAPA_single_index` or `IDT384UMI_dual`)
- *UDI*: Set the UDI of the sample.
- *path_to_R1_R2*: Provide the directory path where raw sequencing files are stored (path to `*R1*.fastq.gz`, `*R2*.fastq.gz`) 
  
In **`config/config_foal_cohort_EquCabAll.yaml`**, ensure to specify the following:

General Settings:
- *units*: Specify name of the samplesheet (for example, `config/units_foal_cohort.txt`). 
- *run_name*: Set a unique name for the run.
- *outdir*: Specify the output directory where results will be stored.

Reference Genome Settings:
- *reference_genome*: Indicate the reference genome to be used.
- *reference_genome_dir*: Provide the directory path where the reference genome is stored.

Kraken2 Classification Settings:
- *database*: Specify the database used for Kraken2 classification.
- *database_dir*: Define the directory path where the Kraken2 database is located.
- *k2_threshold*: Set the threshold value for Kraken2 classification.

### Running the FBI-pipeline on an interactive node
1. Start a screen session. 
2. Request an interactive node for for running the jobs (long enough to finish all jobs of one liquid biopsy sample, e.g. 24 hours), with 450G mem, 16 cores. 
3. Move to the `pipeline/` sub-directory within the cloned Git directory where your workflow resides.
4. Activate your conda environment.
      ```bash
    # Activate the environment "FBI_env"
    conda activate FBI_env
    ```
5. Run the snakemake pipeline.
   ```bash
   snakemake --configfile ./config/config_foal_cohort_EquCabAll.yaml \
   --snakefile workflow/Snakefile_FOALS  --cores all --conda-frontend conda --use-conda
   ```

### Running the FBI-pipeline by submitting jobs via [slurm](https://slurm.schedmd.com/documentation.html) scheduler:
1. Start a screen session. 
2. Request an interactive node for submitting jobs (long enough for all jobs to finish, e.g. 48 hours), with 16G mem, 2 cores.
3. Move to the `pipeline/` sub-directory within the cloned Git directory where your workflow resides.
4. Activate your conda environment.
    ```bash
    # Activate the environment "FBI_env"
    conda activate FBI_env
    ```
5. Run the snakemake pipeline.
   ```bash
   snakemake --configfile ./config/config_foal_cohort_EquCabAll.yaml \
   --snakefile workflow/Snakefile_FOALS  --profile ./profile/slurm --conda-frontend conda --use-conda
   ```
   (Resources defined in each rule will be used. If not defined, default resources defined in the `profile/slurm/config.yaml` will be used.)

More information See [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and page [snakemake slurm](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#executing-on-slurm-clusters). 

### Trouble shooting
- Input fasta folder always should contains following files to start with `*R1*.fastq.gz, *R2*.fastq.gz`
- The current pipeline version includes adapter sequence information utilized by the SRSLY Claret Kit, the KAPA Kit, and 384 IDT UMI's. If another library preparation method is employed, kindly update the `workflow/rules/trim.smk` file and append the adapter index sequences to the `resources/adapter_indexes/` directory.

## FBI-Pipeline Employment for Publication Purposes 
#### Comparing Library Preparations in a Foal's Liquid Biopsy Sample: SRSLY retention protocol (moderate/small) & bead-based size-selection (yes/no)

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
