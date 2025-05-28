# microbial
This folder contains processing scripts for intermediate processing regarding microbial counts between pipeline output and figures. The content of two subfolders are described. The scripts should be executed in the order listed.

## decontam
Decontamination is an important step for microbial cfDNA sequencing in low biomass samples. We adopted method developed and implemented by Davis et al. (https://github.com/benjjneb/decontam). 

In decontam folder, we provide a R script to run_decontam on the pipeline output. 
1. `run_decontam.R`

## preprocessing
Preprocessing folder contains scripts performining tasks to prepare useful data format for plotting. The usage of each script is provided below.

2. `run_combine_kreport.sh`: git clone https://github.com/jenniferlu717/KrakenTools. Modify parameter KRAKEN_TOOL_PATH in the script to where the folder is located. Relative paths are provided. Execute script by `bash run_combine_kreport.sh`.
3. `get_MT_load.sh`: `bash ./get_MT_load.sh <input_folder> <output_file>`.

```
input_folder: <PATH_FBI>/output/01_pipeline/foal_cohort_EquCab3/results/host_mapping 
output_folder: <PATH_FBI>/output/02_tables/00_pipeline_postprocess/MT_DNA_load/mito_minMQ40.txt
```
4. `get_decontam_input.py`: relative paths has been provided within the script. Execute by: `python get_decontam_input.py`

5. `merge_table.py`:
```
python merge_table.py --clinical <path/to/clinical.csv> --kreport <path/to/kreport.txt> --wetlab-isolation <isolation.csv> --wetlab-yield <yield.csv> --stats <output/.../stats> --MT <FBI_output/output/02_tables/00_pipeline_postprocess/MT_DNA_load/> --output </path/to/output>
```

7. `pathogen_table_new.py`: relative paths has been provided within the script. Execute by: `python pathogen_table_new.py`

