#!/bin/bash
# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_folder> <output_file>"
    exit 1
fi

MIN_MQ=10
CHROM_LOC='NC_001640.1:1-16660'
input_folder="$1" #./output/01_pipeline/foals_sepsis_EquCabAll/results/host_mapping
output_file="$2" # ./output/02_tables/00_pipeline_postprocess/MT_DNA_load/mito_minMQ

# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder '$input_folder' not found."
    exit 1
fi

# Loop through each file with the specified suffix in the input folder
for file in "$input_folder"/*_alligned_host_pe_srt.bam; do
    if [ -f "$file" ]; then
        echo "Processing $file..."
        basename=$(basename "$file")
        echo "$basename"  >> "$output_file"
        # Perform the action with samtools stats on each file
        echo "Samtools Coverage"
        samtools coverage -r ${CHROM_LOC} --min-MQ ${MIN_MQ} "$file" >> "$output_file/${MIN_MQ}.txt"
        echo "Done"
    fi
done

echo "Finished processing all files."