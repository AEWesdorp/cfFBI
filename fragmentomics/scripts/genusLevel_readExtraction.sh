## extract reads classified to a Genus (incl. children)
in_dir="../../output/01_pipeline/foals_sepsis_EquCabAll/results/"
out_dir="../../output/04_fragmentomics/genus_level/"

mkdir -p $out_dir

for sample_id in $(cut -f 1 ../../pipeline/config/unit_foal_sepsis_EquCabAll.txt | grep -v "sample_name" ); do 
	echo "sample_id: $sample_id"
	for taxId in 613 590 286 583 745 53335 570 561 547 642 713 469 1301 1279 1350 1386 1578 1637 2735 460517 4930 5206 5820 36862 2836216 1653176 1548510 2742598; do 
		echo $taxId

		../../pipeline/resources/extract_kraken2_taxids_edit.py -r ${in_dir}kraken2_report/after_host_mapping/${sample_id}_EquAllRS_conf0.8.report -t ${taxId} --include-children --output_file ${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_taxIDs.txt
		awk 'NR == FNR { keywords[$1]=1; next; } { if ($3 in keywords) print $2; }' ${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_taxIDs.txt ${in_dir}kraken2_output/after_host_mapping/${sample_id}_EquAllRS_conf0.8.output > ${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_nms.txt
		nr_lines=$( wc -l  ${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_nms.txt | cut -d" " -f 1 )
		echo $nr_lines

		if [ $nr_lines != "0" ]; then
			echo "mapping reads classified to this taxa"
			filterbyname.sh in1=${in_dir}host_mapping/${sample_id}_unmapped_host_r1.fq in2=${in_dir}host_mapping/${sample_id}_unmapped_host_r2.fq out1=${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_R1.fastq out2=${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_R2.fastq names=${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_nms.txt include=TRUE overwrite=TRUE
		fi
	rm ${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_taxIDs.txt
	rm ${out_dir}${sample_id}_EquAllRS_conf0.8_${taxId}_nms.txt

	done
done
