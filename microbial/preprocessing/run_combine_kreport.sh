# First clone KrakenTools from remote repository
# git clone https://github.com/jenniferlu717/KrakenTools.git
CONFIDENCE_SCORE='0.8'
# Combine files from all samples with CONFIDENCE_SCORE Kraken2 threshold
python3 /PATH/TO/KrakenTools/combine_kreports.py \
-r ../../output/01_pipeline/foals_sepsis_EquCabAll/results/kraken2_report/after_host_mapping/*_EquAllRS_conf${CONFIDENCE_SCORE}.report \
-o ../../output/02_tables/00_raw/combined_foals_sepsis_EquAllRS_${CONFIDENCE_SCORE}.txt --display-headers

