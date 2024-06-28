require("Hmisc")
require("dplyr")
require("R.utils")
require("Biostrings")
require("data.table")

DIR_IN_FOALS="../../output/01_pipeline/foal_cohort_EquCabAll/results/host_TLEN_readCount/"
list.files(DIR_IN_FOALS, pattern = "[pos/neg]_nonMT_TLEN_EndMotif.txt") %>% length

if (exists("foal_nonMT_IS_meta")){rm("foal_nonMT_IS_meta")}
for (file in list.files(DIR_IN_FOALS, pattern = "[pos/neg]_nonMT_TLEN_EndMotif.txt")){
    sample_tmp <- unlist(strsplit(file, "_"))[1] 
    side_tmp <- unlist(strsplit(file, "_"))[2]
    
    if(countLines(paste0(DIR_IN_FOALS, file)) > 1){
        pre_tmp_IS <- fread(paste0(DIR_IN_FOALS, file), sep = ",") %>% 
                    dplyr::rename(EndMotif_tmp = EndMotif) %>% 
                    mutate(sample_id = sample_tmp) %>% 
                    mutate(side = side_tmp) %>%
                    filter(TLEN != 0) %>% #remove TLEN not equal to 0
                    filter(abs(TLEN) < 1000) %>% #filter max TLEN of 1000
                    filter(!grepl(EndMotif_tmp, pattern = "N")) %>% #remove reads with end-motif including "N"
                    group_by(TLEN, EndMotif_tmp, side, sample_id) %>% 
                    summarise(Count = sum(Count), .groups = "keep")
        
        if (grepl(side_tmp, pattern = "neg")){
            tmp_IS <- pre_tmp_IS
            tmp_IS$EndMotif <- substring(as.character(reverseComplement(DNAStringSet(pre_tmp_IS$EndMotif_tmp))), 0,1)
        } else {
            tmp_IS <- pre_tmp_IS
            tmp_IS$EndMotif <- substring(pre_tmp_IS$EndMotif_tmp, 0,1)
        }
        tmp_IS_meta <- tmp_IS %>% 
            group_by(side, sample_id, EndMotif) %>% 
            summarise(Count = sum(Count), .groups = "keep")
        
        if (!exists("foal_nonMT_IS_meta")){foal_nonMT_IS_meta <- tmp_IS_meta} else {
            foal_nonMT_IS_meta <- rbind(foal_nonMT_IS_meta, tmp_IS_meta)}
    }
    fwrite(foal_nonMT_IS_meta, "../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif.csv", row.names = FALSE)
}
print("DONE Part1 of 2! You can find the output here: ../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif.csv") 

if (exists("foal_nonMT_IS_meta")){rm("foal_nonMT_IS_meta")}
for (file in list.files(DIR_IN_FOALS, pattern = "[pos/neg]_nonMT_TLEN_EndMotif.txt")){
    sample_tmp <- unlist(strsplit(file, "_"))[1] 
    side_tmp <- unlist(strsplit(file, "_"))[2]
    
    if(countLines(paste0(DIR_IN_FOALS, file)) > 1){
        pre_tmp_IS <- fread(paste0(DIR_IN_FOALS, file), sep = ",") %>% 
                    dplyr::rename(EndMotif_tmp = EndMotif) %>% 
                    mutate(sample_id = sample_tmp) %>% 
                    mutate(side = side_tmp) %>%
                    filter(TLEN != 0) %>% #remove TLEN not equal to 0
                    filter(abs(TLEN) < 1000) %>% #filter max TLEN of 1000
                    filter(!grepl(EndMotif_tmp, pattern = "N")) %>% #remove reads with end-motif including "N"
                    group_by(TLEN, EndMotif_tmp, side, sample_id) %>% 
                    summarise(Count = sum(Count), .groups = "keep")
        
        if (grepl(side_tmp, pattern = "neg")){
            tmp_IS <- pre_tmp_IS
            tmp_IS$EndMotif <- substring(as.character(reverseComplement(DNAStringSet(pre_tmp_IS$EndMotif_tmp))),0,3)
        } else {
            tmp_IS <- pre_tmp_IS
            tmp_IS$EndMotif <- substring(pre_tmp_IS$EndMotif_tmp,0,3)
        }
        tmp_IS_meta <- tmp_IS 
        
        if (!exists("foal_nonMT_IS_meta")){foal_nonMT_IS_meta <- tmp_IS_meta} else {
            foal_nonMT_IS_meta <- rbind(foal_nonMT_IS_meta, tmp_IS_meta)}
    }
    fwrite(foal_nonMT_IS_meta, "../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif3n.csv", row.names = FALSE)
}
print("DONE Part2 of 2! You can find the output here: ../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif3n.csv") 