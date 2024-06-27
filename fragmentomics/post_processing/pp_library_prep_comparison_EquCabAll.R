require("Hmisc")
require("dplyr")
require("R.utils")
require("Biostrings")
require("data.table")

HOST_IN_EquAll="../../output/01_pipeline/library_prep_comparison_EquCabAll/results/host_TLEN_readCount/"
list.files(HOST_IN_EquAll, pattern = "[pos/neg]_nonMT_TLEN_EndMotif.txt") %>% length()

if (exists("foal_nonMT_IS_meta")){rm("foal_nonMT_IS_meta")}
for (file in list.files(HOST_IN_EquAll, pattern = "[pos/neg]_nonMT_TLEN_EndMotif.txt")){
    sample_tmp <- unlist(strsplit(file, "_"))[1] 
    side_tmp <- unlist(strsplit(file, "_"))[2]
    
    if(countLines(paste0(HOST_IN_EquAll, file)) > 1){
        pre_tmp_IS <- fread(paste0(HOST_IN_EquAll, file), sep = ",") %>% 
                    dplyr::rename(EndMotif_tmp = EndMotif) %>% 
                    mutate(sample_id = sample_tmp) %>% 
                    mutate(side = side_tmp) %>%
                    filter(TLEN != 0) %>% #remove TLEN not equal to 0
                    filter(abs(TLEN) < 1000) %>% #filter max TLEN of 1000
                    group_by(TLEN, sample_id, side) %>% 
                    summarise(Count = sum(Count), .groups = "keep")
        
        if (!exists("foal_nonMT_IS_meta")){foal_nonMT_IS_meta <- pre_tmp_IS} else {
            foal_nonMT_IS_meta <- rbind(foal_nonMT_IS_meta, pre_tmp_IS)}
    }
    fwrite(foal_nonMT_IS_meta, "../../output/04_fragmentomics/library_prep_comparison_EquCabAll_nonMT_IS_meta.csv", row.names = FALSE)
}
