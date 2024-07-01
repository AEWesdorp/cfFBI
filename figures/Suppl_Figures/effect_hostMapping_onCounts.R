require("ggpubr")
require("tidyverse")
require("Hmisc")
require("dplyr")
require("purrr")
require("readr")
require("stringr")
require("reshape2")
require("rstatix")

counts_list <- list()
for (d in c("EquCabAll","EquCab3")) {
    if (exists("k2rep_comb")) {rm("k2rep_comb")}
    
    # Loop through files in each directory
    for (f in list.files(paste0("../../output/01_pipeline/foal_cohort_", d, "/results/stats/"),
                         pattern = "R1_05_")) {
        
        sample_nm <- unlist(strsplit(f, '_'))[1]
        QC_c_tmp <- read.csv(file = 
            paste0("../../output/01_pipeline/foal_cohort_", d, "/results/stats/", f), 
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
            as.numeric()
        
        host_c_tmp <- read.csv(file = 
            paste0("../../output/01_pipeline/foal_cohort_", d, "/results/stats/", 
                   sample_nm, "_R1_06_host_mapp_fastq.txt"), 
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
            as.numeric()
        
        k2rep <- read.csv(file = 
            paste0("../../output/01_pipeline/foal_cohort_", d, "/results/kraken2_report/after_host_mapping/", 
                   sample_nm, "_EquAllRS_conf0.8.report"), 
                          header=FALSE, sep = "\t", stringsAsFactors=FALSE)
        k2rep[, ncol(k2rep)] <- str_trim(k2rep[, ncol(k2rep)], side = "left")
        bact_c_tmp <- k2rep[which(k2rep[, ncol(k2rep)] == "Bacteria"),"V2"] 
        
        temp_df <- data.frame(sample_id = sample_nm, 
                              QC_count = QC_c_tmp, 
                              Host_map_count = host_c_tmp, 
                              Bact_k2_count = bact_c_tmp)
        
        if (!exists("k2rep_comb")) {k2rep_comb <- temp_df} else {
            k2rep_comb <- bind_rows(k2rep_comb, temp_df)} 

        counts_list[[d]] <- k2rep_comb
    }
}

options(repr.plot.width=6, repr.plot.height=6)

# Process the EquCabAll data
rel_EquCabAll <- counts_list[["EquCabAll"]] %>%
    mutate(
        rel_EquCabAll_Host_map_count = Host_map_count / QC_count * 100,
        rel_EquCabAll_Bact_k2_count = Bact_k2_count / QC_count * 100
    ) %>% filter(!grepl(sample_id, pattern = "[N/P]")) %>%
    select(sample_id, contains("EquCabAll"))

# Process the EquCab3 data
rel_EquCab3 <- counts_list[["EquCab3"]] %>%
    mutate(
        rel_EquCab3_Host_map_count = Host_map_count / QC_count * 100,
        rel_EquCab3_Bact_k2_count = Bact_k2_count / QC_count * 100
    ) %>% filter(!grepl(sample_id, pattern = "[N/P]")) %>%
    select(sample_id, contains("EquCab3"))

# Merge the two datasets on sample_id
merged_data <- merge(rel_EquCabAll, rel_EquCab3, by = "sample_id")

# Add the relative comparison columns -- Host Count
merged_data_sel <- merged_data %>%
    mutate(
        rel2_Host_EquCabAll = rel_EquCabAll_Host_map_count / rel_EquCab3_Host_map_count,
        rel2_Host_EquCab3 = rel_EquCab3_Host_map_count / rel_EquCab3_Host_map_count
    ) %>%
    select(sample_id, contains("rel2")) 

# Reshape the data for plotting -- Host Count
merged_data_sel <- reshape2::melt(merged_data_sel, id.vars = "sample_id") %>% 
    mutate(
        Host_mapping = gsub("rel2_Host_", "", variable),
        Host_mapping = ifelse(variable == "rel2_Host_EquCab3", "EquCab3.0", Host_mapping),
        Host_mapping = ifelse(variable == "rel2_Host_EquCabAll", "EquCab NCBI all 10", Host_mapping),
        Host_mapping = factor(Host_mapping, levels = c("EquCab3.0", "EquCab NCBI all 10"))
    )

# Plot the data -- Host Count
stat_test_Host <- merged_data_sel %>% 
        t_test(value ~ Host_mapping, paired = TRUE, 
            var.equal = "sample_id", alternative = "greater") 

effect_hostMapping_hostCount <- ggplot(merged_data_sel, aes(x = Host_mapping, y = value)) +
    geom_violin(width = 0.4, size = 1, fill = NA, col = "grey", alpha = 0.5) + 
    geom_dotplot(width = 0.6, binaxis = "y", stackdir = "center", dotsize = 0.42, fill = "grey") + 
    theme_bw() + theme() + ggtitle("Host Mapped (using Bowtie2)") + 
    labs(x = "Genome assembly\nused for host mapping", y = "% of QC reads, normalized", fill = "") + 
    stat_pvalue_manual(size = 5, stat_test_Host, label = "p", y.position = 1.1) 

# Save as pdf & png -- Host Count
ggsave("../../output_figures/effect_hostMapping_hostCount.png", plot = effect_hostMapping_hostCount, 
       width = 6, height = 6, units = "in")
ggsave("../../output_figures/effect_hostMapping_hostCount.pdf", plot = effect_hostMapping_hostCount, 
       width = 6, height = 6, units = "in")
print("Done, part 1 of 2! Figures can be found here: ../../output_figures/effect_hostMapping_hostCount.*")

# Add the relative comparison columns -- Bacterial Count
merged_data_sel <- merged_data %>%
    mutate(
        rel2_Bact_EquCabAll = rel_EquCabAll_Bact_k2_count / rel_EquCab3_Bact_k2_count,
        rel2_Bact_EquCab3 = rel_EquCab3_Bact_k2_count / rel_EquCab3_Bact_k2_count
    ) %>%
    select(sample_id, contains("rel2")) 

# Reshape the data for plotting -- Bacterial Count
merged_data_sel <- reshape2::melt(merged_data_sel, id.vars = "sample_id") %>% 
    mutate(
        Bact_k2 = gsub("rel2_Bact_", "", variable),
        Bact_k2 = ifelse(variable == "rel2_Bact_EquCab3", "EquCab3.0", Bact_k2),
        Bact_k2 = ifelse(variable == "rel2_Bact_EquCabAll", "EquCab NCBI all 10", Bact_k2),
        Bact_k2 = factor(Bact_k2, levels = c("EquCab3.0", "EquCab NCBI all 10"))
    )

# Plot the data -- Bacterial Count
stat_test_Bact <- merged_data_sel %>% 
        t_test(value ~ Bact_k2, paired = TRUE, 
            var.equal = "sample_id", alternative = "greater") 

effect_hostMapping_bactCount <- ggplot(merged_data_sel, aes(x = Bact_k2, y = value)) +
    geom_violin(width = 0.4, size = 1, fill = NA, col = "grey", alpha = 0.5) + 
    geom_dotplot(width = 0.6, binaxis = "y", stackdir = "center", dotsize = 0.42, fill = "grey") + 
    theme_bw() + theme() + ggtitle("Bacterial Classified (using Kraken2)") + 
    labs(x = "Genome assembly\nused for host mapping", y = "% of QC reads, normalized", fill = "") + 
    stat_pvalue_manual(size = 5, stat_test_Bact, label = "p", y.position = 1.001) 

# Save as pdf & png -- Bacterial Count
ggsave("../../output_figures/effect_hostMapping_bactCount.png", plot = effect_hostMapping_bactCount, 
       width = 6, height = 6, units = "in")
ggsave("../../output_figures/effect_hostMapping_bactCount.pdf", plot = effect_hostMapping_bactCount, 
       width = 6, height = 6, units = "in")
print("Done, part 2 of 2! Figures can be found here: ../../output_figures/effect_hostMapping_bactCount.*")
