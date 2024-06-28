require("Hmisc")
require("dplyr")
require("R.utils")
require("Biostrings")
require("data.table")
require("ggplot2")

foal_nonMT_IS <- fread(file = "../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif3n.csv") %>% as.data.frame() #%>% 

metadata <- fread(file = "../../output/02_tables/02_data_merged/summary_basic_stats.csv") %>% 
    dplyr::rename(sample_id = patientID) %>%
    select("sample_id", "FID", "nSIRS_class")
foal_nonMT_IS_meta <- merge(x = metadata, y = foal_nonMT_IS, by = "sample_id")

options(repr.plot.width=7, repr.plot.height=3)
foal_nonMT_1mer_endmotif_heatmap <- foal_nonMT_IS_meta %>% 
    filter(!grepl(sample_id, pattern = "^[P/N]")) %>% 
    mutate(read = substring(side, 0, 2)) %>% 
    mutate(EndMotif = substring(EndMotif, 0, 1)) %>% 
    mutate(nSIRS_class = factor(nSIRS_class, levels = c('H', 'nS-', 'sS-', 'S+'))) %>% 
    group_by(nSIRS_class, FID, read, EndMotif) %>%
    summarise(sCount = sum(Count, na.rm = TRUE), .groups = 'drop') %>%     
    group_by(nSIRS_class, FID, read) %>% 
    mutate(nsCount = log10(sCount / sum(sCount) * 4)) %>% 
    filter(!is.na(nsCount)) %>% 
    ungroup() %>% 
    ggplot(aes(x = FID, y = EndMotif, fill = nsCount)) + 
        geom_tile() + 
        facet_grid(rows = vars(read), cols = vars(nSIRS_class), scales = "free", space = "free") + 
        theme_bw() + 
        scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) + 
        ylab("1-mer end motif") + xlab("foals") + labs(fill = "log10,\no/e-ratio") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background = element_blank()) 
ggsave("../../output_figures/foal_nonMT_1mer_endmotif_heatmap.png", plot = foal_nonMT_1mer_endmotif_heatmap, 
       width = 7, height = 3, units = "in")
ggsave("../../output_figures/foal_nonMT_1mer_endmotif_heatmap.pdf", plot = foal_nonMT_1mer_endmotif_heatmap, 
       width = 7, height = 3, units = "in")

options(repr.plot.width=7, repr.plot.height=7)
foal_nonMT_2mer_endmotif_heatmap <- foal_nonMT_IS_meta %>% 
    filter(!grepl(sample_id, pattern = "^[P/N]")) %>% 
    mutate(read = substring(side, 0, 2)) %>% 
    mutate(EndMotif = substring(EndMotif, 0, 2)) %>% 
    mutate(nSIRS_class = factor(nSIRS_class, levels = c('H', 'nS-', 'sS-', 'S+'))) %>% 
    group_by(nSIRS_class, FID, read, EndMotif) %>%
    summarise(sCount = sum(Count, na.rm = TRUE), .groups = 'drop') %>%     
    group_by(nSIRS_class, FID, read) %>% 
    mutate(nsCount = log10(sCount / sum(sCount) * 16)) %>% 
    filter(!is.na(nsCount)) %>% 
    ungroup() %>% 
    ggplot(aes(x = FID, y = EndMotif, fill = nsCount)) + 
        geom_tile() + 
        facet_grid(rows = vars(read), cols = vars(nSIRS_class), scales = "free", space = "free") + 
        theme_bw() + 
        scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0)+ 
        ylab("2-mer end motif") + xlab("foals") + labs(fill = "log10,\no/e-ratio") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background = element_blank()) 
ggsave("../../output_figures/foal_nonMT_2mer_endmotif_heatmap.png", plot = foal_nonMT_2mer_endmotif_heatmap, 
       width = 7, height = 3, units = "in")
ggsave("../../output_figures/foal_nonMT_2mer_endmotif_heatmap.pdf", plot = foal_nonMT_2mer_endmotif_heatmap, 
       width = 7, height = 3, units = "in")
