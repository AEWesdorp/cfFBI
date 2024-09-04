## -----------------------------------------------------------------------------
# Code accompanying the publication: "Bacterial cell-free DNA profiling reveals co-elevation of multiple bacteria in newborn septic foals"
# Li-Ting Chen, Emmy Wesdorp et al. 

# Please direct any questions or comments to a.e.wesdorp@umcutrecht.nl

## -----------------------------------------------------------------------------
source("../Sourced/Elevated_bacteria.R")

metadata <- all_basic_stats %>%
    filter(grepl(FID, pattern = "^F")) %>%
    mutate("Blood culture" = ifelse(Y.N.Blood.culture.positive == "Y" | 
                                    FID == "F28", yes = "P", no = NA)) %>% 
    mutate("Age (days)"= DescTools::RoundTo(c(Foal.age.at.presentation..hours./24), FUN = "floor")) %>% 
    mutate(Outcome = ifelse(Lived == "Y", yes = "L", no = 
                     ifelse(Lived == "N", yes = "D", no = NA))) %>% 
    select(c("FID", "nSIRS_class", "Age (days)", "Blood culture", "Outcome")) %>%
    reshape2::melt(id.vars = c("FID", "nSIRS_class"), variable.name = "variable", value.name = "label") %>%
    mutate(variable = factor(variable, levels = c("Outcome", "Age (days)", "Blood culture"))) %>% 
    mutate(facets = "metadata") %>% 
    mutate(sum_rel_abundance = NA) %>% 
    mutate(sum_count_abundance = NA) %>% 
    mutate(tmp_col = NA) 

data <- fig3_contaminantFree_sumGenus_maxAbundance_sepsisGenera %>% 
        mutate(tmp_col = ifelse(sum_rel_abundance > max_H & 
                                sum_rel_abundance < max_nSneg &
                               sum_count_abundance >= min_hard_count_genus, yes = "H", no = 
                        ifelse(sum_rel_abundance < max_H & 
                               sum_rel_abundance > max_nSneg &
                               sum_count_abundance >= min_hard_count_genus, yes = "nS-", no =
                        ifelse(sum_rel_abundance > max_H & 
                               sum_rel_abundance > max_nSneg &
                               sum_count_abundance >= min_hard_count_genus, yes = "H & nS-", 
                        ifelse(sum_count_abundance < min_hard_count_genus, yes = "n/a, 'low abundance'", no = "none"))))) %>%
        mutate(tmp_col = factor(tmp_col, levels = c("H", "nS-", "H & nS-", "none", "n/a, 'low abundance'"))) %>%
        mutate(Gram = ifelse(Genus %in% sepsis_gramnegative, yes = "Gram-negative\nbacteria", no = 
                      ifelse(Genus %in% sepsis_grampositive, yes = "Gram-positive\nbacteria", no = NA))) %>% 
        rename(facets = "Gram") %>% 
        rename(variable = "Genus") %>% 
        mutate(label = NA) %>% 
        select(colnames(metadata))

options(repr.plot.width = 20, repr.plot.height = 10)
co_elevation_Splus <- rbind(data, metadata) %>% 
    mutate(tmp_col = factor(tmp_col, levels = c("H", "nS-", "H & nS-", "none", "n/a, 'low abundance'"))) %>%
    filter(nSIRS_class == "S+") %>% 
        ggplot(aes(x = FID, y = variable)) +
            geom_point(col = NA) + 
            geom_point(data = . %>% filter(facets != "metadata") %>% filter(sum_rel_abundance != 0)  %>% 
                           filter(tmp_col %in% c("H & nS-")) , 
                       aes(x = FID, y = variable, size = -log10(sum_rel_abundance), col = tmp_col, shape = tmp_col)) + 
            geom_text(data = . %>% filter(facets == "metadata") %>% filter(!is.na(label)),
                      aes(x = FID, y = variable, label = label), size = 3) + 
            theme_bw() + labs(col = "Abundance exceeds", shape = "Abundance exceeds", size = "Relative abundance,\n-log10", 
                            x = "", y = "") + 
            scale_color_manual(values = c("H" = "#6066B6", 
                                          "H & nS-" = "#A84750", 
                                          "nS-" = "#C7D0D6", 
                                          "n/a, 'low abundance'" = "black",
                                          "none" = "black")) +
            scale_shape_manual(values = c("H" = 16, 
                                          "H & nS-" = 17, 
                                          "nS-" = 15, 
                                          "n/a, 'low abundance'" = 1,
                                          "none" = 3)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                        strip.background = element_rect(color = "white", fill = "white"), 
                        strip.text.x = element_text(size = 10), 
                        panel.spacing = unit(1, "lines"), 
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        strip.text.y = element_text(angle = 0, size = 12, hjust = 0), 
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank()   
                ) + 
            facet_grid(rows = vars(facets), cols = vars(nSIRS_class), scales = "free", space = "free", drop = TRUE) + 
            guides(col = guide_legend(override.aes = list(size = 5))) + 
            scale_size_continuous(range = c(10, 1), breaks = c(1.5, 4.5, 7.5)) 


ggsave(paste0("../../output_figures/Fig4_co-elevation_S+.png"), plot = co_elevation_Splus, 
           width = 20, height = 14, units = "cm")
ggsave(paste0("../../output_figures/Fig4_co-elevation_S+.pdf"), plot = co_elevation_Splus, 
           width = 20, height = 14, units = "cm")
