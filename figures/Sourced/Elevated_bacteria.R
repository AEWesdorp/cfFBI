## -----------------------------------------------------------------------------
# Code accompanying the publication: "Bacterial cell-free DNA profiling reveals co-elevation of multiple bacteria in newborn septic foals"
# Li-Ting Chen, Emmy Wesdorp et al. 

# Please direct any questions or comments to a.e.wesdorp@umcutrecht.nl

## -----------------------------------------------------------------------------

require("dplyr")
require("ggplot2")
require("ggpubr")
require("stringr")
require("ggh4x")
require("cowplot")
require("DescTools")
require("Hmisc")

## -----------------------------------------------------------------------------
min_hard_count_species = 10
min_hard_count_genus = 10

## -----------------------------------------------------------------------------
sepsis_gramnegative <- c("Serratia", "Salmonella", "Pseudomonas", "Proteus", "Pasteurella", "Pantoea", "Klebsiella", "Escherichia", "Enterobacter", "Aeromonas", "Actinobacillus", "Acinetobacter") 
sepsis_grampositive <- c("Streptococcus", "Staphylococcus", "Enterococcus", "Bacillus")
sepsis_genera <- c(sepsis_grampositive, sepsis_gramnegative)

## -----------------------------------------------------------------------------
all_basic_stats <- read.csv(file = "../../output/02_tables/02_data_merged/all_basic_stats.csv") %>% as.data.frame()
fig3_basics = read.csv(file = "../../output/02_tables/04_source_data/contaminantFree_bacteria_long_inc_batched.csv") %>% as.data.frame()

## -----------------------------------------------------------------------------
## Removal of the contaminants
fig3_contaminantFree <- fig3_basics %>% 
    filter(Contaminant == "False") %>% 
    mutate(Contaminant = factor(Contaminant, level = unique(Contaminant))) %>% 
    mutate(Species_short = str_replace(Species, "^[^\\s]+\\s", "")) %>% 

    #removing misclassified species, meaning those placed incorrectly in a higher taxonomic rank, indicated by squared brackets 
    filter(!str_detect(Genus, "\\[")) %>% 
    filter(Genus != "all_filtered_reads_exc_listed_here")   

## -----------------------------------------------------------------------------
## Preprocessing data for genus level assesment
fig3_contaminantFree_sumGenus <- fig3_contaminantFree %>% 
    group_by(FID, nSIRS_class, nSIRS_score, Genus) %>% 
    dplyr::summarize(sum_rel_abundance = sum(Relative.abundance),
                     sum_count_abundance = sum(Exact.abundance), .groups = "keep") %>% 
    ungroup() %>%
    group_by(Genus) %>%
    filter(sum(sum_rel_abundance) != 0) %>%
    ungroup()

max_abundance_H <- fig3_contaminantFree_sumGenus %>%
        filter(nSIRS_class == "H") %>% 
        group_by(Genus) %>%
        dplyr::summarize(max_H = max(sum_rel_abundance), .groups = "keep")
    
max_abundance_nSneg <- fig3_contaminantFree_sumGenus %>%
        filter(nSIRS_class == "nS-") %>% 
        group_by(Genus) %>%
        dplyr::summarize(max_nSneg = max(sum_rel_abundance), .groups = "keep")

max_abundance_Spos <- fig3_contaminantFree_sumGenus %>%
        filter(nSIRS_class == "S+") %>% 
        group_by(Genus) %>%
        dplyr::summarize(max_Spos = max(sum_rel_abundance), .groups = "keep")

fig3_contaminantFree_sumGenus_maxAbundance <- fig3_contaminantFree_sumGenus %>%
        left_join(max_abundance_H, by = c("Genus")) %>% 
        left_join(max_abundance_nSneg, by = c("Genus")) %>% 
        left_join(max_abundance_Spos, by = c("Genus")) %>% 
        mutate(
            max_H = coalesce(max_H, 0),
            max_nSneg = coalesce(max_nSneg, 0),
            max_Spos = coalesce(max_Spos, 0)
        )

## -----------------------------------------------------------------------------
#split table in sepsis causing pathogens and other genera
fig3_contaminantFree_sumGenus_maxAbundance_otherGenera <- 
    fig3_contaminantFree_sumGenus_maxAbundance %>% 
        filter(Genus %nin% sepsis_genera)

fig3_contaminantFree_sumGenus_maxAbundance_sepsisGenera <- 
    fig3_contaminantFree_sumGenus_maxAbundance %>% 
        filter(Genus %in% sepsis_genera) %>% 
        mutate(Genus = factor(Genus, levels = sepsis_genera))
