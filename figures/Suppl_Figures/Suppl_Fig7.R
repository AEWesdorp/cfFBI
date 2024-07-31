## -----------------------------------------------------------------------------
# Code accompanying the publication: "Bacterial cell-free DNA profiling reveals co-elevation of multiple bacteria in newborn septic foals"
# Li-Ting Chen, Emmy Wesdorp et al. 

# Please direct any questions or comments to a.e.wesdorp@umcutrecht.nl

## -----------------------------------------------------------------------------
source("../Sourced/Elevated_bacteria.R")

list_legend <- list()
list_plots <- list()
for (GOI in sepsis_genera){
    df_species_sep <- fig3_contaminantFree %>%    
        filter(Genus == GOI) %>%   

        ungroup() %>%
        group_by(Species_short) %>%
        filter(sum(Exact.abundance,na.rm = TRUE) != 0) %>%
        ungroup() %>% 

        mutate(rel_abundance = ifelse(Relative.abundance != 0, yes = Relative.abundance, no = 10^-8)) %>% 
        mutate(abs_abundance = Exact.abundance) %>% 
        mutate(Taxa = Species_short) %>% 
        select(c(FID, nSIRS_class, nSIRS_score, Taxa, rel_abundance, abs_abundance)) %>%
        mutate(tpe = "Species\nincl. sum species") 

    df_species_sum <- fig3_contaminantFree %>%    
        filter(Genus == GOI) %>%   
        mutate(Taxa = "Sum Species") %>% 
        group_by(FID, nSIRS_class, nSIRS_score, Taxa) %>% 
            dplyr::summarize(rel_abundance = sum(Relative.abundance),
                             abs_abundance = sum(Exact.abundance), .groups = "keep") %>% 
        filter(abs_abundance != 0) %>% 
        mutate(tpe = "Species\nincl. sum species") 

    custom__labels <- function(breaks) {
        labels <- ifelse(breaks <= 1 & breaks !=10^-8, as.character(breaks), "")
        return(labels)
    }

    if(GOI != "Pseudomonas"){rel_widths_tmp = c(1, 1)}
    if(GOI == "Pseudomonas"){rel_widths_tmp = c(1, 3)}
    mainplot_GOI <- 
        rbind(df_species_sep, df_species_sum) %>% 
            mutate(tpe = factor(tpe, levels = c("Species\nincl. sum species", "Sum species\nexcl. 'low abundance'"))) %>% 
            ggplot(aes(y = FID)) + 
                geom_point(data = . %>% filter(tpe == "Species\nincl. sum species")  %>% filter(Taxa != "Sum Species"), 
                           aes(x = rel_abundance, col = Taxa, size = -log10(rel_abundance)), alpha = 0.3) +
                geom_point(data = . %>% filter(tpe == "Species\nincl. sum species")  %>% filter(Taxa == "Sum Species"), 
                           aes(x = rel_abundance, size = -log10(rel_abundance)), 
                               col = "black", alpha = 1, shape = 1, stroke = 0.5) + 
                facet_grid(rows = vars(nSIRS_class), scales = "free", space = "free", drop = TRUE) + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                        strip.background = element_rect(color = "white", fill = "white"), 
                        strip.text.x = element_text(size = 16), 
                        panel.spacing = unit(1, "lines"), 
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        strip.text = element_text(angle = 0, size = 12, hjust = 0), 
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank() 
                        ) + 
                scale_x_log10(labels = custom__labels) + 
                scale_size_continuous(range = c(10, 1)) + 
                ggtitle(GOI) + 
                xlab("Species Relative Abundance") + 
                guides(
                    color = guide_legend(title = GOI, order = 1, ncol = 2), # Color legend first
                    size = guide_legend(title = "Relative abundance,\n-log10", order = 2)   # Size legend second
                ) 

    # Extract the legend using ggplot_gtable
    g <- ggplotGrob(mainplot_GOI + theme(legend.box.margin = margin(0, 0, 0, 0)))

    # Extract legends from the plot
    legend_grobs <- g$grobs[grep("guide-box", g$layout$name)]

    # If multiple legends exist, combine them
    if (length(legend_grobs) > 1) {
      combined_legend <- cowplot::plot_grid(plotlist = legend_grobs, ncol = 1)
    } else {
      combined_legend <- legend_grobs[[1]]
    }

    # Remove the legend from plot
    figure_plot <- mainplot_GOI + theme(legend.position = "none")
    
    # Save in list 
    list_legend[[GOI]] <- combined_legend 
    list_plots[[GOI]] <- figure_plot
}

## -----------------------------------------------------------------------------
options(repr.plot.height = 30, repr.plot.width = 12)
l1.1 <- plot_grid(list_legend[[1]], list_legend[[2]], list_legend[[3]], list_legend[[4]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.1.png"), plot = l1.1, 
           width = 12, height = 30, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.1.pdf"), plot = l1.1, 
           width = 12, height = 30, units = "in")

l1.2 <- plot_grid(list_legend[[5]], list_legend[[6]], list_legend[[7]], list_legend[[8]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.2.png"), plot = l1.2, 
           width = 12, height = 30, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.2.pdf"), plot = l1.2, 
           width = 12, height = 30, units = "in")

l1.3 <- plot_grid(list_legend[[9]], list_legend[[10]], list_legend[[11]], list_legend[[12]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.3.png"), plot = l1.3, 
           width = 12, height = 30, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.3.pdf"), plot = l1.3, 
           width = 12, height = 30, units = "in")

l1.4 <- plot_grid(list_legend[[13]], list_legend[[14]], list_legend[[15]], list_legend[[16]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.4.png"), plot = l1.4, 
           width = 12, height = 30, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_l1.4.pdf"), plot = l1.4, 
           width = 12, height = 30, units = "in")

## -----------------------------------------------------------------------------
p1.1 <- plot_grid(list_plots[[1]], list_plots[[2]], list_plots[[3]], list_plots[[4]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.1.png"), plot = p1.1, 
           width = 12, height = 7, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.1.pdf"), plot = p1.1, 
           width = 12, height = 7, units = "in")

p1.2 <- plot_grid(list_plots[[5]], list_plots[[6]], list_plots[[7]], list_plots[[8]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.2.png"), plot = p1.2, 
           width = 12, height = 7, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.2.pdf"), plot = p1.2, 
           width = 12, height = 7, units = "in")

p1.3 <- plot_grid(list_plots[[9]], list_plots[[10]], list_plots[[11]], list_plots[[12]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.3.png"), plot = p1.3, 
           width = 12, height = 7, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.3.pdf"), plot = p1.3, 
           width = 12, height = 7, units = "in")

p1.4 <- plot_grid(list_plots[[13]], list_plots[[14]], list_plots[[15]], list_plots[[16]], ncol = 4)
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.4.png"), plot = p1.4, 
           width = 12, height = 7, units = "in")
ggsave(paste0("../../output_figures/SupplFig7_abundance_species_incl_sum_vertical_p1.4.pdf"), plot = p1.4, 
           width = 12, height = 7, units = "in")