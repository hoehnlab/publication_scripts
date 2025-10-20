library(ggplot2)
library(viridis)
library(ggpubr)
library(ggExtra)
library(ggtree)
library(cowplot)

# load all the plot objects

font_size <- theme(axis.text = element_text(size=7),
                   axis.title = element_text(size=9),
                   plot.title = element_text(size=9),
                   strip.text = element_text(size=8),
                   legend.text = element_text(size=7),
                   legend.title = element_text(size=8))

legend_size <- theme(legend.key.size = unit(0.25, "cm"), legend.background = element_rect(fill = "transparent"))

selection_strength_plot <- readRDS("~/simble-validation/selection_strength_with_selection_plot.rds")

tweaked_selection <- selection_strength_plot + ggtitle("") + 
  geom_vline(xintercept=0, linetype="dashed", color="black", linewidth=0.25) + theme_bw() + 
  labs(color = "Time (gen)", x = expression(paste("Selection strength (", Sigma, ")"))) +
  scale_color_viridis(discrete = TRUE, direction=-1) +
  theme(panel.grid.major.y = element_blank()) +
  font_size + legend_size + theme(legend.box.margin = margin(t=-7, r=0, b=-5, l=0, unit="mm")) + 
  theme(legend.position = "top") +
  guides(colour = guide_legend(direction="horizontal", reverse=TRUE))
tweaked_selection

neutral_selection_strength_plot <- readRDS("~/simble-validation/selection_strength_neutral_plot.rds")

tweaked_selection_neutral <- neutral_selection_strength_plot + ggtitle("") + 
  geom_vline(xintercept=0, linetype="dashed", color="black", linewidth=0.25) + theme_bw() + 
  # theme(panel.grid.major = element_line(colour="darkgrey", size=0.5)) +
  labs(color = "Time (gen)", x = expression(paste("Selection strength (", Sigma, ")"))) +
  scale_color_viridis(discrete = TRUE, direction=-1) +
  theme(panel.grid.major.y = element_blank()) +
  font_size + legend_size + theme(legend.box.margin = margin(t=-7, r=0, b=-5, l=0, unit="mm")) + 
  theme(legend.position = "top") +
  # ylim(layer_scales(tweaked_selection)$y$range$range) +
  guides(colour = guide_legend(direction="horizontal", reverse=TRUE))
tweaked_selection_neutral


######## Migration and cell types plots

celltype_levels = c("Memory B Cell", "Plasma Cell", "Default")
celltype_palette = c("#FFB000", "#648FFF", "#B0B0B0")

names(celltype_palette) = celltype_levels

angle_x_text <- theme(axis.text.x = element_text(angle = 25))

migration_plot <- readRDS("~/simble-validation/simble_migration_plot.rds")
migration_plot <- migration_plot + font_size + 
  scale_fill_manual(values=celltype_palette) + 
  xlab("Time of Differentiation (gen)") +
  angle_x_text +
  theme(axis.title.x = element_text(margin = margin(t = -5))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
  theme(axis.text.x=element_text(size=7)) +
  theme(plot.margin = margin(t=0, r=1, b=0, l=1, unit="mm"))

affinity_plot <- readRDS("~/simble-validation/affinity_plot.rds")
cross_reactivity_plot <- readRDS("~/simble-validation/cross_reactivity_plot.rds")
comps <- list(c("Memory B Cell", "Plasma Cell"))
affinity_plot <- affinity_plot + font_size + 
  scale_color_manual(values=celltype_palette) + 
  ylab("Affinity Score") + 
  stat_compare_means(method = "wilcox.test", label="p", paired = TRUE, size=8/.pt) +
  theme(legend.position = "none") +
  theme(plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text.x = element_text(size=9, color="black")) +
  theme(axis.title.x=element_blank())
cross_reactivity_plot <- cross_reactivity_plot + font_size + 
  scale_color_manual(values=celltype_palette) + 
  ylab("") +
  stat_compare_means(method = "wilcox.test", label="p", paired = TRUE, size=8/.pt) +
  theme(legend.position = "none") +
  theme(plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text.x = element_text(size=9, color="black")) +
  theme(axis.title.x=element_blank())

affinity_crossreactivity_plot <- plot_grid(affinity_plot, NULL, cross_reactivity_plot, nrow=1, rel_widths=c(1, -0.15, 1))
affinity_crossreactivity_plot



############# SHM Plots

flu_data_means <- readRDS(flu_means_file)
means_flu <- ggplot(flu_data_means, aes(x = means_heavy, y = means_light)) + theme_bw() + 
  geom_point(aes(colour = NA), alpha = 0.7) +
  # geom_point(colour = "darkgrey", alpha = 0.7) + 
  geom_abline(colour = "black") + xlab("Heavy Chain SHM") +
  xlim(0, 0.146) + 
  ylim(0,0.146) +
  labs(colour = "Time\n(gen)") +
  ylab("") + ggtitle("Influenza")+ 
  theme_bw() +
  theme(legend.position= "hidden") + font_size + theme(
    panel.background = element_rect(fill = "transparent", color = NA), # Background of the main panel
    plot.background = element_rect(fill = "transparent", color = NA)    # Background of the entire plot area
  ) +
  # theme(plot.title = element_text(vjust=-2.5)) 
  theme(plot.title = element_blank()) +
  annotate("text", x = -Inf, y = Inf, label = "Influenza",
           hjust = -0.1, vjust = 1.2, size = 3)

means_flu<- ggMarginal(means_flu, colour = "darkgrey")


sim_data_means <- readRDS("~/simble-validation/sim_data_means.rds")

means_sim <- ggplot(sim_data_means, aes(x = means_heavy, y = means_light)) + theme_bw() + 
  geom_point(aes(colour = as.numeric(sample_time.x)), alpha = 0.7) + geom_abline(colour = "black") + 
  xlab("Heavy Chain SHM") +
  xlim(0, 0.146) + 
  ylim(0,0.146) +
  ylab("Light Chain SHM") +
  ggtitle("SimBLE")+
  labs(colour = "Time \n(gen)") +
  scale_colour_viridis(direction=1, limits=c(0, 500), breaks=c(100, 300, 500)) +
  theme_bw() + font_size + legend_size + theme(
    panel.background = element_rect(fill = "transparent", color = NA), # Background of the main panel
    plot.background = element_rect(fill = "transparent", color = NA)    # Background of the entire plot area
  ) +
  theme(legend.text = element_text(angle = 45, hjust=1, vjust=1.2),
        legend.title.position = "left") +
  guides(colour = guide_colorbar(direction="horizontal", reverse=FALSE)) +
  theme(plot.title = element_blank()) +
  annotate("text", x = -Inf, y = Inf, label = "SimBLE",
             hjust = -0.1, vjust = 1.2, size = 3)

legend <- get_legend(means_sim)

means_sim <- means_sim + theme(legend.position = "hidden")

means_sim<- ggMarginal(means_sim, colour = "darkgrey")
shm_plot <- plot_grid(means_sim, NULL, means_flu, rel_widths=c(1, -0.1, 1), nrow=1)

shm_plot_legend <- plot_grid(legend, NULL, shm_plot, ncol=1, rel_heights = c(0.1, -0.1, 1))
shm_plot_legend


####### Trees

neutral_tree <- readRDS("~/simble-validation/neutral_tree.rds")
hiv_tree <- readRDS("~/simble-validation/hiv_tree.rds")
selection_tree <- readRDS("~/simble-validation/selection_tree.rds")

scale <- geom_treescale(fontsize = 7/.pt, width=0.02, offset=1.5)

neutral_tree <- neutral_tree + font_size + legend_size + 
  coord_cartesian(clip = "off") + theme(plot.title = element_text(hjust=1)) + 
  scale
selection_tree <- selection_tree + font_size + legend_size + 
  coord_cartesian(clip = "off") + theme(plot.title = element_text(hjust=1)) + 
  theme(legend.position = c(0.2, 0.8)) + 
  scale
hiv_tree <- hiv_tree + font_size + legend_size + 
  coord_cartesian(clip = "off") + theme(plot.title = element_text(hjust=1)) + 
  theme(legend.position = c(0.86, 0.35)) + 
  scale



###### Imbalance stat
neutral_imbalance <- readRDS("~/simble-validation/imbalance_metric_neutral_trees.rds")
selection_imbalance <- readRDS("~/simble-validation/imbalance_metric_selection_trees.rds")
hiv_imbalance <- readRDS("~/simble-validation/imbalance_metric_hiv_trees.rds")

highlight_ids <- readRDS("~/simble-validation/example_tree_clone_ids.rds")

neutral_imbalance$highlight <- ifelse(neutral_imbalance$clone_id == highlight_ids["neutral"], "TRUE", "FALSE")
selection_imbalance$highlight <- ifelse(selection_imbalance$clone_id == highlight_ids["selection"], "TRUE", "FALSE")
hiv_imbalance$highlight <- ifelse(hiv_imbalance$clone_id == highlight_ids["hiv"], "TRUE", "FALSE")

imbalance_colors <- c("neutral" = "#FE6100", "selection" = "#DC267F", "hiv" = "#785EF0")
all_trees <- neutral_imbalance[,c("clone_id", "imbalance_metric", "model", "highlight")]
all_trees <- rbind(all_trees, selection_imbalance[,c("clone_id", "imbalance_metric", "model", "highlight")])
all_trees <- rbind(all_trees, hiv_imbalance[,c("clone_id", "imbalance_metric", "model", "highlight")])
comps <- list(c("neutral", "selection"), c("neutral", "hiv"), c("selection", "hiv"))
imbalance_plot <- ggplot(all_trees, aes(y=imbalance_metric, x=model)) + 
  geom_boxplot(aes(fill=model), outlier.shape=NA) + geom_jitter(data = filter(all_trees, highlight == "FALSE"), color="black", width=0.15, height = 0, size=0.7, alpha=0.7) +
  geom_jitter(data = filter(all_trees, highlight == "TRUE"), color="#00FF00", width=0.15, height = 0, size=0.7, alpha=1) +
  scale_fill_manual(values=imbalance_colors)
dataorder <- c( "neutral", "selection", "hiv")
imbalance_plot <- imbalance_plot + 
  labs(y="Tree Balance") + 
  ggtitle("") +
  scale_x_discrete(limits=dataorder,labels=c("Neutral","Selection","HIV\nInfection")) +
  stat_compare_means(comparison = comps, method = "wilcox.test", size=8/.pt, label.y = 3, step.increase = 0.19) +
  theme_bw() + 
  font_size +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=9, color="black")) +
  theme(axis.title.x = element_blank()) +
  # ylim(NA, 0.33) +
  font_size  +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))


####### Combine them into a single figure with space for the simble outline

row1 <- plot_grid(NULL, labels=c("A)"), nrow = 1, label_fontfamily="Helvetica", label_fontface = "plain")
row2 <- plot_grid(neutral_tree + theme(legend.position = "none"), 
                  selection_tree, 
                  hiv_tree,
                  imbalance_plot,
                  labels=c("B)", "", "", "C)"), nrow = 1, rel_widths=c(0.8, 1, 1, 1.1), label_fontfamily="Helvetica", label_fontface = "plain")
row3 <- plot_grid(tweaked_selection, shm_plot_legend,
                  nrow=1, labels=c("D)", "E)"),rel_widths = c(1.1, 1),
                  label_fontfamily="Helvetica", label_fontface = "plain")
row4 <- plot_grid(migration_plot, 
                  affinity_crossreactivity_plot, 
                  labels=c("F)", "G)"), nrow=1, rel_widths = c(1.5, 1),
                  label_fontfamily="Helvetica", label_fontface = "plain")


figure2 <- plot_grid(row1, row2, row3, row4, ncol=1, rel_heights=c(0.7, 1.2, 1, 1.1), rel_widths=c(1, 1, 1, 1))
figure2

ggsave(plot=figure2, filename="~/simble-validation/figure2.pdf", width=7, height=7.5, units="in")
ggsave(plot=figure2, filename="~/simble-validation/figure2cowplot.svg", width=7, height=7.5, units="in")

ggsave(plot=tweaked_selection_neutral, filename="~/simble-validation/SupplementalFig1.pdf", width=4.5, height=3.5, units="in")

###### mutations per site (supp 2)
mutations_per_site <- readRDS("~/simble-validation/mutations_per_site_plot.rds")
mutations_per_site <- mutations_per_site + font_size + legend_size
mutations_per_site

ggsave(plot=mutations_per_site, filename="~/simble-validation/SupplementalFig2.pdf", width=6, height=5.5, units="in")
