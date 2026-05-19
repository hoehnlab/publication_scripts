library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

main_summary <- read.csv("/Volumes/HoehnK/Sherry/beast_workspace/TyCHE/figures/dec2025/main_analysis/plot_data/summary_data_main.csv")

main_summary <- main_summary %>%
  mutate(model_display = recode(model_display, "TyCHE\nest clocks"="TyCHE"))
main_summary$model_display <- factor(main_summary$model_display, levels=c("TyCHE", "SC", "UCLD"))
main_summary$tree_height_prop_error <- abs(main_summary$tree_height_prop_error)
main_summary$tree_length_prop_error <- abs(main_summary$tree_length_prop_error)

font_size <- theme(axis.text = element_text(size=7),
                   axis.title = element_text(size=9),
                   plot.title = element_text(size=9),
                   strip.text = element_text(size=8),
                   legend.text = element_text(size=7),
                   legend.title = element_text(size=8))


standard_additions <- function() {
  my_comparisons <- list( c("TyCHE", "SC"), c("TyCHE", "UCLD"))
  
  list(
    scale_fill_manual(values=c("TyCHE" = "darkorange", "SC"="green", "UCLD"="purple")),
    
    theme_bw(),
    
    theme(legend.position = "none"),
    
    stat_compare_means(
      label="p.signif", 
      comparisons = my_comparisons, 
      method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                         symbols = c("***", "**", "*", "ns")), 
      paired=TRUE, 
      size=8/.pt,
      vjust = 0.5),
    
    font_size
  )
}
height <- ggplot(main_summary, aes(x = model_display, y=tree_height_prop_error, fill = model_display)) + 
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) + 
  geom_jitter(shape=21,  # control dodge grouping
               size = 1,
               position = position_jitterdodge(
    jitter.height = 0,
    jitter.width = 0.5,
    dodge.width = 0.85)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_wrap(~facet_label_with_sim, scales="free_y", labeller=as_labeller(c("Selective Evolution" = "Primary GC", "Selective Evolution\n(GC re-entry)"= "Recall GC"))) +
  standard_additions() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(x="", y="Tree Height Proportional Error") 

length <- ggplot(main_summary, aes(x = model_display, y=tree_length_prop_error, fill = model_display)) + 
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) + 
  geom_jitter(shape=21,  # control dodge grouping
               size = 1,
               position = position_jitterdodge(
    jitter.height = 0,
    jitter.width = 0.5,
    dodge.width = 0.85)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_wrap(~facet_label_with_sim, scales="free_y") +
  standard_additions() +
  labs(x="", y="Tree Length Proportional Error") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) 
length

angle_x_text <- theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"))
rf <- ggplot(main_summary, aes(x = model_display, y=rf_distance, fill = model_display)) + 
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) + 
  geom_jitter(shape=21,  # control dodge grouping
              size = 1, position = position_jitterdodge(
    jitter.height = 0,
    jitter.width = 0.5,
    dodge.width = 0.85)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_wrap(~facet_label_with_sim, scales="free_y") +
  standard_additions() +
  labs(y="RF distance") + 
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  angle_x_text
rf

combined_plot <- plot_grid(height, length, rf, ncol=1, align="v", rel_heights=c(0.92, 0.87, 1))
combined_plot
ggsave(combined_plot, width=2.6, height=4.8, filename="recall/Fig3C_selection.pdf")
