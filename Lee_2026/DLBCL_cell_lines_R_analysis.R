library(dplyr)
library(ggplot2)
library(ggtree)
library(cowplot)


samples <- c(
  "Bcl2",
  "CBP1",
  "CBP2",
  "CMBP3",
  "CMBP7",
  "EBMP2",
  "EBMP3",
  "EZBP1",
  "EZBP6",
  "MCDP2",
  "MCDP4",
  "MCDP6",
  "SN2F",
  "SN2M",
  "WP2",
  "WP6"
)

samples_to_display_names <- c(
  "Bcl2"="B2",
  "CBP1"="C11-P1",
  "CBP2"="C11-P2",
  "CMBP3"="MC11-P3",
  "CMBP7"="MC11-P7",
  "EBMP2"="EZB-MYC-P2",
  "EBMP3"="EZB-MYC-P3",
  "EZBP1"="EZB-P1",
  "EZBP6"="EZB-P6",
  "MCDP2"="MCD-P2",
  "MCDP4"="MCD-P4",
  "MCDP6"="MCD-P6",
  "SN2F"="SN2-F-P1",
  "SN2M"="SN2-M-P7",
  "WP2"="SB2-P2",
  "WP6"="SB2-P6"
)

main_figure_display_names <- c(
  "B2",
  "C11-P1",
  "C11-P2",
  "EZB-P1",
  "EZB-P6",
  "EZB-MYC-P2",
  "SN2-F-P1",
  "SN2-M-P7"
)

supp_figure_display_names <- c(
  "SB2-P2",
  "SB2-P6",
  "EZB-MYC-P3",
  "MC11-P3",
  "MC11-P7",
  "MCD-P2",
  "MCD-P4",
  "MCD-P6"
)



data <- data.frame()
folder <- "define_clones_combined/"
for (i in seq(1, length(samples))) {
  filename = paste(folder, samples[i], "-final_collapse-unique_atleast-2_igblast_clone-pass.tsv", sep="")
  curr_data <- airr::read_rearrangement(filename)
  curr_data$sample_id <- samples_to_display_names[samples[i]]
  data <- rbind(data, curr_data)
}

define_clones_heavy <- filter(data, locus == "IGH")
define_clones_light <- filter(data, locus %in% c("IGL", "IGK"))

clonality_define_clones_heavy <- define_clones_heavy %>% group_by(sample_id, clone_id) %>% summarise(n=sum(duplicate_count)) %>% 
  mutate(freq=n/sum(n)) %>% group_by(sample_id) 

clonality_define_clones_heavy <- clonality_define_clones_heavy %>% group_by(sample_id) %>% 
  mutate(rank = rank(-freq, ties.method = "first"))

clonality_define_clones_heavy$chain <- "heavy"

clonality_define_clones_light <- define_clones_light %>% group_by(sample_id, clone_id) %>% summarise(n=sum(duplicate_count)) %>% 
  mutate(freq=n/sum(n)) %>% group_by(sample_id) %>% group_by(sample_id) %>% 
  mutate(rank = rank(-freq, ties.method = "first"))

clonality_define_clones_light$chain <- "light"

clonality_polished <- rbind(clonality_define_clones_heavy, clonality_define_clones_light)

clonality_polished$rank <- unlist(lapply(clonality_polished$rank, function(x) min(x, 30)))

clonality_polished$rank <- factor(clonality_polished$rank, levels=1:30)



make_plot <- function(data) {
  clonality_plot <- ggplot(data, aes(x = 3, y = freq, fill = rank, group = freq)) +
    geom_bar(stat = "identity", width = 1.5, colour="black") +
    geom_text(data=subset(data, freq > 0.14), aes(label = paste0(round(freq * 100, digits=1), "%")),
              position=position_stack(vjust=0.5), size=3) +
    coord_radial(theta="y", inner.radius= 0.6, expand=FALSE, rotate.angle = TRUE) +
    geom_text(
      data = reframe(group_by(data, sample_id, chain), n=sum(n), freq=freq, rank=rank, chain=chain),
      aes(label = paste0(formatC(n, big.mark = ",")), group = sample_id),
      stat = 'summary', fun = sum,
      fontface = "bold",
      size = 3.8,
      vjust = 4.5
    ) +
    scale_fill_brewer(palette = "Accent") +
    theme_void() +
    facet_grid(vars(chain), vars(sample_id), switch = "y") + 
    theme(legend.position = "none") +
    theme(panel.spacing.x=unit(0, "lines")) + 
    theme(panel.spacing.y=unit(-0.75, "lines")) +
    theme(panel.background = element_rect(fill = "lightgrey")) +
    theme(
      strip.text.x = element_text(
        size = 12, margin=margin(b=0.75)
      ),
      strip.text.y = element_text(
        size = 11, face = "italic"
      )
    )
  return(clonality_plot)
}


heavy_light_combined <- plot_grid(
  make_plot(filter(clonality_polished, sample_id %in% main_figure_display_names)),
  NULL,
  make_plot(filter(clonality_polished, sample_id %in% supp_figure_display_names)),
  ncol=1,
  rel_heights = c(1, 0.02, 1)
)

ggsave("combined_heavy_light_clonality_plot.pdf", plot=heavy_light_combined, width = 9, height = 6, units = "in", dpi = 300)

