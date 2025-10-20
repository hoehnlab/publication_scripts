library(alakazam)
library(data.table)
library(dowser)
library(dplyr)
library(ggplot2)
library(scoper)
library(shazam)
library(ggtree)
library(ggExtra)
library(ggpubr)
library(viridis)


references <- readIMGT(dir = "~/vdj")

flu_means_file <- "~/simble-validation/flu_data_means.rds"

if (file.exists(flu_means_file)) {
  flu_data_means <- readRDS(flu_means_file)
} else {
  flu_data <- readChangeoDb("~/docker-share/bcr-db_3_airr.tsv")
  
  flu_data <- flu_data %>% filter(productive)

  multi_heavy <- table(filter(flu_data, locus == "IGH")$cell_id)
  multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]
  
  flu_data <- filter(flu_data, !cell_id %in% multi_heavy_cells)
  
  # split cells by heavy and light chains
  heavy_cells <- filter(flu_data, locus == "IGH")$cell_id
  light_cells <- filter(flu_data, locus == "IGK" | locus == "IGL")$cell_id
  no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]
  
  flu_data <- filter(flu_data, !cell_id %in% no_heavy_cells)
  
  dist_nearest <- distToNearest(filter(flu_data, locus == "IGH"), nproc = 1)
  p1 <- ggplot(subset(dist_nearest, !is.na(dist_nearest)),
               aes(x = dist_nearest)) +
    theme_bw() +
    xlab("Hamming distance") + ylab("Count") +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    geom_histogram(color = "white", binwidth = 0.02) +
    theme(axis.title = element_text(size = 18))
  plot(p1)
  
  # find threshold for cloning automatically
  threshold_output <- shazam::findThreshold(dist_nearest$dist_nearest,
                                            method = "gmm", model = "gamma-norm",
                                            cutoff = "user", spc = 0.995)
  threshold <- threshold_output@threshold
  print(threshold)
  # TODO: make a note of this threshold
  
  plot(threshold_output, binwidth = 0.02, silent = TRUE) +
    theme(axis.title = element_text(size = 18))
  
  results <- hierarchicalClones(flu_data, cell_id = 'cell_id',
                                threshold = threshold, only_heavy = TRUE,
                                split_light = FALSE, summarize_clones = FALSE, nproc=3)
  
  
  # data <- filter(data, nchar(SEQUENCE_IMGT) == nchar(GERMLINE_IMGT))
  
  
  example <- resolveLightChains(results)
  
  resolved <- createGermlines(example, references = references, clone = "clone_subgroup_id", nproc = 1)
  
  hlclones = formatClones(resolved,chain="HL",
                          split_light=TRUE,heavy="IGH",cell="cell_id",
                          trait=c("vj_gene","c_call"), minseq = 5, collapse = TRUE, nproc = 3)
  
  # TODO make a note that we're only using clones that had at least 5 distinct seqs
  # TODO note that we are processing this data in the same way as Jensen et al. ....
  keep_clone_ids <- unique(hlclones$clone_id)
  results <- filter(resolved, clone_subgroup_id %in% keep_clone_ids)
  shmu_h <- filter(results, locus == "IGH")
  shmu_l <- filter(results, locus != "IGH")
  heavy_shmu <- observedMutations(shmu_h,
                                  germlineColumn="germline_alignment",
                                  regionDefinition=IMGT_V,
                                  frequency=TRUE, 
                                  combine=TRUE,
                                  nproc=3)
  means <- c()
  clone_id <- c()
  for(clone in unique(heavy_shmu$clone_id)){
    sub_obs <- filter(heavy_shmu, clone_id == clone)
    means <- append(means, mean(sub_obs$mu_freq))
    clone_id <- append(clone_id, clone)
  }
  mean_heavy_shmu <- data.frame(means_heavy = means, 
                                clone_id = clone_id)
  
  light_shmu <- observedMutations(shmu_l, sequenceColumn="sequence_alignment",
                                  germlineColumn="germline_alignment",
                                  regionDefinition=IMGT_V,
                                  frequency=TRUE, 
                                  combine=TRUE,
                                  nproc=3)
  means <- c()
  clone_id <- c()
  for(clone in unique(light_shmu$clone_id)){
    sub_obs <- filter(light_shmu, clone_id == clone)
    means <- append(means, mean(sub_obs$mu_freq))
    clone_id <- append(clone_id, clone)
  }
  mean_light_shmu <- data.frame(means_light = means,
                                clone_id = clone_id)
  
  flu_data_means <- merge(mean_heavy_shmu, mean_light_shmu, by = "clone_id")
  
  saveRDS(flu_data_means, file = "~/simble-validation/flu_data_means.rds")
}

means_flu <- ggplot(flu_data_means, aes(x = means_heavy, y = means_light)) + theme_bw() + 
  geom_point(aes(colour = NA), alpha = 0.7) +
  # geom_point(colour = "darkgrey", alpha = 0.7) + 
  geom_abline(colour = "black") + xlab("Heavy Chain SHM") +
  xlim(0, 0.146) + 
  ylim(0,0.146) +
  labs(colour = "Sample time \n(generations)") +
  ylab("") + ggtitle("Influenza") + theme_bw() +
  theme(legend.position= "hidden")

means_flu<- ggMarginal(means_flu, colour = "darkgrey")
means_flu

################################################
file = "/Volumes/HoehnK/jessie/simble-validation/500gen_selection_150clones_seed253437/all_samples_airr.tsv"

shm_full_data <- airr::read_rearrangement(file)
times <- c()
max <- 10
for (i in 1:max) {
  times <- append(times, rep(i*500/max, 150/max))
}
# times <- c(rep(100, 150/6), rep(150, 150/6), rep(200, 150/6), rep(250, 150/6), rep(300, 150/6), rep(350, 150/6))
names(times) <- 1:150
shm_full_data <- shm_full_data %>%
  split(shm_full_data$clone_id) %>%
  lapply(function(x) {filter(x, sample_time==times[x$clone_id[1]])}) %>%
  bind_rows()


shmu_h <- filter(shm_full_data, locus == "IGH")
shmu_l <- filter(shm_full_data, locus != "IGH")
heavy_shmu <- observedMutations(shmu_h,
                                germlineColumn="germline_alignment",
                                regionDefinition=IMGT_V,
                                frequency=TRUE, 
                                combine=TRUE,
                                nproc=3)
means <- c()
medians <- c()
clone_id <- c()
sample_times <- c()
for(clone in unique(heavy_shmu$clone_id)){
  sub_obs <- filter(heavy_shmu, clone_id == clone)
  means <- append(means, mean(sub_obs$mu_freq))
  medians <- append(medians, median(sub_obs$mu_freq))
  clone_id <- append(clone_id, clone)
  sample_times <- append(sample_times, times[clone][[1]])
}
mean_heavy_shmu <- data.frame(means_heavy = means, medians_heavy = medians, 
                              clone_id = clone_id, sample_time = sample_times)

light_shmu <- observedMutations(shmu_l, sequenceColumn="sequence_alignment",
                                germlineColumn="germline_alignment",
                                regionDefinition=IMGT_V,
                                frequency=TRUE, 
                                combine=TRUE,
                                nproc=3)
means <- c()
medians <- c()
clone_id <- c()
sample_times <- c()
for(clone in unique(light_shmu$clone_id)){
  sub_obs <- filter(light_shmu, clone_id == clone)
  means <- append(means, mean(sub_obs$mu_freq))
  medians <- append(medians, median(sub_obs$mu_freq))
  clone_id <- append(clone_id, clone)
  sample_times <- append(sample_times, times[clone][[1]])
}
mean_light_shmu <- data.frame(means_light = means, medians_light = medians, 
                              clone_id = clone_id, sample_time = sample_times)

heavy_shmu <- heavy_shmu[order(heavy_shmu$sequence_id),]
light_shmu <- light_shmu[order(light_shmu$sequence_id),]
mean_heavy_shmu <- mean_heavy_shmu[order(mean_heavy_shmu$clone_id),]
mean_light_shmu <- mean_light_shmu[order(mean_light_shmu$clone_id),]


sim_data <- data.frame(heavy = heavy_shmu$mu_freq, clone = heavy_shmu$clone_id, sample_time = heavy_shmu$sample_time)
sim_data <- merge(sim_data, data.frame(light= light_shmu$mu_freq, clone= light_shmu$clone_id), by = "clone")
sim_data_means <- merge(mean_heavy_shmu, mean_light_shmu, by = "clone_id")

saveRDS(sim_data_means, file = "~/simble-validation/sim_data_means.rds")
sim_data_means <- readRDS("~/simble-validation/sim_data_means.rds")

library(cowplot)
means_sim <- ggplot(sim_data_means, aes(x = means_heavy, y = means_light)) + theme_bw() + 
  geom_point(aes(colour = as.numeric(sample_time.x)), alpha = 0.7) + geom_abline(colour = "black") + xlab("Heavy Chain SHM") +
  xlim(0, 0.146) + 
  ylim(0,0.146) +
  ylab("Light Chain SHM") + ggtitle("Simble") +
  labs(colour = "Sample time \n(generations)") +
  scale_colour_viridis(direction=1, limits=c(0, 500)) +
  theme_bw() +
  theme(legend.key.size=unit(0.3, 'cm'))

legend <- get_legend(means_sim)
cowplot_legend <- cowplot::get_legend((means_sim))

means_sim <- means_sim + theme(legend.position = "hidden")

means_sim<- ggMarginal(means_sim, colour = "darkgrey")
means_sim

# shm_plot_legend <- plot_grid(means_sim, means_flu, nrow=1, legend=legend, rel_widths = c(1, 1, 0.4))
# shm_plot_legend
# shm_plot_legend <- ggarrange(means_sim, means_flu, ggplotify::as.ggplot(legend), nrow=1, widths=c(1, 1, 0.3))
shm_plot_legend <- cowplot::plot_grid(means_sim, means_flu, legend, rel_widths=c(1, 1, 0.5), rel_heights=c(1, 1, 0.5), nrow=1)
shm_plot_legend

ggsave(plot=shm_plot_legend,
       height=1.7,width=3.4,filename=paste0("~/simble-validation/heavy_light_shm_legend.pdf"))
ggsave(plot=shm_plot_legend,
       height=5,width=10.82,filename=paste0("~/simble-validation/heavy_light_shm_legend.pdf"))
ggsave(plot=shm_plot_legend,
       height=1.7,width=3.4,filename=paste0("~/simble-validation/heavy_light_shm_legend.svg"))

shm_plot=ggarrange(means_sim,means_flu,ncol=2)
shm_plot

ggsave(plot=shm_plot,
       height=5.26,width=10.82,filename=paste0("~/simble-validation/heavy_light_shm.pdf"))
saveRDS(shm_plot, file = "~/simble-validation/shm_plot.rds")
saveRDS(means_sim, file = "~/simble-validation/means_sim.rds")
saveRDS(means_flu, file = "~/simble-validation/means_flu.rds")
