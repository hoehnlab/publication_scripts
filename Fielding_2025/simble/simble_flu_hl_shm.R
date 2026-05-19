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
source("./flu_data_processing.R")


references <- readIMGT(dir = "~/vdj")

flu_means_file <- "~/simble-validation/flu_data_means.rds"

if (file.exists(flu_means_file)) {
  flu_data_means <- readRDS(flu_means_file)
} else {
  results <- getProcessedFluResults()
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
  geom_abline(colour = "black") + xlab("Heavy Chain SHM") +
  xlim(0, 0.146) + 
  ylim(0,0.146) +
  labs(colour = "Sample time \n(generations)") +
  ylab("") + ggtitle("Influenza") + theme_bw() +
  theme(legend.position= "hidden")

means_flu<- ggMarginal(means_flu, colour = "darkgrey")
means_flu

################################################
file = "200gen_selection_150clones_seed253437/all_samples_airr.tsv"

shm_full_data <- airr::read_rearrangement(file)
times <- c()
max <- 10
for (i in 1:max) {
  times <- append(times, rep(i*200/max, 150/max))
}

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
  scale_colour_viridis(direction=1, limits=c(0, 200)) +
  theme_bw() +
  theme(legend.key.size=unit(0.3, 'cm'))

legend <- get_legend(means_sim)
cowplot_legend <- cowplot::get_legend((means_sim))

means_sim <- means_sim + theme(legend.position = "hidden")

means_sim<- ggMarginal(means_sim, colour = "darkgrey")
means_sim

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

