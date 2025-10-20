library(dplyr)
library(airr)

folder <- "/Volumes/HoehnK/jessie/simble-validation/150gen_selection_100clones_differentiation_study_seed344276/"
simble_data <- read_rearrangement(paste0(folder, "all_samples_airr.tsv", sep=""))

simble_trees <- treeio::read.beast(paste0(folder, "all_simplified_trees.nex"))
simble_tip_data <- fortify(simble_trees) %>% filter(isTip)
simble_tip_data$time_of_differentiation <- as.numeric(simble_tip_data$time_of_differentiation)

all_other_samples <- simble_tip_data %>%
  mutate(celltype = recode(celltype, "plasma_cell" = "Plasma Cell", "memory_b_cell" = "Memory B Cell", "gc_b_cell" = "GC B Cell")) %>%
  filter(celltype != "GC B Cell")

simble_other_binned <- all_other_samples %>%
  mutate(differentiation_time_binned = cut(time_of_differentiation,
                               breaks = c(0:15*10, Inf), # Define bin boundaries
                               labels = paste0(0:15*10, "-", c(1:16*10)), # Define bin labels
                               include.lowest = TRUE))

simble_other_binned_proportions <- simble_other_binned %>%
  group_by(celltype, differentiation_time_binned) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))


saveRDS(simble_other_binned_proportions, file="~/simble-validation/simble_celltype_proportions_binned.rds")
simble_other_binned_proportions <- readRDS("~/simble-validation/simble_celltype_proportions_binned.rds")
distribution_celltypes_by_differentiation_time <- ggplot(simble_other_binned_proportions, aes(x=differentiation_time_binned)) +
  facet_wrap(~ celltype, nrow=2, scales = "free_y") +
  geom_bar(stat="identity", aes(y=proportion*100, fill=celltype), width = 0.75, color="black") +
  labs(x = "Time of Differentiation (generations)",
       y = "Proportion of cell type (%)") +
  theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.minor.y = element_blank()) +
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=10))

distribution_celltypes_by_differentiation_time
saveRDS(distribution_celltypes_by_differentiation_time, file="~/simble-validation/simble_migration_plot.rds")

