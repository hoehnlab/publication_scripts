library(ggtree)
library(treeio)
library(dowser)
library(dplyr)
library(viridis)


neutral_imbalance <- readRDS("~/simble-validation/imbalance_metric_neutral_trees.rds")
selection_imbalance <- readRDS("~/simble-validation/imbalance_metric_selection_trees.rds")
hiv_imbalance <- readRDS("~/simble-validation/imbalance_metric_hiv_trees.rds")

# to pick example trees, we will select trees with imbalance metric close to median

# calculate median
hiv_median <- median(hiv_imbalance$imbalance_metric)

# calculate difference from median for each tree
hiv_imbalance$diff_from_median <- abs(hiv_imbalance$imbalance_metric - hiv_median)

# select the five trees closest to the median
hiv_closest_to_median <- hiv_imbalance %>% slice_min(diff_from_median, n=5)

hiv_closest_to_median <- scaleBranches(hiv_closest_to_median, edge_type = "mutations")
# plot and format a representative tree close to median
hiv_plots <- plotTrees(hiv_closest_to_median, scale=10)
g1 <- ggtree(hiv_plots[[2]]$data, ladderize = TRUE, linewidth=0.4) +
  geom_tippoint(aes(fill=time), pch=21) +
  scale_fill_viridis() + labs(fill="Time\n(weeks)") + 
  geom_treescale(width=10)
hiv_tree <- g1 + ggtitle("Chronic HIV\nInfection")
hiv_tree

# save the clone_id of the example tree to plot later
hiv_example_clone_id <- hiv_closest_to_median$clone_id[1]

# save the tree plot
saveRDS(hiv_tree, file = "~/simble-validation/hiv_tree.rds")


# as above, calculate median and difference from median for simble data
selection_median <- median(selection_imbalance$imbalance_metric)
selection_imbalance$diff_from_median <- abs(selection_imbalance$imbalance_metric - selection_median)

# simble data has the same number of sequences per tree, so we can just take the 5 trees closest to median
selection_closest_to_median <- selection_imbalance %>% slice_min(diff_from_median, n=5)

# plot, format and save a representative tree from those closest to median
selection_closest_to_median <- scaleBranches(selection_closest_to_median, edge_type = "mutations")
selection_plots <- plotTrees(selection_closest_to_median, scale=10)
g1 <- ggtree(selection_plots[[5]]$data, ladderize = TRUE, linewidth=0.4) +
  geom_tippoint(aes(fill=sample_time), pch=21) +
  scale_fill_viridis() + labs(fill="Time\n(gen)") + 
  geom_treescale(width=10)
selection_tree <- g1 + ggtitle("SimBLE\nSelection")
selection_tree
saveRDS(selection_tree, file = "~/simble-validation/selection_tree.rds")

# save the clone_id of the example tree to plot later
selection_example_clone_id <- selection_closest_to_median$clone_id[2]


# as above, calculate median and difference from median for simble data, and select trees closest to median
neutral_median <- median(neutral_imbalance$imbalance_metric)
neutral_imbalance$diff_from_median <- abs(neutral_imbalance$imbalance_metric - neutral_median)
neutral_closest_to_median <- neutral_imbalance %>% slice_min(diff_from_median, n=5)

# plot, format and save a representative tree from those closest to median
neutral_closest_to_median <- scaleBranches(neutral_closest_to_median, edge_type = "mutations")
neutral_plots <- plotTrees(neutral_closest_to_median, scale=10)

g1 <- ggtree(neutral_plots[[3]]$data, ladderize = TRUE, linewidth=0.4) +
  geom_tippoint(aes(fill=sample_time), pch=21) +
  scale_fill_viridis() + labs(fill="Time\n(gen)") + 
  geom_treescale(width=10)
neutral_tree <- g1 + ggtitle("SimBLE\nNeutral")
neutral_tree
saveRDS(neutral_tree, file = "~/simble-validation/neutral_tree.rds")

# save the clone_id of the example tree to plot later
neutral_example_clone_id <- neutral_closest_to_median$clone_id[1]


# process and save the clone ids that were picked as example trees
example_tree_clone_ids <- c(hiv_example_clone_id, selection_example_clone_id, neutral_example_clone_id)
names(example_tree_clone_ids) <- c("hiv", "selection", "neutral")
saveRDS(example_tree_clone_ids, file = "~/simble-validation/example_tree_clone_ids.rds")


