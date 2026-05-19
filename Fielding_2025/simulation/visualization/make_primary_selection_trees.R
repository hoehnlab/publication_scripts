########################
## Hunter J. Melton   ##
## Jessie J. Fielding ##
###### 4/28/2026 #######

# Makes tree plots for Fig 3A.

library(airr)
library(dplyr)
library(dowser)
library(ggtree)
library(purrr)
library(treeio)
library(data.table)
library(patchwork)
source("./plot_extras.R")

simulation_id = "sel"
folder_name = "primary"

path_to_results <- paste0("./output_", folder_name, "/beast_output/", simulation_id, "/")

model_names <- c(
  "ExpectedOccupancy_EstClockRates",
  "StrictClock_AncestralReconstruction",
  "UCRelaxedClock_AncestralReconstruction"
)

true_tree_height <- 200
ignore_ess_params <- c("freqParameter", "traitfrequencies", "rateIndicator", "rootType", "rateCategories")


#####################################################################
# Load in the log.tsv data, record tree heights, record convergence #
#####################################################################
all_results <- list()

for (i in 1:length(model_names)){
  all_results[[model_names[i]]] <- data.frame(clone_id = 1:20, tree_height = NA, tree_height_prop_error = NA, converged = NA)
  
  for (j in 1:20){
    file_name <- paste0(path_to_results, model_names[i], "_", j, "_log.tsv" )
    if (file.exists(file_name)) {
      
      # Load in log.tsv
      log <- fread(file_name, header=TRUE) %>% as.data.frame()
      
      # Check for convergence
      check_ess_params <- unique(log$item)
      for(regex in ignore_ess_params){
        check_ess_params <- check_ess_params[!grepl(regex, check_ess_params)]
      }
      
      ess_fail_check <- log %>% filter(item %in% check_ess_params) %>% filter(ESS < 100) %>% nrow()
      ess_fail_200 <- log %>% filter(item %in% check_ess_params) %>% filter(ESS < 200) %>% nrow()
      if(ess_fail_200 == 0){
        all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "converged"] <- "Yes"
      } else if(ess_fail_check == 0) {
        all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "converged"] <- "100"
      } else {
        all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "converged"] <- "No"
      }
      
      # Get the tree height
      tree_height <- log %>% filter(item == "TreeHeight") %>% pull(mean)
      all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "tree_height"] <- tree_height
      all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "tree_height_prop_error"] <- abs(tree_height - true_tree_height) / true_tree_height
      
    } else {
      all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "converged"] <- "No"
    }
    
  }
  cat("Finished processing model:", model_names[i], "\n")
}


# Bind together the results with a model name column and filter for convergence

final_results <- bind_rows(lapply(names(all_results), function(model_name) {
  cbind(model = model_name, all_results[[model_name]])
})) %>%
  mutate(model = recode(model,
                        "ExpectedOccupancy_EstClockRates" = "EO Est",
                        "StrictClock_AncestralReconstruction" = "SC",
                        "UCRelaxedClock_AncestralReconstruction" = "UCLD"))
final_results$model <- factor(final_results$model, levels = c("EO Est", "SC", "UCLD"))




all_trees <- list()

colors = RColorBrewer::brewer.pal(2,"Set1")
names(colors) = c("GC", "other")

for (i in 1:length(model_names)){
  all_trees[[model_names[i]]] <- list()
  all_trees_plots[[model_names[i]]] <- list()
  branches_list[[model_names[i]]] <- list()
  
  for (j in 1:20){
    file_name <- paste0(path_to_results, model_names[i], "_", j, "_tree_with_trait.tree" )
    if (file.exists(file_name) & all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "converged"] == "Yes") {
      
      # Load in tree
      tree <- read.beast(file_name)
      
      # Get correct relative times for first and second rxn
      tree_height <- all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "tree_height"]
      
      # Store tree
      all_trees[[model_names[i]]][[as.character(j)]] <- tree
      
    } else {
      all_trees[[model_names[i]]][[as.character(j)]] <- NULL
    }
    
  }
  cat("Finished processing trees for model:", model_names[i], "\n")
}




true_data_folder <- "/Volumes/HoehnK/Sherry/beast_workspace/TyCHE/tltt_12_19/data/raw/config_ratio_1to1_sel/"
true_simple_trees <- read.beast(file.path(true_data_folder, "all_simplified_time_trees.nex"))
true_full_trees <- read.beast(file.path(true_data_folder, "all_pruned_time_trees.nex"))



names(colors) = c("germinal_center", "other", "Ambig.")

node_size <- 1.0

create_plots <- function(clone_id){
  temp_sc <- all_trees[["StrictClock_AncestralReconstruction"]][[as.character(clone_id)]]
  if (temp_sc %>% is.null()) {
    return(NULL)
  }
  temp_sc@data$prob_gc <- NA
  for (i in 1:nrow(temp_sc@data)) {
    celltype_set <- temp_sc@data$location.set[i][[1]]
    celltype_set_prob <- temp_sc@data$location.set.prob[i][[1]]
    names(celltype_set_prob) <- celltype_set
    temp_sc@data$prob_gc[i] <- ifelse("germinal_center" %in% celltype_set, celltype_set_prob["germinal_center"], 0.0) %>% as.numeric()
  }
  temp_tyche <- all_trees[["ExpectedOccupancy_EstClockRates"]][[as.character(clone_id)]]
  temp_tyche@data$prob_gc <- NA
  for (i in 1:nrow(temp_tyche@data)) {
    celltype_set <- temp_tyche@data$location.set[i][[1]]
    celltype_set_prob <- temp_tyche@data$location.set.prob[i][[1]]
    names(celltype_set_prob) <- celltype_set
    temp_tyche@data$prob_gc[i] <- ifelse("germinal_center" %in% celltype_set, celltype_set_prob["germinal_center"], 0.0) %>% as.numeric()
  }
  
  tyche_plot <- ggtree(temp_tyche, aes(color = as.numeric(expectedOccupancies)), linewidth = 0.3) +
    geom_nodepoint(aes(fill = prob_gc), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    geom_tippoint(aes(fill = prob_gc), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    ggtitle("TyCHE") +
    eo_node_fill_scale +
    eo_line_color_scale +
    theme(legend.position = "none") +
    coord_cartesian(clip="off")
  tyche_plot$data$x = tyche_plot$data$x - max(tyche_plot$data$x) + 200
  
  temp_simple_true <- true_simple_trees[[clone_id]]
  simple_ids <- temp_simple_true@data$cell_id
  temp_full_true <- true_full_trees[[clone_id]]
  
  temp_full_true@data$to_show <- FALSE
  temp_full_true@data$to_show[temp_full_true@data$cell_id %in% simple_ids] <- TRUE
  temp_full_true@data$branch_location <- temp_full_true@data$location
  temp_full_true@data$branch_location[temp_full_true@data$antigen == "1000" & temp_full_true@data$generation == "1100"] <- "other"
  
  
  true_plot <- ggtree(temp_full_true, aes(color = branch_location), linewidth = 0.3) +
    geom_nodepoint(aes(subset = to_show, fill = location), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    geom_tippoint(aes(fill = location), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    ggtitle("True") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    theme(legend.position = "none") +
    coord_cartesian(clip="off")
  
  strict_plot <- ggtree(temp_sc, linewidth = 0.25) +
    geom_nodepoint(aes(fill = prob_gc), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    geom_tippoint(aes(fill = prob_gc), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    ggtitle("SC") +
    eo_node_fill_scale +
    theme(legend.position = "none") +
    coord_cartesian(clip="off")
  strict_plot$data$x = strict_plot$data$x - max(strict_plot$data$x) + 200
  strict_plot
  
  align_plots_by_xscale(true_plot, tyche_plot, strict_plot, extra = 50)
}




for (i in 1:20) {
  pdf(paste0("./tree_plots/gc_primary_trees_clone_", i, ".pdf"), width = 5, height = 2.7)
  print(create_plots(i))
  dev.off()
  cat("Plotted trees for clone:", i, "\n")
  
}

pdf(paste0("./tree_plots/gc_primary_trees_clone_", 2, ".pdf"), width = 5, height = 2.7)
print(create_plots(2))
dev.off()

