########################
## Hunter J. Melton   ##
## Jessie J. Fielding ##
###### 4/13/2026 #######

# Exploring output from gc_reentry simulations
# Makes tree plots for Fig 3B.
# Makes recall start plots for Fig 3E.

############################
# Load packages and set up #
############################

library(airr)
library(dplyr)
library(dowser)
library(ggtree)
library(purrr)
library(treeio)
library(data.table)

simulation_id = "sel"
folder_name = "4_7"

# read in an argument for which simulation to use
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No simulation argument provided.")
}

if (length(args) == 1) {
  stop("No folder argument provided.")
}
simulation_id <- args[1]
folder_name <- args[2]

if (simulation_id == "sel") {
  # Use selection simulation data
  bcrs <- read_rearrangement("/Volumes/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/gc_recall_sims/simble_sims_gc_reentry_12_18/all_samples_airr.tsv")
} else if (simulation_id == "un") {
  # Use neutral simulation data
  bcrs <- read_rearrangement("/Volumes/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/gc_recall_sims/simble_sims_gc_reentry_uniform_neutral_12_18/all_samples_airr.tsv")
} else {
  stop("Invalid simulation argument provided. Use 'selection' or 'neutral'.")
}


path_to_results <- paste0("./output_", folder_name, "/beast_output/", simulation_id, "/")

model_names <- c(
  "expectedOccupancy_EstTraitClockRates_EmpFreq",
  "strictClock_AncestralReconstruction_EmpFreq",
  "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq"
)

true_tree_height <- 1200
ignore_ess_params <- c("freqParameter", "traitfrequencies", "rateIndicator", "rootType", "rateCategories")

# Get the tips that were sampled after the first and second reactions
bcrs_heavy <- bcrs %>% 
  filter(locus == "IGH") %>% 
  mutate(celltype = ifelse(celltype == "default", "GC", "other")) %>%
  mutate(sample_time = as.numeric(sample_time))
first_rxn_tips <- list()
second_rxn_tips <- list()
for (i in 1:20){
  first_rxn_tips[[i]] <- bcrs_heavy %>% filter(sample_time %in% c(50,100) & clone_id == i) %>% pull(sequence_id)
  second_rxn_tips[[i]] <- bcrs_heavy %>% filter(sample_time %in% c(1150, 1200) & clone_id == i) %>% pull(sequence_id)
}

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

final_results <- bind_rows(lapply(names(all_results), function(model_name) {
  cbind(model = model_name, all_results[[model_name]])
})) %>%
  mutate(model = recode(model,
                        "expectedOccupancy_EstTraitClockRates_EmpFreq" = "EO Est",
                        "strictClock_AncestralReconstruction_EmpFreq" = "SC",
                        "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq" = "UCLD"))
final_results$model <- factor(final_results$model, levels = c("EO Est", "SC", "UCLD"))


##################
# Plot the trees #
##################

all_trees <- list()
branches_list <- list()

colors = RColorBrewer::brewer.pal(2,"Set1")
names(colors) = c("GC", "other")

for (i in 1:length(model_names)){
  all_trees[[model_names[i]]] <- list()
  branches_list[[model_names[i]]] <- list()
  
  for (j in 1:20){
    file_name <- paste0(path_to_results, model_names[i], "_", j, "_tree_with_trait.tree" )
    if (file.exists(file_name) & all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "converged"] == "Yes") {
      
      # Load in tree
      tree <- read.beast(file_name)
      
      # Get correct relative times for first and second rxn
      tree_height <- all_results[[model_names[i]]][all_results[[model_names[i]]]$clone_id == j, "tree_height"]
      first_rxn_tree_time <- tree_height - 1200
      second_rxn_tree_time <- tree_height - 100
      
      # Store tree
      all_trees[[model_names[i]]][[as.character(j)]] <- tree
      
      # Get branches dataframe for later
      tree_data <- tree@data %>% select(node, height, celltype) %>%
        mutate(height = as.numeric(height),
               time = max(height) - height) %>%
        as.data.frame()
      rownames(tree_data) <- tree_data$node
      
      branches <- as_tibble(tree@phylo$edge) %>%
        rename(parent = V1, child = V2) %>%
        mutate(
          start_time = tree_data[as.character(parent), "time"],
          end_time = tree_data[as.character(child), "time"],
          type_parent = tree_data[as.character(parent), "celltype"],
          type_child = tree_data[as.character(child), "celltype"],
          type_branch = if_else(type_parent == type_child, type_parent, "mixed"),
        ) %>% 
        arrange(start_time, end_time) 
      
      # Adjust the start and end times to be relative to the first vaccine
      branches$start_time <- branches$start_time - first_rxn_tree_time
      branches$end_time <- branches$end_time - first_rxn_tree_time
      
      branches_list[[model_names[i]]][[as.character(j)]] <- branches
      
      
    } else {
      all_trees[[model_names[i]]][[as.character(j)]] <- NULL
      branches_list[[model_names[i]]][[as.character(j)]] <- NULL
    }
    
  }
  cat("Finished processing trees for model:", model_names[i], "\n")
}



true_simple_trees <- read.beast("/Volumes/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/gc_recall_sims/simble_sims_gc_reentry_12_18/all_simplified_time_trees.nex")
true_full_trees <- read.beast("/Volumes/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/gc_recall_sims/simble_sims_gc_reentry_12_18/all_pruned_time_trees.nex")

source("./plot_extras.R")

node_size <- 1.0

create_plots <- function(clone_id){
  temp_sc <- all_trees[["strictClock_AncestralReconstruction_EmpFreq"]][[clone_id]]
  temp_sc@data$prob_gc <- NA
  for (i in 1:nrow(temp_sc@data)) {
    celltype_set <- temp_sc@data$celltype.set[i][[1]]
    celltype_set_prob <- temp_sc@data$celltype.set.prob[i][[1]]
    names(celltype_set_prob) <- celltype_set
    temp_sc@data$prob_gc[i] <- ifelse("GC" %in% celltype_set, celltype_set_prob["GC"], 0.0) %>% as.numeric()
  }
  temp_tyche <- all_trees[["expectedOccupancy_EstTraitClockRates_EmpFreq"]][[clone_id]]
  temp_tyche@data$prob_gc <- NA
  for (i in 1:nrow(temp_tyche@data)) {
    celltype_set <- temp_tyche@data$celltype.set[i][[1]]
    celltype_set_prob <- temp_tyche@data$celltype.set.prob[i][[1]]
    names(celltype_set_prob) <- celltype_set
    temp_tyche@data$prob_gc[i] <- ifelse("GC" %in% celltype_set, celltype_set_prob["GC"], 0.0) %>% as.numeric()
  }
  temp_simple_true <- true_simple_trees[[clone_id]]
  simple_ids <- temp_simple_true@data$cell_id
  temp_full_true <- true_full_trees[[clone_id]]
  
  temp_full_true@data$to_show <- FALSE
  temp_full_true@data$to_show[temp_full_true@data$cell_id %in% simple_ids] <- TRUE
  
  tyche_plot <- ggtree(temp_tyche, aes(color = as.numeric(expectedOccupancies)), linewidth = 0.3) +
    geom_nodepoint(aes(fill = prob_gc), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    geom_tippoint(aes(fill = prob_gc), pch = 21, size = node_size, stroke = 0.5, color = "black") +
    ggtitle("TyCHE") +
    eo_node_fill_scale +
    eo_line_color_scale +
    theme(legend.position = "none") +
    coord_cartesian(clip="off")
  tyche_plot$data$x = tyche_plot$data$x - max(tyche_plot$data$x) + 1200
  
  names(colors) = c("germinal_center", "other", "Ambig.")
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
  strict_plot$data$x = strict_plot$data$x - max(strict_plot$data$x) + 1200
  strict_plot
  
  align_plots_by_xscale(true_plot, tyche_plot, strict_plot)
}




for (i in 1:20) {
  pdf(paste0("./tree_plots/gc_recall_trees_clone_", i, ".pdf"), width = 5, height = 2.7)
  print(create_plots(i))
  dev.off()
  cat("Plotted trees for clone:", i, "\n")

}

####################################
# Inferred GC reaction start times #
####################################

getDiffPoint = function(tree, node){
  type = filter(tree@data, !!node==node)$celltype
  edge = tree@phylo$edge[tree@phylo$edge[,2] == node,]
  if(length(edge) == 0){
    return(tibble(diffnode=node, type="root", height=filter(tree@data, !!node==node)$height))
  }
  if(!is.null(nrow(edge))){
    stop("weird")
  }
  parent = as.character(edge[1])
  parent_type = filter(tree@data, node==parent)$celltype
  parent_height = filter(tree@data, node==parent)$celltype
  child = as.character(edge[2])
  child_type = filter(tree@data, node==child)$celltype
  child_height = filter(tree@data, node==child)$height
  
  if(parent_type == type){
    return(getDiffPoint(tree, parent))
  }else{
    return(tibble(diffnode=child, type=child_type, height=child_height))
  }
}

getDiffPoints = function(tree){
  diffpoints = tibble()
  for(l in tree@phylo$tip.label){
    #print(l)
    d = filter(tree@data, node == which(tree@phylo$tip.label == l))
    df = getDiffPoint(tree, which(tree@phylo$tip.label == l))
    temp = tibble(tip=l, tip_type=d$celltype, tip_height=d$height, tip_node=d$node)
    diffpoints = bind_rows(diffpoints, bind_cols(temp, df))
  }
  diffpoints$height = as.numeric(diffpoints$height)
  diffpoints$tip_height = as.numeric(diffpoints$tip_height)
  diffpoints
}


alldiffs = tibble()
meandiffs = tibble()
for(id in model_names){
  id_trees = all_trees[[id]]
  cat(length(id_trees), " trees\n")
  
  diffs = bind_rows(Map(function(x, nm) {
    d = getDiffPoints(x)
    d$clone_id = nm
    d
  }, id_trees, names(id_trees))) %>% 
    filter(tip_type == "GC")
  
  diffs$difftime <- NA
  diffs$sampletime <- NA
  diffs$rxn <- NA
  for (i in 1:nrow(diffs)){
    branches <- branches_list[[id]][[diffs$clone_id[i]]]
    diffs$difftime[i] <- branches[branches$parent == diffs$diffnode[i], ]$start_time[1] 
    diffs$rxn[i] <- ifelse(diffs$tip[i] %in% first_rxn_tips[[as.numeric(diffs$clone_id[i])]], "first", NA)
    diffs$rxn[i] <- ifelse(diffs$tip[i] %in% second_rxn_tips[[as.numeric(diffs$clone_id[i])]], "second", diffs$rxn[i])
  }
  diffs <- diffs %>% filter(!(diffnode == tip_node)) # filter out any GC B cell tips that came directly from a memory B cell
  diffs$model <- id
  alldiffs <- bind_rows(alldiffs, diffs)
  cat("Finished processing model:", id, "\n")
}

alldiffs <- alldiffs %>%
  mutate(first_rxn = 0) %>%
  mutate(second_rxn = 1100) %>%
  mutate(model = recode(model, 
                        "expectedOccupancy_EstTraitClockRates_EmpFreq" = "EO Est",
                        "strictClock_AncestralReconstruction_EmpFreq" = "SC",
                        "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq" = "UCLD"))

# Get all of the unique GC explosion times for all clones and models
gc_explosions <- alldiffs %>% 
  select(clone_id, model, diffnode, type, difftime, first_rxn, second_rxn, rxn) %>%
  filter(type == "GC" | type == "root") %>%
  distinct() %>%
  mutate(model = factor(model)) 
# If there isn't a reaction time for the model clone combo, it means that the model predicts one chronic GC reaction
# So we'll add a line with the germline because that's the GC MRCA
long_model_names <- c("expectedOccupancy_EstTraitClockRates_EmpFreq",
                      "strictClock_AncestralReconstruction_EmpFreq", "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq")
names(long_model_names) <- c("EO Est", "SC", "UCLD")
for (id in gc_explosions$model %>% unique()){
  for (cloneid in 1:20){
    filtered_gc_explosions <- gc_explosions %>% filter(model == id, clone_id == cloneid)
    if (nrow(filtered_gc_explosions) == 0){
      if (!(cloneid) %in% names(branches_list[[long_model_names[id]]])){
        cat("Skipping", cloneid, "in", long_model_names[id], "\n")
        next
      }
      branches <- branches_list[[long_model_names[id]]][[cloneid]]
      gc_rxn_start_time <- branches[branches$type_branch == "GC", ]$start_time %>% min()
      gc_explosions <- bind_rows(gc_explosions, tibble(clone_id = as.character(cloneid),
                                                       model = id,
                                                       diffnode = "germline",
                                                       type = "GC",
                                                       difftime = gc_rxn_start_time,
                                                       first_rxn = 0,
                                                       second_rxn = 1100,
                                                       rxn = c("first", "second")))
    }
    
  }
}


# Set the model names for consistency
gc_explosions <- gc_explosions %>%
  mutate(figure_names = recode(model,
                               "EO Est" = "TyCHE",
                               "SC" = "SC",
                               "UCLD" = "UCLD")) %>%
  mutate(figure_names = factor(figure_names, levels = c("TyCHE", "SC", "UCLD")))

# Set the colors to match the other simulation metric plots
# colors <- RColorBrewer::brewer.pal(n = 3, name = "Pastel1")
colors <- c("darkorange", "green", "purple")
names(colors) <- c("TyCHE", "SC", "UCLD")


font_size <- theme(axis.text = element_text(size=7),
                   axis.title = element_text(size=9),
                   plot.title = element_text(size=9),
                   strip.text = element_text(size=8),
                   legend.text = element_text(size=7),
                   legend.title = element_text(size=8))
angle_x_text <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

my_comparisons <- list( c("TyCHE", "SC"), c("TyCHE", "UCLD"))
recall_start_plot <- ggplot(filter(gc_explosions, rxn == "second"), aes(y = difftime, x = figure_names, fill = figure_names)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(aes(group = interaction(model, rxn), fill=figure_names), shape=21,  # control dodge grouping
              size = 1,
              position = position_jitterdodge(
                jitter.height = 0,
                jitter.width = 0.5,
                dodge.width = 0.85)) +
  geom_hline(yintercept = 1100, linetype = "dashed", color = "black") +
  labs(y = "Recall Start Date", x=NULL, fill = "Reaction") +
  theme_bw() +
  font_size + angle_x_text + 
  scale_fill_manual(values = colors) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size=8/.pt)
recall_start_plot
pdf("./gc_recall_times.pdf", width = 1.5, height = 1.9)
print(recall_start_plot)
dev.off()
