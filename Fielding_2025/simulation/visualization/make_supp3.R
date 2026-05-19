library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ape)
library(data.table)
library(airr)
library(dowser)
library(purrr)
library(treeio)


model_names_recall <- c(
  "Unconstrained"="TyCHE_EO_Est_Emp_no_constraints",
  "MSD"="TyCHE_EO_Est_Emp_only_max_start",
  "Root"="TyCHE_EO_Est_Emp_only_root_prior"
)


true_tree_height <- 1200
ignore_ess_params <- c("freqParameter", "traitfrequencies", "rateIndicator", "rootType", "rateCategories")

simulation_id = "sel"
folder_name = "5_12"
path_to_results <- paste0("./output_", folder_name, "/beast_output/", simulation_id, "/")


#####################################################################
# Load in the log.tsv data, record tree heights, record convergence #
#####################################################################
all_results <- list()

for (i in 1:length(model_names_recall)){
  all_results[[model_names_recall[i]]] <- data.frame(clone_id = 1:20, tree_height = NA, tree_height_prop_error = NA, converged = NA)
  
  for (j in 1:20){
    file_name <- paste0(path_to_results, model_names_recall[i], "_", j, "_log.tsv" )
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
        all_results[[model_names_recall[i]]][all_results[[model_names_recall[i]]]$clone_id == j, "converged"] <- "Yes"
      } else if(ess_fail_check == 0) {
        all_results[[model_names_recall[i]]][all_results[[model_names_recall[i]]]$clone_id == j, "converged"] <- "100"
      } else {
        all_results[[model_names_recall[i]]][all_results[[model_names_recall[i]]]$clone_id == j, "converged"] <- "No"
      }
      
      # Get the tree height
      tree_height <- log %>% filter(item == "TreeHeight") %>% pull(mean)
      all_results[[model_names_recall[i]]][all_results[[model_names_recall[i]]]$clone_id == j, "tree_height"] <- tree_height
      all_results[[model_names_recall[i]]][all_results[[model_names_recall[i]]]$clone_id == j, "tree_height_prop_error"] <- abs(tree_height - true_tree_height) / true_tree_height
      
    } else {
      all_results[[model_names_recall[i]]][all_results[[model_names_recall[i]]]$clone_id == j, "converged"] <- "No"
    }
  }
  
  cat("Finished processing model:", model_names_recall[i], "\n")
}


final_results <- bind_rows(lapply(names(all_results), function(model_name) {
  cbind(model = model_name, all_results[[model_name]])
})) %>%
  mutate(model = recode(model,
                        "TyCHE_EO_Est_Emp_no_constraints"="Unconstrained",
                        "TyCHE_EO_Est_Emp_only_max_start"="MSD",
                        "TyCHE_EO_Est_Emp_only_root_prior"="Root"))
final_results$model <- factor(final_results$model, levels = c("Unconstrained", "MSD", "Root"))


##################
# Plot the trees #
##################

all_trees <- list()

colors = RColorBrewer::brewer.pal(2,"Set1")
names(colors) = c("GC", "other")

for (i in 1:length(model_names_recall)){
  all_trees[[model_names_recall[i]]] <- list()
  
  for (j in 1:20){
    file_name <- paste0(path_to_results, model_names_recall[i], "_", j, "_tree_with_trait.tree" )
    if (file.exists(file_name)) {
      
      # Load in tree
      tree <- read.beast(file_name)
      
      # Store tree
      all_trees[[model_names_recall[i]]][[as.character(j)]] <- tree
      
      
    } else {
      all_trees[[model_names_recall[i]]][[as.character(j)]] <- NULL
    }
    
  }
  cat("Finished processing trees for model:", model_names_recall[i], "\n")
}


true_trees <- read.beast("/Volumes/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/gc_recall_sims/simble_sims_gc_reentry_12_18/all_simplified_time_trees.nex")
get_summary_stats <- function(clone_id_input, model_input) {
  row = final_results %>% filter(clone_id == clone_id_input, model == model_input)
  if (nrow(row) > 1) {
    cat("Warning: More than one row found for clone", clone_id_input, "model", model_input, "\n")
    print(row)
  }
  curr_tree <- all_trees[[model_names_recall[model]]][[as.character(clone_id_input)]]
  curr_true_phylo <- true_trees[[clone_id_input]]@phylo
  
  beast_tips <- gsub("_heavy$", "", curr_tree@phylo$tip.label)
  curr_tree@phylo$tip.label <- beast_tips
  true_tips <- curr_true_phylo$tip.label
  common_tips <- intersect(beast_tips, true_tips)
  
  # Prune trees to common tips
  beast_keep <- curr_tree@phylo$tip.label[beast_tips %in% common_tips]
  true_keep <- curr_true_phylo$tip.label[true_tips %in% common_tips]
  
  curr_tree_phylo <- ape::keep.tip(curr_tree@phylo, beast_keep)
  true_pruned_phylo <- ape::keep.tip(curr_true_phylo, true_keep)
  
  rf_dist <- tryCatch({
    dowser::calcRF(di2multi(curr_tree_phylo, tol = 5),
                   di2multi(true_pruned_phylo, tol = 5))
  }, error = function(e) {
    cat("  Warning: RF calculation failed for clone", clone_id_input, "and model", model_input, ":", e$message, "\n")
    NA
  })
  curr_tree_length <- sum(curr_tree_phylo$edge.length)
  true_tree_length <- sum(true_pruned_phylo$edge.length)
  tree_length_prop_error <- abs(curr_tree_length - true_tree_length) / true_tree_length
  summary_row <- data.frame(
    clone_id = clone_id_input,
    model_display = model_input,
    tree_height_prop_error = row$tree_height_prop_error,
    tree_length_prop_error = tree_length_prop_error,
    rf_distance = rf_dist,
    converged = row$converged
  )
  return(summary_row)
}

all_summary <- data.frame(clone_id = integer(),
                          model_display = character(),
                          tree_height_prop_error = numeric(),
                          tree_length_prop_error = numeric(),
                          rf_distance = numeric(),
                          converged = character(),
                          stringsAsFactors = FALSE
)
for (model in levels(final_results$model)) {
  for (clone_id in 1:20) {
    summary_row <- get_summary_stats(clone_id, model)
    all_summary <- rbind(all_summary, summary_row)
  }
}

all_summary$model_display <- factor(all_summary$model_display, levels=c("Unconstrained", "MSD", "Root"))
all_summary$clone_id <- factor(all_summary$clone_id, levels=1:20)
all_summary$tree_height_prop_error <- abs(all_summary$tree_height_prop_error)
all_summary$tree_length_prop_error <- abs(all_summary$tree_length_prop_error)

saveRDS(all_summary, "./supp_4_summary.rds")

all_summary <- readRDS("./supp_4_summary.rds")

main_summary <- read.csv("/Volumes/HoehnK/Sherry/beast_workspace/TyCHE/figures/dec2025/main_analysis/plot_data/summary_data_main.csv")

main_summary <- main_summary %>%
  mutate(model_display = recode(model_display, "TyCHE\nest clocks"="TyCHE"))
main_summary$model_display <- factor(main_summary$model_display, levels=c("TyCHE", "SC", "UCLD"))
main_summary$tree_height_prop_error <- abs(main_summary$tree_height_prop_error)
main_summary$tree_length_prop_error <- abs(main_summary$tree_length_prop_error)

main_summary$converged <- "Yes"
relevant <- main_summary %>% filter(model_display == "TyCHE", simulation_name == "gc_reentry_12_18") %>% select(colnames(all_summary))
all_summary <- rbind(all_summary, relevant)
all_summary$model_display <- factor(all_summary$model_display, levels=c("Unconstrained", "MSD", "Root", "TyCHE"))

font_size <- theme(axis.text = element_text(size=7),
                   axis.title = element_text(size=9),
                   plot.title = element_text(size=9),
                   strip.text = element_text(size=8),
                   legend.text = element_text(size=7),
                   legend.title = element_text(size=8))

standard_additions <- function() {
  my_comparisons <- list( c("TyCHE", "SC"), c("TyCHE", "UCLD"))

  list(
    scale_fill_manual(values=c("Unconstrained"="#E41A1C", "MSD"="#377EB8", "Root"="#FFFF33", "TyCHE"="darkorange")),
  

    scale_shape_manual(values=c("Yes"=21, "100"=23, "No"=23)),

    theme_bw(),

    theme(legend.position = "none"),

    font_size
  )
}

point_size <- 2
height <- ggplot(all_summary, aes(x = model_display, y=tree_height_prop_error, fill = model_display)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(aes(shape=converged, group = model_display),  # control dodge grouping
              size = point_size,
              position = position_jitterdodge(
                jitter.height = 0,
                jitter.width = 0.5,
                dodge.width = 0.85)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  standard_additions() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(x="", y="Tree Height\nProportional Error")
height

length <- ggplot(all_summary, aes(x = model_display, y=tree_length_prop_error, fill = model_display)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(aes(shape=converged, group = model_display),
              size = point_size,
              position = position_jitterdodge(
                jitter.height = 0,
                jitter.width = 0.5,
                dodge.width = 0.85)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  standard_additions() +
  labs(x="", y="Tree Length\nProportional Error") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
length

angle_x_text <- theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"))
rf <- ggplot(all_summary, aes(x = model_display, y=rf_distance, fill = model_display)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(aes(shape=converged, group = model_display),
              size = point_size, position = position_jitterdodge(
                jitter.height = 0,
                jitter.width = 0.5,
                dodge.width = 0.85)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  standard_additions() +
  labs(y="RF distance") +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + font_size +
  angle_x_text
rf

combined_plot <- plot_grid(height, length, rf, ncol=1, align="v", rel_heights=c(0.92, 0.87, 1))
combined_plot
ggsave(combined_plot, width=4, height=6, filename="recall/supp3_plot.pdf")

all_summary$converged <- factor(all_summary$converged, levels=c("No", "100", "Yes"))
convergence_summary_plot <- ggplot(all_summary, aes(model_display, fill=converged)) +
  geom_bar() + theme_bw() +
  scale_fill_manual(values=c("Yes"="#2c7bb6", "100"="#fdae61", "No"="#d7191c")) +
  font_size

convergence_summary_plot

ggsave(convergence_summary_plot, width=2, height=3, filename="recall/supp3_convergence_plot.pdf")
