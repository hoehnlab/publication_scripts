######################
## Hunter J. Melton ##
###### 4/20/2026 #####


# Supplementary flu figures for the TyCHE manuscript
library(airr)
library(dplyr)
library(purrr)
library(ggpubr)
library(dowser)
library(ggtree)
library(tidyr)
library(stringr)
library(viridis)
library(tidytree)




# Times of vaccination
# Year 2 refers to after the second vaccine, p04 was week 35 (7*35 = 245), p05 was week 38 (7*38 = 266), p11 was week 17 (7*17 = 119)
vax_times <- matrix(c(0, 245, 0, 266), byrow = TRUE, ncol = 2, dimnames = list(c("p04", "p05"), c("first", "second")))


# Figure out counts of each type at each time point 
comp_at_time <- function(t, branches, occupancy = FALSE) {
    active <- branches %>%
        filter(start_time <= t & end_time >= t)
    if (occupancy) {
        out <- tibble(
            time = t,
            GC = sum(active$occupancy),
            other = sum(1 - active$occupancy)
        )
    } else {
        out <- tibble(
            time = t,
            GC = sum(active$type_branch == "GC") + 0.5 * sum(active$type_branch == "mixed"),
            other = sum(active$type_branch == "Other") + 0.5 * sum(active$type_branch == "mixed")
        )
    }
    return(out)
}


# Utility function to construct occupancy proportions on each branch
construct_occupancy_prop_GC <- function(branches, alpha, beta) {
    branches <- branches %>%
        mutate(occupancy = mapply(get_occupancy_prop_GC, type_parent, type_child, (end_time - start_time), 
                                  alpha, beta))
    return(branches)
}


# Read in clones
f_clones <- readRDS("./clones_not_bulk_downsampled.rds")


# Read in the output and attach to clones
ids <- c("EO_Est_", "SC_", "UCLD_")


beast_path <- "~/software/beast/bin/"
trees <- list()
for (id in ids) {
    trees[[id]] <- f_clones
    trees[[id]] <- readBEAST(f_clones, beast_path, dir = "./flu_output_4_13/beast_output/", trait = "gex_annotation", id = id, nproc = 4)
}


#########################################################
# Figures (a), (b) - trees and GC proportion line plots #
#########################################################
# Settings for the plot
shapes = c("BM" = 24, "LN" = 22, "PBMC" = 23, "unsampled" = 21)
colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC", "Other", "Ambig.")
eo_line_color_scale <- scale_color_gradient(
    low = colors["Other"],
    high = colors["GC"],
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)
)
eo_node_fill_scale <- scale_fill_gradient(
    low = colors["Other"],
    high = colors["GC"],
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)
)


tl_mrcas <- list()
tl_plots <- list()
tl_data <- list()
logs <- list()
time_prop_lines <- list()
branches_list <- list()
for (id in ids) {
	print(id)
	tl_plots[[id]] <- list()
	tl_data[[id]] <- list()
	logs[[id]] <- tibble()
	tl_mrcas[[id]] <- tibble()
    time_prop_lines[[id]] <- list()
    branches_list[[id]] <- list()


	id_trees <- trees[[id]]
    for (cloneid in id_trees$clone_id) {
		print(cloneid)
		
		t <- filter(id_trees, clone_id == cloneid)$trees[[1]]
        clone_data <- filter(id_trees, clone_id == cloneid)$data[[1]]@data
		edges <- t@phylo$edge
		data <- t@data
		node <- ape::getMRCA(t@phylo, tip=t@phylo$tip.label) # mrca (of all tips, i.e., root) cell type


        # Append the root cell type
		tl_mrcas[[id]] <- bind_rows(tl_mrcas[[id]], tibble(clone_id=cloneid, 
			gex_annotation=filter(t@data, !!node==node)$gex_annotation,
			height=filter(t@data, !!node==node)$height))


        # Keep the tree
		tl_data[[id]][[cloneid]] <- t


        # Keep the dowser log object
		log <- filter(id_trees, clone_id == cloneid)$parameters[[1]]
		log$model <- id
		log$clone_id <- cloneid
		logs[[id]] <- bind_rows(logs[[id]], log)


        # Get the timepoints of the vaccines relative to the tree height
        # Treeheight is the time of the last sample, need to go back to timepoint = 0 for first vax time
        tree_height <- filter(log, item == "TreeHeight")$mean
        max_timepoint <- max(filter(id_trees, clone_id == cloneid)$data[[1]]@data$timepoint) # latest observed tip
        # Set the time between first and second vaccine depending on the donor (only clone 70058 comes from p04)
        # ggtree sets root time to 0, latest tip to tree_height, so we adjust vaccine times accordingly
        if (cloneid == "70058") {
            vax_time_1 <- vax_times["p04", "first"]
            vax_time_2 <- vax_times["p04", "second"]
        } else {
            vax_time_1 <- vax_times["p05", "first"]
            vax_time_2 <- vax_times["p05", "second"]
        }
        diff_timepoint_vax_second <- max_timepoint - vax_time_2
        diff_timepoint_vax_first <- max_timepoint - 0
        plot_timepoint_vax_second <- tree_height - diff_timepoint_vax_second
        plot_timepoint_vax_first <- tree_height - diff_timepoint_vax_first


        # Convert tree height into time, then put 0 at the time of the first vaccine
        t@data$raw_time <- 0 - as.numeric(t@data$height)
        t@data$time <- t@data$raw_time + diff_timepoint_vax_first
        
        # Convert location names into GC and Other and ambig.
        t@data$cell_type <- recode(t@data$gex_annotation, "germinal_center" = "GC", "other" = "Other", "germinal_center+other" = "Ambig.")


        # Note the sequences that bind to flu antigen
        fortified_t <- fortify(t)
        fortified_t$node <- as.character(fortified_t$node)
        t@data <- left_join(t@data, fortified_t[, c("node", "label", "isTip")], by = "node")
        t@data <- left_join(t@data, clone_data[, c("sequence_id", "is.agseq", "tissue")], by = c("label" = "sequence_id"))
        t@data$is.agseq <- ifelse(is.na(t@data$is.agseq), "Unknown", t@data$is.agseq)
        agseq_nodes <- t@data %>% filter(is.agseq == "Yes") %>% pull(node)
        t@data$tissue <- ifelse(is.na(t@data$tissue), "unsampled", t@data$tissue)


		if (id == "EO_Est_") {
             t@data$expectedOccupancies <- as.numeric(t@data$expectedOccupancies)
            
             t@data$prob_gc <- NA
             for (i in 1:nrow(t@data)) {
                 gex_annot_set <- t@data$gex_annotation.set[i][[1]]
                 gex_annot_set_prob <- t@data$gex_annotation.set.prob[i][[1]]
                 names(gex_annot_set_prob) <- gex_annot_set
                 t@data$prob_gc[i] <- ifelse("germinal_center" %in% gex_annot_set, gex_annot_set_prob["germinal_center"], 0) %>% as.numeric()
             }


            # Make an EO tree
            tl_plots[[id]][[cloneid]] <- ggtree(t, aes(color = expectedOccupancies),linewidth = 0.75) + 
                                            geom_nodepoint(aes(fill = prob_gc), pch = 21, size = 1.7, stroke = 0.5, color = "black") + 
                                            geom_tippoint(aes(fill = prob_gc, pch = tissue), size = 1.7, stroke = 0.5, color = "black") +
                                            ggtitle(paste0(cloneid)) +
                                            geom_tippoint(aes(subset = node %in% agseq_nodes, fill = prob_gc, pch = tissue), size = 1.5, stroke = 0.75, color = "black", show.legend = FALSE) + # Add a slightly thicker stroke to the flu antigen-binding sequences
                                            # scale_fill_manual(values = colors) +
                                            scale_shape_manual(values = shapes) +
                                            eo_line_color_scale +
                                            eo_node_fill_scale +
                                            geom_treescale(width = 365, y = -5) + 
                                            geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black")  +
                                            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") +
                                            theme(legend.position = "none") +
                                            coord_cartesian(clip="off")


        } else {
            tl_plots[[id]][[cloneid]] <- ggtree(t,linewidth = 0.4) + 
                                            geom_nodepoint(aes(fill = cell_type), pch = 21, size = 1.7, stroke = 0.5, color = "black") + 
                                            geom_tippoint(aes(fill = cell_type, pch = tissue), size = 1.7, stroke = 0.5, color = "black") +
                                            ggtitle(paste0(cloneid)) +
                                            geom_tippoint(aes(subset = node %in% agseq_nodes, fill = cell_type, pch = tissue), size = 1.5, stroke = 0.75, show.legend = FALSE) + # Add a slightly thicker stroke to the flu antigen-binding sequences
                                            scale_fill_manual(values = colors) +
                                            scale_shape_manual(values = shapes) +
                                            geom_treescale(width = 365, y = -5) + 
                                            geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black")  +
                                            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") +
                                            theme(legend.position = "none") +
                                            coord_cartesian(clip="off")
        }
        
        


        # Make a node data frame with occupancy column
        if ("expectedOccupancies" %in% colnames(t@data)) {
            tree_data <- t@data %>% select(node, height, cell_type, expectedOccupancies) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
        } else {
            tree_data <- t@data %>% select(node, height, cell_type) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
            tree_data$expectedOccupancies <- NA  # If no occupancies, set to NA
        }
        
        # Make a data frame for the branches
        # ape has height of branches starting at 0 at the tip, but we want time so we need the root to be 0
        branches <- as_tibble(t@phylo$edge) %>%
        rename(parent = V1, child = V2) %>%
        mutate(
            start_time = tree_data[as.character(parent), "time"],
            end_time = tree_data[as.character(child), "time"],
            type_parent = tree_data[as.character(parent), "cell_type"],
            type_child = tree_data[as.character(child), "cell_type"],
            type_branch = if_else(type_parent == type_child, type_parent, "mixed"),
            occupancy = as.numeric(tree_data[as.character(child), "expectedOccupancies"])
        ) %>% 
        arrange(start_time, end_time) 
        branches_list[[id]][[cloneid]] <- branches


        # Discretize the time
        time_grid <- seq(from = min(branches$start_time), to = max(branches$end_time), length.out = 2000) 


        # Apply to all timepoints and get proportions of each cell type
        do_occupancy <- !all(is.na(tree_data$expectedOccupancies))
        time_comps <- map_dfr(time_grid, comp_at_time, branches = branches, occupancy = do_occupancy)


        # Set time = 0 to the first vaccine
        time_comps$time <- time_comps$time - plot_timepoint_vax_first
        plot_timepoint_vax_second <- plot_timepoint_vax_second - plot_timepoint_vax_first
        plot_timepoint_vax_first <- 0


        # Make long format for plotting
        time_comps_long <- time_comps %>%
            pivot_longer(cols = c(GC, other), names_to = "type", values_to = "count") 


        # Calculate proportions
        time_props_GC <- time_comps_long %>%
        group_by(time) %>%
        mutate(proportion = count / sum(count)) %>%
        filter(type == "GC") %>%
        ungroup()


        # Line plot with proportions
        time_prop_lines[[id]][[cloneid]] <- ggplot(time_props_GC, aes(x = time, y = proportion)) +
            geom_line(linewidth = 1, color = "#E41A1C") +
            labs(title = "",
                x = "", y = "") +
            theme_bw() +
            geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black") +
            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") + 
            theme(legend.position = "none") + 
            theme(axis.text = element_text(color = "black"), panel.border = element_rect(color = "black"), axis.ticks = element_line(color = "black")) 
	}
  
    # Save the trees
    pdf(paste0("./figures_4_20/other_figures_4_21/", id, "_time_trees.pdf"), width = 2.5, height = 5)
    for (cloneid in id_trees$clone_id) {
        print(tl_plots[[id]][[cloneid]])
    }
    dev.off()




    # Save the population plots
    pdf(paste0("./figures_4_20/other_figures_4_21/", id, "_time_populations.pdf"), width = 5, height = 2.5)
    for (cloneid in id_trees$clone_id) {
        p4 <- time_prop_lines[[id]][[cloneid]]
        print(p4)
    }
    dev.off()
}




#####################################
# Figure C - Bayesian skyline plots #
#####################################
mediansky_tyche <- getSkylines(trees$EO_Est_, 
   time = "timepoint",
   dir = "./flu_output_4_13/beast_output/",
   id = "EO_Est_",
   max_height = "median",
   exclude_germline = TRUE,
   nproc = 8
)


# set up the data
skyline_data <- list()
skyline_data <- mediansky_tyche$skyline
names(skyline_data) <- c(trees$EO_Est_$clone_id)
skyline_plots <- list()


for (i in 1:8){
    if (trees$EO_Est_$clone_id[i] == "70058") {
        vax_time_1 <- vax_times["p04", "first"]
        vax_time_2 <- vax_times["p04", "second"]
    } else {
        vax_time_1 <- vax_times["p05", "first"]
        vax_time_2 <- vax_times["p05", "second"]
    }


    skyline_plot <- ggplot(skyline_data[[i]], aes(x = bin)) + 
        geom_line(aes(y = median), color = "#730028", size = 1) +
        geom_line(aes(y = lci), color = "#730028", size = 0.5) +
        geom_line(aes(y = uci), color = "#730028", size = 0.5) +
        geom_vline(xintercept = vax_time_1, linetype = "dashed", color = "black") +
        geom_vline(xintercept = vax_time_2, linetype = "dashed", color = "black") +
        labs(title = names(skyline_data)[i], x = "", y = "Pop. Size") + 
        theme_bw() + 
        theme(axis.text = element_text(color = "black"), panel.border = element_rect(color = "black"), axis.ticks = element_line(color = "black")) +
        # scale_x_continuous(breaks = seq(-1500, 400, 500)) + 
        scale_y_log10()
    skyline_plots[[i]] <- skyline_plot        
}


pdf(paste0("./figures_4_20/other_figures_4_21/tyche_skylines.pdf"), width = 5, height = 2.5)
    for (i in 1:8) {
        print(skyline_plots[[i]])
    }
dev.off()
