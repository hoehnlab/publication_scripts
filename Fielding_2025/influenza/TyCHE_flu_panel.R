######################
## Hunter J. Melton ##
###### 4/20/2026 #####

# Figure 5 for the TyCHE manuscript

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
    cloneid <- "70058"
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


        # Get the cell type probabilities
        t@data$prob_gc <- NA
             for (i in 1:nrow(t@data)) {
                 gex_annot_set <- t@data$gex_annotation.set[i][[1]]
                 gex_annot_set_prob <- t@data$gex_annotation.set.prob[i][[1]]
                 names(gex_annot_set_prob) <- gex_annot_set
                 t@data$prob_gc[i] <- ifelse("germinal_center" %in% gex_annot_set, gex_annot_set_prob["germinal_center"], 0) %>% as.numeric()
             }


		if (id == "EO_Est_") {
             t@data$expectedOccupancies <- as.numeric(t@data$expectedOccupancies)




            # Make an EO tree
            tl_plots[[id]][[cloneid]] <- ggtree(t, aes(color = expectedOccupancies),linewidth = 0.75) + 
                                            geom_nodepoint(aes(fill = prob_gc), pch = 21, size = 1.7, stroke = 0.5, color = "black") + 
                                            geom_tippoint(aes(fill = prob_gc, pch = tissue), size = 1.7, stroke = 0.5, color = "black") +
                                            ggtitle(paste0(cloneid)) +
                                            geom_tippoint(aes(subset = node %in% agseq_nodes, fill = prob_gc, pch = tissue), size = 1.5, stroke = 0.75, color = "black", show.legend = FALSE) + # Add a slightly thicker stroke to the flu antigen-binding sequences
                                            eo_node_fill_scale +
                                            scale_shape_manual(values = shapes) +
                                            eo_line_color_scale +
                                            geom_treescale(width = 365, y = -5) + 
                                            geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black")  +
                                            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") +
                                            theme(legend.position = "none") +
                                            coord_cartesian(clip="off")


        } else {
            tl_plots[[id]][[cloneid]] <- ggtree(t,linewidth = 0.4) + 
                                            geom_nodepoint(aes(fill = prob_gc), pch = 21, size = 1.7, stroke = 0.5, color = "black") + 
                                            geom_tippoint(aes(fill = prob_gc, pch = tissue), size = 1.7, stroke = 0.5, color = "black") +
                                            ggtitle(paste0(cloneid)) +
                                            geom_tippoint(aes(subset = node %in% agseq_nodes, fill = prob_gc, pch = tissue), size = 1.5, stroke = 0.75, color = "black", show.legend = FALSE) + # Add a slightly thicker stroke to the flu antigen-binding sequences
                                            eo_node_fill_scale +
                                            scale_shape_manual(values = shapes) +
                                            eo_line_color_scale +
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
	# }
  
    # Save the trees
    pdf(paste0("./figures_4_20/", id, "_p04_p05_time_trees.pdf"), width = 2.5, height = 5)
    print(tl_plots[[id]][[cloneid]])
    dev.off()




    # Save the population plots
    pdf(paste0("./figures_4_20/", id, "_p04_p05_time_population.pdf"), width = 5, height = 2.5)
        p4 <- time_prop_lines[[id]][[cloneid]]
        print(p4)
    dev.off()
}


#######################################
# Figure (c) - Bayesian Skyline Plots #
#######################################


mediansky_tyche <- getSkylines(trees$EO_Est_[1,], 
   time = "timepoint",
   dir = "./flu_output_4_13/beast_output/",
   id = "EO_Est_",
   max_height = "median",
   exclude_germline = TRUE
)


mediansky_sc <- getSkylines(trees$SC_[1,], 
   time = "timepoint",
   dir = "./flu_output_4_13/beast_output/",
   id = "SC_",
   max_height = "median",
   exclude_germline = FALSE
)


mediansky_ucld <- getSkylines(trees$UCLD_[1,], 
   time = "timepoint",
   dir = "./flu_output_4_13/beast_output/",
   id = "UCLD_",
   max_height = "median",
   exclude_germline = FALSE
)




# set up the data
skyline_data <- list()
skyline_data[[1]] <- mediansky_tyche$skyline[[1]]
skyline_data[[2]] <- mediansky_sc$skyline[[1]]
skyline_data[[3]] <- mediansky_ucld$skyline[[1]]
names(skyline_data) <- c("TyCHE", "SC", "UCLD")


for (i in 1:3){
    vax_time_1 <- vax_times["p04", "first"]
    vax_time_2 <- vax_times["p04", "second"]


    skyline_plot <- ggplot(skyline_data[[i]], aes(x = bin)) + 
        geom_line(aes(y = median), color = "#730028", size = 1) +
        geom_line(aes(y = lci), color = "#730028", size = 0.5) +
        geom_line(aes(y = uci), color = "#730028", size = 0.5) +
        geom_vline(xintercept = vax_time_1, linetype = "dashed", color = "black") +
        geom_vline(xintercept = vax_time_2, linetype = "dashed", color = "black") +
        labs(title = "", x = "", y = "Pop. Size") + 
        theme_bw() + 
        theme(axis.text = element_text(color = "black"), panel.border = element_rect(color = "black"), axis.ticks = element_line(color = "black")) +
        scale_y_log10()
    pdf(paste0("./figures_4_20/", names(skyline_data)[i], "_skyline.pdf"), width = 5, height = 2.5)
    print(skyline_plot)
    dev.off()
}






################################################
# Figure D - Finding time of GC reaction start #
################################################
# So for each clone, we want to find the child of the most recent "other" ancestor of each GC tip, 


getDiffPoint2 = function(tree, node){
	type = filter(tree@data, !!node==node)$gex_annotation
	edge = tree@phylo$edge[tree@phylo$edge[,2] == node,]
	if(length(edge) == 0){
		return(tibble(diffnode=node, type="root", height=filter(tree@data, !!node==node)$height))
	}
	if(!is.null(nrow(edge))){
		stop("weird")
	}
    parent = as.character(edge[1])
	parent_type = filter(tree@data, node==parent)$gex_annotation
	parent_height = filter(tree@data, node==parent)$height
	child = as.character(edge[2])
	child_type = filter(tree@data, node==child)$gex_annotation
	child_height = filter(tree@data, node==child)$height


	if(parent_type == type){
		return(getDiffPoint2(tree, parent))
	}else{
		return(tibble(diffnode=child, type=child_type, height=child_height))
	}
}


getDiffPoints2 = function(tree){
	diffpoints = tibble()
	for(l in tree@phylo$tip.label){
		d = filter(tree@data, node == which(tree@phylo$tip.label == l))
		df = getDiffPoint2(tree, which(tree@phylo$tip.label == l))
		temp = tibble(tip=l, tip_type=d$gex_annotation, tip_height=d$height, tip_node=d$node)
		diffpoints = bind_rows(diffpoints, bind_cols(temp, df))
	}
	diffpoints$height = as.numeric(diffpoints$height)
	diffpoints$tip_height = as.numeric(diffpoints$tip_height)
	diffpoints
}




# Need branches list for each clone


# Get all of the root cell types
tl_data <- list()
branches_list <- list()
for (id in ids) {
	print(id)
	tl_data[[id]] <- list()
    branches_list[[id]] <- list()


 	id_trees <- trees[[id]]
	for (cloneid in id_trees$clone_id) {
		print(cloneid)
		
		t <- filter(id_trees, clone_id == cloneid)$trees[[1]]
        clone_data <- filter(id_trees, clone_id == cloneid)$data[[1]]@data
		edges <- t@phylo$edge
		data <- t@data
		node <- ape::getMRCA(t@phylo, tip=t@phylo$tip.label)


		tl_data[[id]][[cloneid]] <- t


        # Treeheight is the time of the last sample, need to go back to timepoint = 0 for first vax time
        tree_height <- max(as.numeric(t@data$height))
        max_timepoint <- max(filter(id_trees, clone_id == cloneid)$data[[1]]@data$timepoint)
        # Set the time between first and second vaccine depending on the donor (only clone 70058 comes from p04)
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


        # Note the sequences that bind to flu antigen
        fortified_t <- fortify(t)
        fortified_t$node <- as.character(fortified_t$node)
        t@data <- left_join(t@data, fortified_t[, c("node", "label", "isTip")], by = "node")
        t@data <- left_join(t@data, clone_data[, c("sequence_id", "is.agseq", "tissue")], by = c("label" = "sequence_id"))
        t@data$is.agseq <- ifelse(is.na(t@data$is.agseq), "Unknown", t@data$is.agseq)
        agseq_nodes <- t@data %>% filter(is.agseq == "Yes") %>% pull(node)
        t@data$tissue <- ifelse(is.na(t@data$tissue), "unsampled", t@data$tissue)


        if ("expectedOccupancies" %in% colnames(t@data)) {
            tree_data <- t@data %>% select(node, height, gex_annotation, expectedOccupancies) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
        } else {
            tree_data <- t@data %>% select(node, height, gex_annotation) %>%
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
            type_parent = tree_data[as.character(parent), "gex_annotation"],
            type_child = tree_data[as.character(child), "gex_annotation"],
            type_branch = if_else(type_parent == type_child, type_parent, "mixed"),
            occupancy = as.numeric(tree_data[as.character(child), "expectedOccupancies"])
        ) %>% 
        arrange(start_time, end_time) 


        # Adjust the start and end times to be relative to the first vaccine
        branches$start_time <- branches$start_time - plot_timepoint_vax_first
        branches$end_time <- branches$end_time - plot_timepoint_vax_first


        branches_list[[id]][[cloneid]] <- branches


    }
}




alldiffs = tibble()
meandiffs = tibble()
for(id in ids){
	id_trees = trees[[id]]


	diffs = bind_rows(lapply(id_trees$trees, function(x){
		d = getDiffPoints2(x)
		d$clone_id = x@info$name
		d
	})) %>% 
    filter(tip_type == "germinal_center") %>%
    filter(!(tip == "Germline")) # filter out the germline tip


    # Map GC ancestor time into teh same time space as the vaccines
    diffs$difftime <- NA
    diffs$sampletime <- NA
    diffs$sampleyear <- NA
    for (i in 1:nrow(diffs)){
        branches <- branches_list[[id]][[diffs$clone_id[i]]]
        diffs$difftime[i] <- branches[branches$parent == diffs$diffnode[i], ]$start_time[1] 
        diffs$sampletime[i] <- id_trees[id_trees$clone_id == diffs$clone_id[i], ]$data[[1]]@data %>% filter(sequence_id == diffs$tip[i]) %>% pull(timepoint)
        diffs$sampleyear[i] <- id_trees[id_trees$clone_id == diffs$clone_id[i], ]$data[[1]]@data %>% filter(sequence_id == diffs$tip[i]) %>% pull(year)
        cat(i)
    }
    diffs <- diffs %>% filter(!(diffnode == tip_node)) # filter out any GC B cells that came directly from a memory B cell
    diffs$model <- id
    alldiffs <- bind_rows(alldiffs, diffs)
}


alldiffs <- alldiffs %>%
    mutate(vax_first = ifelse(clone_id == "70058", vax_times["p04", "first"], vax_times["p05", "first"]),
           vax_second = ifelse(clone_id == "70058", vax_times["p04", "second"], vax_times["p05", "second"])) %>%
    mutate(model = recode(model, 
        "EO_Est_" = "TyCHE",
        "SC_" = "SC",
        "UCLD_" = "UCLD"))


# Get all of the unique inferred GC reaction start times for all clones and models
gc_explosions <- alldiffs %>% 
    select(clone_id, model, diffnode, type, difftime, vax_first, vax_second, sampleyear) %>%
    filter(type == "germinal_center") %>%
    distinct() %>%
    mutate(patient = ifelse(clone_id == "70058", "p04", "p05")) %>% 
    mutate(model = factor(model, levels = c("UCLD", "SC", "TyCHE"))) %>%
    mutate(year = ifelse(sampleyear == 1, "2018", "2019")) 


shapes = c("p04" = 23, "p05" = 21)
gc_explosion_plot <- ggplot(gc_explosions, aes(x = difftime, y = model, fill = year)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_jitter(aes(shape = patient, group = interaction(model, year)),  # control dodge grouping
                size = 2,     
                position = position_jitterdodge(
                    jitter.height = 0,
                    jitter.width = 0.25,
                    dodge.width = 0.85)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = vax_times["p05", "second"], linetype = "dashed", color = "black") +
    labs(x = "Inferred GC Reaction Start Date", y = "", fill = "Year", shape = "Patient") +
    theme_bw() +
    scale_fill_manual(values = c("2018" = "#FDB863", "2019" = "#B2ABD2")) +
    scale_shape_manual(values = shapes) +
    theme(legend.position = "none") +
    # theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(color = "black"), panel.border = element_rect(color = "black"), axis.ticks = element_line(color = "black")) 


pdf("./figures_4_20/gc_explosion_times.pdf", width = 4.5, height = 2.5)
print(gc_explosion_plot)
dev.off()


# Check proportion of inferred GC reactions that start after the vaccine time points
gc_explosions_check_start <- gc_explosions
gc_explosions_check_start$after_vax <- NA
for (i in 1:nrow(gc_explosions_check_start)){
    vax_time <- ifelse(gc_explosions_check_start$sampleyear[i] == 1, gc_explosions_check_start$vax_first[i], gc_explosions_check_start$vax_second[i])
    gc_explosions_check_start$after_vax[i] <- ifelse(gc_explosions_check_start$difftime[i] > vax_time, "Yes", "No")
}
gc_explosions_prop_after <- gc_explosions_check_start %>% 
                                group_by(model) %>%
                                summarize(prop_after = mean(after_vax == "Yes"))
#   model prop_after
#   <fct>      <dbl>
# 1 UCLD       0.5  
# 2 SC         0.421
# 3 TyCHE      0.792


# Check proportion of GCs sampled after 2019/2020 vax that are estimated to start before the 2018/2019 vacccine
gc_explosions_check_start2 <- gc_explosions %>% filter(sampleyear == 2)
gc_explosions_check_start2$before_first_vax <- ifelse(gc_explosions_check_start2$difftime < gc_explosions_check_start2$vax_first, "Yes", "No")
gc_explosions_prop_before <- gc_explosions_check_start2 %>%
                                group_by(model) %>%
                                summarize(prop_before = mean(before_first_vax == "Yes"))
#   model prop_before
#   <fct>       <dbl>
# 1 UCLD        0
# 2 SC          0.286
# 3 TyCHE       0  




############################################
# Figure F - Differentiation time analysis #
############################################


# Get the differentation time from GC for all different cell types
# Need to append cell types to the trees
p04 <- read_rearrangement("data/P04_BCR_03182022_airr.tsv") 
p05 <- read_rearrangement("data/P05_BCR_03182022_airr.tsv") 


trees_w_ct <- trees


for (id in ids) {
    print(id)
    id_trees <- trees_w_ct[[id]]
    for (cloneid in id_trees$clone_id) {
        print(cloneid)
        clone_data <- filter(id_trees, clone_id == cloneid)$data[[1]]@data
        if (cloneid == "70058") {
            rearrangement <- p04
        } else {
            rearrangement <- p05
        }
        clone_data <- clone_data %>% left_join(rearrangement %>% select(sequence_id, gex_annotation) %>% mutate(cell_type = gex_annotation) %>% select(sequence_id, cell_type), by = c("sequence_id"))
        clone_data <- clone_data %>% mutate(ct_use = case_when(
            cell_type == "GC" ~ "GC",
            cell_type == "ABC" & tissue == "LN" ~ "MBC_LN",
            cell_type == "RMB" & tissue == "LN" ~ "MBC_LN",
            cell_type == "PB" & tissue == "LN" ~ "PC_LN",
            cell_type == "ABC" & tissue == "PBMC" ~ "MBC_PBMC",
            cell_type == "RMB" & tissue == "PBMC" ~ "MBC_PBMC",
            cell_type == "PB" & tissue == "PBMC" ~ "PC_PBMC",
            cell_type == "ABC" & tissue == "BM" ~ "MBC_BM",
            cell_type == "RMB" & tissue == "BM" ~ "MBC_BM",
            cell_type == "PB" & tissue == "BM" ~ "PC_BM"
        )) 
        trees_w_ct[[id]][trees_w_ct[[id]]$clone_id == cloneid, ]$data[[1]]@data <- clone_data
        cat(table(clone_data$ct_use))
    }
}


# Get the differentiation time from GC for all different cell types
diff_times <- getDiffPoints(trees_w_ct$EO_Est_, trait = "gex_annotation", tip_traits = c("ct_use", "cell_type", "tissue"), height = "height", eo_adjust = TRUE, eo_type = "germinal_center")
diff_times_other <- diff_times %>% filter(tip_type == "other", !(is.na(ct_use))) %>%
    mutate(raw_time = NA, time = NA, tip_raw_time = NA, tip_time = NA, normalized_time = NA) 


# Convert tree height to time
clones <- c("70058", "26299", "121056", "129238", "19989", "111394", "101085", "43876")
tree_heights <- lapply(trees, function(x) {
    x$parameters %>% bind_rows(.id = "clone_id") %>% mutate(clone_id = clones[as.numeric(clone_id)])
}) %>%
    bind_rows(.id = "model") %>%
    filter(item %in% c("TreeHeight") & model == "EO_Est_") %>%
    select(model, clone_id, mean) %>%
    rename(tree_height = mean)






for (i in 1:nrow(diff_times_other)){
    cloneid <- diff_times_other$clone_id[i]
    tree_height <- tree_heights %>% filter(clone_id == cloneid) %>% pull(tree_height)
    if (cloneid == "70058") {
        vax_time_1 <- vax_times["p04", "first"]
        vax_time_2 <- vax_times["p04", "second"]
    } else {
        vax_time_1 <- vax_times["p05", "first"]
        vax_time_2 <- vax_times["p05", "second"]
    }


    max_timepoint <- max(filter(trees_w_ct$EO_Est_, clone_id == cloneid)$data[[1]]@data$timepoint) # latest observed tip
    # Set the time between first and second vaccine depending on the donor (only clone 70058 comes from p04)
    # ggtree sets root time to 0, latest tip to tree_height, so we adjust vaccine times accordingly
    
    diff_timepoint_vax_second <- max_timepoint - vax_time_2
    diff_timepoint_vax_first <- max_timepoint - 0
    plot_timepoint_vax_second <- tree_height - diff_timepoint_vax_second
    plot_timepoint_vax_first <- tree_height - diff_timepoint_vax_first


    # Convert tree height into time, then put 0 at the time of the first vaccine
    diff_times_other$raw_time[i] <- 0 - diff_times_other$node_height[i]
    diff_times_other$time[i] <- diff_times_other$raw_time[i] + diff_timepoint_vax_first
    diff_times_other$normalized_time[i] <- 1 - (diff_times_other$node_height[i] / tree_height)
    diff_times_other$tip_raw_time[i] <- 0 - diff_times_other$tip_height[i]
    diff_times_other$tip_time[i] <- diff_times_other$tip_raw_time[i] + diff_timepoint_vax_first
}


# Now make boxplots like the gc_explosion plots but with each category of ct_use on the y-axis and the differentiation time on the x-axis
diff_times_other <- diff_times_other %>% 
    mutate(patient = ifelse(clone_id == "70058", "p04", "p05")) %>%
    mutate(ct_use = factor(ct_use, levels = c("MBC_LN", "PC_LN", "MBC_PBMC", "PC_PBMC", "MBC_BM", "PC_BM"))) %>%
    mutate(tip_epoch = case_when(
        tip_time < 0 ~ "Pre-Vax",
        tip_time >= 0 & patient == "p05" & tip_time < vax_times["p05", "second"] ~ "2018",
        tip_time >= 0 & patient == "p04" & tip_time < vax_times["p04", "second"] ~ "2018",
        tip_time >= vax_times["p04", "second"] & patient == "p04" ~ "2019",
        tip_time >= vax_times["p05", "second"] & patient == "p05" ~ "2019"
    )) %>%
    mutate(tip_epoch = factor(tip_epoch, levels = c("Pre-Vax", "2018", "2019"))) %>%
    mutate(cell_type = ifelse(cell_type == "ABC", "MBC", ifelse(cell_type == "RMB", "MBC", ifelse(cell_type == "PB", "PC", cell_type)))) %>%
    mutate(cell_type = factor(cell_type, levels = c("MBC", "PC"))) %>% 
    mutate(tissue_plot = recode(tissue, "LN" = "Lymph Node", "PBMC" = "Blood", "BM" = "Bone Marrow")) %>%
    mutate(tissue_plot = factor(tissue_plot, levels = c("Lymph Node", "Blood", "Bone Marrow")))




shapes = c("p04" = 23, "p05" = 21)


diff_plot <- ggplot(rbind(diff_times_other ), aes(x = time, y = cell_type, group = cell_type)) +
    geom_boxplot(aes(fill = factor(cell_type, levels = c("MBC", "PC"))), position = position_dodge(width = 0.85, preserve = "single"), outlier.shape = NA, color = "black") +
    geom_jitter(aes(shape = patient, group = ct_use, fill = cell_type),  # control dodge grouping
                size = 2,    
                position = position_jitterdodge(
                    jitter.height = 0,
                    jitter.width = 0.35,
                    dodge.width = 0.85)) +
    scale_fill_manual(values = c("MBC" = "#377EB8", "PC" = "#4DAF4A", "Pre-Vax" = "black", "2018" = "#FDB863", "2019" = "#B2ABD2")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = vax_times["p05", "second"], linetype = "dashed", color = "black") +
    labs(x = "Est. date of differentiation from GC (days)", y = "", shape = "Patient") +
    theme_bw() +
    scale_shape_manual(values = shapes) +
    facet_wrap(~tissue_plot, ncol = 1) +
    theme(legend.position = "none") +
    # theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
    theme(axis.text = element_text(color = "black"), panel.border = element_rect(color = "black"), axis.ticks = element_line(color = "black")) 


pdf("./figures_4_20/diff_times.pdf", width = 4.5, height = 3)
print(diff_plot)
dev.off()




#########################################################################
# Figure E - Clone 123238 to show differentiation from recall reactions #
#########################################################################


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
id = ids[1]
	print(id)
	tl_plots[[id]] <- list()
	tl_data[[id]] <- list()
	logs[[id]] <- tibble()
	tl_mrcas[[id]] <- tibble()
    time_prop_lines[[id]] <- list()
    branches_list[[id]] <- list()


	id_trees <- trees_w_ct[[id]]
    cloneid <- "129238"
		print(cloneid)
		
		t <- filter(id_trees, clone_id == cloneid)$trees[[1]]
        clone_data <- filter(id_trees, clone_id == cloneid)$data[[1]]@data
        clone_data <- clone_data %>% mutate(ct = case_when(
            cell_type == "GC" ~ "GC",
            cell_type == "ABC" ~ "MBC",
            cell_type == "RMB" ~ "MBC",
            cell_type == "PB" ~ "PC"
        )) 


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
        t@data <- left_join(t@data, clone_data[, c("sequence_id", "is.agseq", "tissue", "ct")], by = c("label" = "sequence_id"))
        t@data$is.pc <- ifelse(t@data$ct == "PC", "Yes", "No")
        pc_nodes <- t@data %>% filter(is.pc == "Yes") %>% pull(node)
        t@data$tissue <- ifelse(is.na(t@data$tissue), "unsampled", t@data$tissue)


        # Get the cell type probabilities
        t@data$prob_gc <- NA
             for (i in 1:nrow(t@data)) {
                 gex_annot_set <- t@data$gex_annotation.set[i][[1]]
                 gex_annot_set_prob <- t@data$gex_annotation.set.prob[i][[1]]
                 names(gex_annot_set_prob) <- gex_annot_set
                 t@data$prob_gc[i] <- ifelse("germinal_center" %in% gex_annot_set, gex_annot_set_prob["germinal_center"], 0) %>% as.numeric()
             }


		
             t@data$expectedOccupancies <- as.numeric(t@data$expectedOccupancies)




            # Make an EO tree
            tl_plots[[id]][[cloneid]] <- ggtree(t, aes(color = expectedOccupancies),linewidth = 0.5) + 
                                            geom_nodepoint(aes(fill = prob_gc), pch = 21, size = 0.75, stroke = 0.25, color = "black") + 
                                            geom_tippoint(aes(fill = prob_gc, pch = tissue), size = 0.75, stroke = 0.25, color = "black") +
                                            ggtitle(paste0(cloneid)) +
                                            eo_node_fill_scale +
                                            scale_shape_manual(values = shapes) +
                                            eo_line_color_scale +
                                            geom_treescale(width = 365, y = -5) + 
                                            geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black")  +
                                            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") +
                                            theme(legend.position = "none") +
                                            coord_cartesian(clip="off")
        
        


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
        
       
  
    # Save the trees
    pdf(paste0("./figures_4_20/", cloneid, "_p04_p05_time_trees.pdf"), width = 1.5, height = 2)
    print(tl_plots[[id]][[cloneid]])
    dev.off()
