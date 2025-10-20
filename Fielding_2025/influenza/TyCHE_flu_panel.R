######################
## Hunter J. Melton ##
###### 8/9/2025 ######

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
            GC = sum(active$type_branch == "germinal_center") + 0.5 * sum(active$type_branch == "mixed"),
            other = sum(active$type_branch == "other") + 0.5 * sum(active$type_branch == "mixed")
        )
    }
    return(out)
}

# Get occupancy of each branch (same calculation as in TyCHE EO)
get_occupancy_prop_GC <- function(parent, child, time, rate_GC2other, rate_other2GC) {
    
    # Set up the Q matrix as BEAST does
    Q <- matrix(0, nrow = 2, ncol = 2)
    # Start with the rates
    Q[1, 2] <- rate_GC2other
    Q[2, 1] <- rate_other2GC
    Q[1, 1] <- -rate_GC2other
    Q[2, 2] <- -rate_other2GC
    # Multiply by the frequency of the types
    Q[1, 2] <- Q[1, 2] * 0.5
    Q[2, 1] <- Q[2, 1] * 0.5
    Q[1, 1] <- -Q[1, 2]
    Q[2, 2] <- -Q[2, 1]
    # Normalize the rate matrix to one expected substitution per unit time
    subst <- sum(diag(-0.5*Q))
    Q <- Q / subst
    alpha <- Q[1, 2]  # GC to other rate
    beta <- Q[2, 1]   # Other to GC rate


    k <- sum(alpha, beta)
    expmkt <- exp(-k * time)

    # Catch mixed nodes
    if (parent == "other+germinal_center" || parent == "germinal_center+other") {
        parent <- "mixed"
    }
    if (child == "other+germinal_center" || child == "germinal_center+other") {
        child <- "mixed"
    }

    if (parent == "germinal_center" && child == "germinal_center") {
        exp_occupancy <- 1/k * ((beta^2*time + ((2*alpha*beta)/(k))*(1 - expmkt) + alpha^2*time*expmkt) / (beta + alpha*expmkt))
    } else if (parent == "other" && child == "other") {
        exp_occupancy <- 1/k * ((alpha*beta*time - ((2*alpha*beta)/(k))*(1 - expmkt) + alpha*beta*time*expmkt) / (alpha + beta*expmkt))
    } else if (parent == "germinal_center" && child == "other") {
        exp_occupancy <- 1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k)
    } else if (parent == "other" && child == "germinal_center") {
        exp_occupancy <- 1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k)
    } # Catch any mixed nodes
    else if (parent == "mixed" && child == "mixed") 
    {
        exp_occupancy <- mean(1/k * ((beta^2*time + ((2*alpha*beta)/(k))*(1 - expmkt) + alpha^2*time*expmkt) / (beta + alpha*expmkt)), 
                              1/k * ((alpha*beta*time - ((2*alpha*beta)/(k))*(1 - expmkt) + alpha*beta*time*expmkt) / (alpha + beta*expmkt)),
                              1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k),
                              1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k))
    } else if (parent == "mixed" && child == "germinal_center") {
        exp_occupancy <- mean(1/k * ((beta^2*time + ((2*alpha*beta)/(k))*(1 - expmkt) + alpha^2*time*expmkt) / (beta + alpha*expmkt)),
                              1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k))
    } else if (parent == "germinal_center" && child == "mixed") {
        exp_occupancy <- mean(1/k * ((beta^2*time + ((2*alpha*beta)/(k))*(1 - expmkt) + alpha^2*time*expmkt) / (beta + alpha*expmkt)),
                              1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k))
    } else if (parent == "mixed" && child == "other") {
        exp_occupancy <- mean(1/k * ((alpha*beta*time - ((2*alpha*beta)/(k))*(1 - expmkt) + alpha*beta*time*expmkt) / (alpha + beta*expmkt)),
                              1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k))
    } else if (parent == "other" && child == "mixed") {
        exp_occupancy <- mean(1/k * ((alpha*beta*time - ((2*alpha*beta)/(k))*(1 - expmkt) + alpha*beta*time*expmkt) / (alpha + beta*expmkt)),
                              1/k * ( (beta*time - alpha*time*expmkt) / (1 - expmkt) + (alpha - beta) / k)) 
    } else {
        stop("Unknown parent-child type combination")
    }
    return(exp_occupancy / time)
}

# Utility function to construct occupancy proportions on each branch
construct_occupancy_prop_GC <- function(branches, alpha, beta) {
    branches <- branches %>%
        mutate(occupancy = mapply(get_occupancy_prop_GC, type_parent, type_child, (end_time - start_time), 
                                  alpha, beta))
    return(branches)
}

# Read in the output
ids <- c("expectedOccupancy_FixedTraitClockRates_EmpFreq", "expectedOccupancy_EstTraitClockRates_EmpFreq", "instantSwitch_EstTraitClockRates_EmpFreq", "strictClock_AncestralReconstruction_EmpFreq", "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq")

trees <- list()
for (id in ids) {
    trees[[id]] <- readRDS(paste0("./flu_beast_8_25/", id, "_p04_p05_type-linked_time_tree_no_bulk.rds"))
}

#########################################################
# Figures (a), (b) - trees and GC proportion line plots #
#########################################################
# Settings for the plot
shapes = c("BM" = 24, "LN" = 22, "PBMC" = 23, "unsampled" = 21)
colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC", "Other", "Ambig.")

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
			location=filter(t@data, !!node==node)$location,
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
        t@data$cell_type <- recode(t@data$location, "germinal_center" = "GC", "other" = "Other", "germinal_center+other" = "Ambig.")

        # Note the sequences that bind to flu antigen
        fortified_t <- fortify(t)
        fortified_t$node <- as.character(fortified_t$node)
        t@data <- left_join(t@data, fortified_t[, c("node", "label", "isTip")], by = "node")
        t@data <- left_join(t@data, clone_data[, c("sequence_id", "is.agseq", "tissue")], by = c("label" = "sequence_id"))
        t@data$is.agseq <- ifelse(is.na(t@data$is.agseq), "Unknown", t@data$is.agseq)
        agseq_nodes <- t@data %>% filter(is.agseq == "Yes") %>% pull(node)
        t@data$tissue <- ifelse(is.na(t@data$tissue), "unsampled", t@data$tissue)

		tl_plots[[id]][[cloneid]] <- ggtree(t, linewidth = 0.15) + 
			geom_nodepoint(aes(fill = cell_type), pch = 21, size = 1.5, stroke = 0.5) + 
			geom_tippoint(aes(fill = cell_type, pch = tissue), size = 1.5, stroke = 0.5) +
            ggtitle(paste0(cloneid)) +
            geom_tippoint(aes(subset = node %in% agseq_nodes, fill = cell_type, pch = tissue), size = 1.5, stroke = 0.75, show.legend = FALSE) + # Add a slightly thicker stroke to the flu antigen-binding sequences
            scale_fill_manual(values = colors) +
            scale_shape_manual(values = shapes) +
			geom_treescale(width = 365) + 
			geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black")  +
            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") +
            theme(legend.position = "bottom") 

        # Make a node data frame with occupancy column
        if ("occupancies" %in% colnames(t@data)) {
            tree_data <- t@data %>% select(node, height, location, occupancies) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
        } else {
            tree_data <- t@data %>% select(node, height, location) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
            tree_data$occupancies <- NA  # If no occupancies, set to NA
        }
        
        # Make a data frame for the branches
        # ape has height of branches starting at 0 at the tip, but we want time so we need the root to be 0
        branches <- as_tibble(t@phylo$edge) %>%
        rename(parent = V1, child = V2) %>%
        mutate(
            start_time = tree_data[as.character(parent), "time"],
            end_time = tree_data[as.character(child), "time"],
            type_parent = tree_data[as.character(parent), "location"],
            type_child = tree_data[as.character(child), "location"],
            type_branch = if_else(type_parent == type_child, type_parent, "mixed"),
            occupancy = as.numeric(tree_data[as.character(child), "occupancies"])
        ) %>% 
        arrange(start_time, end_time) 
        branches_list[[id]][[cloneid]] <- branches

        # If occupancy isn't available, we need to construct it
        if (sum(is.na(branches$occupancy)) > 0 & (startsWith(id, "expectedOccupancy") || startsWith(id, "TraitLinkedExpectedOccupancy"))) {
            # Get the alpha and beta values from the log
            traitClockRate <- filter(log, item == "typeSwitchClockRate")$mean
            rate_GC2other <- filter(log, item == "relativeGeoRates.type.1")$mean 
            rate_other2GC <- filter(log, item == "relativeGeoRates.type.2")$mean 
            branches <- construct_occupancy_prop_GC(branches, rate_GC2other, rate_other2GC)
        } else if (sum(is.na(branches$occupancy)) > 0 & !(startsWith(id, "expectedOccupancy") || startsWith(id, "TraitLinkedExpectedOccupancy"))) {
           branches$occupancy <- ifelse(branches$type_child == "germinal_center", 1, 0)
        }

        # Discretize the time
        time_grid <- seq(from = min(branches$end_time), to = max(branches$end_time), length.out = 500) 

        # Apply to all timepoints and get proportions of each cell type
        time_comps <- map_dfr(time_grid, comp_at_time, branches = branches, occupancy = TRUE)

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
            geom_line(size = 1, color = "#E41A1C") +
            labs(title = "",
                x = "", y = "GC Proportion") +
            theme_bw() +
            geom_vline(xintercept = plot_timepoint_vax_first, linetype = "dashed", color = "black") +
            geom_vline(xintercept = plot_timepoint_vax_second, linetype = "dashed", color = "black") + 
            theme(legend.position = "none")
	# }
  
    # Save the trees
    pdf(paste0("./TyCHE_figures_9_2/", id, "_p04_p05_time_trees.pdf"), width = 2.5, height = 5)
    print(tl_plots[[id]][[cloneid]])
    dev.off()


    # Save the population plots
    pdf(paste0("./TyCHE_figures_9_2/", id, "_p04_p05_time_population.pdf"), width = 5, height = 2.5)
        p4 <- time_prop_lines[[id]][[cloneid]]
        print(p4)
    dev.off()
}

#######################################
# Figure (c) - Bayesian Skyline Plots #
#######################################
colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC", "non-GC", "Ambig.")

mediansky_tyche <- getSkylines(trees$expectedOccupancy_EstTraitClockRates_EmpFreq[1,], 
   time = "timepoint",
   dir = "/dartfs-hpc/scratch/f007p0j/flu_beast_8_25/",
   id = "expectedOccupancy_EstTraitClockRates_EmpFreq",
   max_height = "median"
)

mediansky_sc <- getSkylines(trees$strictClock_AncestralReconstruction_EmpFreq[1,], 
   time = "timepoint",
   dir = "/dartfs-hpc/scratch/f007p0j/flu_beast_8_25/",
   id = "strictClock_AncestralReconstruction_EmpFreq",
   max_height = "median"
)

mediansky_ucld <- getSkylines(trees$UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq[1,], 
   time = "timepoint",
   dir = "/dartfs-hpc/scratch/f007p0j/flu_beast_8_25/",
   id = "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq",
   max_height = "median"
)


# Read in the data
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
        scale_x_continuous(breaks = seq(-1500, 400, 500)) + 
        scale_y_log10()
    pdf(paste0("./TyCHE_figures_9_2/", names(skyline_data)[i], "_skyline.pdf"), width = 5, height = 2.5)
    print(skyline_plot)
    dev.off()
}

################################
# Figure (e) - Root cell types #
################################

# Get all of the root cell types
tl_data <- list()
tl_mrcas <- list()
logs <- list()
branches_list <- list()
for (id in ids) {
	print(id)
	tl_data[[id]] <- list()
	logs[[id]] <- tibble()
	tl_mrcas[[id]] <- tibble()
    branches_list[[id]] <- list()

	id_trees <- trees[[id]]
	for (cloneid in id_trees$clone_id) {
		print(cloneid)
		
		t <- filter(id_trees, clone_id == cloneid)$trees[[1]]
        clone_data <- filter(id_trees, clone_id == cloneid)$data[[1]]@data
		edges <- t@phylo$edge
		data <- t@data
		node <- ape::getMRCA(t@phylo, tip=t@phylo$tip.label)
		

		tl_mrcas[[id]] <- bind_rows(tl_mrcas[[id]], tibble(clone_id=cloneid, 
			location=filter(t@data, !!node==node)$location,
			height=filter(t@data, !!node==node)$height))

		tl_data[[id]][[cloneid]] <- t

		log <- filter(id_trees, clone_id == cloneid)$parameters[[1]]
		log$model <- id
		log$clone_id <- cloneid
		logs[[id]] <- bind_rows(logs[[id]], log)

        # Treeheight is the time of the last sample, need to go back to timepoint = 0 for first vax time
        tree_height <- filter(log, item == "TreeHeight")$mean
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

        colors = RColorBrewer::brewer.pal(3,"Set1")
        names(colors) = c("GC", "non-GC", "Ambig.")

        # Note the sequences that bind to flu antigen
        fortified_t <- fortify(t)
        fortified_t$node <- as.character(fortified_t$node)
        t@data <- left_join(t@data, fortified_t[, c("node", "label", "isTip")], by = "node")
        t@data <- left_join(t@data, clone_data[, c("sequence_id", "is.agseq", "tissue")], by = c("label" = "sequence_id"))
        t@data$is.agseq <- ifelse(is.na(t@data$is.agseq), "Unknown", t@data$is.agseq)
        agseq_nodes <- t@data %>% filter(is.agseq == "Yes") %>% pull(node)
        t@data$tissue <- ifelse(is.na(t@data$tissue), "unsampled", t@data$tissue)

        if ("occupancies" %in% colnames(t@data)) {
            tree_data <- t@data %>% select(node, height, location, occupancies) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
        } else {
            tree_data <- t@data %>% select(node, height, location) %>%
            mutate(height = as.numeric(height),
                time = max(height) - height) %>%
            as.data.frame()
            rownames(tree_data) <- tree_data$node
            tree_data$occupancies <- NA  # If no occupancies, set to NA
        }
        
        # Make a data frame for the branches
        # ape has height of branches starting at 0 at the tip, but we want time so we need the root to be 0
        branches <- as_tibble(t@phylo$edge) %>%
        rename(parent = V1, child = V2) %>%
        mutate(
            start_time = tree_data[as.character(parent), "time"],
            end_time = tree_data[as.character(child), "time"],
            type_parent = tree_data[as.character(parent), "location"],
            type_child = tree_data[as.character(child), "location"],
            type_branch = if_else(type_parent == type_child, type_parent, "mixed"),
            occupancy = as.numeric(tree_data[as.character(child), "occupancies"])
        ) %>% 
        arrange(start_time, end_time) 

        # Adjust the start and end times to be relative to the first vaccine
        branches$start_time <- branches$start_time - plot_timepoint_vax_first
        branches$end_time <- branches$end_time - plot_timepoint_vax_first

        branches_list[[id]][[cloneid]] <- branches

    }
}


tl_mrcas_combined <- bind_rows(tl_mrcas, .id = "model")
tl_mrcas_plot_df <- tl_mrcas_combined %>%
  group_by(model, location) %>%
  summarise(freq = n()) %>%
  mutate(prop = freq / sum(freq)) %>%
  mutate(model_abbrev = recode(model, "expectedOccupancy_FixedTraitClockRates_EmpFreq" = "TyCHE EO Fixed",
                               "expectedOccupancy_EstTraitClockRates_EmpFreq" = "TyCHE",
                               "instantSwitch_EstTraitClockRates_EmpFreq" = "TyCHE IS Est",
                               "strictClock_AncestralReconstruction_EmpFreq" = "SC",
                               "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq" = "UCLD"),
        cell_type = recode(location, "germinal_center" = "GC", "other" = "Other", "other+germinal_center" = "Ambig.")) %>%
  filter(model_abbrev %in% c("TyCHE", "SC", "UCLD")) %>%
  mutate(model_abbrev = factor(model_abbrev, levels = c("TyCHE", "SC", "UCLD"))) %>%
  ungroup()

# Settings for the plots
colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC", "Other", "Ambig.")

mrcas_loc_plot <- ggplot(tl_mrcas_plot_df, aes(x = model_abbrev, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", col = "black") +
  labs(x = "", y = "Prop. Clones", fill = "Root\nCell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_fill_manual(values = colors) 
ggsave("./TyCHE_figures_9_2/mrcas_location_plot.pdf", plot = mrcas_loc_plot, width = 3.5, height = 3)


##################################
# Figure (f) - Type Switch Rates #
##################################
# Attach the clone_id to a column in the logs
clones <- c("70058", "26299", "121056", "129238", "19989", "111394", "101085", "43876")
models <- ids[1:2]
names(clones) <- 1:8
names(models) <- 1:2
params_EO <- trees[[1]]$parameters %>% bind_rows(.id = "clone_id") %>% mutate(clone_id = clones[as.numeric(clone_id)])

params <- lapply(trees, function(x) {
    x$parameters %>% bind_rows(.id = "clone_id") %>% mutate(clone_id = clones[as.numeric(clone_id)])
}) %>%
    bind_rows(.id = "model") %>%
    filter(item %in% c("relativeGeoRates.type.1", "relativeGeoRates.type.2", "relativeGeoRates.newTrait.1", "relativeGeoRates.newTrait.2")) %>%
    mutate(item = recode(item, 
        "relativeGeoRates.type.1" = "GC->Other",
        "relativeGeoRates.type.2" = "Other->GC",
        "relativeGeoRates.newTrait.1" = "GC->Other",
        "relativeGeoRates.newTrait.2" = "Other->GC")) %>%
    filter(model != "TraitLinkedExpectedOccupancy_EstimatedTraitClockRates_EmpFreq") %>%
    # Normalize the rates for each model
    group_by(model, clone_id) %>%
    mutate(normalized = mean / mean(mean)) %>%
    ungroup()


params <- params %>% 
    mutate(model = recode(model, 
        "expectedOccupancy_FixedTraitClockRates_EmpFreq" = "TyCHE EO Fixed",
        "expectedOccupancy_EstTraitClockRates_EmpFreq" = "TyCHE",
        "instantSwitch_EstTraitClockRates_EmpFreq" = "TyCHE IS Est",
        "strictClock_AncestralReconstruction_EmpFreq" = "SC",
        "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq" = "UCLD")) %>%
    filter(model %in% c("TyCHE", "SC", "UCLD")) %>%
    mutate(model = factor(model, levels = c("TyCHE", "SC", "UCLD"))) %>%
    mutate(paired_id = paste(model, clone_id, sep = "_"))

colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC->Other", "Other->GC", "Ambig.")
rates_plot <- ggplot(params, aes(x = item, y = mean, fill = item)) +
    geom_boxplot() +
    geom_point(color = "black", size = 0.5) +
    geom_line(aes(group = paired_id), size=0.5, color="black", alpha = 0.4) +
    facet_wrap(~model) +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p", size = 1.5, label.x = 1.8) + 
    labs(x = "", y = "Transition Rate") +
    theme_bw() +
    scale_fill_manual(values = colors) +
    ylim(0, 0.22) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", panel.spacing = unit(0.1, "lines"))
pdf("./TyCHE_figures_9_2/trait_clock_rates.pdf", width = 3, height = 3)
print(rates_plot)
dev.off()


################################################
# Figure E - Finding time of GC reaction start #
################################################
# So for each clone, we want to find the child of the most recent "other" ancestor of each GC tip, 

getDiffPoint = function(tree, node){
	type = filter(tree@data, !!node==node)$location
	edge = tree@phylo$edge[tree@phylo$edge[,2] == node,]
	if(length(edge) == 0){
		return(tibble(diffnode=node, type="root", height=filter(tree@data, !!node==node)$height))
	}
	if(!is.null(nrow(edge))){
		stop("weird")
	}
    parent = as.character(edge[1])
	parent_type = filter(tree@data, node==parent)$location
	parent_height = filter(tree@data, node==parent)$height
	child = as.character(edge[2])
	child_type = filter(tree@data, node==child)$location
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
		temp = tibble(tip=l, tip_type=d$location, tip_height=d$height, tip_node=d$node)
		diffpoints = bind_rows(diffpoints, bind_cols(temp, df))
	}
	diffpoints$height = as.numeric(diffpoints$height)
	diffpoints$tip_height = as.numeric(diffpoints$tip_height)
	diffpoints
}


alldiffs = tibble()
meandiffs = tibble()
for(id in ids){
	id_trees = trees[[id]]

	diffs = bind_rows(lapply(id_trees$trees, function(x){
		d = getDiffPoints(x)
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
        "expectedOccupancy_FixedTraitClockRates_EmpFreq" = "TyCHE EO Fixed",
        "expectedOccupancy_EstTraitClockRates_EmpFreq" = "TyCHE",
        "instantSwitch_EstTraitClockRates_EmpFreq" = "TyCHE IS Est",
        "strictClock_AncestralReconstruction_EmpFreq" = "SC",
        "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq" = "UCLD")) %>%
    filter(model %in% c("TyCHE", "SC", "UCLD"))

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
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(pch = patient), size = 2, width = 0, height = 0.25) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = vax_times["p05", "second"], linetype = "dashed", color = "black") +
    labs(x = "Inferred GC Reaction Start Date", y = "", fill = "Year", shape = "Patient") +
    theme_bw() +
    scale_fill_manual(values = c("2018" = "skyblue3", "2019" = "palegreen3")) +
    scale_shape_manual(values = shapes) +
    theme(legend.position = "right") + 
    theme(axis.text.y = element_text(angle = 45, hjust = 1))

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
    scale_fill_manual(values = c("2018" = "skyblue3", "2019" = "palegreen3")) +
    scale_shape_manual(values = shapes) +
    theme(legend.position = "right") +
    theme(axis.text.y = element_text(angle = 45, hjust = 1))

pdf("./TyCHE_figures_9_2/gc_explosion_times.pdf", width = 4.5, height = 2.5)
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
# 1 UCLD       0.462
# 2 SC         0.292
# 3 TyCHE      0.852

# Check proportion of GCs sampled after 2019/2020 vax that are estimated to start before the 2018/2019 vacccine
gc_explosions_check_start2 <- gc_explosions %>% filter(sampleyear == 2)
gc_explosions_check_start2$before_first_vax <- ifelse(gc_explosions_check_start2$difftime < gc_explosions_check_start2$vax_first, "Yes", "No")
gc_explosions_prop_before <- gc_explosions_check_start2 %>%
                                group_by(model) %>%
                                summarize(prop_before = mean(before_first_vax == "Yes"))
#   model prop_before
#   <fct>       <dbl>
# 1 UCLD        0.143
# 2 SC          0.375
# 3 TyCHE       0  

