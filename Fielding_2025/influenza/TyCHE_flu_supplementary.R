######################
## Hunter J. Melton ##
##### 08/12/2025 #####

# Supplementary Figure(s)?

library(airr)
library(dplyr)
library(purrr)
library(ggpubr)
library(dowser)
library(ggtree)
library(tidyr)
library(patchwork)
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

# Get occupancy at each branch
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
construct_occupancy_prop_GC <- function(branches, alpha, beta) {
    # Calculate the occupancy for each branch
    branches <- branches %>%
        mutate(occupancy = mapply(get_occupancy_prop_GC, type_parent, type_child, (end_time - start_time), 
                                  alpha, beta))
    return(branches)
}

# Read in the output
ids <- c("expectedOccupancy_EstTraitClockRates_EmpFreq")

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
time_pop_bars <- list()
time_pop_lines <- list()
time_prop_bars <- list()
time_prop_lines <- list()
branches_list <- list()
for (id in ids) {
	print(id)
	tl_plots[[id]] <- list()
	tl_data[[id]] <- list()
	logs[[id]] <- tibble()
	tl_mrcas[[id]] <- tibble()
    time_pop_bars[[id]] <- list()
    time_pop_lines[[id]] <- list()
    time_prop_bars[[id]] <- list()
    time_prop_lines[[id]] <- list()
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
        
        # Convert location names into GC and Other
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

        # Should we use the occupancy or the type of the branch?
        use_occupancy <- ifelse(sum(is.na(branches$occupancy)) == 0, TRUE, FALSE)
        cat("Using occupancy:", use_occupancy, "\n")

        # Apply to all timepoints
        time_comps <- map_dfr(time_grid, comp_at_time, branches = branches, occupancy = use_occupancy)

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
	
  
        # Save the plots
        pdf(paste0("./TyCHE_supplementary_figures/", id, "_p04_p05_time_trees_", cloneid, ".pdf"), width = 2.5, height = 5)
            print(tl_plots[[id]][[cloneid]])
        dev.off()

        pdf(paste0("./TyCHE_supplementary_figures/", id, "_p04_p05_time_population_", cloneid, ".pdf"), width = 5, height = 2.5)
            print(time_prop_lines[[id]][[cloneid]])
        dev.off()
    }

}

#######################################
# Figure (c) - Bayesian Skyline Plots #
#######################################
colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC", "non-GC", "Ambig.")

mediansky_tyche_supp <- getSkylines(trees$expectedOccupancy_EstTraitClockRates_EmpFreq[2:8,], 
   time = "timepoint",
   dir = "/dartfs-hpc/scratch/f007p0j/flu_beast_8_25/",
   id = "expectedOccupancy_EstTraitClockRates_EmpFreq",
   max_height = "median",
   nproc = 4
)

# Read in the skyline data
for (i in 1:nrow(mediansky_tyche_supp)){
    vax_time_1 <- vax_times["p05", "first"]
    vax_time_2 <- vax_times["p05", "second"]

    skyline_plot <- ggplot(mediansky_tyche_supp$skyline[[i]], aes(x = bin)) + 
        geom_line(aes(y = median), color = "#730028", size = 1) +
        geom_line(aes(y = lci), color = "#730028", size = 0.5) +
        geom_line(aes(y = uci), color = "#730028", size = 0.5) +
        geom_vline(xintercept = vax_time_1, linetype = "dashed", color = "black") +
        geom_vline(xintercept = vax_time_2, linetype = "dashed", color = "black") +
        labs(title = "", x = "", y = "Pop. Size") + 
        theme_bw() + 
        scale_x_continuous(breaks = seq(-1500, 400, 500)) + 
        scale_y_log10()
    pdf(paste0("./TyCHE_supplementary_figures/", mediansky_tyche_supp$clone_id[i], "_skyline.pdf"), width = 5, height = 2.5)
    print(skyline_plot)
    dev.off()
}

# How many sequences bind to flu in each clone - need for annotating in illustrator
clones_supp <- trees$expectedOccupancy_EstTraitClockRates_EmpFreq[2:8,]
for (i in 1:nrow(clones_supp)){
    cat(clones_supp$clone_id[i],"\n")
    num_binding <- sum(clones_supp$data[[i]]@data$is.agseq == "Yes")
    cat("Number of binding sequences:", num_binding, "\n\n\n")
}

