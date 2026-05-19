######################
## Hunter J. Melton ##
###### 4/13/2026 #####


# Analyzing recall GC reactions after influenza vaccination
# Data from McIntire, K. M. et al. Maturation of germinal center B cells after influenza virus vaccination in humans. Journal of Experimental Medicine 221, e20240668 (2024)


##########################
# Load packages and data #
##########################


library(airr)
library(dplyr)
library(dowser)
library(ggtree)
library(purrr)


# Get array ID to run different models in parallel
job.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


# Set the path to the xml-writer directory
xml_writer_path <- "./xml-writer/"


# Beast path
beast_path = "~/software/beast/bin/"


# Load clones
clones_not_bulk_downsampled <- readRDS("./data/clones_not_bulk_downsampled.rds")


# Prep input for TyCHE
p04_clones <- c("70058")
p05_clones <- c("26299", "121056", "129238", "19989", "111394", "101085", "43876")
clone_data <- data.frame(
	sample = rep(c("P04", "P05"), c(length(p04_clones), length(p05_clones))),
	clone_id = c("70058", "26299", "121056", "129238","19989", "111394", "101085", "43876"))
input_data <- tibble(
    model = c(
		"EO_Est",
		"SC",
		"UCLD"
	),
		template_path = c(
			paste0(xml_writer_path, "TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml"),
			paste0(xml_writer_path, "StrictClock/StrictClock_AncestralReconstruction_EmpFreq.xml"),
			paste0(xml_writer_path, "UCLD/UCRelaxedClock_AncestralReconstruction_EmpFreq.xml")
			),
	clone_data = list(clone_data),
	ignore = list(c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator", "rateCategories")),
	TYPE_SWITCH_CLOCK_RATE = 1/365, # 1 transition per year as init value when estimating
	RATE_INDICATORS = "1 1",
	# Expected occupancy inputs
	# RELATIVE RATE PRIORS
	TRANSITION_RATE_ALPHA_1 = c("0.1"),
	TRANSITION_RATE_BETA_1 = c("1.0"),
	TRANSITION_RATE_ALPHA_2 = c("0.1"),
	TRANSITION_RATE_BETA_2 = c("1.0"),
	# RELATIVE RATE INIT VALUES
	TRANSITION_RATE_1_INIT = c("1"),
	TRANSITION_RATE_2_INIT = c("1"),
	TRAIT_RATE_MEAN_1 = 4.9E-3/7, 
	TRAIT_RATE_MEAN_2 = "0.000001",
	TRAIT_RATE_SIGMA_1 = 0.01 * (4.9E-3/7),
	TRAIT_RATE_SIGMA_2 = c("0.001"), 
	KAPPA_PRIOR_M = "0.67",
	KAPPA_PRIOR_S = "0.2",
	# TYPE SWITCH RATE PRIORS
	TYPE_SWITCH_INIT = c("0.005"), #c("0.002"),
	TYPE_SWITCH_RATE_ALPHA = c("0.001"), #c("0.001"), #c("10"),
	TYPE_SWITCH_RATE_BETA = c("5"), #c("5"), #c("0.00005"),
	# Strict clock/UCLD inputs
	CLOCK_RATE_INIT = 4.9E-3/7,
	TRANSITION_RATE_ALPHA = "0.1",
	TRANSITION_RATE_BETA = "1.0",
	UCLD_SIGMA_INIT = "0.5",
	# MAX START DATE
	max_start_date = list(-1460), # 4 years before the first time point
	# ROOT PRIOR
	root_trait=list("germinal_center")
)


max_iter <- 10
trees <- clones_not_bulk_downsampled


##################
# Get Time Trees #
##################
trees <- getTimeTreesIterate(trees,
	iterations = max_iter,
	ignore = input_data$ignore[[job.id]],
	beast = beast_path, 
	trait = "gex_annotation",
	time = "timepoint",
	dir = "./flu_output_4_13/beast_output/", 
	id = paste0(input_data$model[job.id], "_"),
	template = input_data$template_path[job.id],
	nproc = 8,
	include_germline = TRUE,
	mcmc_length = 500000000,
    log_every = "auto",
	log_target = 2000,
	# max_start_date = -1460, # 4 years before the first time point
	RATE_INDICATORS = input_data$RATE_INDICATORS[job.id],
	TRANSITION_RATE_ALPHA_1 = input_data$TRANSITION_RATE_ALPHA_1[job.id],
	TRANSITION_RATE_BETA_1 = input_data$TRANSITION_RATE_BETA_1[job.id],
	TRANSITION_RATE_ALPHA_2 = input_data$TRANSITION_RATE_ALPHA_2[job.id],
	TRANSITION_RATE_BETA_2 = input_data$TRANSITION_RATE_BETA_2[job.id],
	TRANSITION_RATE_1_INIT = input_data$TRANSITION_RATE_1_INIT[job.id],
	TRANSITION_RATE_2_INIT = input_data$TRANSITION_RATE_2_INIT[job.id],
	TRAIT_RATE_MEAN_1 = input_data$TRAIT_RATE_MEAN_1[job.id],
	TRAIT_RATE_MEAN_2 = input_data$TRAIT_RATE_MEAN_2[job.id],
	TRAIT_RATE_SIGMA_1 = input_data$TRAIT_RATE_SIGMA_1[job.id],
	TRAIT_RATE_SIGMA_2 = input_data$TRAIT_RATE_SIGMA_2[job.id],
	KAPPA_PRIOR_M = input_data$KAPPA_PRIOR_M[job.id],
	KAPPA_PRIOR_S = input_data$KAPPA_PRIOR_S[job.id],
	# Strict clock specific inputs
	CLOCK_RATE_INIT = input_data$CLOCK_RATE_INIT[job.id],
	TRANSITION_RATE_ALPHA = input_data$TRANSITION_RATE_ALPHA[job.id],
	TRANSITION_RATE_BETA = input_data$TRANSITION_RATE_BETA[job.id],
	# UCLD specific inputs
	UCLD_SIGMA_INIT = input_data$UCLD_SIGMA_INIT[job.id],
	seed = 89243 + job.id,
	TYPE_SWITCH_INIT = input_data$TYPE_SWITCH_INIT[job.id],
	TYPE_SWITCH_CLOCK_RATE = input_data$TYPE_SWITCH_CLOCK_RATE[job.id],
	TYPE_SWITCH_ALPHA = input_data$TYPE_SWITCH_RATE_ALPHA[job.id],
	TYPE_SWITCH_BETA = input_data$TYPE_SWITCH_RATE_BETA[job.id],
	max_start_date = input_data$max_start_date[[job.id]],
	root_trait = input_data$root_trait[[job.id]]
)


# Save the output
saveRDS(trees, file = paste0("./flu_output_4_13/", input_data$model[job.id], "_p04_p05_type-linked_time_tree_no_bulk.rds"))
