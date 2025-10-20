######################
## Hunter J. Melton ##
###### 8/20/2025 #####

# Analyzing simulated GC reentry data with TyCHE - uniform neutral data
##########################
# Load packages and data #
##########################

library(airr)
library(dplyr)
library(dowser)
library(ggtree)
library(purrr)

bcrs <- read_rearrangement("./simble_sims_gc_reentry_uniform_neutral_8_29/all_samples_airr.tsv")
bcrs_heavy <- bcrs %>% 
  filter(locus == "IGH") %>% 
  mutate(celltype = ifelse(celltype == "default", "GC", "other")) %>%
  mutate(sample_time = as.numeric(sample_time))

# Use all clones
clones_using <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
				  "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
bcrs_heavy <- bcrs_heavy %>%
	filter(clone_id %in% clones_using)

# Get array ID to run different models in parallel
job.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set the path to the xml-writer directory
xml_writer_path <- "../../xml-writer/"

# Set up the path to the beast executable
beast_path = "~/software/beast/bin/"

# Other variables
first_gc_end_time <- 100

# Function to fix empty fields in uniform neutral data - required for dowser analysis
ensure_neutral_fields <- function(df) {
  if (!"productive" %in% names(df))                df$productive <- TRUE
  if (!"v_call" %in% names(df))                    df$v_call <- "IGHV1-1*01"
  if (!"d_call" %in% names(df))                    df$d_call <- "IGHD1-1*01"
  if (!"j_call" %in% names(df))                    df$j_call <- "IGHJ1*01"
  if (!"rev_comp" %in% names(df))                  df$rev_comp <- FALSE
  if (!"germline_alignment_d_mask" %in% names(df)) df$germline_alignment_d_mask <- df$germline_alignment
  df
}
bcrs_heavy <- ensure_neutral_fields(bcrs_heavy)

##############################
# Find GC B cells clock rate #
##############################

# We need to figure out reasonable estimates for the clock rate of GC B cells
# We'll do this by fitting a strict clock model to the GC B cells sampled during the first GC reaction at generations 50 and 100

# Only fit strict clock once - job.id == 1 only, other arrays will be submitted dependent on array=1
if (job.id == 1){

	# Filter to the first GC reaction (sample_time = 100 and earlier)
	bcrs_heavy_gc_first_rxn <- bcrs_heavy %>% filter(sample_time <= first_gc_end_time, celltype == "GC")

	# Format the clones
	f_clones_gc_first_rxn <- formatClones(bcrs_heavy_gc_first_rxn, traits = c("celltype", "sample_time", "location"),
										  germ = "germline_alignment", filterstop=FALSE)

	# Build the original trees
	trees_gc_first_rxn <- getTrees(f_clones_gc_first_rxn, build = "pml", sub_model= "HKY", nproc = 10)

	# Run the strict clock
	trees <- trees_gc_first_rxn
	trees <- getTimeTreesIterate(trees,
		iterations = 10,
		ignore = c("freqParameter"),
		beast = beast_path, 
		time = "sample_time",
		dir = "./get_GC_rates_SC_uniform_neutral/",
		id = "get_GC_rates_SC", 
		template = paste0(xml_writer_path, "templates/custom/StrictClock/StrictClock_Standard_EmpFreq.xml"),
		nproc = 20,
		include_germline = TRUE,
		mcmc_length = 100000000,
    	log_every = "auto",
		log_target = 2000,
		KAPPA_PRIOR_M = "0.67",
		KAPPA_PRIOR_S = "0.2",
		CLOCK_RATE_INIT = "0.001",
		seed = 89243 + job.id
	)

	saveRDS(trees, "./get_GC_rates_SC_uniform_neutral/get_GC_rates_SC_trees.rds")

} else { # Do the actual analysis

	# Figure out the clock rates for the GC B cells for each clone
	# Make sure the file exists first
	if(!file.exists("./get_GC_rates_SC_uniform_neutral/get_GC_rates_SC_trees.rds")) {
		stop("Where's the SC GC rates file?")
	}
	SC_trees <- readRDS("./get_GC_rates_SC_uniform_neutral/get_GC_rates_SC_trees.rds")

	# Filter to whatever clones we're currently using only
	SC_trees <- SC_trees %>% filter(clone_id %in% clones_using)

	# Make sure (some of) the SC trees have converged
	SC_trees_conv <- SC_trees %>% filter(below_ESS == 0)
	if (nrow(SC_trees_conv) == 0) {
		stop("No converged trees found in the SC GC rates file. You need to check this.")
	}

	gc_clock_rates <- sapply(SC_trees_conv$parameters, function(x) {
  		x$mean[x$item == "geneticClockRate"]
	})
	# Flag if there are any NAs
	if (sum(is.na(gc_clock_rates)) > 0) {
		cat("Warning: There are NA clock rates calculated for the GC B cells.")
		num_nas <- sum(is.na(gc_clock_rates))
		cat("Number of NA clock rates:", num_nas)
	}
	gc_clock_rate_mean <- mean(gc_clock_rates, na.rm = TRUE)
	
	# Format clones
	f_clones <- formatClones(bcrs_heavy, traits = c("celltype", "sample_time"), germ = "germline_alignment", filterstop = FALSE)

	# Build initial genetic distance trees
	trees <- getTrees(f_clones, build = "pml", sub_model= "HKY", nproc = 10)

	# Set up to run analysis
	input_data <- tibble(
		model = c("expectedOccupancy_FixedTraitClockRates_EmpFreq", "expectedOccupancy_EstTraitClockRates_EmpFreq",
					"instantSwitch_EstTraitClockRates_EmpFreq",
					"strictClock_AncestralReconstruction_EmpFreq", "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq"),
		template_path = c(paste0(xml_writer_path, "templates/custom/TypeLinked/TraitLinkedExpectedOccupancy_FixedTraitClockRates_EmpFreq.xml"),
						  paste0(xml_writer_path, "templates/custom/TypeLinked/TraitLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml"),
						  paste0(xml_writer_path, "templates/custom/TypeLinked/TraitLinkedInstantSwitch_EstTraitClockRates_EmpFreq.xml"),
						  paste0(xml_writer_path, "templates/custom/StrictClock/StrictClock_AncestralReconstruction_EmpFreq.xml"),
						  paste0(xml_writer_path, "templates/custom/UCLD/UCRelaxedClock_AncestralReconstruction_EmpFreq.xml")),
		ignore = list(c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator"),
					  c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator"),
					  c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator"),
					  c("freqParameter", "traitfrequencies", "rateIndicator"),
					  c("freqParameter", "traitfrequencies", "rateIndicator", "rateCategories")),
		RATE_INDICATORS = "1 1",
		# Expected occupancy inputs
		TRANSITION_RATE_ALPHA_1 = "0.1",
		TRANSITION_RATE_BETA_1 = "1.0",
		TRANSITION_RATE_ALPHA_2 = "0.1",
		TRANSITION_RATE_BETA_2 = "1.0", 
		TRAIT_RATE_MEAN_1 = gc_clock_rate_mean, 
		TRAIT_RATE_MEAN_2 = "0.000001",
		TRAIT_RATE_SIGMA_1 = 0.01 * gc_clock_rate_mean,
		TRAIT_RATE_SIGMA_2 = "0.001",
		KAPPA_PRIOR_M = "0.67",
		KAPPA_PRIOR_S = "0.2",
		# Strict clock/UCLD inputs
		CLOCK_RATE_INIT = gc_clock_rate_mean/3,
		TRANSITION_RATE_ALPHA = "0.1",
		TRANSITION_RATE_BETA = "1.0",
		UCLD_SIGMA_INIT = "0.5"
	)

	##################
	# Get Time Trees #
	##################

	# Set up and check index (job.id - 1)
	index <- job.id - 1
	if(index < 1 || index > nrow(input_data)) stop(paste0("Index ", index, " is out of bounds"))

	trees <- getTimeTreesIterate(trees,
		iterations = 10,
		ignore = input_data$ignore[[index]],
		beast = beast_path, 
		trait = "celltype",
		time = "sample_time",
		dir = "/dartfs-hpc/scratch/f007p0j/gc_reentry_sims_uniform_neutral/9_9/", 
		id = input_data$model[index], 
		template = input_data$template_path[index],
		nproc = 20,
		include_germline = TRUE,
		mcmc_length = 1000000000,
		log_every = "auto",
		log_target = 2000,
		RATE_INDICATORS = input_data$RATE_INDICATORS[index],
		TRANSITION_RATE_ALPHA_1 = input_data$TRANSITION_RATE_ALPHA_1[index],
		TRANSITION_RATE_BETA_1 = input_data$TRANSITION_RATE_BETA_1[index],
		TRANSITION_RATE_ALPHA_2 = input_data$TRANSITION_RATE_ALPHA_2[index],
		TRANSITION_RATE_BETA_2 = input_data$TRANSITION_RATE_BETA_2[index],
		TRAIT_RATE_MEAN_1 = input_data$TRAIT_RATE_MEAN_1[index],
		TRAIT_RATE_MEAN_2 = input_data$TRAIT_RATE_MEAN_2[index],
		TRAIT_RATE_SIGMA_1 = input_data$TRAIT_RATE_SIGMA_1[index],
		TRAIT_RATE_SIGMA_2 = input_data$TRAIT_RATE_SIGMA_2[index],
		KAPPA_PRIOR_M = input_data$KAPPA_PRIOR_M[index],
		KAPPA_PRIOR_S = input_data$KAPPA_PRIOR_S[index],
		# Strict clock specific inputs
		CLOCK_RATE_INIT = input_data$CLOCK_RATE_INIT[index],
		TRANSITION_RATE_ALPHA = input_data$TRANSITION_RATE_ALPHA[index],
		TRANSITION_RATE_BETA = input_data$TRANSITION_RATE_BETA[index],
		# UCLD specific inputs
		UCLD_SIGMA_INIT = input_data$UCLD_SIGMA_INIT[index],
		seed = 89243 + index
	)

	# Save the output
	saveRDS(trees, paste0("./output_uniform_neutral_9_9/gc_reentry_sims_uniform_neutral", input_data$model[index], "_", paste(clones_using, collapse = "_"), "_trees.rds"))


}
