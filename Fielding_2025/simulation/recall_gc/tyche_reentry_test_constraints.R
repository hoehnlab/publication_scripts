######################
## Hunter J. Melton ##
## Jessie J. Fielding ##
##### 4/28/2025 #####

# Analyzing simulated GC reentry data with TyCHE
##########################
# Load packages and data #
##########################

library(airr)
library(dplyr)
library(dowser)
library(ggtree)
library(purrr)


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
    bcrs <- read_rearrangement("./simble_sims_gc_reentry_12_18/all_samples_airr.tsv")
} else if (simulation_id == "un") {
    # Use neutral simulation data
    bcrs <- read_rearrangement("./simble_sims_gc_reentry_uniform_neutral_12_18/all_samples_airr.tsv")
} else {
    stop("Invalid simulation argument provided. Use 'selection' or 'neutral'.")
}


bcrs_heavy <- bcrs %>% 
  filter(locus == "IGH") %>% 
  mutate(celltype = ifelse(celltype == "gc_b_cell", "GC", "other")) %>%
  mutate(sample_time = as.numeric(sample_time))

# Get array ID to run different models in parallel
job.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set the path to the xml-writer directory
xml_writer_path <- "~/xml-templates/"

beast_path = "~/beast/bin/"

sc_dir <- "./get_GC_rates_SC"
if (simulation_id == "un") {
    sc_dir <- paste0(sc_dir, "_uniform_neutral")
}
sc_dir <- paste0(sc_dir, "/")


# Other variables
first_gc_end_time <- 100

##############################
# Find GC B cells clock rate #
##############################

# Use the estimate from previous analysis

# Make sure the file exists first
if(!file.exists(paste0(sc_dir, "get_GC_rates_SC_trees.rds"))) {
    print(paste0("Looking for file at: ", paste0(sc_dir, "get_GC_rates_SC_trees.rds")))
    stop("Where's the SC GC rates file?")
}
SC_trees <- readRDS(paste0(sc_dir, "get_GC_rates_SC_trees.rds"))

# # Filter to whatever clones we're currently using only
# SC_trees <- SC_trees %>% filter(clone_id %in% clones_using)

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
if (simulation_id == "un") {
    f_clones <- formatClones(bcrs_heavy, filterstop=FALSE, traits = c("celltype", "sample_time"), germ = "germline_alignment")
} else {
    f_clones <- formatClones(bcrs_heavy, traits = c("celltype", "sample_time"), germ = "germline_alignment")
}

# Build initial genetic distance trees
trees <- getTrees(f_clones, build = "pml", sub_model= "HKY", nproc = 10)


TRANSITION_RATE_ALPHA <- "0.1"
TRANSITION_RATE_BETA <- "1.0"


# Set up to run analysis
input_data <- tibble(
    model = c(
        "TyCHE_EO_Est_Emp_no_constraints",
        "TyCHE_EO_Est_Emp_only_max_start",
        "TyCHE_EO_Est_Emp_only_root_prior"
        ),
    template_path = paste0(xml_writer_path, "TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml"),
    ignore = list(c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator", "rootType"),
                    c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator", "rootType"),
                    c("freqParameter", "traitfrequencies", "typeLinkedRates", "rateIndicator", "rootType")),
    RATE_INDICATORS = "1 1",
    # Expected occupancy inputs
    TRANSITION_RATE_ALPHA_1 = TRANSITION_RATE_ALPHA,
    TRANSITION_RATE_BETA_1 = TRANSITION_RATE_BETA,
    TRANSITION_RATE_ALPHA_2 = TRANSITION_RATE_ALPHA,
    TRANSITION_RATE_BETA_2 = TRANSITION_RATE_BETA,
    TRANSITION_RATE_1_INIT = "0.2",
    TRANSITION_RATE_2_INIT = "0.0004",
    TRAIT_RATE_MEAN_1 = gc_clock_rate_mean, 
    TRAIT_RATE_MEAN_2 = "0.000001",
    TRAIT_RATE_SIGMA_1 = 0.01 * gc_clock_rate_mean,
    TRAIT_RATE_SIGMA_2 = "0.001",
    KAPPA_PRIOR_M = "0.67",
    KAPPA_PRIOR_S = "0.2",
    TYPE_SWITCH_INIT = "0.002",
    TYPE_SWITCH_ALPHA = 0.001,
    TYPE_SWITCH_BETA = 5.0,
)

##################
# Get Time Trees #
##################

index <- job.id
if(index < 1 || index > nrow(input_data)) stop(paste0("Index ", index, " is out of bounds"))

if (any(grep("max_start", input_data$model[index]))) {
    max_start_date <- -1000
} else {
    max_start_date <- NULL
}

if (any(grep("root_prior", input_data$model[index]))) {
    root_prior <- 0
} else {
    root_prior <- NULL
}


trees <- getTimeTreesIterate(trees,
    iterations = 10,
    ignore = input_data$ignore[[index]],
    beast = beast_path, 
    trait = "celltype",
    time = "sample_time",
    dir = paste0("/scratch/", folder_name, "_", simulation_id, "/"), 
    id = input_data$model[index], 
    template = input_data$template_path[index],
    nproc = 20,
    include_germline = TRUE,
    mcmc_length = 500000000,
    # mcmc_length = 5000,
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
    TYPE_SWITCH_INIT = input_data$TYPE_SWITCH_INIT[index],
    TYPE_SWITCH_ALPHA = input_data$TYPE_SWITCH_ALPHA[index],
    TYPE_SWITCH_BETA = input_data$TYPE_SWITCH_BETA[index],
    TRANSITION_RATE_1_INIT = input_data$TRANSITION_RATE_1_INIT[index],
    TRANSITION_RATE_2_INIT = input_data$TRANSITION_RATE_2_INIT[index],
    # Strict clock specific inputs
    CLOCK_RATE_INIT = input_data$CLOCK_RATE_INIT[index],
    TRANSITION_RATE_ALPHA = input_data$TRANSITION_RATE_ALPHA[index],
    TRANSITION_RATE_BETA = input_data$TRANSITION_RATE_BETA[index],
    # UCLD specific inputs
    UCLD_SIGMA_INIT = input_data$UCLD_SIGMA_INIT[index],
    NODES_TYPE_INIT=0,
    seed = 89243 + index,
    max_start_date = max_start_date,
    root_trait=root_prior,
    continue=TRUE
    )

# Save the output
saveRDS(trees, paste0("test_constraints/output_", folder_name, "/gc_reentry_sims_", input_data$model[index], "_", simulation_id, "_trees.rds"))


