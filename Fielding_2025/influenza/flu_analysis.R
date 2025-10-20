######################
## Hunter J. Melton ##
###### 8/19/2025 #####

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

# Get the multi-group balancing version of sampleClones
# source("../../dowser/R/Clones.R")

# Load data
p04 <- read_rearrangement("data/P04_BCR_03182022_airr.tsv") 
p05 <- read_rearrangement("data/P05_BCR_03182022_airr.tsv") 

# Get array ID to run different models in parallel
job.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set the path to the xml-writer directory
xml_writer_path <- "../../xml-writer/"

#################
# Filter clones #
#################

# We want to analyze clones that:
# 1) Have GC B cells at both time points
# 2) Are measurably evolving over time
# 3) Contain at least one sequence that binds to flu antigen

# Set the time(point) of samples
# Year 2 refers to after the second vaccine, p04 was week 35 (7*35 = 245), p05 was week 38 (7*38 = 266)
p04 <- p04 %>% mutate(day = as.numeric(day), year = as.numeric(year), timepoint = ifelse(year == 2, day + 245, day))
p05 <- p05 %>% mutate(day = as.numeric(day), year = as.numeric(year), timepoint = ifelse(year == 2, day + 266, day))

# Standardize the annotated cell types column
p04 <- p04 %>% mutate(gex_annotation = ifelse(is.na(gex_annotation), "missing", gex_annotation))
p05 <- p05 %>% mutate(gex_annotation = ifelse(is.na(gex_annotation), "missing", gex_annotation))

# Find all clones that have GC B cells at both time points
p04_multi_year_GC_clones <- p04 %>% group_by(clone_id) %>%
	filter(gex_annotation == "GC") %>%
	summarize(y1 = sum(year==1), y2=sum(year==2)) %>%
	filter(y1 > 1, y2 > 1) %>%
	arrange(desc(y1+y2))
#   clone_id    y1    y2
#   <chr>    <int> <int>
# 1 150569     150     2
# 2 97116       76     5
# 3 70058       27    10

p05_multi_year_GC_clones <- p05 %>% group_by(clone_id) %>%
	filter(gex_annotation == "GC") %>%
	summarize(y1 = sum(year==1), y2=sum(year==2)) %>%
	filter(y1 > 1, y2 > 1) %>%
	arrange(desc(y1+y2))

# Filter out sequences from bulk data (we don't know the cell type), standardize the sampling locations, swap out NAs for "Unknown" in sequence/clone binding to flu antigen column
p04_not_bulk <- p04 %>% filter(datatype != "bulk") %>%
	mutate(is.agseq = ifelse(is.na(is.agseq), "Unknown", is.agseq),
		   is.agclone = ifelse(is.na(is.agclone), "Unknown", is.agclone),
		   tissue = recode(sampletype, "GC" = "LN", "BMPC" = "BM", "FNA" = "LN", "IgD-BCell"="PBMC", "PBMC-ASC"="PBMC"))
p05_not_bulk <- p05 %>% filter(datatype != "bulk") %>%
	mutate(is.agseq = ifelse(is.na(is.agseq), "Unknown", is.agseq),
		   is.agclone = ifelse(is.na(is.agclone), "Unknown", is.agclone),
		   tissue = recode(sampletype, "GC" = "LN", "BMPC" = "BM", "FNA" = "LN", "IgD-BCell"="PBMC", "PBMC-ASC"="PBMC"))

# Format clones
f04_not_bulk <- formatClones(filter(p04_not_bulk, clone_id %in% p04_multi_year_GC_clones$clone_id),
	traits=c("gex_annotation", "timepoint", "year" ,"day", "is.agseq", "is.agclone", "tissue")) 
f05_not_bulk <- formatClones(filter(p05_not_bulk, clone_id %in% p05_multi_year_GC_clones$clone_id),
	traits=c("gex_annotation", "timepoint", "year" ,"day", "is.agseq", "is.agclone", "tissue")) 

# Build trees for testing for measurable evolution
trees04_not_bulk <- getTrees(f04_not_bulk, build = "pml", sub_model= "HKY", nproc = 4)
trees05_not_bulk <- getTrees(f05_not_bulk, build = "pml", sub_model= "HKY", nproc = 4)

# Get correlation of SHM and time 
set.seed(892431) # TyCHE1
ct04_not_bulk <- correlationTest(trees04_not_bulk, time="timepoint", permutations=10000, nproc = 4)
set.seed(892432) # TyCHE2
ct05_not_bulk <- correlationTest(trees05_not_bulk, time="timepoint", permutations=10000, nproc = 4)

# Select clones that have a significant positive correlation between SHM and time
select_clones04 <- ct04_not_bulk %>% filter(p < 0.05 & slope > 0) %>% pull(clone_id)
select_clones05 <- ct05_not_bulk %>% filter(p < 0.05 & slope > 0) %>% pull(clone_id)
# Clone 70058 for P04
# Clones 26299, 69914, 121056, 129238, 19989, 111394, 101085, 43876

# Filter clones for evolution over time
clones04_not_bulk <- filter(trees04_not_bulk, clone_id %in% select_clones04)
clones05_not_bulk <- filter(trees05_not_bulk, clone_id %in% select_clones05)
# Filter clones based on binding to flu antigen
include_p4 <- c()
include_p5 <- c()
include_p4 <- lapply(clones04_not_bulk$data, function(x) {
  include_p4 = c(include_p4, ifelse(x@data$is.agclone[1] == "Yes", TRUE, FALSE))
  include_p4
}) %>% unlist()
include_p5 <- lapply(clones05_not_bulk$data, function(x) {
  include_p5 = c(include_p5, ifelse(x@data$is.agclone[1] == "Yes", TRUE, FALSE))
  include_p5
}) %>% unlist()
p4_include_clones_agclone <- clones04_not_bulk$clone_id[include_p4]
p5_include_clones_agclone <- clones05_not_bulk$clone_id[include_p5] # lose 69914 in P05
clones04_not_bulk <- filter(clones04_not_bulk, clone_id %in% p4_include_clones_agclone)
clones05_not_bulk <- filter(clones05_not_bulk, clone_id %in% p5_include_clones_agclone)

##########################################
# Prep clones and other inputs for TyCHE #
##########################################

# Rename the GC and other B cells
clones04_not_bulk <- clones04_not_bulk %>%
  mutate(data = map(data, function(x) {
    x@data <- x@data %>%
      mutate(gex_annotation = if_else(gex_annotation == "GC", "germinal_center", "other"))
    x
  }))
clones05_not_bulk <- clones05_not_bulk %>%
  mutate(data = map(data, function(x) {
	x@data <- x@data %>%
	  mutate(gex_annotation = if_else(gex_annotation == "GC", "germinal_center", "other"))
	x
  }))

# Combine data from p04 and p05
clones_not_bulk <- bind_rows(clones04_not_bulk, clones05_not_bulk)

# Downsample to a maximum of 100 sequences per clone (only affects 26299), balance over gex_annotation and tissue
set.seed(8924337)
clones_not_bulk_downsampled <- sampleClones(clones_not_bulk, size = 100, group = c("gex_annotation", "tissue"))

# Save the intermediate step
saveRDS(clones_not_bulk_downsampled, "./clones_not_bulk_downsampled.rds")

# Prep input for TyCHE
p04_clones <- c("70058")
p05_clones <- c("26299", "121056", "129238", "19989", "111394", "101085", "43876")
clone_data <- data.frame(sample = rep(c("P04", "P05"), c(length(p04_clones), length(p05_clones))),
						 clone_id = c("70058", "26299", "121056", "129238", 
									  "19989", "111394", "101085", "43876"))
input_data <- tibble(
	model = c("expectedOccupancy_FixedTraitClockRates_EmpFreq", "expectedOccupancy_EstTraitClockRates_EmpFreq", "instantSwitch_EstTraitClockRates_EmpFreq",
				"strictClock_AncestralReconstruction_EmpFreq", "UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq"),
	clone_data = list(clone_data),
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
	TRAIT_RATE_MEAN_1 = 4.9E-3/7, 
	TRAIT_RATE_MEAN_2 = "0.000001",
	TRAIT_RATE_SIGMA_1 = 0.01 * (4.9E-3/7),
	TRAIT_RATE_SIGMA_2 = "0.001",
	KAPPA_PRIOR_M = "0.67",
	KAPPA_PRIOR_S = "0.2",
	# Strict clock/UCLD inputs
	CLOCK_RATE_INIT = 4.9E-3/20,
	TRANSITION_RATE_ALPHA = "0.1",
	TRANSITION_RATE_BETA = "1.0",
	UCLD_SIGMA_INIT = "0.5",
)

# Other input
beast_path = "~/software/beast/bin/"
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
	dir = "/dartfs-hpc/scratch/f007p0j/flu_beast_8_25/", 
	id = input_data$model[job.id], 
	template = input_data$template_path[job.id],
	nproc = 8,
	include_germline = TRUE,
	mcmc_length = 1000000000,
    log_every = "auto",
	log_target = 2000,
	max_start_date = -1460, # 4 years before the first time point
	RATE_INDICATORS = input_data$RATE_INDICATORS[job.id],
	TRANSITION_RATE_ALPHA_1 = input_data$TRANSITION_RATE_ALPHA_1[job.id],
	TRANSITION_RATE_BETA_1 = input_data$TRANSITION_RATE_BETA_1[job.id],
	TRANSITION_RATE_ALPHA_2 = input_data$TRANSITION_RATE_ALPHA_2[job.id],
	TRANSITION_RATE_BETA_2 = input_data$TRANSITION_RATE_BETA_2[job.id],
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
	seed = 89243 + job.id
)

# Save the output
saveRDS(trees, file = paste0("./flu_beast_8_25/", input_data$model[job.id], "_p04_p05_type-linked_time_tree_no_bulk.rds"))
