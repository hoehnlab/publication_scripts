library(airr)
library(dplyr)
library(treestats)
library(dowser)

set.seed(847445)
downsample_size <- 40
min_seqs <- downsample_size

# function to calculate imbalance metric
get_imbalance_metric <- function(tree) {
  # use the maximum width over the maximum depth of the tree 
  # (depth is essentially divergence of tree with all branch lengths = 1)
  mw_over_md(tree)
}


# process HIV data, slightly different pre-processing required
# read in pre-processed trees from external data source
hiv_trees <- readRDS("~/Downloads/Cluster Trees Gen.rds")

# filter to only include trees with enough sequences
hiv_trees <- filter(hiv_trees, seqs > min_seqs)

# ignore previously built trees so we can rebuild from clones object
hiv_clones <- hiv_trees[,1:6]

# for fair comparison with simble data, sample a set number of sequences per clone, sampling evenly across timepoints
hiv_clones <- sampleClones(hiv_clones, size = downsample_size, group="time")

number_of_hiv_clones <- length(unique(hiv_clones$clone_id))

# build trees using "pml" method to match simble processing
hiv_trees <- getTrees(hiv_clones, build = "pml", nproc = 3)

# measure imbalance and write to RDS
hiv_trees$imbalance_metric <- unlist(lapply(hiv_trees$trees, function(x) get_imbalance_metric(x)))
hiv_trees$model <- "hiv"
saveRDS(hiv_trees, file = "~/simble-validation/imbalance_metric_hiv_trees.rds")


make_imbalance_rds <- function(input_file, output_file, model) {
  # read in data, format data, format clones
  data <- read_rearrangement(input_file)
  filtered_sim <- transform(data, sample_time = as.numeric(sample_time))
  clones <- formatClones(filtered_sim, minseq = 3, filterstop = FALSE, traits = c("sample_time"), germ = "germline_alignment", nproc=3)
  clones <- sampleClones(clones, size = downsample_size, group="sample_time")
  
  # make sure to use the correct number of clones (the same number as hiv data)
  clones <- clones[1:number_of_hiv_clones, ]
  
  # build trees
  trees <- getTrees(clones, build = "pml", nproc = 3)
  
  # measure imbalance metric
  trees$imbalance_metric <- unlist(lapply(trees$trees, function(x) get_imbalance_metric(x)))
  trees$model <- model
  
  # saveRDS
  saveRDS(trees, file = output_file)
}

# process selection simble data
selection_folder = "/Volumes/HoehnK/jessie/simble-validation/200gen_selection233853/"
file <- paste0(selection_folder, "all_samples_airr.tsv")
outfile <- "~/simble-validation/imbalance_metric_selection_trees.rds"
make_imbalance_rds(file, outfile, "selection")

# process neutral simble data
neutral_folder = "/Volumes/HoehnK/jessie/simble-validation/200gen_neutral866228/"
file <- paste0(neutral_folder, "all_samples_airr.tsv")
outfile <- "~/simble-validation/imbalance_metric_neutral_trees.rds"
make_imbalance_rds(file, outfile, "neutral")

