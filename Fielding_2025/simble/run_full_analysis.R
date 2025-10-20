library(devtools)

# uninstall any existing dowser version
remove.packages("dowser")

# install dev version of dowser with sampleClones function
devtools::install_github("immcantation/dowser", ref="3c7a91f25930536e9319e8a052b3e213cf3d3419")

# source all the individual scripts to process data for each figure

# example trees and tree imbalance (must calculate tree imbalance first)
source("imbalance_metric.R")
source("simble-vs-hiv-trees.R")

# baseline plots (for selection in fig 2 and neutral for supplemental)
source("simble-validation-baseline.R")

# shm plots
source("simble_flu_hl_shm.R")

# time of differentiation mbc vs pc
source("migration_differentiation.R")

# affinity and cross-reactivity plots
source("simble_affinity_cross_reactivity.R")

# supplemental plot for mutations per site
source("mutations_per_site.R")

# source the script to format all the plots to match the theme of the paper and
# compile them
source("compile_plots.R")

print_loaded_packages <- function() {
  si <- sessionInfo()
  
  # attached packages
  attached <- sub("^package:", "", si$otherPkgs |> names())
  attached_pkgs <- si$otherPkgs
  
  # versions
  attached_info <- vapply(attached, function(pkg) {
    ver <- attached_pkgs[[pkg]]$Version
    paste0(pkg, " v", ver)
  }, character(1))
  
  cat(attached_info, sep = "\n")
}
print_loaded_packages()
