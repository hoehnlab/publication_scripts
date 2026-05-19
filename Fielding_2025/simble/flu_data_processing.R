library(airr)
library(dplyr)
library(shazam)

flu_processed_file <- "~/simble-validation/flu_data_processed.rds"

getProcessedFluResults <- function(use_rds=TRUE) {
  if (file.exists(flu_processed_file) && use_rds) {
    flu_results <- readRDS(flu_processed_file)
  } else {
    flu_data <- alakazam::readChangeoDb("~/docker-share/bcr-db_3_airr.tsv")
    
    flu_data <- flu_data %>% filter(productive)
    
    multi_heavy <- table(filter(flu_data, locus == "IGH")$cell_id)
    multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]
    
    flu_data <- filter(flu_data, !cell_id %in% multi_heavy_cells)
    
    # split cells by heavy and light chains
    heavy_cells <- filter(flu_data, locus == "IGH")$cell_id
    light_cells <- filter(flu_data, locus == "IGK" | locus == "IGL")$cell_id
    no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]
    
    flu_data <- filter(flu_data, !cell_id %in% no_heavy_cells)
    
    dist_nearest <- distToNearest(filter(flu_data, locus == "IGH"), nproc = 1, cellIdColumn = "cell_id")
    p1 <- ggplot(subset(dist_nearest, !is.na(dist_nearest)),
                 aes(x = dist_nearest)) +
      theme_bw() +
      xlab("Hamming distance") + ylab("Count") +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      geom_histogram(color = "white", binwidth = 0.02) +
      theme(axis.title = element_text(size = 18))
    plot(p1)
    
    # find threshold for cloning automatically
    threshold_output <- shazam::findThreshold(dist_nearest$dist_nearest,
                                              method = "gmm", model = "gamma-norm",
                                              cutoff = "user", spc = 0.995)
    threshold <- threshold_output@threshold
    cat(threshold, file = "flu_threshold.txt")
    
    
    plot(threshold_output, binwidth = 0.02, silent = TRUE) +
      theme(axis.title = element_text(size = 18))
    
    results <- scoper::hierarchicalClones(flu_data, cell_id = 'cell_id',
                                          threshold = threshold, only_heavy = TRUE,
                                          split_light = FALSE, summarize_clones = FALSE, nproc=3)
    
    
    
    
    example <- resolveLightChains(results)
    
    resolved <- createGermlines(example, references = references, clone = "clone_subgroup_id", nproc = 1)
    
    hlclones = formatClones(resolved,chain="HL",
                            split_light=TRUE,heavy="IGH",cell="cell_id",
                            trait=c("vj_gene","c_call"), minseq = 5, collapse = TRUE, nproc = 3)
    
    keep_clone_ids <- unique(hlclones$clone_id)
    flu_results <- filter(resolved, clone_subgroup_id %in% keep_clone_ids)
    saveRDS(flu_results, file = flu_processed_file)
  }
  return(flu_results)
}