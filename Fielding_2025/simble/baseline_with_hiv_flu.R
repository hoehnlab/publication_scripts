library(airr)
library(dplyr)
library(shazam)
library(viridis)
source("./flu_data_processing.R")

file = "200gen_selection233853/all_samples_airr.tsv"
sim_full_data <- read_rearrangement(file)
sim_full_data$sample_time <- factor(sim_full_data$sample_time, levels = c("200", "150", "100", "50"))
sim_h <- filter(sim_full_data, locus == "IGH")

##########################################################
hiv_clones <- readRDS("~/simble-validation/imbalance_metric_hiv_trees.rds")
simble_selection_clones <- readRDS("~/simble-validation/imbalance_metric_selection_trees.rds")
simble_neutral_clones <- readRDS("~/simble-validation/imbalance_metric_neutral_trees.rds")

hiv_tsv <- read_rearrangement("Liao_2013_filtered_genotyped_germ-pass.tsv")
hiv_tsv$new_sequence_id <- paste0(hiv_tsv$sequence_id,"-", rownames(hiv_tsv))
used_seq_ids <- c()
clone_ids <- c()

for (i in 1:nrow(hiv_clones)) {
  used_seq_ids <- c(used_seq_ids, hiv_clones$data[[i]]@data$sequence_id)
  curr_clone_ids <- rep(hiv_clones$clone_id[i], length(used_seq_ids))
  names(curr_clone_ids) <- used_seq_ids
  clone_ids <- c(clone_ids, curr_clone_ids)
}
hiv_filtered <- filter(hiv_tsv, new_sequence_id %in% used_seq_ids)
hiv_filtered$clone_id <- clone_ids[hiv_filtered$new_sequence_id]

# simble selection
used_seq_ids <- c()
for (i in 1:nrow(simble_selection_clones)) {
  used_seq_ids <- c(used_seq_ids, simble_selection_clones$data[[i]]@data$sequence_id)
}

sim_selection_filtered <- filter(sim_h, sequence_id %in% used_seq_ids)

hiv_filtered$source <- "hiv"


z_hiv <- collapseClones(hiv_filtered, fields=c("source"))
b_hiv <- calcBaseline(z_hiv, regionDefinition = IMGT_V)
g_hiv <- groupBaseline(b_hiv, groupBy = "source")

z_sim <- collapseClones(sim_selection_filtered, germlineColumn = "germline_alignment", fields=c("sample_time"))
b_sim <- calcBaseline(z_sim, regionDefinition = IMGT_V)
g_sim <- groupBaseline(b_sim, groupBy = "sample_time")


### Potential flu baseline
flu_results = getProcessedFluResults()

flu_h <- filter(flu_results, locus == "IGH")
flu_l <- filter(flu_results, locus != "IGH")

flu_h$source <- "flu"
flu_h$timepoint <- factor(flu_h$timepoint, levels = c("d0", "d5", "d12", "d28", "d60"))

z_flu <- collapseClones(flu_h, fields=c("source"))
b_flu <- calcBaseline(z_flu, regionDefinition = IMGT_V)
g_flu <- groupBaseline(b_flu, groupBy = "source")


labeller <- c("fwr" = "FWR", "cdr" = "CDR")
all_baselines <- ggplot(plot(g_sim, "sample_time", sigmaLimits = c(-1.5, 1.5), silent = F)$data) +
  geom_line(data=plot(g_flu, "source", sigmaLimits = c(-1.5, 1.5), silent = F)$data, aes(x=SIGMA, y=DENSITY, linetype=source), color="red", linewidth=0.5) +
  geom_line(aes(y=DENSITY, x=SIGMA, color=sample_time), linewidth=0.8, linetype="solid") +
  geom_line(data=plot(g_hiv, "source", sigmaLimits = c(-1.5, 1.5), silent = F)$data, aes(x=SIGMA, y=DENSITY, linetype=source), color="red", linewidth=0.5) +
  # geom_line(data = temp_hiv$data, aes(x=SIGMA, y=DENSITY), color="grey", linewidth=1, linetype=2)
  facet_wrap(~region, ncol=1, labeller = as_labeller(labeller), strip.position = "right") + 
  scale_linetype_manual(values = c("hiv" = "21", "flu" = "4121"), labels=c("flu"="Flu", "hiv"="HIV")) +
  scale_color_viridis(discrete = TRUE, direction=-1, labels=time_labs) + theme_bw()

all_baselines

saveRDS(g_sim, file="~/simble-validation/sim_selection_baseline.rds")
saveRDS(g_flu, file="~/simble-validation/flu_baseline.rds")
saveRDS(g_hiv, file="~/simble-validation/hiv_baseline.rds")
saveRDS(all_baselines, file="~/simble-validation/all_baselines_plot.rds")
