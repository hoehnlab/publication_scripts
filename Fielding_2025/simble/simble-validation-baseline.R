library(airr)
library(dplyr)
library(shazam)
library(viridis)
library(dowser)

# baseline_with_hiv_flu.R will do selection, so only do neutral here

file_neutral = "200gen_neutral866228/all_samples_airr.tsv"
sim_full_data_neutral <- read_rearrangement(file_neutral)
sim_h_neutral <- filter(sim_full_data_neutral, locus == "IGH")

# Filter to only include sequences without early stop codons
formatted_clones_neutral <- formatClones(sim_h_neutral, minseq = 3, filterstop = TRUE, traits = c("sample_time"), germ = "germline_alignment", nproc=3)
productive_seq_ids <- c()
for (i in 1:nrow(formatted_clones_neutral)) {
  productive_seq_ids <- c(productive_seq_ids, formatted_clones_neutral$data[[i]]@data$sequence_id)
}

sim_h_neutral_filtered <- filter(sim_h_neutral, sequence_id %in% productive_seq_ids)
fileConn<-file("~/simble-validation/filtered_stops_baseline_neutral.txt")
writeLines(
  c(paste0(
    "Filtered ", 
    nrow(sim_h_neutral) - nrow(sim_h_neutral_filtered), 
    " sequences with early stop codons, out of ", 
    nrow(sim_h_neutral), 
    " sequences."
    )), 
  fileConn)
close(fileConn)

# Calculate clonal consensus and selection using the BASELINe method
sim_h_neutral_filtered$sample_time <- factor(sim_h_neutral_filtered$sample_time, levels = c("200", "150", "100", "50"))
z_neutral <- collapseClones(sim_h_neutral_filtered, germlineColumn = "germline_alignment", fields=c("sample_time"))
b_neutral <- calcBaseline(z_neutral, regionDefinition = IMGT_V)

## calcBaseline will calculate observed and expected mutations for clonal_sequence using clonal_germline as a reference.

# Combine selection scores for all clones in each group
g_neutral <- groupBaseline(b_neutral, groupBy = "sample_time")

g_neutral@regions <- toupper(g_neutral@regions)
g_neutral@stats$region <- toupper(g_neutral@stats$region) 
g_neutral@pdfs$FWR <- g_neutral@pdfs$fwr
g_neutral@pdfs$CDR <- g_neutral@pdfs$cdr

## Grouping BASELINe probability density functions...
## Calculating BASELINe statistics...

# Plot probability densities for the selection pressure
neutral <- plot(g_neutral, "sample_time", sigmaLimits = c(-3, 2), silent = F) + ggtitle("Simble 200 generations, neutral (sampled at 50, 100, 150, 200), 100 clones")

# Save the plots
saveRDS(neutral, file="~/simble-validation/selection_strength_neutral_plot.rds")
