library(airr)
library(dplyr)

folder = "/Volumes/HoehnK/jessie/simble-validation/mutations_per_site_754659"
file = paste0(folder, "/selection/all_samples_airr.tsv")
file_neutral = paste0(folder, "/neutral/all_samples_airr.tsv")
file_uniform = paste0(folder, "/uniform/all_samples_airr.tsv")

uniform <- read_rearrangement(file_uniform)
neutral <- read_rearrangement(file_neutral)
selection <- read_rearrangement(file)

# Let's only look at heavy chains
uniform <- uniform %>%
  filter(locus == "IGH")
neutral <- neutral %>%
  filter(locus == "IGH")
selection <- selection %>%
  filter(locus == "IGH")

seqs <- list(uniform = uniform, neutral = neutral, selection = selection)
plot_data <- list()
plot_data_single <- list()

# For one clone per model, count the number of differences at each position in 
# the sequence_alignment vs. the germline alignment

models <- c("uniform", "neutral", "selection")
for (model in models) {
  cloneid <- seqs[[model]]$clone_id[1]
  bcrs_temp <- seqs[[model]] %>%
    filter(clone_id == cloneid)
  
  # Count the number of differences at each position from each time point
  plot_df <- data.frame(
    position = 1:nchar(bcrs_temp$sequence_alignment[1]),
    differences_per_seq = NA,
    clone_id = cloneid,
    model = model
  )
  
  seq <- bcrs_temp$sequence_alignment
  germline <- bcrs_temp$germline_alignment
  
  for (i in 1:nchar(bcrs_temp$sequence_alignment[1])) {
    # Count the number of differences at this position
    differences <- sum(substr(seq, i, i) != substr(germline, i, i))
    plot_df$differences_per_seq[i] <- differences / nrow(bcrs_temp)
  }
  
  plot_data_single[[paste0(model, "_", cloneid)]] <- plot_df
}

# get the CDR3 length -- this is the same for selection and neutral so just use one
cdr3_length <- as.numeric(selection$cdr3_aa_length[[1]])

get_region <- function(position, model) {
  if (model == "uniform") {
    return(NA)
  }
  if (position < (27*3)) {
    return("FWR")
  } else if (position >= 27*3 & position < 39*3) {
    return("CDR")
  } else if (position >= 39*3 & position < 56*3) {
    return("FWR")
  } else if (position >= 56*3 & position < 66*3) {
    return("CDR")
  } else if (position >= 66*3 & position < 105*3) {
    return("FWR")
  } else if (position >= 105*3 & position < (105+cdr3_length)*3) {
    return("CDR")
  } else {
    return("FWR")
  }
}

# Plot the data
plot_df <- bind_rows(plot_data_single)

plot_df$region <- unlist(mapply(get_region, plot_df$position, plot_df$model))
plot_df$model <- factor(plot_df$model, levels = c("selection", "neutral", "uniform"))

labeller <- c("uniform" = "Uniform Neutral", "neutral" = "Neutral (BCR + S5F)", "selection" = "Selection (BCR + S5F + Selection)")
plot <- ggplot(plot_df, aes(x = position, y = differences_per_seq, fill = region)) +
  geom_bar(stat = "identity") +
  facet_wrap(~model, nrow=3, labeller = as_labeller(labeller)) +
  labs(
    x = "Position in nucleotide sequence",
    y = "Proportion of sequences with a difference"
  ) +
  guides(fill=guide_legend(title="Region")) +
  theme(legend.position = "top") +
  theme_bw()

plot
ggsave("mutations_per_site.png", plot, width = 8, height = 10, dpi = 300)
saveRDS(plot, file="~/simble-validation/mutations_per_site_plot.rds")




