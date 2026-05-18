library(dplyr)


samples <- c(
  "Bcl2",
  "CBP1",
  "CBP2",
  "CMBP3",
  "CMBP7",
  "EBMP2",
  "EBMP3",
  "EZBP1",
  "EZBP6",
  "MCDP2",
  "MCDP4",
  "MCDP6",
  "SN2F",
  "SN2M",
  "WP2",
  "WP6",
)

data <- data.frame()
neb_igblast <- data.frame()
folder <- "NEB_output/igblast_output_assembled_first"
for (i in seq(1, length(samples))) {
  filename = paste(folder, samples[i], "-final_collapse-unique_atleast-2_igblast_db-pass.tsv", sep="")
  if (!file.exists(filename)) {
    next

  }
  curr_data <- airr::read_rearrangement(filename)
  curr_data$sample_id <- samples[i]
  neb_igblast <- rbind(neb_igblast, curr_data)
}
neb_igblast$data_source <- "NEB"

miseq_igblast <- data.frame()
folder <- "MiSeq_output/igblast_output_assembled_first/"
for (i in seq(1, length(samples))) {
  filename = paste(folder, samples[i], "_atleast-2_igblast_db-pass.tsv", sep="")
  if (!file.exists(filename)) {
    next
  }
  curr_data <- airr::read_rearrangement(filename)
  curr_data$sample_id <- samples[i]
  miseq_igblast <- rbind(miseq_igblast, curr_data)
}

miseq_igblast$data_source <- "miseq"

combined_igblast <- bind_rows(neb_igblast, miseq_igblast)
for (i in seq(1, length(samples))) {
  curr_sample_combined_igblast <- combined_igblast %>% filter(sample_id == samples[i])
  write.table(curr_sample_combined_igblast,
            file = paste0("combined_igblast/", samples[i], "-final_collapse-unique_atleast-2_igblast.tsv"),
            sep="\t", row.names = FALSE)
}
