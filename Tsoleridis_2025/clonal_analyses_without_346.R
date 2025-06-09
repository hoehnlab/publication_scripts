# characterizing wuhan vs omicron binding b cell clones 
# R script for analyses
# ishita singh
# 04/28/2024

# only need to run once 
install.packages("tidyverse")
install.packages("scoper")
install.packages("BiocManager")
BiocManager::install("dowser")
install.packages("shazam")
install.packages("ggpubr")
install.packages("alakazam")
install.packages("gridExtra")
install.packages("BSDA")

# load packages 
library(tidyverse)
library(dplyr)
library(scoper)
library(BiocManager)
library(dowser)
library(shazam)
library(ggpubr)
library(alakazam)
library(gridExtra)
library(parallel)
library(stringr)
library(BSDA)

# setting and checking the working directory
# if this doesn't work, try session -> set wd -> choose wd 
setwd("/Users/ishita/Desktop/Research/hoehn_lab/COVID_data/")
getwd()

# reading tsv files from igb output to individual data frames 
d003O <- read.delim("../../igb_output/003O_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d003W <- read.delim("../../igb_output/003W_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d009O <- read.delim("../../igb_output/009O_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d009W <- read.delim("../../igb_output/009W_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d025O <- read.delim("../../igb_output/025O_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d025W <- read.delim("../../igb_output/025W_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d049O <- read.delim("../../igb_output/049O_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d049W <- read.delim("../../igb_output/049W_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d418O <- read.delim("../../igb_output/418O_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)
d418W <- read.delim("../../igb_output/418W_atleast-2_reheader_igblast_db-pass.tsv", header = TRUE)

# compiling to one data frame 
db_compiled <- bind_rows(d003O, d003W, d009O, d009W, d025O, d025W, d049O, d049W, d418O, d418W)

# add donor and sort as the first columns
db_compiled <- db_compiled %>%
  mutate(donor = substr(as.character(sample), 1, 3),
         sort = str_sub(as.character(sample), -1)) %>%
  select(donor, sort, duplicate_count, everything()) 

# creating a summary table of duplicate counts
db_dup_counts <- db_compiled %>%
  filter(productive == T) %>%
  group_by(donor, sort) %>%
  summarize(
    n = n(), 
    n2 = sum(duplicate_count > 2),
    n3 = sum(duplicate_count > 3),
    n4 = sum(duplicate_count > 4),
    n5 = sum(duplicate_count > 5), 
    n6 = sum(duplicate_count > 6),
    n7 = sum(duplicate_count > 7),
    n8 = sum(duplicate_count > 8),
    n9 = sum(duplicate_count > 9),
    n10 = sum(duplicate_count > 10),
    n11 = sum(duplicate_count > 11),
    n12 = sum(duplicate_count > 12),
    n13 = sum(duplicate_count > 13),
    n14 = sum(duplicate_count > 14),
    n15 = sum(duplicate_count > 15),
    n16 = sum(duplicate_count > 16),
    n17 = sum(duplicate_count > 17),
    n18 = sum(duplicate_count > 18),
    n19 = sum(duplicate_count > 19),
    n20 = sum(duplicate_count > 20),
    n21 = sum(duplicate_count > 21),
  )

# writing the summary table to a tsv file 
write.table(db_dup_counts, file = "dup_count_summary_w346.tsv", sep = "\t", row.names = FALSE)

# filtering to get cell counts equal to flow counts from Table I
db_compiled <- db_compiled %>%
  filter(productive == T) %>%
  group_by(donor) %>%
  filter(
    (donor == "025" & sort == "W" & duplicate_count > 19) |
      (donor == "025" & sort == "O" & duplicate_count > 20) |
      (donor == "418" & sort == "W" & duplicate_count > 18) |
      (donor == "418" & sort == "O" & duplicate_count > 11) |
      (donor == "049" & sort == "W" & duplicate_count > 6) |
      (donor == "049" & sort == "O" & duplicate_count > 5) |
      (donor == "009" & sort == "W" & duplicate_count > 2) |
      (donor == "009" & sort == "O" & duplicate_count > 2) |
      (donor == "003" & sort == "W" & duplicate_count > 14) |
      (donor == "003" & sort == "O" & duplicate_count > 4)
  ) 

# using SCOPer to get clonal clusters
obj_clusters <- hierarchicalClones(db_compiled, threshold=0.1, fields='donor', summarize_clones = FALSE)
db_clusters <- as.data.frame(obj_clusters)
db_clusters <- db_clusters %>%
  arrange(donor, sort)

# writing the clonal clustering results and summary to a tsv file 
write.table(db_clusters, file = "clonal_clustering_results_w346.tsv", sep = "\t", row.names = FALSE)
write.table(summary(obj_clusters), file = "clonal_clustering_summary_w346.tsv", sep = "\t", row.names = FALSE)

# getting reference files copied from Docker container's /usr/local/share/germlines/imgt/human/vdj
references <- readIMGT(dir = "imgt_human_vdj")

# using Dowser to create germlines
obj_germlines = createGermlines(db_clusters, references, nproc=4)

# writing the reconstructed germlines to a tsv file 
write.table(obj_germlines, file = "reconstructed_germlines_w346.tsv", sep = "\t", row.names = FALSE)

# using SHazaM to get somatic hypermutation frequency
db_mutations <- observedMutations(obj_germlines, sequenceColumn="sequence_alignment",
                                  germlineColumn="germline_alignment_d_mask",
                                  regionDefinition=IMGT_V,
                                  frequency=TRUE, combine=TRUE,
                                  nproc=4)

# z-test for whether the mutation frequencies are significant 
# single sample test against null hypothesis = 0 with significance threshold alpha
z_test <- function(mu_freq, alpha = 0.05) {
  z_test_result <- z.test(mu_freq, mu = 0, sigma.x = sd(mu_freq))
  z_score <- z_test_result$statistic
  p_value <- z_test_result$p.value
  significant <- ifelse(p_value < alpha, TRUE, FALSE)
  return(data.frame(z_score = z_score, p_value = p_value, significant = significant))
}

db_mu_sig <- db_mutations %>%
  group_by(donor, sort, cprimer) %>%
  summarize(
    # the conditional here is because z.test from BSDA required a min of 3 sequences
    test_results = if(n() > 2) list(z_test(mu_freq)) else list(NULL),
    z_score = if(!is.null(test_results[[1]])) signif(test_results[[1]]$z_score, 2) else NA,
    p_value = if(!is.null(test_results[[1]])) signif(test_results[[1]]$p_value, 2) else NA,
    significant = if(!is.null(test_results[[1]])) test_results[[1]]$significant else NA,
    avg_mu_freq = mean(mu_freq),
    sequence_count = n(),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  select(-test_results)

# writing the reconstructed germlines to a tsv file 
write.table(db_mu_sig, file = "shm-significance-summary-w346.tsv", sep = "\t", row.names = FALSE)

# filtering to just igg for plotting and further analysis
db_mu_sig_igg <- db_mu_sig %>%
  filter (cprimer == "IgG")

# max p value - used in results section
max_p_value <- max(db_mu_sig_igg$p_value, na.rm = TRUE)
print(max_p_value)

# custom sort labels for plot 
sort_labels <- c(
  "O" = "Omicron",
  "W" = "Wuhan"
)

# shm plot - Figure I
ggplot(db_mutations %>% filter(cprimer == "IgG"), aes(x = donor, y = mu_freq)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, size = 1) + 
  theme_bw() + 
  facet_grid(sort~., labeller = labeller(sort = sort_labels)) + 
  xlab("Donor") + 
  ylab("Mutation Frequency")

ggsave("shm-significant-w346.pdf", width = 6, height = 7)

# function to return sequences that closely match to a particular sequence in the CoV-AbDab
match_aa = function(i, subject, query){
  if(i %% 10 == 0){print(paste(i, nrow(query)))}
  temp = query[i,]
  v = getGene(temp$Heavy.V.Gene)
  j = getGene(temp$Heavy.J.Gene)
  cdrh3 = gsub("\\s+","",temp$CDRH3)
  pmatch = filter(subject, locus == "IGH" &
                    grepl(paste0(v,"\\*"), v_call) & 
                    grepl(paste0(j,"\\*"), j_call) &
                    cdr3_aal == nchar(cdrh3))
  if(nrow(pmatch) > 0){
    pmatch$aa_dist = unlist(lapply(pmatch$cdr3_aa, function(x)
      seqDist(cdrh3, x, getAAMatrix())))/nchar(cdrh3)
    pmatch$match = temp$Name
  }
  return(pmatch)
}

# map to publicly characterized covid antibodies
db_ref <- db_mutations
db_ref$cdr3_aa = substr(db_ref$junction_aa, 2, nchar(db_ref$junction_aa) - 1)
db_ref$cdr3_aal = nchar(db_ref$cdr3_aa)
covid = read.csv("new-cov.csv") # the most up to date version of the database 
covid = filter(covid, grepl("uman",Heavy.V.Gene))
names(covid)[1] = "Name"

# find matches
matches = bind_rows(parallel::mclapply(1:nrow(covid),
                                       function(x)match_aa(x, db_ref, covid),
                                       mc.cores=4))

# writing the summary table to a tsv file 
write.table(matches, file = "CoV-AbDab-matches-w346.tsv", sep = "\t", row.names = FALSE)

# processing matches
epitope = covid$Protein...Epitope
names(epitope) = covid$Name

neut = covid$Neutralising.Vs
names(neut) = covid$Name

clone_size = table(filter(db_ref, locus=="IGH")$clone_id)

matches_long = matches %>%
  filter(aa_dist <= 0.2) %>%
  select(sequence_id, donor, sort, c_call, clone_id, match, aa_dist) %>%
  mutate(clone_size=as.numeric(clone_size[clone_id]), 
         epitope = epitope[match], neutralizing=neut[match]) %>%
  arrange(desc(clone_size), clone_id)

matches_seq = matches_long %>%
  group_by(sequence_id, donor, sort, c_call, clone_id, clone_size) %>%
  summarize(
    epitopes = paste(unique(epitope),collapse=","),
    matches = paste(unique(match),collapse=","),
    neutralizing = paste(unique(neutralizing),collapse=",")) %>%
  arrange(desc(clone_size), clone_id)

write.csv(matches_long, file="matches_long_w346.csv",row.names=FALSE)
write.csv(matches_seq,  file="matches_w346.csv",row.names=FALSE)
matches_long = read.csv("matches_long_w346.csv")
matches_seq =  read.csv("matches_w346.csv")

m = match(db_ref$sequence_id, matches_seq$sequence_id)
db_ref$public = !is.na(matches_seq$matches[m])
db_ref$matches = matches_seq$matches[m]
db_ref$epitopes = matches_seq$epitopes[m]
db_ref$neutralizing = matches_seq$neutralizing[m]
db_ref_public = db_ref %>% filter(public)

write.table(db_ref_public, file = "public-matches_w346.tsv", sep = "\t", row.names = FALSE)

# summarizing matches by clone
db_matches_summary_clone <- db_ref_public %>%
  group_by(clone_id) %>%
  summarise(
    donor = paste(unique(donor)),
    sort = paste(unique(sort)),
    # this data cleaning pipeline is required to keep only unique v genes/matches etc. across match sequences
    v_genes = paste(unique(trimws(unlist(str_split(v_call, ",|;")))), collapse = ", "),
    matches = paste(unique(trimws(unlist(str_split(matches, ",|;")))), collapse = ", "),
    epitopes = paste(unique(trimws(unlist(str_split(epitopes, ",|;")))), collapse = ", "),
    neutralizing = paste(unique(trimws(unlist(str_split(neutralizing, ",|;")))), collapse = ", "),
    isotypes = paste(unique(cprimer), collapse = ", "),
    avg_mu_freq = signif(mean(mu_freq, na.rm = TRUE), 2),
    sequence_count = n(), 
    match_count = str_count(matches, ",") + 1 
  )


write.table(db_matches_summary_clone, file = "matches_summary_clone_w346.tsv", sep = "\t", row.names = FALSE)

# msa
# get sequences of chosen clone, chosen because it is the largest that is 0, W neutralizing
# recognizing it by the v-genes instead of clone id which is randomized
chosen_clone = db_matches_summary_clone %>%
  filter(grepl("Omicron", neutralizing)) %>%
  arrange(desc(sequence_count)) %>%
  slice(1) %>%
  pull(clone_id)

print(chosen_clone)
c_seqs = filter(db_ref_public, clone_id == chosen_clone)

# get matches corresponding to the sequences of the chosen clone 
chosen_matches = c_seqs$matches

# function to replace non-leading dots in germline and sequence alignment
replace_dots <- function(seq) {
  # temp 
  leading_dots <- regmatches(seq, regexpr("^\\.+", seq))
  # remove the leading dots from the sequence
  remaining_seq <- sub("^\\.+", "", seq)
  # remove all the dots in the remaining sequence
  cleaned_seq <- gsub("\\.", "", remaining_seq)
  # add back the leading dots to the cleaned sequence
  result <- paste0(leading_dots, cleaned_seq)
  return(result)
}

# germline AA
glines = table(substr(c_seqs$germline_alignment,1,312))
gline_max = names(glines[which.max(glines)])
gline_max_mod <- replace_dots(gline_max)
AA = strsplit(translateDNA(gline_max_mod),split="") # remove the replace dot with dot
# add a if else to remove X's not in the beginning - "\\.",""
# for germline and sequence
heavy_aa = bind_rows(lapply(1:length(AA),function(x){
  tibble(sequence_id="GermlineV",position=1:length(AA[[x]]),AA=AA[[x]])
}))

# sequence AA
seq_mod <- replace_dots(c_seqs$sequence_alignment)
# remove the first X that causes a frameshift in the alignment 
seq_mod <- sub("^\\.\\.\\.", "", seq_mod)
AA = strsplit(translateDNA(seq_mod),split="")
heavy_aa = bind_rows(heavy_aa,lapply(1:length(AA),function(x){
  tibble(sequence_id=c_seqs$sequence_id[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

# matches AA
ids = strsplit(chosen_matches,split=",")[[1]]
fc = filter(covid, Name %in% ids & VHorVHH != "ND")
AA = strsplit(fc$VHorVHH,split="")
heavy_aa = bind_rows(heavy_aa,lapply(1:length(AA),function(x){
  tibble(sequence_id=fc$Name[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

heavy_aa$sequence_id = factor(heavy_aa$sequence_id,
                              levels=unique(heavy_aa$sequence_id))

# AA palette for plotting
aapalette = c(
  "H" = "#a6cee3",
  "K" = "#a6cee3",
  "R" = "#a6cee3",
  "D" = "#fb9a99",
  "E" = "#fb9a99",
  "S" = "#1f78b4",
  "T" = "#1f78b4",
  "N" = "#1f78b4",
  "Q" = "#1f78b4",
  "A" = "white",
  "V" = "white",
  "L" = "white",
  "I" = "white",
  "M" = "white",
  "F" = "#cab2d6",
  "Y" = "#cab2d6",
  "W" = "#cab2d6",
  "P" = "#fb9a99",
  "G" = "#fb9a99",
  "C" = "#fdbf6f",
  "B" = "grey",
  "Z" = "grey",
  "X" = "grey")

# save Figure II
pdf("msa-w346.pdf",width=9.5,height=5,useDingbats=FALSE)
print(
  ggplot(heavy_aa, aes(x=position,y=sequence_id, fill=AA, label=AA)) + geom_tile(width=1.1, height=1.1) + 
    geom_text(size=1.7) + 
    theme_bw() + theme(legend.position = "none") + xlab("Position") + ylab("") +
    scale_x_continuous(expand=c(0,0.1)) +
    theme(text = element_text(size = 9), axis.text =element_text(size = 9)) +
    scale_fill_manual(values=aapalette) +
    scale_y_discrete(limits=rev)+
    ggtitle(paste0("Clone ",chosen_clone,
                   " heavy chain matches"))
)
dev.off()

# IgA-IgG ratio for BCR data 
db_ratio <- db_mutations

db_ratio <- db_ratio %>%
  group_by(donor, sort, cprimer) %>%
  summarize(cprimer_count = n(), .groups = 'drop')

db_ratio <- db_ratio %>%
  pivot_wider(names_from = cprimer, values_from = cprimer_count, values_fill = list(count = 0))

db_ratio <- db_ratio %>%
  mutate(IgA_IgG = IgA / IgG)  %>%
  replace_na(list(IgA_IgG = 0))  %>%
  replace_na(list(IgA = 0))

# IgA-IgG ratio for flow data
flow_stats <- read.csv("flow-stats.csv", header = TRUE)

patients <- c("003", "009", "25", "49", "418")
sorts <- c("wuhan", "Wuhan", "Omicron", "omicron")

flow_stats_filtered <- flow_stats %>%
  filter(str_detect(Data.Set, paste(patients, collapse="|")),
         Input.Gate == "[Spike++]",
         Gate %in% c("IgA+IgD-IgG-IgM-", "IgA-IgD-IgG+IgM-"))

flow_stats_filtered <- flow_stats_filtered %>%
  mutate(
    donor = str_extract(Data.Set, paste(patients, collapse = "|")),
    sort = str_extract(Data.Set, paste(sorts, collapse = "|"))
  )

flow_stats_filtered <- flow_stats_filtered %>%
  mutate(
    sort = case_when(
      str_detect(sort, "wuhan|Wuhan") ~ "W",
      str_detect(sort, "Omicron|omicron") ~ "O"
    )
  )

flow_stats_filtered <- flow_stats_filtered %>%
  mutate(donor = str_pad(donor, 3, pad = "0"))

flow_stats_filtered <- flow_stats_filtered %>%
  select(donor, sort, Input.Gate, Gate, X.Gated)

flow_stats_filtered <- flow_stats_filtered %>%
  mutate(X.Gated = as.numeric(X.Gated))

flow_stats_filtered <- flow_stats_filtered %>%
  group_by(donor, sort) %>%
  summarise(
    IgA = sum(X.Gated[Gate == "IgA+IgD-IgG-IgM-"], na.rm = TRUE),
    IgG = sum(X.Gated[Gate == "IgA-IgD-IgG+IgM-"], na.rm = TRUE)
  ) %>%
  mutate(
    IgA_IgG = IgA / IgG
  )

db_ratio <- db_ratio %>%
  left_join(flow_stats_filtered %>% 
              select(donor, sort, IgA_IgG) %>%
              rename(IgA_IgG_flow = IgA_IgG),
            by = c("donor", "sort"))

# create a new row for totals and mean - values mentioned in results section
totals <- db_ratio %>%
  summarize(IgA = sum(IgA, na.rm = TRUE),
            IgG = sum(IgG, na.rm = TRUE),
            IgA_IgG = mean(IgA_IgG, na.rm = TRUE), 
            IgA_IgG_flow = mean(IgA_IgG_flow, na.rm = TRUE))

totals_row <- data.frame(
  IgA = totals$IgA,
  IgG = totals$IgG,
  IgA_IgG = totals$IgA_IgG, 
  IgA_IgG_flow = totals$IgA_IgG_flow
)
db_ratio <- bind_rows(db_ratio, totals_row)

db_ratio <- db_ratio %>%
  mutate(IgA_IgG = signif(IgA_IgG, 2),
         IgA_IgG_flow = signif(IgA_IgG_flow, 2))

write.table(db_ratio, file = "isotype-ratio-w346.tsv", sep = "\t", row.names = FALSE)

# end of script