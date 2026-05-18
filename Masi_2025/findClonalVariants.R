# Kenneth B. Hoehn
# 11/19/2024
# Combine all BCR data sources, identify clonal clusterss

library(alakazam)
library(scoper)
library(dowser)
library(dplyr)
library(tigger)
library(shazam)
library(tidyr)
library(stringr)
library(airr)
library(ggtree)
library(RColorBrewer)

print(sessionInfo())

nproc=5

#Fitchner et al 2022 data
fitchner_data = read_rearrangement("data/combined_data_genotyped.tsv")

# Gianvito's new IgA mAb
mabs = read_rearrangement("data/mabs_igblast_db-pass.tsv")
gv = filter(mabs, sequence_id == "Gianvito_iga")
gv$sample = "MYG385"
gv$patient = "MuSK-1"
gv$type = "mAb"
gv$c_call = "IGHA"

# new patient data
new_samples = read_rearrangement("data/Patient_A.tsv")
new_samples$sample = new_samples$myg_id
new_samples$patient = "MuSK-1"

# combine, filter to MuSK-1 (A)
comb = bind_rows(fitchner_data, gv)
comb = bind_rows(comb, new_samples)
m1 = filter(comb, patient == "MuSK-1")

m1$sequence_id = paste0(m1$sequence_id, "-", 1:nrow(m1))

#set patient id to A
m1$patient = "A"

# correct sample times
sample_time = c(
"MYG037"=0,
"MYG069"=17,
"MYG079"=18,
"MYG108"=26,
"MYG122"=31,
"MYG135"=36,
"MYG139"=38,
"MYG166"=44,
"MYG190"=50,
"MYG245"=61,
"MYG265"=68,
"MYG274"=69,
"MYG297"=78,
"MYG385"=127,
"MYG444"=138)

table(m1$sample)

m1$old_timepoint = m1$timepoint
m1$timepoint = sample_time[m1$sample]
sum(is.na(m1$timepoint))

# remove sequences with low ATCG
m1$atcg = stringr::str_count(m1$sequence_alignment,"[ATCG]")
m1 = filter(m1, atcg >= 300)


# infer individual genotype
novel <- findNovelAlleles(m1, SampleGermlineIGHV, nproc=nproc)
novel_rows <- selectNovel(novel)
geno <- inferGenotype(m1, germline_db=SampleGermlineIGHV, novel=novel,
                    find_unmutated=TRUE)

# Save the genotype sequences to a vector
genotype_db <- genotypeFasta(geno, SampleGermlineIGHV, novel)
# Visualize the genotype and sequence counts
print(geno)
m1_genotyped <- reassignAlleles(m1, genotype_db)


# Find threshold using density method
m1_genotyped_dist <- distToNearest(m1_genotyped, 
                          sequenceColumn="junction", 
                          vCallColumn="v_call_genotyped", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=nproc)
output <- findThreshold(m1_genotyped_dist$dist_nearest, method="density")
threshold <- output@threshold
print(threshold)

# Plot distance histogram, density estimate and optimum threshold
# sanity check ATCG threshold and new timepoints
pdf("results/sanity_checks.pdf")
plot(output, title="Density Method")

ggplot(m1, aes(as.numeric(old_timepoint),timepoint)) + 
geom_point() + geom_abline() 

ggplot(m1, aes(x=atcg)) + geom_histogram() + 
geom_vline(xintercept=300)
dev.off()

# identify clones
m1s = hierarchicalClones(m1_genotyped_dist, 
  v_call="v_call_genotyped", threshold=threshold, nproc=nproc, summarize=FALSE)
write_rearrangement(m1s, "results/MuSK-A_combined_clones.tsv")

# reconstruct germlines
imgt = readIMGT("germlines/imgt/human/vdj")
m1sg = createGermlines(m1s, imgt, v_call="v_call_genotyped", nproc=nproc,
  trim_lengths=TRUE)

write_rearrangement(m1sg, "results/MuSK-A_combined_clones_germlines.tsv")
m1sg = read_rearrangement("results/MuSK-A_combined_clones_germlines.tsv")

# get mab clone
gc = filter(m1sg, grepl("Gianvito_iga",sequence_id))$clone_id
gva = filter(m1sg, clone_id %in% gc)

write_rearrangement(gva, "results/MuSK-A_mAb_clone_uncollapsed.tsv")

gvac = gva %>%
  group_by(sample, type, clone_id, c_call) %>%
  do(collapseDuplicates(., text_fields=c("sample", "type", "clone_id", "c_call"), 
    add_count=TRUE))

write_rearrangement(gvac, "results/MuSK-A_mAb_clone_collapsed.tsv")
gvac = read_rearrangement("results/MuSK-A_mAb_clone_collapsed.tsv")

# inspect alignment of clonal sequences
gvac = gvac[order(as.numeric(gvac$timepoint), decreasing=TRUE),]
fa = paste0(">",gvac$sample,"_",gvac$timepoint,"\n",gvac$sequence_alignment)
fa = c(fa, paste0(">Germline\n",gvac$germline_alignment_d_mask[1]))
writeLines(fa, con="results/MuSK-A_mAb_clone_collapsed.fa")

gvacAll = gva %>%
  do(collapseDuplicates(., text_fields=c("sample", "type", "clone_id", "c_call")))
  

gvacSeqs = gva %>%
  group_by(sample, type, clone_id, c_call) %>%
  do(collapseDuplicates(., text_fields=c("sequence_id","sample", "type", "clone_id", "c_call"), 
    add_count=TRUE))

write_rearrangement(gvacSeqs, "results/MuSK-A_mAb_clone_collapsed_seqids.tsv")

gvac <- observedMutations(gvac,
            sequenceColumn="sequence_alignment",
            germlineColumn="germline_alignment_d_mask",
            regionDefinition=IMGT_V,
            frequency=TRUE, 
            combine=TRUE,
            nproc=1)

quantile(gvac$mu_freq)
#       0%       25%       50%       75%      100% 
#0.1263158 0.1438596 0.1754386 0.2350877 0.2421053 

# confirm similar SHM values before genotyping and clonal clustering
m1s = read_rearrangement("results/MuSK-A_combined_clones.tsv")

m1s_gva <- observedMutations(filter(m1s, sequence_id %in% gvac$sequence_id),
            sequenceColumn="sequence_alignment",
            germlineColumn="germline_alignment",
            regionDefinition=IMGT_V,
            frequency=TRUE, 
            combine=TRUE,
            nproc=1)

# very high SHM even without clonal germline
hist(m1s_gva$mu_freq)
quantile(m1s_gva$mu_freq)
#       0%       25%       50%       75%      100% 
#0.1263158 0.1438596 0.1719298 0.2315789 0.2385965 


clones = formatClones(gvac, traits=c("c_call","timepoint","sample", "collapse_count"), 
  columns="type", nproc=nproc, minseq=1, collapse=FALSE)
saveRDS(clones, "results/mab_clone.rds")

trees = getTrees(clones, build="igphyml", exec="~/igphyml/src/igphyml")
saveRDS(trees, "results/mab_tree_igphyml.rds")

trees = readRDS("results/mab_tree_igphyml.rds")
trees = scaleBranches(trees)

trees$data[[1]]@data$c_call[trees$data[[1]]@data$c_call == "IGHA"] = "IgA"
trees$data[[1]]@data$c_call[grepl("Gianvito_iga",trees$data[[1]]@data$sequence_id)] = "aMu1 IgA"
timepoints = unique(trees$data[[1]]@data$timepoint)
time_scale = brewer.pal(n=length(timepoints), "RdYlBu")

timepoints = as.character(sort(as.numeric(timepoints)))
names(time_scale) = rev(timepoints)
div = max(getDivergence(trees$trees[[1]]))


pdf("results/mab_tree.pdf")
plotTrees(trees)[[1]] + geom_tippoint(aes(fill=as.numeric(timepoint), size=collapse_count), pch=21) +
scale_fill_distiller(palette="RdYlBu", name="Timepoint") + geom_tiplab() +
xlim(0, div*1.5)

plotTrees(trees, scale=10)[[1]] + geom_tippoint(aes(fill=(timepoint), size=collapse_count), pch=21) +
scale_fill_manual(values=time_scale, name="Time (months)", limits=names(time_scale)) + 
geom_tiplab(aes(label=c_call), offset=2) +
labs(size="Duplicates") +
xlim(0, div*1.5)

plotTrees(trees)[[1]] + geom_tippoint(aes(fill=as.numeric(timepoint), size=collapse_count), pch=21) +
scale_fill_distiller(palette="RdYlBu", name="Timepoint") + geom_tiplab(aes(label=sample)) +
xlim(0, div*1.5)
dev.off()

pdf("results/mab_tree_figure.pdf",width=6,height=5)
plotTrees(trees, scale=10)[[1]] + geom_tippoint(aes(fill=(timepoint), size=collapse_count), pch=21) +
scale_fill_manual(values=time_scale, name="Time (months)", limits=names(time_scale)) + 
geom_tiplab(aes(label=c_call), offset=2) +
guides(fill = guide_legend(override.aes = list(size = 3) ) ) +
labs(size="Duplicates") +
xlim(0, div*1.5)
dev.off()

trees$data[[1]]@data$time = as.numeric(trees$data[[1]]@data$timepoint)

# perform correlation test for measurable evolution
ct = correlationTest(trees, perm_type="uniform", permutations=100000, verbose=TRUE)
saveRDS(ct, "results/correlationTest.rds")

print(ct)


# compare to MFI
divs = getDivergence(trees$trees[[1]])

mfi = read.csv("data/MFI.csv")

m = match(gvac$sequence_id, names(divs))
gvac$divergence = divs[m]

m = match(gvac$sample, mfi[,1])
gvac$mfi = mfi[m,4]

mean_mu_mfi = gvac %>%
  group_by(sample, timepoint) %>%
  summarize(mu_freq = mean(divergence), mfi=unique(mfi))

reg = summary(lm(mfi ~ mu_freq, mean_mu_mfi))
print(reg)

pdf("results/MFI_vs_SHM.pdf",width=4,height=3)
ggplot(mean_mu_mfi, aes(x=mu_freq, y=mfi)) + geom_smooth(method="lm") + 
geom_point(pch=21, size=3, aes(fill=timepoint)) +
theme_bw() + 
scale_fill_manual(values=time_scale, name="Time (months)", limits=names(time_scale)) +
guides(fill = guide_legend(override.aes = list(size = 3)))+
  xlab("Mean divergence from Germline") +
  ylab("IgA binding (MFI ratio)")
dev.off()


