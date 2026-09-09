# Kenneth B. Hoehn
# 9/19/25
# Make figures from pseudomonas data

library(dowser)
library(dplyr)
library(tidyr)
library(ggtree)
library(ape)
library(phangorn)
library(ggtree)
library(treeio)
library(viridis)

print(sessionInfo())


# read in results
t = readRDS("intermediates/pseudomonas_all.rds")
s = readRDS("intermediates/pseudomonas_sc_all.rds")
u = readRDS("intermediates/pseudomonas_ucld_all.rds")

dates = c(0, 9, 50, 66, 76,101)
labs = c(
"Admission",
"Intubation",
"Sepsis",
"Ab completion",
"Hypermutator detect",
"last"
	)

shapes = c("H"=24, "N"=21, "Both"=23)

types = RColorBrewer::brewer.pal(3,"Set1")
names(types) = c("H", "N", "Both")


pdf("results/hypermutator_densitree.pdf",width=5,height=2)
plotTrees(t, densitree=TRUE, layout="slanted", alpha=0.05,
    show_occupancy=TRUE, palette=c("H"="red","H-N"="purple","N"="blue"), 
    tips="location", scale=FALSE)[[1]] +
geom_vline(xintercept=dates - 101, linetype="dashed") +
geom_text(label=labs[1],x=dates[1]-101, y=10, angle=90,check_overlap=TRUE, color="black")+
geom_text(label=labs[2],x=dates[2]-101, y=10, angle=90,check_overlap=TRUE, color="black")+
geom_text(label=labs[3],x=dates[3]-101, y=10, angle=90,check_overlap=TRUE, color="black") +
geom_text(label=labs[4],x=dates[4]-101, y=10, angle=90,check_overlap=TRUE, color="black") +
geom_text(label=labs[5],x=dates[5]-101, y=10, angle=90,check_overlap=TRUE, color="black")
dev.off()


hmrca = getMRCA(t$trees[[1]]@phylo, 
  tip=t$trees[[1]]@phylo$tip.label[grepl("_H_", t$trees[[1]]@phylo$tip.label)])
101 - as.numeric(filter(t$trees[[1]]@data, node==hmrca)$height)
# H MRCA date
#[1] 46.52153

101 - max(as.numeric(t$trees[[1]]@data$height))
# MRCA date
#[1] 14.06611

# get genetic distance tree
fasta = readFasta("data/hypermutator_paper_alignment_L.fasta")
seqs = strsplit(unlist(fasta), split="")
data = phangorn::phyDat(ape::as.DNAbin(t(as.matrix(dplyr::bind_rows(seqs)))))
dm <- phangorn::dist.ml(data)
treeNJ <- ape::multi2di(phangorn::NJ(dm), random = FALSE)
treeNJ$edge.length[treeNJ$edge.length < 0] <- 0
pml <- phangorn::pml(ape::unroot(treeNJ), data = data)
fit <- tryCatch(phangorn::optim.pml(pml, model = "HKY", 
    optNni = TRUE, optQ = TRUE, optGamma = FALSE, rearrangement = "NNI", 
    control = phangorn::pml.control(epsilon = 1e-08, 
        maxit = 10, trace = 0)), error = function(e) e)

namesplit = strsplit(names(seqs), split="_")

data = tibble(taxa=names(seqs))
data$timepoint = as.numeric(sapply(namesplit, function(x)x[4])) + 46
data$hypermutator = sapply(namesplit, function(x)x[2])

pdf("results/distance_tree_daylight_HKY.pdf", width=6,height=2)
stree = fit$tree
stree$edge.length = stree$edge.length*length(seqs[[1]])
ggtree(stree, layout="daylight") %<+% data +
geom_tippoint(aes(fill=timepoint,pch=hypermutator),size=2) +
scale_fill_viridis(direction=1, limits=c(0, end=101)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Strain") +
geom_treescale(width=25)
dev.off()

pdf("results/distance_tree_daylight_big_HKY.pdf", width=60,height=20)
stree = fit$tree
stree$edge.length = stree$edge.length*length(seqs[[1]])
ggtree(stree, layout="daylight") %<+% data +
geom_tippoint(aes(fill=timepoint,pch=hypermutator),size=2) +
scale_fill_viridis(direction=1, limits=c(0, end=101)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Strain") +
geom_treescale(width=1)
dev.off()

# parameter plots
hl = readRDS("intermediates/htree.rds")$parameters[[1]]
hl$run = "HL"
nl = readRDS("intermediates/ntree.rds")$parameters[[1]]
nl$run = "NL"
tl = t$parameters[[1]]
tl$run = "TL"
sc = s$parameters[[1]]
sc$run = "SC"
uc = u$parameters[[1]]
uc$run = "UC"
comb = bind_rows(hl,nl,tl,sc,uc)

heights = filter(comb, item %in% c("Tree.height", "TreeHeight"))

heights$strain = "Both"
heights$strain[heights$run == "HL"] = "H"
heights$strain[heights$run == "NL"] = "N"

rates = filter(comb, grepl("typeLinkedRates", item))
rates = bind_rows(rates, filter(comb, grepl("geneticClockRate", item)))

rates$strain = "Both"
rates$strain[rates$run == "HL"] = "H"
rates$strain[rates$run == "NL"] = "N"
rates$strain[rates$run == "TL" & rates$item == "typeLinkedRates.1"] = "H"
rates$strain[rates$run == "TL" & rates$item == "typeLinkedRates.2"] = "N"

hr = bind_rows(heights, rates)

rates$run = factor(rates$run, 
    levels=rev(c("NL","HL","TL","SC","UC")))
heights$run = factor(heights$run, 
    levels=rev(c("NL","HL","TL","SC","UC")))

size = 2
pdf("results/treeheights.pdf", width=3,height=2.5)
g1 = ggplot(heights, aes(y=run, x=mean, xmax=X95.HPDup, xmin=X95.HPDlo,
	fill=strain, shape=strain)) +
geom_vline(xintercept=101, linetype="dashed") + 
geom_errorbar(width=0.1) +
geom_point(size=size) + 
scale_x_log10() +
theme_bw() +
scale_shape_manual(values=shapes) +
xlab("Tree height (time to MRCA)")+
scale_fill_manual(values=types)

g2 = ggplot(rates, aes(y=run, x=mean, xmax=X95.HPDup, xmin=X95.HPDlo,
	fill=strain, pch=strain)) +
geom_errorbar(width=0.1) +
geom_point(size=size) + 
scale_x_log10() +
theme_bw() +
xlab("Clock rate (mu/site/day)") +
scale_shape_manual(values=shapes)+
scale_fill_manual(values=types)

gridExtra::grid.arrange(g1, g2, ncol=1)
dev.off()

heights %>% 
    group_by(run) %>%
    summarize(101-mean, (101-mean)/365)
#1 UC         -4270.            -11.7   
#2 SC        -15110.            -41.4   
#3 TL            14.1             0.0387
#4 HL            40.3             0.110 
#5 NL            35.4             0.0970

101 - as.numeric(filter(tl, item=="TreeHeight")$X95.HPDlo)
# 24.1708
101 - as.numeric(filter(tl, item=="TreeHeight")$X95.HPDup)
# 3.0795

rates %>% 
    group_by(run) %>%
    select(mean)
#1 TL    0.00369  
#2 TL    0.000238 
#3 HL    0.00369  
#4 NL    0.000238 
#5 SC    0.0000546
#6 UC    0.00325  

hmrca = ape::getMRCA(t$trees[[1]]@phylo, 
    tip=filter(t$data[[1]]@data, location=="H")$sequence_id)

hnode = filter(t$trees[[1]]@data, node==hmrca)
101 - as.numeric(hnode$height)
# 46.52153
101 - as.numeric(hnode$height_0.95_HPD[[1]])
# 50.49623 42.52766


