# Kenneth B. Hoehn
# 9/19/25
# Make figure panels from HIV analysis

library(dowser)
library(airr)
library(ggtree)
library(dplyr)
library(shazam)
library(tidyr)
library(treeio)
library(viridis)
library(ggpubr)

print(sessionInfo())

getDiffPoint = function(tree, node){
	type = filter(tree@data, !!node==node)$location
	edge = tree@phylo$edge[tree@phylo$edge[,2] == node,]
	if(length(edge) == 0){
		return(tibble(diffnode=node, type="root", height=filter(tree@data, !!node==node)$height))
	}
	if(!is.null(nrow(edge))){
		stop("weird")
	}
	parent = as.character(edge[1])
	parent_type = filter(tree@data, node==parent)$location
	parent_height = filter(tree@data, node==parent)$height

	if(parent_type == type){
		return(getDiffPoint(tree, parent))
	}else{
		return(tibble(diffnode=parent, type=parent_type, height=parent_height))
	}
}

getDiffPoints = function(tree){
	diffpoints = tibble()
	for(l in tree@phylo$tip.label){
		#print(l)
		d = filter(tree@data, node == which(tree@phylo$tip.label == l))
		df = getDiffPoint(tree, which(tree@phylo$tip.label == l))
		temp = tibble(tip=l, tip_type=d$location, tip_height=d$height)
		diffpoints = bind_rows(diffpoints, bind_cols(temp, df))
	}
	diffpoints$height = as.numeric(diffpoints$height)
	diffpoints$tip_height = as.numeric(diffpoints$tip_height)
	diffpoints
}

patients = c("HIV1", "HIV2", "HIV3")

patient_palette = c("HIV3"="orange", "HIV1"="darkgreen",
	"HIV2"="navyblue")

# make MRCA timing plots
models = c(
	"strict",
	"typelinked-irrev",
	"typelinked-eo-est",
	"ucld"
		)

types = RColorBrewer::brewer.pal(3,"Set1")
names(types) = c("GC", "MBC", "Ambig.")

t = read.csv("data/allclock_results_clustered.csv")
me = filter(t, p < 0.05)
liao = filter(me, study=="Liao_2013")
gclock = mean(liao$slope)

pdf("figure/clock_rate.pdf",width=2, height=1.5)
ggplot(liao, aes(x=slope)) + geom_histogram() +
theme_bw() + geom_vline(xintercept=gclock, linetype="dashed", color="red") +
xlab("Slope (mu/site/week)") + ylab("Clones") +
theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1))
dev.off()



allheights = tibble()
alltrees = tibble()
for(model in models){
	for(patient in patients){
		trees = readRDS(paste0("intermediates/",patient,"_v003_",model,".trees.rds"))
		heights = bind_rows(lapply(trees$parameters, function(x){
			f = filter(x, item=="TreeHeight")
			f[,!names(f) == "item"] = apply(f[,!names(f) == "item"],2,as.numeric)
			f
		}))
		heights$clone_id = trees$clone_id
		heights$patient = patient
		heights$model = model

		mrca_locations = sapply(trees$trees, function(x){
			node = ape::getMRCA(x@phylo, tip=x@phylo$tip.label)
			filter(x@data, !!node==node)$location
		})
		mrca_probs = sapply(trees$trees, function(x){
			node = ape::getMRCA(x@phylo, tip=x@phylo$tip.label)
			as.numeric(filter(x@data, !!node==node)$location.prob)
		})

		heights$location = mrca_locations
		heights$location.prob = mrca_probs

		allheights = bind_rows(allheights, heights)
		trees$model=model
		trees$patient = patient
		alltrees = bind_rows(alltrees, trees)
	}
}

ignore = c("traitRates", "typeLinkedRates", "freqParameter", "clockRate", "traitfrequencies", "geneticClockRate", "rateCategories")

params = alltrees$parameters
for (regex in ignore) {
    params = lapply(params, function(x) {
        filter(x, !grepl(regex, item))
    })
}
alltrees$converged = sapply(params, function(x) sum(x$ESS[!x$item %in% 
        ignore] < 200, na.rm = TRUE) == 0)

table(alltrees$converged, alltrees$model)

alltrees$patient_clone = paste0(alltrees$patient,"_",alltrees$clone_id)
nc_clones = unique(filter(alltrees, !converged)$patient_clone)

# remove nonconverged clones
allheights$patient_clone = paste0(allheights$patient, "_",allheights$clone_id)
allheights = filter(allheights, !patient_clone %in% nc_clones)
alltrees = filter(alltrees, !patient_clone %in% nc_clones)

# label which model is TyCHE
allheights$model[allheights$model == "typelinked-eo-est"] = "ztyche"
alltrees$model[alltrees$model == "typelinked-eo-est"] = "ztyche"


mrca_locations = select(allheights, patient, model, clone_id, location)
for(patient in patients){
	gdt = readRDS(paste0("intermediates/",patient,"_gd_parsimony_trees.rds"))
	for(cloneid in gdt$clone_id){
		t = filter(gdt, clone_id == cloneid)$trees[[1]]
		node = ape::getMRCA(t, tip=t$tip.label)
		mrca_locations = bind_rows(mrca_locations, tibble(patient=patient, clone_id=cloneid, 
			location=t$state[node], model="parsimony"))
	}
}
mrca_locations$patient_clone = paste0(mrca_locations$patient, "_",mrca_locations$clone_id)
mrca_locations = filter(mrca_locations, !patient_clone %in% nc_clones)



pdf("results/UCA_dates_TL.pdf", width=3.0,height=2)
ggplot(filter(allheights, model=="typelinked-irrev"), 
	aes(x=0-mean, y=patient)) + 
annotate("rect", xmin = -52/12*5, xmax = -52/12*7, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["HIV3"]) +
annotate("rect", xmin = -52, xmax = -78, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["HIV1"]) +
annotate("rect", xmin = -104, xmax = -156, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["HIV2"]) +
geom_boxplot(outliers=FALSE, aes(fill=patient), alpha=0.3) +
#geom_jitter(height=0.1, width=0,size=0.1)+
facet_wrap(.~model, scales="free_x") +
theme_bw() +
xlab("UCA date (weeks before sample)") +
ylab("Patient") +
scale_fill_manual(values=patient_palette) +
guides(fill="none") +
scale_y_discrete(limits=rev)

ggplot(filter(allheights, model=="ztyche"), 
	aes(x=0-mean, y=patient)) + 
annotate("rect", xmin = -52/12*5, xmax = -52/12*7, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["HIV3"]) +
annotate("rect", xmin = -52, xmax = -78, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["HIV1"]) +
annotate("rect", xmin = -104, xmax = -156, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["HIV2"]) +
geom_boxplot(outliers=FALSE, aes(fill=patient), alpha=0.3) +
#geom_jitter(height=0.1, width=0,size=0.1)+
facet_wrap(.~model, scales="free_x") +
theme_bw() +
xlab("UCA date (weeks before sample)") +
ylab("Patient") +
scale_fill_manual(values=patient_palette) +
guides(fill="none") +
scale_y_discrete(limits=rev)
dev.off()


# UCA locations
# Probabilities
allheights$gcprob = allheights$location.prob
allheights$gcprob[allheights$location == "other"] = 1 - allheights$gcprob[allheights$location == "other"]
pdf("results/UCA_probabilities.pdf", width=4,height=5)
ggplot(filter(allheights, model!="typelinked-irrev"), aes(x=gcprob, fill=model)) + 
geom_histogram(position='identity', alpha=0.8) + 
geom_vline(xintercept=0.5) +
theme_bw() +
facet_grid(patient~., scales="free_y") +
xlab("Prob(UCA = GC)") + ylab("Clones")
dev.off()

comb_count = mrca_locations %>%
		group_by(patient, model, location) %>%
		summarize(n=n()) %>%
		mutate(freq = n/sum(n))

comb_count$location = gsub("other", "MBC", comb_count$location)
comb_count$location = gsub("germinal_center", "GC", comb_count$location)
comb_count$location = gsub("germinal-center", "GC", comb_count$location)
comb_count$location = gsub("GC,MBC", "Ambig.", comb_count$location)

# discrete locations
pdf(paste0("results/UCA_locations.pdf"),width=4, height=2)
print(
ggplot(filter(comb_count, model!="typelinked-irrev"), aes(y=freq, fill=location, 
	x=patient)) + geom_bar(stat="identity", color="black") +
theme_bw() +
xlab("") +
ylab("Proportion of clones") +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
facet_grid(.~model) +
scale_fill_manual(values=types)
)
print(
ggplot(comb_count, aes(y=n, fill=location, 
	x=model)) + geom_bar(stat="identity", color="black") +
theme_bw() +
xlab("") +
ylab("Clones") +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
facet_grid(.~patient) +
scale_fill_manual(values=types)
)
dev.off()

# compare relative rate params
alltrees$GC_MBC = sapply(alltrees$parameters, function(x)filter(x, grepl("relativeGeoRates.*1",item))$mean)
alltrees$MBC_GC = sapply(alltrees$parameters, function(x)filter(x, grepl("relativeGeoRates.*2",item))$mean)
alltrees$GM_MG = alltrees$GC_MBC/alltrees$MBC_GC

g = alltrees$GC_MBC
m = alltrees$MBC_GC

relativeGeos = alltrees %>%
	filter(model %in% c("strict", "ztyche","ucld")) %>%
	select(clone_id, patient, model, GC_MBC, MBC_GC)

gc_mbc = select(relativeGeos, clone_id, patient, model, GC_MBC)
gc_mbc$Class = "GC->MBC"
gc_mbc = mutate(gc_mbc,rate=GC_MBC)

mbc_gc = select(relativeGeos, clone_id, patient, model, MBC_GC)
mbc_gc$Class = "MBC->GC"
mbc_gc = mutate(mbc_gc,rate=MBC_GC)

rates = bind_rows(gc_mbc, mbc_gc)

ratepal = types
ratepal["GC->MBC"] = types["GC"]
ratepal["MBC->GC"] = types["MBC"]

pdf("results/relativeRates.pdf", width=6.0,height=2.25)
ggplot(rates, aes(x=Class, y=rate, fill=Class)) + 
geom_boxplot(outlier.shape=NA) +
geom_point(size=0.5) +
geom_line(aes(group=clone_id), size=0.5,alpha=0.5) +
stat_compare_means(method = "wilcox.test", paired=TRUE, label="p.format", size=2.5) +
facet_grid(.~paste0(model,patient)) +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
ylab("Relative rate") + xlab("") +
ylim(0,0.23) +
scale_fill_manual(values=ratepal)
dev.off()

# compare trees across programs
patient = "HIV1"
strict = readRDS(paste0("intermediates/",patient,"_v003_strict.trees.rds"))
tl = readRDS(paste0("intermediates/",patient,"_v003_typelinked-eo-est.trees.rds"))
gd = readRDS(paste0("intermediates/",patient,"_gd_parsimony_trees.rds"))
uc = readRDS(paste0("intermediates/",patient,"_v003_ucld.trees.rds"))

clone = "107368"
print(paste0(patient,"_",clone) %in% nc_clones)

strict = filter(strict, clone_id == clone)
tl = filter(tl, clone_id == clone)
gd = filter(gd, clone_id == clone)
uc = filter(uc, clone_id == clone)

gd$data[[1]]@data$location = gsub("-","_",gd$data[[1]]@data$location)
gd$data[[1]]@data$time = 0
gd$trees[[1]]$state = gsub("-","_",gd$trees[[1]]$state)

table(tl$trees[[1]]@data$location)
table(gd$trees[[1]]$state)

tl$trees[[1]]@data$location[tl$trees[[1]]@data$location =="germinal_center"] = "GC"
tl$trees[[1]]@data$location[tl$trees[[1]]@data$location =="other"] = "MBC"
tl$trees[[1]]@data$location[tl$trees[[1]]@data$location =="other+germinal_center"] = "Ambig."

strict$trees[[1]]@data$location[strict$trees[[1]]@data$location =="germinal_center"] = "GC"
strict$trees[[1]]@data$location[strict$trees[[1]]@data$location =="other"] = "MBC"
strict$trees[[1]]@data$location[strict$trees[[1]]@data$location =="germinal_center+other"] = "Ambig."

uc$trees[[1]]@data$location[uc$trees[[1]]@data$location =="germinal_center"] = "GC"
uc$trees[[1]]@data$location[uc$trees[[1]]@data$location =="other"] = "MBC"
uc$trees[[1]]@data$location[uc$trees[[1]]@data$location =="germinal_center+other"] = "Ambig."

gd$data[[1]]@data$location[gd$data[[1]]@data$location == "germinal_center"] = "GC"
gd$data[[1]]@data$location[gd$data[[1]]@data$location == "other"] = "MBC"
gd$trees[[1]]$state[gd$trees[[1]]$state == "germinal_center"] = "GC"
gd$trees[[1]]$state[gd$trees[[1]]$state == "other"] = "MBC"
gd$trees[[1]]$state[gd$trees[[1]]$state == "germinal_center,other"] = "Ambig."

tree_height = max(as.numeric(tl$trees[[1]]@data$height))
tree_height = filter(tl$parameters[[1]], item=="TreeHeight")$mean

tlp = ggtree(tl$trees[[1]]) + 
geom_nodepoint(aes(fill=location),size=2, pch=21) +
geom_tippoint(aes(fill=location),size=2, pch=21) +
scale_fill_manual(values=types) +
labs(fill="Day", shape="Type") +
geom_treescale(width=10)

slp = ggtree(strict$trees[[1]]) + 
geom_nodepoint(aes(fill=location),size=2, pch=21) +
geom_tippoint(aes(fill=location),size=2, pch=21) +
scale_fill_manual(values=types) +
labs(fill="Day", shape="Type") +
geom_treescale(width=10)

ucp = ggtree(uc$trees[[1]]) + 
geom_nodepoint(aes(fill=location),size=2, pch=21) +
geom_tippoint(aes(fill=location),size=2, pch=21) +
scale_fill_manual(values=types) +
labs(fill="Day", shape="Type") +
geom_treescale(width=10)

gdp = plotTrees(gd, title=FALSE)[[1]] + 
geom_nodepoint(aes(fill=gd$trees[[1]]$state),pch=21,size=2) +
geom_tippoint(aes(fill=location),size=2,pch=21) +
scale_fill_manual(values=types) +
labs(fill="Day", shape="Type") 

pdf("figure/hiv_trees_color.pdf", width=7.5,height=2.5,useDingbats=FALSE)
gridExtra::grid.arrange(gdp, slp, ucp, tlp, ncol=4)
dev.off()


# MBC differentiation timing
alldiffs = tibble()
meandiffs = tibble()
ignore = c("traitRates", "typeLinkedRates", "freqParameter", "clockRate", "traitfrequencies", "geneticClockRate")
for(patient in patients){
	trees = readRDS(paste0("intermediates/",patient,"_v003_gc-origin_subset3_trees.rds"))

	params = trees$parameters
	for (regex in ignore) {
	    params = lapply(params, function(x) {
	        filter(x, !grepl(regex, item))
	    })
	}
	trees$converged = sapply(params, function(x) sum(x$ESS[!x$item %in% 
	        ignore] < 200, na.rm = TRUE) == 0)

	print(sum(!trees$converged))
	trees = filter(trees, converged)

	diffs = bind_rows(lapply(trees$trees, function(x){
		d = getDiffPoints(x)
		d$clone_id = x@info$name
		d
	}))

	heights = sapply(trees$parameters, function(x)filter(x, item=="TreeHeight")$mean)
	names(heights) = trees$clone_id

	diffs$relative_height = diffs$height/heights[diffs$clone_id]
	diffs$tree_height = heights[diffs$clone_id]
	diffs$patient = patient
	alldiffs = bind_rows(alldiffs, diffs)

	mean_relative = diffs %>%
		group_by(clone_id, tip_type) %>%
		summarize(relative_height=mean(relative_height))
	mean_relative$patient = patient
	meandiffs = bind_rows(meandiffs, mean_relative)	
}

pdf("results/diff_times_pretty_all.pdf", width=3.5,height=1.75)
ggplot(filter(meandiffs,tip_type != "GC"), aes(x=1-relative_height, y=tip_type)) + 
geom_boxplot(outlier.shape=NA) + 
geom_jitter(aes(color=patient), height=0.1,width=0) +
theme_bw() +
xlab("Normalized time from UCA") + ylab("Tip celltype") +
scale_color_manual(values=patient_palette)

ggplot(alldiffs, aes(x=patient, y=tree_height)) + geom_boxplot()

ggplot(filter(meandiffs,tip_type != "GC"), aes(x=1-relative_height, y=patient, color=tip_type)) + 
geom_boxplot(outlier.shape=NA) + 
theme_bw() +
xlab("Normalized time from UCA") + ylab("Tip celltype") 

ggplot(filter(meandiffs,tip_type != "GC"), aes(y=1-relative_height, x=patient, color=tip_type)) + 
geom_boxplot(outlier.shape=NA) + 
theme_bw() +
ylab("Normalized time from UCA") + xlab("Tip celltype") 
dev.off()

