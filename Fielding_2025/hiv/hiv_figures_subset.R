# Kenneth B. Hoehn
# 5/18/26
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

patients = c("HIV1", "HIV2", "HIV3")

patient_palette = c("H3"="orange", "H1"="darkgreen",
	"H2"="navyblue")

models = c(
	"typelinked-irrev"
		)

types = RColorBrewer::brewer.pal(3,"Set1")
names(types) = c("GC", "MBC", "Ambig.")

types["Ambig."] = "grey"
types = types[c("GC","Ambig.","MBC")]

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

# read in output files
allheights = tibble()
alltrees = tibble()
for(model in models){
	for(patient in patients){
		if(model %in% c("typelinked-irrev","typelinked-eo-est","typelinked-eo-estflip")){
			trees = readRDS(paste0("intermediates/",patient,"_v009op4_",model,"_fullposterior.trees.rds"))
		}else{
			trees = readRDS(paste0("intermediates/",patient,"_v009op4_",model,"_fullposterior.trees.rds"))
		}
		print(patient)
		print(table(unlist(sapply(trees$data, function(x)x@data$subset))))

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

		if(model %in% c("typelinked-irrev","typelinked-eo-est")){
			plots = plotTrees(trees, tips="location", nodes=TRUE, scale=52, show_occupancy=TRUE,
			palette=c("germinal_center"="red","other+germinal_center"="purple","other"="blue"))
			treesToPDF(plots,paste0("results/",patient,"_v009_",model,".trees.pdf"))
		}else{
			plots = plotTrees(trees, tips="location", nodes=TRUE, scale=52, 
			palette=c("other"="blue","other+germinal_center"="purple","germinal_center"="red"))
			treesToPDF(plots,paste0("results/",patient,"_v009_",model,".trees.pdf"))
		}
	}
}

ignore = c("traitRates", "typeSwitchClockRate", "typeLinkedRates", "freqParameter", 
	"clockRate", "traitfrequencies", "geneticClockRate", "rateCategories")

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

# compare tyche to starting things with other and diff order
alltrees$posterior = sapply(alltrees$parameters, function(x)filter(x, item=="posterior")$mean)
alltrees$likelihood = sapply(alltrees$parameters, function(x)filter(x, item=="likelihood")$mean)
allheights$patient = gsub("IV","",allheights$patient)

saveRDS(alltrees,"results/alltrees_fullposterior.rds")
alltrees = readRDS("results/alltrees_fullposterior.rds")

for(i in 1:nrow(alltrees)){
	alltrees$trees[[i]]@info$trees_with_traits_posterior = alltrees$trees[[i]]@info$tree_posterior
}

# make UCA timing plots
pdf("results/UCA_dates_TL_subset.pdf", width=3.0,height=2)
ggplot(filter(allheights, model=="typelinked-irrev"), 
	aes(x=0-mean, y=patient)) + 
annotate("rect", xmin = -52/12*4, xmax = -52/12*8, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["H3"]) +
annotate("rect", xmin = -52, xmax = -52/12*18, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["H1"]) +
annotate("rect", xmin = -52/12*24, xmax = -52/12*36, ymin = 0, ymax = 3.5, alpha = .5, fill=patient_palette["H2"]) +
geom_boxplot(outliers=FALSE, aes(fill=patient), alpha=0.3) +
geom_jitter(height=0.1, width=0,size=0.1)+
facet_wrap(.~model, scales="free_x") +
theme_bw() +
xlab("UCA date (weeks before sample)") +
ylab("Patient") +
scale_fill_manual(values=patient_palette) +
guides(fill="none") +
scale_y_discrete(limits=rev)
dev.off()

pdf("results/UCA_dates_TL_norect_subset.pdf", width=3.0,height=2)
ggplot(filter(allheights, model=="typelinked-irrev"), 
	aes(x=0-mean, y=patient)) + 
geom_boxplot(outliers=FALSE, aes(fill=patient), alpha=0.3) +
geom_jitter(height=0.1, width=0,size=0.1)+
theme_bw() +
xlab("UCA date (weeks before sample)") +
ylab("Patient") +
scale_fill_manual(values=patient_palette) +
guides(fill="none") +
scale_y_discrete(limits=rev)
dev.off()


## get differentiation points
alltrees$clone_id = paste0(alltrees$patient, "-",alltrees$clone_id)
diff_points = dowser::getDiffPoints(alltrees, trait="location",
	tip_traits=c("subset","collapse_count"),
	eo_adjust=TRUE, eo_type="germinal_center", full_posterior=TRUE,
	nproc=16, verbose=TRUE)

saveRDS(diff_points, "results/diffpoints.rds")
diff_points = readRDS("results/diffpoints.rds")

diff_points$patient = sapply(strsplit(diff_points$clone_id, split="-"), function(x)x[1])
diff_points$date = 0 - diff_points$node_height_mean
diff_points$patient = gsub("IV","",diff_points$patient)

dclones = diff_points %>%
	group_by(patient, clone_id, subset) %>%
	summarize(mean_date = mean(date), tree_height=paste0(unique(node_height_mean), collapse=","))

my_comparisons <- list(c("UnMem", "MemHi"), c("UnMem", "MemLo"), c("MemHi", "MemLo"))

pdf("results/TLdiff_fullposterior.pdf",width=3,height=3)
ggplot(filter(diff_points, tip != "Germline" & subset != "GC"), aes(y = date, x = subset)) +
geom_boxplot(outlier.shape=NA,aes(fill=patient), alpha=0.25) +
stat_compare_means(comparisons=my_comparisons,
	method="wilcox.test",label="p.format", size=2.5,fill="patient") +
geom_jitter(width=0.1,size=0.25) +
facet_wrap(patient~., scales="free",ncol=1,strip.position="right")+
coord_flip() +
scale_fill_manual(values=patient_palette) +
geom_hline(yintercept=0, linetype="dashed") +
theme_bw() +
theme(legend.position = "none")+
scale_y_continuous(expand = expansion(mult = 0.06))
dev.off()

# find good example tree
# add subset to tips
for(i in 1:nrow(alltrees)){
	pd = alltrees$trees[[i]]
	tips = pd@phylo$tip.label
	nodes = as.numeric(pd@data$node)
	pd@data$label = NA
	pd@data$label[nodes <= length(tips)] = tips[nodes[nodes <= length(tips)]]
	m1  = match(pd@data$label,alltrees$data[[i]]@data$sequence_id)

	pd@data$subset = alltrees$data[[i]]@data$subset[m1]
	alltrees$trees[[i]] = pd
}
alltrees$types = sapply(alltrees$data, function(x)n_distinct(x@data$subset))
a4 = filter(alltrees, types == 4 & model=="typelinked-irrev")

shapes = c("GC"=21,"germinal_center"=21,"other"=22,"UnMem"=23,"MemHi"=24,"MemLo"=25)
p4 = plotTrees(a4)
p4s = lapply(p4, function(x){
	x + geom_tippoint(aes(shape=subset, fill=location))+
	geom_nodepoint(aes(fill=location, shape=location)) +
	scale_shape_manual(values=shapes)
})
treesToPDF(p4s, "results/p4s.pdf",nrow=1,ncol=2)

# filter down to example clone
palette = c("germinal_center" = "red", "germinal_center, other"="purple", "other"="blue")
tree = filter(alltrees, clone_id == "HIV2-59968")$trees[[1]]
dps = getDiffPoints(a4,trait="location",
	tip_traits=c("subset","collapse_count"),
	eo_adjust=FALSE, eo_type="germinal_center",
	height="CAheight_mean")
dps = filter(dps, clone_id == "HIV2-59968")
dphs = unique(filter(dps, tip_type == "other")$parent_height)
mh = max(as.numeric(tree@data$CAheight_mean),na.rm=TRUE)

# add in GC probabilities
tree@data$prob_gc <- NA
 for (i in 1:nrow(tree@data)) {
       tree@data$prob_gc[i] <- ifelse(grepl("germinal_center", tree@data$location[i]), 
       	as.numeric(tree@data$location.prob[i]), 
     	1-as.numeric(tree@data$location.prob[i]))
 }

pdf("results/extree_wider.pdf",width=7,height=6)
ggtree(tree, aes(color=as.numeric(expectedOccupancies)),size=1) + 
geom_tippoint(aes(shape=subset, fill=prob_gc),color="black",size=3) +
geom_nodepoint(aes(fill=prob_gc, shape=location),color="black",size=3) +
scale_shape_manual(values=shapes) +
scale_color_gradient2(
            low = palette[3],
            mid = palette[2],
            high = palette[1],
            midpoint = 0.5,
            breaks = c(0, 0.25, 0.5, 0.75, 1),
            limits = c(0, 1)) +
scale_fill_gradient2(
            low = palette[3],
            mid = palette[2],
            high = palette[1],
            midpoint = 0.5,
            breaks = c(0, 0.25, 0.5, 0.75, 1),
            limits = c(0, 1)) +
labs(color=paste("Occupancy\nin GC"), fill="Cell type\nProb. GC") +
geom_treescale(width=mh, fontsize=1)+
geom_vline(xintercept=mh-dphs, linetype="dashed",linewidth=0.5,color="darkgrey")
dev.off()

pdf("results/extree_wider_densitree.pdf",width=7,height=6)
plotTrees(filter(alltrees, clone_id == "HIV2-59968"), 
	densitree=TRUE, layout="slanted", alpha=0.1,
	show_occupancy=TRUE, palette=c("germinal_center"="red",
		"germinal_center-other"="purple","other"="blue"), 
	tips="location",
	scale=FALSE)[[1]] +
geom_treescale(width=25, x=0,y=0)
dev.off()

