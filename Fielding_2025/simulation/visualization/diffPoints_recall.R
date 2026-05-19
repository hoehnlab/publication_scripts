
library(dowser)
library(dplyr)
library(ggplot2)
library(treeio)
library(ggtree)

tt = treeio::read.beast("/Volumes/HoehnK/Hunter/Type_Linked_Clock/gc_reentry/gc_recall_sims/simble_sims_gc_reentry_12_18//all_simplified_time_trees.nex")

true = tibble()
for(i in 1:length(tt)){
	tips = tt[[i]]@phylo$tip.label
	data = tt[[i]]@data
	data$node = as.numeric(data$node)
	dtips = filter(data, node <= length(tips))
	dtips$sequence_id = paste0(dtips$cell_id,"_heavy")
	dtips$clone_id = i
	for (m in tips) {
	  parentM = parent(tt[[i]], tt[[i]]@data$node[m])
	  diff = as.numeric(tt[[i]]@data$time_of_differentiation[m])
	  if(tt[[i]]@data$location[m] == "other" & as.numeric(tt[[i]]@data$generation[m]) > 1100 & as.numeric(tt[[i]]@data$generation[node == parentM]) != diff) {
	    print("in here, diff was: ")
	    print(diff)
	    diff = as.numeric(tt[[i]]@data$generation[node == parentM])
	    print("diff is: ")
	    print(diff)
	  }
	}
	# true = bind_rows(true, dtips)
}

models = c("TyCHE", "TyCHE-nonEO", "SC", "UCLD")
dps = tibble()
for(model in models){
	print(model)
	if(grepl("TyCHE", model)){
		t = readRDS("recall/expectedOccupancy_EstTraitClockRates_EmpFreq_trees.rds")
	}
	if(model == "SC"){
		t = readRDS("recall/strictClock_AncestralReconstruction_EmpFreq_trees.rds")
	}
	if(model == "UCLD"){
		t = readRDS("recall/UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq_trees.rds")
	}

	if(model == "TyCHE"){
		dp = getDiffPoints(t, "celltype", eo_adjust=TRUE,
			eo_type="germinal_center")
	}else{
		dp = getDiffPoints(t, "celltype")
	}
	dp$model = model
	m = match(dp$tip, true$sequence_id)
	parentM = parent(true, true$node[m])
	diff = as.numeric(true$time_of_differentiation[m])
	if(true$location[m] == "other" & as.numeric(true$generation[m]) > 1100 & as.numeric(true$generation[node == parentM]) != diff) {
	  
	  diff = as.numeric(true$generation[node == parentM])
	}
	dp$true_date = diff
	dps = bind_rows(dps, dp)
}


write.csv(dps, "recall/dps.csv")


odps = filter(dps, tip_type != "GC")

odps$model = factor(odps$model, levels=models)

table(odps$model)

palette = c("other" = "#0173B2", "PC" = "#E69F00", "GC"="red", "germinal_center"="red")

pdf("recall/diff_time.pdf", width=6,height=4)
ggplot(odps,aes(y=node_height, x=1200-true_date, color=tip_type)) + 
geom_point() + geom_abline() + theme_bw() +
xlab("True GC exit date (generations)") + ylab("Estimated GC exit date (generations)") +
facet_wrap(.~model, scales="free")
dev.off()

pdf("recall/diff_time_fig.pdf", width=5.5,height=1.75)
ggplot(filter(odps, model != "TyCHE-nonEO"), 
	aes(y=1200-node_height, x=true_date, color=tip_type)) + 
geom_point(size=1,alpha=0.75) + geom_abline() + theme_bw() +
xlab("True differentiation date (generation)") + 
ylab("Est. differentiation date") +
facet_wrap(.~model, scales="free",nrow=1) +
scale_color_manual(values=palette)
dev.off()


tp = treeio::read.beast("recall/all_pruned_time_trees.nex")

pdf("recall/truetrees.pdf", height=2, width=3)
for(i in 1:length(tp)){
	tr = tp[[i]]
	print(ggtree(tr, aes(color=location)) +
	scale_color_manual(values=palette)) +
	geom_treescale(width=200)
}
dev.off()

pdf("recall/truetree.pdf", height=2, width=3)
ggtree(tr, aes(color=location)) +
scale_color_manual(values=palette) +
geom_treescale(width=200)
dev.off()
