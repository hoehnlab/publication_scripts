# Kenneth B. Hoehn
# 5/19/2026
# Generate Fig 3D showing differentiation times under different models

library(dowser)
library(dplyr)
library(ggplot2)
library(treeio)
library(ggtree)

tt = treeio::read.beast("all_simplified_time_trees.nex")

true = tibble()
for(i in 1:length(tt)){
	tips = tt[[i]]@phylo$tip.label
	data = tt[[i]]@data
	data$node = as.numeric(data$node)
	dtips = filter(data, node <= length(tips))
	dtips$sequence_id = paste0(dtips$cell_id,"_heavy")
	dtips$clone_id = i
	true = bind_rows(true, dtips)
}

models = c("TyCHE", "TyCHE-nonEO", "SC", "UCLD")
dps = tibble()
for(model in models){
	print(model)
	if(grepl("TyCHE", model)){
		t = readRDS("config_ratio_1to1_sel_ExpectedOccupancy_EstClockRates_all_trees.rds")
	}
	if(model == "SC"){
		t = readRDS("config_ratio_1to1_sel_StrictClock_AncestralReconstruction_all_trees.rds")
	}
	if(model == "UCLD"){
		t = readRDS("config_ratio_1to1_sel_UCRelaxedClock_AncestralReconstruction_all_trees.rds")
	}

	if(model == "TyCHE"){
		dp = getDiffPoints(t, "location", tip_trait="celltype", eo_adjust=TRUE,
			eo_type="germinal_center")
	}else{
		dp = getDiffPoints(t, "location", tip_trait="celltype")
	}
	dp$model = model
	m = match(dp$tip, true$sequence_id)
	dp$true_date = as.numeric(true$time_of_differentiation[m])
	dps = bind_rows(dps, dp)
}
write.csv(dps, "dps.csv")

odps = filter(dps, tip_type != "germinal_center")
odps$model = factor(odps$model, levels=models)
table(odps$model)

palette = c("MBC" = "#0173B2", "PC" = "#E69F00", "GC"="red")
odps$celltype[odps$celltype=="memory_b_cell"] = "MBC"
odps$celltype[odps$celltype=="plasma_cell"] = "PC"

pdf("diff_time.pdf", width=6,height=4)
ggplot(odps,aes(y=node_height, x=200-true_date, color=celltype)) + 
geom_point() + geom_abline() + theme_bw() +
xlab("True GC exit date (generations)") + ylab("Estimated GC exit date (generations)") +
facet_wrap(.~model, scales="free")
dev.off()

pdf("diff_time_fig.pdf", width=5.5,height=1.75)
ggplot(filter(odps, model != "TyCHE-nonEO"), 
	aes(y=200-node_height, x=true_date, color=celltype)) + 
geom_point(size=1,alpha=0.75) + geom_abline() + theme_bw() +
xlab("True differentiation date (generation)") + 
ylab("Est. differentiation date") +
facet_wrap(.~model, scales="free",nrow=1) +
scale_color_manual(values=palette)
dev.off()

tp = treeio::read.beast("all_pruned_time_trees.nex")

pdf("truetrees.pdf", height=2, width=2)
for(i in 1:length(tp)){
	tr = tp[[i]]
	tr@data$celltype[tr@data$celltype=="memory_b_cell"] = "MBC"
	tr@data$celltype[tr@data$celltype=="plasma_cell"] = "PC"
	tr@data$celltype[tr@data$celltype=="gc_b_cell"] = "GC"
	print(ggtree(tr, aes(color=celltype)) +
	scale_color_manual(values=palette)) +
	geom_treescale(width=200)
}
dev.off()

pdf("truetree.pdf", height=2, width=3)
ggtree(tr, aes(color=celltype)) +
scale_color_manual(values=palette) +
geom_treescale(width=200)
dev.off()












