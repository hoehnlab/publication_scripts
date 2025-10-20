# Kenneth B. Hoehn
# 9/19/25
# Run TyCHE analysis on HIV data

library(dowser)
library(airr)
library(ggtree)
library(dplyr)
library(shazam)
library(tidyr)
library(treeio)

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


t = read.csv("data/allclock_results_clustered.csv")
me = filter(t, p < 0.05)
liao = filter(me, study=="Liao_2013")
gclock = mean(liao$slope)

beast = "~/Programs/beast/bin"
igphyml = "~/Programs/igphyml/src/igphyml"
trait_list <- c("germinal_center", "other")
patient = "HIV1"
runid = "v003"
minseqs = 0
cores = 8
max_group = 50
resample = FALSE
regd_trees = FALSE
mcmc_length = 1e+08
iterations = 15
seed = 12345
patients = c("HIV1", "HIV2", "HIV3")

patients = commandArgs(trailingOnly=TRUE)[1]

for(patient in patients){
	print(patient)
	if(resample){
		# Read the input
		data <- read_rearrangement(paste0("data/",patient,"_germ-pass_parse-delete_airr.tsv"))

		data$time = 0
		data$location = "germinal_center"
		data$location[grepl("Mem", data$subset)] = "other"

		clones = formatClones(data, 
			traits=c("subset", "time", "location"), columns=c("subset","location"),
		    nproc=cores, minseq=10)
		saveRDS(clones,paste0("intermediates/",patient,"_mbc_clones.rds"))

		clones = sampleClones(clones, 100, group="subset")

		clones$MemHi = sapply(clones$data, function(x)sum(x@data$subset == "MemHi"))
		clones$MemLo = sapply(clones$data, function(x)sum(x@data$subset == "MemLo"))
		clones$UnMem = sapply(clones$data, function(x)sum(x@data$subset == "UnMem"))
		clones$GC = sapply(clones$data, function(x)sum(x@data$subset == "GC"))
		saveRDS(clones, paste0("intermediates/",patient,"_clones_subsetsampled.rds"))
	}

	clones = readRDS(paste0("intermediates/",patient,"_clones_subsetsampled.rds"))
	print(clones)
	clones$types = (clones$GC > 0) + (clones$MemHi > 0) + (clones$MemLo > 0) + (clones$UnMem > 0)

	clones = filter(clones, types == 4)
	clones = filter(clones, seqs >= minseqs)
	print(nrow(clones))

	if(regd_trees){
		gtrees = getTrees(clones, nproc=cores, build="pml")
		gdt = getTrees(gtrees, igphyml=igphyml,
			trait="subset", fixtrees=TRUE, dir=paste0("~/Documents/hiv_beast/igphyml_",patient),
			check_divergence=FALSE, nproc=cores)
		saveRDS(gdt, paste0("intermediates/",patient,"subset3_gd_parsimony_trees.rds"))

		# SP test and GD trees
		p = plotTrees(gdt, tips="subset", nodes=TRUE, common_scale=TRUE)
		treesToPDF(p, file=paste0("results/",patient,"subset3_gd_parsimony_trees.pdf"), ncol=1,nrow=2)
	}

	runs = c(
	"gc-origin"
		)

	templates = c(
	"gc-origin"=paste0("templates/TraitLinkedInstantSwitch_EstTraitClockRates_EmpFreq_4state.xml")
		)

	ignore = c("traitRates", "typeLinkedRates", "freqParameter", "clockRate", "traitfrequencies", "geneticClockRate")

	for(run in runs){
		xtemplate = templates[[run]]
		print(paste(patient, xtemplate))

		TRAIT_MEAN_1 = gclock
		TRAIT_MEAN_2 = 0.000001
		TRAIT_MEAN_3 = 0.000001
		TRAIT_MEAN_4 = 0.000001
		RATE_INDICATORS = "1 1 1 0 0 0 0 0 0 0 0 0"
		TRAIT_RATE_SIGMA_1 = TRAIT_MEAN_1 * 0.001
		TRAIT_RATE_SIGMA_2 = 0.001
		TRAIT_RATE_SIGMA_3 = 0.001
		TRAIT_RATE_SIGMA_4 = 0.001
		TRANSITION_RATE_ALPHA_1 = 0.1
		TRANSITION_RATE_BETA_1 = 1.0
		
		trees = getTimeTreesIterate(clones, beast=beast, trait="subset", time="time",
			dir=paste0("~/Documents/hiv_beast/", patient), id=paste0(run, "_",runid), 
			template=xtemplate, nproc=cores, log_every="auto",
			INITIAL_STATE=0,
			KAPPA_PRIOR_M=0.67,
			KAPPA_PRIOR_S = 0.2,
			log_target=2000,
			CLOCK_RATE_INIT=TRAIT_MEAN_1,
			TRAIT_RATE_MEAN_1=TRAIT_MEAN_1,
			TRAIT_RATE_MEAN_2=TRAIT_MEAN_2,
			TRAIT_RATE_MEAN_3=TRAIT_MEAN_3,
			TRAIT_RATE_MEAN_4=TRAIT_MEAN_4,
			TRAIT_RATE_SIGMA_1=TRAIT_RATE_SIGMA_1,
			TRAIT_RATE_SIGMA_2=TRAIT_RATE_SIGMA_2,
			TRAIT_RATE_SIGMA_3=TRAIT_RATE_SIGMA_3,
			TRAIT_RATE_SIGMA_4=TRAIT_RATE_SIGMA_4,
			RATE_INDICATORS=RATE_INDICATORS,
			TRANSITION_RATE_ALPHA=TRANSITION_RATE_ALPHA_1,
			TRANSITION_RATE_BETA=TRANSITION_RATE_BETA_1,
			mcmc_length=mcmc_length, iterations=iterations,
			ignore=ignore, seed=seed)

		saveRDS(trees, paste0("intermediates/",patient,"_",runid,"_",run,"_subset3_trees.rds"))
		p = plotTrees(trees, tips="location", nodes=TRUE, scale=100)
		treesToPDF(p, paste0("results/mbc_trees_",runid,"_",patient,"_subset3.pdf"))
	}
}















