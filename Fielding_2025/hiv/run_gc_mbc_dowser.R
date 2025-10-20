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

beast = "~/Programs/beast/bin"
trait_list <- c("germinal_center", "other")
patient = "HIV1"
minseqs = 20
cores = 8
max_group = 25
resample = FALSE
regd_trees = FALSE
runid = "v003"
patients = c("HIV1", "HIV3", "HIV2")
mcmc_length = 1e+08
iterations = 15
seed = 12345

t = read.csv("data/allclock_results_clustered.csv")
me = filter(t, p < 0.05)
liao = filter(me, study=="Liao_2013")
gclock = mean(liao$slope)

patients = commandArgs(trailingOnly=TRUE)[1]
print(patients)

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
		saveRDS(clones,paste0("intermediates/",patient,"_clones.rds"))

		clones$other = sapply(clones$data, function(x)sum(x@data$location == "other"))
		clones$gc = sapply(clones$data, function(x)sum(x@data$location == "germinal_center"))

		clones = clones[order(clones$seqs, decreasing=TRUE),]
		clones = filter(clones, gc > 0 & other > 0)

		for(i in 1:nrow(clones)){
			ct = table(clones$data[[i]]@data$location)
			size = min(c(ct, max_group))
			gc = sample(which(clones$data[[i]]@data$location == "germinal_center"), size=size, replace=FALSE)
			other = sample(which(clones$data[[i]]@data$location == "other"), size=size, replace=FALSE)

			clones$data[[i]]@data = clones$data[[i]]@data[c(gc, other),]
			clones$seqs[[i]] = nrow(clones$data[[i]]@data)
		}
		clones = clones[order(clones$seqs, decreasing=TRUE),]
		clones$other = sapply(clones$data, function(x)sum(x@data$location == "other"))
		clones$gc = sapply(clones$data, function(x)sum(x@data$location == "germinal_center"))
		saveRDS(clones, paste0("intermediates/",patient,"_clones_sampled.rds"))
	}

	clones = readRDS(paste0("intermediates/",patient,"_clones_sampled.rds"))
	print(clones)
	print(clones$seqs)

	clones = filter(clones, seqs >= minseqs)

	if(regd_trees){
		gtrees = getTrees(clones, nproc=cores, build="pml")
		gdt = getTrees(gtrees, igphyml="~/Programs/igphyml/src/igphyml",
			trait="location", fixtrees=TRUE, dir=paste0("~/Documents/hiv_beast/igphyml_",patient),
			check_divergence=FALSE, nproc=cores)
		saveRDS(gdt, paste0("intermediates/",patient,"_gd_parsimony_trees.rds"))

		# SP test and GD trees
		p = plotTrees(gdt, tips="location", nodes=TRUE, common_scale=TRUE)
		treesToPDF(p, file=paste0("results/",patient,"_gd_parsimony_trees.pdf"), ncol=1,nrow=2)
	}

	runs = c(
	"strict",
	"typelinked-irrev",
	"typelinked-eo-est",
	"ucld"
		)

	templates = c(
	"strict"=paste0("templates/StrictClock_AncestralReconstruction_FixedClockRate_EmpFreq.xml"),
	"typelinked-irrev"=paste0("templates/TraitLinkedExpectedOccupancy_FixedTraitClockRates_EmpFreq.xml"),
	"typelinked-eo-est"="templates/TraitLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml",
	"ucld"="templates/UCRelaxedClock_AncestralReconstruction_FixedTraitClockRates_EmpFreq.xml"
		)

	ignore = c("traitRates", "typeLinkedRates", "freqParameter", "clockRate", "traitfrequencies", "geneticClockRate",
		"rateCategories")

	for(run in runs){
		xtemplate = templates[[run]]
		print(paste(patient, xtemplate))

		# iterate until all parameters in all clones ESS > 100 (or max_iter reached)
		TRAIT_MEAN_1 = gclock
		TRAIT_MEAN_2 = 0.000001
		RATE_INDICATORS = "1 1"
		TRAIT_RATE_SIGMA_1 = TRAIT_MEAN_1 * 0.01
		TRAIT_RATE_SIGMA_2 = 0.001
		TRANSITION_RATE_ALPHA_1 = 0.1
		TRANSITION_RATE_ALPHA_2 = 0.1
		TRANSITION_RATE_BETA_1 = 1.0
		TRANSITION_RATE_BETA_2 = 1.0
		if(run == "typelinked-irrev"){
			RATE_INDICATORS = "1 0"
		}
			
		trees = getTimeTreesIterate(clones, beast=beast, trait="location", time="time",
			dir=paste0("~/Documents/hiv_beast/", patient), id=paste0(run,"_",runid), 
			template=xtemplate, nproc=cores, 
			log_every="auto",
			INITIAL_STATE=0,
			KAPPA_PRIOR_M=0.67,
			KAPPA_PRIOR_S = 0.2,
			CLOCK_RATE_INIT=TRAIT_MEAN_1,
			TRAIT_RATE_MEAN_1=TRAIT_MEAN_1,
			TRAIT_RATE_MEAN_2=TRAIT_MEAN_2,
			TRAIT_RATE_SIGMA_1=TRAIT_RATE_SIGMA_1,
			TRAIT_RATE_SIGMA_2=TRAIT_RATE_SIGMA_2,
			RATE_INDICATORS=RATE_INDICATORS,
			TRANSITION_RATE_ALPHA=TRANSITION_RATE_ALPHA_1,
			TRANSITION_RATE_BETA=TRANSITION_RATE_BETA_1,
			TRANSITION_RATE_ALPHA_1=TRANSITION_RATE_ALPHA_1,
			TRANSITION_RATE_ALPHA_2=TRANSITION_RATE_ALPHA_2,
			TRANSITION_RATE_BETA_1=TRANSITION_RATE_BETA_1,
			TRANSITION_RATE_BETA_2=TRANSITION_RATE_BETA_2,
			UCLD_SIGMA_INIT=0.5,
			mcmc_length=mcmc_length,
			log_target=2000,
			ess_cutoff=200, ignore=ignore,
			iterations=iterations, seed=seed)

		saveRDS(trees, paste0("intermediates/",patient,"_",runid,"_",run,".trees.rds"))
		plots = plotTrees(trees, tips="location", nodes=TRUE, scale=52)
		treesToPDF(plots,paste0("results/",patient,"_",runid,"_",run,".trees.pdf"))
	}
}















