# Kenneth B. Hoehn
# 7/24/26
# Re-read in HIV results as new Dowser objects

library(dowser)
library(airr)
library(ggtree)
library(dplyr)
library(shazam)
library(tidyr)
library(treeio)

print(sessionInfo())

beast = "~/Programs/beast/bin"
cores = 16
runid = "v009op4"
patients = c("HIV2", "HIV3", "HIV1")

for(patient in patients){
	print(patient)

	runs = c(
	"typelinked-irrev"
		)

	for(run in runs){

		trees = readRDS(paste0("intermediates/",patient,"_",runid,"_",run,".trees.rds"))

		treesr = readBEAST(trees, beast=beast, trait="location", 
			dir=paste0("~/Documents/hiv_beast/", patient), id=paste0(run,"_",runid), 
			nproc=cores, posterior="all")

		saveRDS(treesr, paste0("intermediates/",patient,"_",runid,"_",run,"_allposterior.trees.rds"))
	}
}


t = readRDS("intermediates/HIV1_v009op4_typelinked-irrev_allposterior.trees.rds")

