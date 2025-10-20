# Kenneth B. Hoehn
# 9/19/25
# Run TyCHE analysis on Pseudomonas data

library(dowser)
library(airr)
library(dplyr)

sessionInfo()

beast = "/Applications/BEAST\ 2.7.7/bin"
cores =1
seed = 12345
mcmc_length = 1e+09
iterations = 10

if(TRUE){

    # read in data
    fasta = readFasta("data/hypermutator_paper_alignment_L.fasta")
    data = tibble(sequence_alignment = unlist(fasta))
    data$sequence_id = names(fasta)

    # parse metadata
    str = strsplit(names(fasta), split="_")
    data$sample_time = as.numeric(sapply(str, function(x)x[4])) + 46
    data$location = sapply(str, function(x)x[2])
    data$BL = sapply(str, function(x)x[3])

    # All lung data
    data = filter(data, BL == "L")
    data$germline_alignment = paste0(rep("N", 
        length=nchar(fasta[[1]])), collapse="")
    data$clone_id = "mixed_lung"
    data$v_call = "N"
    data$j_call = "N"
    data$junction_length = 0
    data$locus="IGH"
    clones = formatClones(data, filterstop=FALSE,
        traits=c("location", "sample_time"), chain="H", 
        nproc=1, collapse = FALSE, 
        germ="germline_alignment", use_regions=FALSE)
    saveRDS(clones, "intermediates/clones.rds")

    # split by H and H
    ndata = filter(data, location=="N")
    hdata = filter(data, location=="H")
    nclones = formatClones(ndata, filterstop=FALSE,
        traits=c("location", "sample_time"), chain="H", 
        nproc=1, collapse = FALSE, 
        germ="germline_alignment", use_regions=FALSE)

    hclones = formatClones(hdata, filterstop=FALSE,
        traits=c("location", "sample_time"), chain="H", 
        nproc=1, collapse = FALSE, 
        germ="germline_alignment", use_regions=FALSE)

    # estimate clock rates for H and N separately
    ntree = getTimeTreesIterate(nclones, beast=beast, time="sample_time",
        dir=paste0("~/Documents/pseudomonas/"), id=paste0("N_strict_v003"), 
        template="templates/StrictClock_Standard_EmpFreq.xml", nproc=cores, 
        mcmc_length=mcmc_length, iterations=iterations,
        log_every="auto", include_germline=FALSE,
        CLOCK_RATE_INIT=0.001, seed=seed,
        KAPPA_PRIOR_M=1.25,
        KAPPA_PRIOR_S = 0.5, ignore="freqParameter")
    saveRDS(ntree, "intermediates/ntree.rds")

    htree = getTimeTreesIterate(hclones, beast=beast, time="sample_time",
        dir=paste0("~/Documents/pseudomonas/"), id=paste0("H_strict_v003"), 
        template="templates/StrictClock_Standard_EmpFreq.xml", nproc=cores, 
        mcmc_length=mcmc_length, iterations=iterations,
        log_every="auto", include_germline=FALSE,
        CLOCK_RATE_INIT=0.001, seed=seed,
        KAPPA_PRIOR_M=1.25,
        KAPPA_PRIOR_S = 0.5, ignore="freqParameter")
    saveRDS(htree, "intermediates/htree.rds")
}

clones = readRDS("intermediates/clones.rds")
ntree = readRDS("intermediates/ntree.rds")
htree = readRDS("intermediates/htree.rds")

runs = c(
    "ucld",
    "strict",
    "typelinked-est-irrev",
    "typelinked-est-irrev-relaxed"
        )

templates = c(
    "typelinked-est-irrev"=paste0("templates/TraitLinkedExpectedOccupancy_EstimatedTraitClockRates_EmpFreq.xml"),
    "typelinked-est-irrev-relaxed"=paste0("templates/TraitLinkedExpectedOccupancy_EstimatedTraitClockRates_EmpFreq.xml"),
    "strict"=paste0("templates/StrictClock_AncestralReconstruction_EmpFreq.xml"),
    "ucld"="templates/UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq.xml"
        )
ignore = c("traitfrequencies" ,"freqParameter", "rateCategories")

for(run in runs){
    xtemplate = templates[[run]]
    print(paste(xtemplate))

    TRAIT_MEAN_1 = filter(htree$parameters[[1]], item=='geneticClockRate')$mean
    TRAIT_MEAN_2 = filter(ntree$parameters[[1]], item=='geneticClockRate')$mean
    RATE_INDICATORS = "0 1"
    TRAIT_RATE_SIGMA_1 = TRAIT_MEAN_1 * 0.001
    TRAIT_RATE_SIGMA_2 = TRAIT_MEAN_2 * 0.001
    TRANSITION_RATE_ALPHA_1 = 0.1
    TRANSITION_RATE_ALPHA_2 = 0.1
    TRANSITION_RATE_BETA_1 = 1.0
    TRANSITION_RATE_BETA_2 = 1.0

    if(run=="typelinked-est-irrev-relaxed"){
        TRAIT_RATE_SIGMA_1 = 0.001
        TRAIT_RATE_SIGMA_2 = TRAIT_MEAN_2 * 0.1
    }
        
    trees = getTimeTreesIterate(clones, beast=beast, 
        trait="location", time="sample_time",
        dir=paste0("~/Documents/pseudomonas/"), 
        id=paste0(run,"_v003"), 
        template=xtemplate, nproc=cores, 
        mcmc_length=mcmc_length,
        iterations=iterations,
        INITIAL_STATE=1,
        UCLD_SIGMA_INIT=1.0,
        KAPPA_PRIOR_M=1.25,
        KAPPA_PRIOR_S = 0.5,
        log_target=2000,
        log_every="auto", include_germline=FALSE,
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
        seed=seed, ignore=ignore)

    saveRDS(trees, paste0("intermediates/",run,"_v003_trees.rds"))
    plots = plotTrees(trees, tips="location", nodes=TRUE, scale=108)
    treesToPDF(plots,paste0("results/",run,"_v003_trees.pdf"))
}

