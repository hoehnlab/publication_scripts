# Kenneth B. Hoehn
# 9/19/25
# Run TyCHE analysis on cancer data

library(dplyr)
library(dowser)
library(ggtree)

sessionInfo()

seed = 12345
mcmc_length = 1e+09
iterations = 20
beast = "/Applications/BEAST\ 2.7.7/bin"

# data for patient 4F0A
patient = "4F0A"
fasta = readFasta("data/GLSS-HF-4F0A/GLSS-HF-4F0A_snp_aln_times.fasta")

if(FALSE){
    data = tibble(sequence_alignment = unlist(fasta))
    data$sequence_id = names(fasta)

    data$germline_alignment = filter(data, grepl("-PN",sequence_id))$sequence_alignment
    data = filter(data, !grepl("-PN",sequence_id))

    str = strsplit(data$sequence_id, split="-")
    data$sample_time = as.numeric(sapply(str, function(x)x[6]))
    data$location = sapply(str, function(x)x[5])

    # data and distance tree for all sequences
    data$clone_id = patient
    data$v_call = "N"
    data$j_call = "N"
    data$junction_length = 0
    data$locus="IGH"
    clones = formatClones(data, filterstop=FALSE,
        traits=c("location", "sample_time"), chain="H", 
        nproc=1, collapse = FALSE, 
        germ="germline_alignment", use_regions=FALSE,
        pad_end=FALSE, mod3=FALSE)

    trees = getTrees(clones, build="dnaml",
    	exec="~/Programs/phylip-3.695/exe/dnaml.app/Contents/MacOS/dnaml")
    saveRDS(trees, paste0("intermediates/",patient,"_gd_HN_tree.rds"))

    pdf(paste0("results/",patient,"_gd_tree.pdf"), width=3,height=3)
    plotTrees(trees, scale=0.1)[[1]] + geom_tippoint(aes(color=sample_time, shape=location)) 
    dev.off()

    # distance tree for N only
    N = formatClones(filter(data, location != "H"), filterstop=FALSE,
        traits=c("location", "sample_time"), chain="H", 
        nproc=1, collapse = FALSE, 
        germ="germline_alignment", use_regions=FALSE, mod3=FALSE)

    Ntrees = getTrees(N, build="dnaml",
    	exec="~/Programs/phylip-3.695/exe/dnaml.app/Contents/MacOS/dnaml")
    Ntrees = correlationTest(Ntrees, time="sample_time")

    saveRDS(Ntrees, paste0("intermediates/",patient,"_gd_N_tree.rds"))

    pdf(paste0("results/",patient,"_gd_N_tree.pdf"), width=5,height=2)
    plotTrees(Ntrees, tips="location")[[1]]
    dev.off()

    ntree = getTimeTreesIterate(N, beast=beast, time="sample_time",
        dir=paste0("~/Documents/beast_cancer/"), id=paste0("N_strict_constant"), 
        template="templates/StrictClock_Standard_EmpFreq_Constant.xml", nproc=1, 
        mcmc_length=mcmc_length, iterations=iterations,
        log_every="auto", include_germline=TRUE,
        CLOCK_RATE_INIT=0.001, seed=seed,
        KAPPA_PRIOR_M=1.25,
        KAPPA_PRIOR_S = 0.5, ignore="freqParameter")
    saveRDS(ntree, paste0("intermediates/strict_N_v003_germline_trees_constant.rds"))

    # distance tree for H only
    H = formatClones(filter(data, sample_time > 0), filterstop=FALSE,
        traits=c("location", "sample_time"), chain="H", 
        nproc=1, collapse = FALSE, 
        germ="germline_alignment", use_regions=FALSE)

    Htrees = getTrees(H, build="dnaml",
    	exec="~/Programs/phylip-3.695/exe/dnaml.app/Contents/MacOS/dnaml")
    Htrees = correlationTest(Htrees, time="sample_time")
    saveRDS(Htrees, paste0("intermediates/",patient,"_gd_H_tree.rds"))

    pdf(paste0("results/",patient,"_gd_H_tree.pdf"), width=5,height=2)
    plotTrees(Htrees, tips="location")[[1]]
    dev.off()
}

# build time trees
Nsc = readRDS("intermediates/strict_N_v003_germline_trees_constant.rds")
Hgd = readRDS(paste0("intermediates/",patient,"_gd_H_tree.rds"))
allgd = readRDS(paste0("intermediates/",patient,"_gd_HN_tree.rds"))

runs = c(
    "ucld",
    "typelinked-est-irrev",
    "strict"
        )

# TODO make sure consistnet with new templates
templates = c(
    "typelinked-est-irrev"=paste0("templates/TraitLinkedExpectedOccupancy_EstimatedTraitClockRates_EmpFreq_Constant.xml"),
    "strict"=paste0("templates/StrictClock_AncestralReconstruction_EmpFreq_Constant.xml"),
    "ucld"="templates/UCRelaxedClock_AncestralReconstruction_EstTraitClockRates_EmpFreq_Constant.xml")

ignore = c("traitfrequencies" ,"freqParameter", "rateCategories")

for(run in runs){
    xtemplate = templates[[run]]
    print(paste(xtemplate))

    TRAIT_MEAN_1 = Hgd$slope
    TRAIT_MEAN_2 = filter(Nsc$parameters[[1]], item=="geneticClockRate")$mean
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
        
    trees = getTimeTreesIterate(allgd, beast=beast, 
        trait="location", time="sample_time",
        dir=paste0("~/Documents/beast_cancer/"), 
        id=paste0(run,"_v003_germline"), 
        template=xtemplate, nproc=1, 
        mcmc_length=mcmc_length,
        iterations=iterations,
        INITIAL_STATE=1,
        KAPPA_PRIOR_M=1.25,
        KAPPA_PRIOR_S = 0.5,
        log_every="auto", include_germline=TRUE,
        UCLD_SIGMA_INIT=1.0,
        log_target=2000,
        CLOCK_RATE_INIT=TRAIT_MEAN_2,
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
        seed=seed, ignore=ignore, germline_range = c(-100000, 100000))

    saveRDS(trees, paste0("intermediates/",run,"_v003_germline_trees.rds"))
    plots = plotTrees(trees, tips="location", nodes=TRUE, scale=108)
    treesToPDF(plots,paste0("results/",run,"_v003_germline_trees.pdf"))
}


