library(tibble)
library(dplyr)


sampleCloneEvenly = function(clone, size, weight=NULL, group=NULL){
  if(size > nrow(clone@data)){
    return(clone)
  }
  if(sum(!group %in% names(clone@data)) != 0){
    stop(paste("One or more grouping columns not found in clone data"))
  }
  
  if(is.null(group)){
    if(is.null(weight)){
      sd = sample(clone@data$sequence_id, size=size)
    } else {
      sd = sample(clone@data$sequence_id, size=size, prob=clone@data[[weight]])
    }
  } else {
    
    combos = lapply(1:length(group), function(x){
      unique(clone@data[[group[x]]])
    })
    
    if (length(group) > 1){
      combo_sets = expand.grid(combos) 
    } else {
      combo_sets = data.frame(combos[[1]])
    }
    colnames(combo_sets) = group
    combo_sets = combo_sets[sample(1:nrow(combo_sets)), , drop = FALSE]
    ncombo_sets = nrow(combo_sets)
    
    group_list = vector("list", ncombo_sets)
    group_sizes = numeric(ncombo_sets)
    
    # Get all sequences in each combination of group(s)
    for (i in seq_len(ncombo_sets)) {
      subset_idx = Reduce(`&`, lapply(group, function(g) clone@data[[g]] == combo_sets[[g]][i]))
      temp = clone@data[subset_idx, ]
      if (nrow(temp) == 0) next
      if (is.null(weight)) {
        group_list[[i]] = sample(temp$sequence_id, size=nrow(temp))
      } else {
        group_list[[i]] = sample(temp$sequence_id, size=nrow(temp), prob=temp[[weight]])
      }
      group_sizes[i] = length(group_list[[i]])
    }
    
    # Sample sequences evenly (as much as possible) from each group
    sd <- c()
    counters <- rep(1, ncombo_sets)
    group_index <- 1
    sd_index <- 1
    
    target_size_per_group <- floor(size / ncombo_sets)
    group_targets <- rep(target_size_per_group, ncombo_sets)
    deficit <- size - sum(group_targets)
    maxed_out_groups <- c()
    for (i in seq_len(ncombo_sets)) {
      if (group_sizes[i] < target_size_per_group) {
        deficit <- deficit + target_size_per_group - group_sizes[i]
        group_targets[i] <- group_sizes[i]
        maxed_out_groups <- c(maxed_out_groups, i)
      }
    }
    
    groups_that_can_still_give <- setdiff(1:ncombo_sets, maxed_out_groups)
    while (deficit > 0) {
      # redistribute the deficit as evenly as possible among groups that have enough sequences, randomly change the order of groups that can still give after each round of redistribution to avoid biasing towards earlier groups
      for (i in sample(groups_that_can_still_give)) {
        if (group_sizes[i] > group_targets[i]) {
          group_targets[i] = group_targets[i] + 1
          deficit = deficit - 1
          if (deficit == 0) break
        }
        else {
          groups_that_can_still_give = setdiff(groups_that_can_still_give, i)
        }
      }
    }
    
    for (i in seq_len(ncombo_sets)) {
      for (j in seq_len(group_targets[i])) {
        sd[sd_index] <- group_list[[i]][j]
        sd_index <- sd_index + 1
        if (sd_index > size) break
      }
      if (sd_index > size) break
    }
    
    
    clone@data = clone@data[clone@data$sequence_id %in% sd,]
    return(clone)
  }
}

sampleClonesEvenly = function(clones, size, weight=NULL, group=NULL){
  clones$data <- lapply(clones$data, function(x) sampleCloneEvenly(x, size, weight, group))
  clones$seqs <- sapply(clones$data, function(x)nrow(x@data))
  return(clones)
}


