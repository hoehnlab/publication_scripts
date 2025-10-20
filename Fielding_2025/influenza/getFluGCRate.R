# Kenneth B. Hoehn
# 9/22/25
# get clock rate estimate for GC B cells in influenza

library(dplyr)

# main results file from https://doi.org/10.7554/eLife.70873
t = read.csv("allclock_results_clustered.csv")

# filter to Turner et al 2020 and measurably evolving lineages
print(median(filter(t, study=="turner_2020_all" & p < 0.05)$slope))
# 0.004924567 mut/site/week