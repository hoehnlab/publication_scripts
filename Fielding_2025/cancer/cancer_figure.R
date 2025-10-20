# Kenneth B. Hoehn
# 9/19/25
# Make figures from cancer analysis

library(dowser)
library(dplyr)
library(tidyr)
library(ggtree)
library(ape)
library(phangorn)
library(ggtree)
library(treeio)
library(viridis)


#geom_range doesn't center HPD height intervals properly
#fix based on this code
#https://github.com/YuLab-SMU/ggtree/issues/306
ggtree_height_bars = function(tr, color="blue", alpha=0.4, size=1.5){
	#make sure everything is numeric
	for(i in 1:nrow(tr@data)){
		if(is.na(tr@data$height_0.95_HPD[[i]][1])){
			tr@data$height_0.95_HPD[i] = 
				list(c(tr@data$height[[i]],
					tr@data$height[[i]]))
		}
		tr@data$height_0.95_HPD[[i]] = as.numeric(tr@data$height_0.95_HPD[[i]])
	}
	
	tree_img = revts(ggtree::ggtree(tr))
	minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
	
	# get X, Y, and error bar mins and maxes
	bar_df = as.data.frame(minmax) %>%
	  rename(min = 2,
	         max = 1) %>%
	  dplyr::mutate_all(~-.x) %>%
	  dplyr::bind_cols(dplyr::select(tree_img$data, y))
	
	# add error bars as line segmnents
	tree_img + 
	  ggplot2::geom_segment(aes(x=min, y=y, xend=max, yend=y), 
	               data=bar_df, 
	               color=color,
	               alpha = alpha,
	               size = size) 
}


patient = "4F0A"
gd = readRDS(paste0("intermediates/",patient,"_gd_HN_tree.rds"))
ngd = readRDS(paste0("intermediates/",patient,"_gd_N_tree.rds"))
tlt = readRDS(paste0("intermediates/typelinked-est-irrev_v003_germline_trees.rds"))
st = readRDS(paste0("intermediates/strict_v003_germline_trees.rds"))
uc = readRDS(paste0("intermediates/ucld_v003_germline_trees.rds"))
no = readRDS(paste0("intermediates/strict_N_v003_germline_trees_constant.rds"))

dates = c(0, 3042, 3772)
labs = c(
"Primary tumor",
"1st recurrence",
"2nd recurrence"
	)

shapes = c("H"=24, "N"=21, "Both"=23)

types = RColorBrewer::brewer.pal(3,"Set1")
names(types) = c("H", "N", "Both")

xmax = max(c(unlist(tlt$trees[[1]]@data$height_0.95_HPD),
	unlist(uc$trees[[1]]@data$height_0.95_HPD)), na.rm=TRUE)

minh = 3372 - max(c(unlist(as.numeric(tlt$trees[[1]]@data$height)),
	unlist(as.numeric(uc$trees[[1]]@data$height))), na.rm=TRUE)

print(tlt$below_ESS)
print(st$below_ESS)
print(uc$below_ESS)
print(st$below_ESS)
print(no$below_ESS)

bind_cols(uc$parameters[[1]]$item, uc$parameters[[1]]$ESS)

pdf(paste0("results/",patient,"_hypermutator_tree_v003.pdf"),width=5,height=1.25)
ggtree_height_bars(tlt$trees[[1]]) +
geom_nodepoint(aes(fill= 3772-(as.numeric(height)),pch=location),size=2) +
geom_tippoint(aes(fill=3772-(as.numeric(height)),pch=location),size=2) +
geom_vline(xintercept=dates-3772, linetype="dashed") +
scale_fill_viridis(direction=1, limits=c(minh, end=3772)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Clone") +
xlim(-xmax, 0)

ggtree_height_bars(uc$trees[[1]]) +
geom_nodepoint(aes(fill= 3772-(as.numeric(height)),pch=location),size=2) +
geom_tippoint(aes(fill=3772-(as.numeric(height)),pch=location),size=2) +
geom_vline(xintercept=dates-3772, linetype="dashed") +
scale_fill_viridis(direction=1, limits=c(minh, end=3772)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Clone") +
xlim(-xmax, 0)

ggtree_height_bars(st$trees[[1]]) +
geom_nodepoint(aes(fill= 3772-(as.numeric(height)),pch=location),size=2) +
geom_tippoint(aes(fill=3772-(as.numeric(height)),pch=location),size=2) +
geom_vline(xintercept=dates-3772, linetype="dashed") +
scale_fill_viridis(direction=1, limits=c(min(3772-as.numeric(st$trees[[1]]@data$height)), end=3772)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Clone")
dev.off()

pdf(paste0("results/",patient,"_hypermutator_tree_v003_TL.pdf"),width=5,height=1.25)
ggtree_height_bars(tlt$trees[[1]]) +
geom_nodepoint(aes(fill= 3772-(as.numeric(height)),pch=location),size=2) +
geom_tippoint(aes(fill=3772-(as.numeric(height)),pch=location),size=2) +
geom_vline(xintercept=dates-3772, linetype="dashed") +
scale_fill_viridis(direction=1) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Clone") 
dev.off()


mintl = min(3772-(as.numeric(tlt$trees[[1]]@data$height)))

gd = scaleBranches(gd)

pdf("results/distance_tree_gd.pdf", width=4,height=2)
plotTrees(gd, scale=100)[[1]] +
geom_tippoint(aes(fill=sample_time,pch=location),size=2) +
scale_fill_viridis(direction=1, limits=c(mintl, end=3772)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Strain")
dev.off()

pdf("results/distance_tree_gd_big.pdf", width=60,height=2)
plotTrees(gd, scale=10)[[1]] +
geom_tippoint(aes(fill=sample_time,pch=location),size=2) +
scale_fill_viridis(direction=1, limits=c(mintl, end=3772)) +
scale_shape_manual(values=shapes) +
labs(fill="Day", shape="Strain")
dev.off()



tl = tlt$parameters[[1]]
tl$mean = as.numeric(tl$mean)
tl$run = "TL"
sc = st$parameters[[1]]
sc$run = "SC"
up = uc$parameters[[1]]
up$run = "UC"
np = no$parameters[[1]]
np$run = "NS"
comb = bind_rows(tl, sc, up, np)

heights = filter(comb, item %in% c("Tree.height", "TreeHeight"))

heights$strain = "Both"
heights$strain[heights$run == "NS"] = "N"

rates = filter(comb, grepl("typeLinkedRates", item))
rates = bind_rows(rates, filter(comb, grepl("geneticClockRate", item)))

rates$strain = "Both"
rates$strain[rates$run == "NS"] = "N"
rates$strain[rates$run == "TL" & rates$item == "typeLinkedRates.1"] = "H"
rates$strain[rates$run == "TL" & rates$item == "typeLinkedRates.2"] = "N"

hr = bind_rows(heights, rates)

rates$run = factor(rates$run, 
    levels=rev(c("NS","TL","SC","UC")))
heights$run = factor(heights$run, 
    levels=rev(c("NS","TL","SC","UC")))

heights$earliest_tip = 3772
heights$earliest_tip[heights$run=="NS"] = 3024

dates = heights
dates$mean = heights$earliest_tip - heights$mean
dates$X95.HPDup = heights$earliest_tip - heights$X95.HPDup
dates$X95.HPDlo = heights$earliest_tip - heights$X95.HPDlo

size = 2

pdf("results/treeheights.pdf", width=3,height=2.5)
g1 = ggplot(dates, aes(y=run, x=-mean, xmax=-X95.HPDup, xmin=-X95.HPDlo,
	fill=strain, shape=strain)) +
geom_vline(xintercept=0, linetype="dashed") + 
geom_errorbar(width=0.1) +
geom_point(size=size) + 
scale_x_log10() +
theme_bw() +
scale_shape_manual(values=shapes) +
xlab("Root date (days before primary surgery")+
scale_fill_manual(values=types)

g2 = ggplot(rates, aes(y=run, x=mean, xmax=X95.HPDup, xmin=X95.HPDlo,
	fill=strain, pch=strain)) +
geom_errorbar(width=0.1) +
geom_point(size=size) + 
scale_x_log10() +
theme_bw() +
xlab("Clock rate (mu/site/day)") +
scale_shape_manual(values=shapes)+
scale_fill_manual(values=types)

gridExtra::grid.arrange(g1, g2, ncol=1)
dev.off()

select(dates, run, mean, X95.HPDlo, X95.HPDup)

dates %>%
	group_by(run) %>%
	summarize(mean/365, X95.HPDlo/365, X95.HPDup/365)
#  run   `mean/365` `X95.HPDlo/365` `X95.HPDup/365`
#  <fct>      <dbl>           <dbl>           <dbl>
#1 UC       -58.9           -0.0382        -226.   
#2 SC        -0.147         -0.0885          -0.207
#3 TL       -14.8           -9.79           -20.6  
#4 NS       -15.1           -5.93           -26.7  

dates %>%
	group_by(run) %>%
	summarize(mean, X95.HPDlo, X95.HPDup)
#1 UC    -21505.      -13.9  -82505. 
#2 SC       -53.5     -32.3     -75.4
#3 TL     -5420.    -3574.    -7515. 
#4 NS     -5508.    -2164.    -9763.


plotTrees(tlt, node_nums=TRUE)

# do more programatically
hmrca = 7
(3772 - as.numeric(filter(tlt$trees[[1]]@data, node==hmrca)$height))/365
# 2.724447

(3772 - as.numeric(filter(tlt$trees[[1]]@data, node==hmrca)$height_0.95_HPD[[1]]))/365
# 4.1302828 0.8940102






