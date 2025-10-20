library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)

file = "/Volumes/HoehnK/jessie/simble-validation/cross_reactivity_226348/affinity_data.csv"
data_affinity = read.csv(file)
# filter to 10 clones
clone_ids = unique(data_affinity$clone_id)
data_affinity = filter(data_affinity, clone_id %in% clone_ids[1:10])

affinity_avg <- data_affinity %>%
  mutate(celltype = recode(celltype, "plasma_cell" = "Plasma Cell", "memory_b_cell" = "Memory B Cell", "gc_b_cell" = "GC B Cell")) %>%
  filter(celltype %in% c("Plasma Cell", "Memory B Cell")) %>%
  group_by(clone_id, celltype, .drop=FALSE) %>%
  summarize(mean_affinity = mean(affinity), mean_cross_reactivity=mean(cross_reactivity))

affinity_avg$clone_id <- as.factor(affinity_avg$clone_id)

celltype_levels = c("Memory B Cell", "Plasma Cell", "GC B Cell")
celltype_palette = c("#648FFF", "#FFB000", "#B0B0B0")

names(celltype_palette) = celltype_levels

affinity_avg$celltype <- factor(affinity_avg$celltype, levels=celltype_levels)

format_plot = function(x, accuracy) {
  x <- x +
    scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format(
                         "log2", 
                         math_format(
                           2^.x, 
                           format=function(x){number(x, accuracy = accuracy)}
                         )
                       ),
                       expand = expansion(mult = c(0.1, 0.15))
    )+
    geom_line(aes(group = clone_id), color="grey") + 
    geom_point(aes(color=celltype))+
    scale_color_manual(values=celltype_palette)+
    theme_bw()+xlab("") +
    theme(panel.grid.minor.y = element_blank())+
    scale_x_discrete(labels= c("Memory B Cell" = "Memory",
                                          "Plasma Cell" = "Plasma"))
  
  return(x)
}
a <- ggplot(affinity_avg, aes(y=mean_affinity, x=celltype))
a <- format_plot(a, 1) + ggtitle("Affinity")

b <- ggplot(affinity_avg, aes(y=mean_cross_reactivity, x=celltype))
b <- format_plot(b, 1) + ggtitle("Cross-reactivity")

saveRDS(a, "~/simble-validation/affinity_plot.rds")
saveRDS(b, "~/simble-validation/cross_reactivity_plot.rds")

g=ggarrange(a, b, common.legend = TRUE, legend="none")
g=annotate_figure(g,
                  left = text_grob("affinity score", size=17, rot = 90),
)
g
