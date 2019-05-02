library(tidyverse)
library(cowplot)

setwd("~/Box/GHdata/50HGI/ChemoR/phylo/ML/")

#########################################################################################################
######################                                                             ######################
######################                        Figure 1B                            ######################
######################                                                             ######################
#########################################################################################################

tree <- read_csv("family_clades.csv")

# load in reference file matching species <-> clade
reference <- read.csv("~/Box/ZamanianLab/Data/Genomics/Phylogenetics/ChemoR/tree/clade_species_phylum.csv", header = FALSE) %>%
  rename(Species = V1, Clade = V2, Phylum = V3)
species_info <- read.csv("~/GitHub/50HGI/auxillary/species_info.csv", header = TRUE) %>%
  select(-BioProject, -Classification)

blast <- read_tsv("../ChemoR_parse.blastout", col_names = FALSE) %>%
  distinct(X2, .keep_all = TRUE) %>%
  select(Species = X1, Gene_ID = X2, Best_Hit = X3) %>%
  left_join(., reference, by = "Species") %>%
  left_join(., tree, by = c("Best_Hit" = "Gene_ID")) %>%
  select(-Phylum, -Best_Hit, -Species.y, -Seq_ID, -Clade.y) %>%
  rename(Species = Species.x, Clade = Clade.x)

df <- tree %>%
  bind_rows(., blast) %>%
  distinct(Gene_ID, .keep_all = TRUE) %>%
  mutate(Gene_ID = ifelse(Species == "celeg", paste0(Gene_ID, "-", Seq_ID), Gene_ID)) %>%
  select(Species, Gene_ID, Clade, Family, -Seq_ID)

df.p <- df %>%
  group_by(Family, Species, Clade) %>%
  summarise(n()) %>%
  rename(Count = "n()") %>%
  ungroup() %>%
  filter(!grepl("pep", Family)) %>%
  group_by(Species) %>%
  mutate(Count_Norm = Count / sum(Count)) %>%
  ungroup %>%
  left_join(., species_info)

chemo_families <- data.frame(Family = c("srh", "str", "sri", "srd", "srj", "srm", "srn", 
                                        "sre", "sra", "srab", "srb",
                                        "srx", "srt", "srg", "sru", "srv", "srxa", 
                                        "srw", "srz", "srbc", "srsx", "srr"),
                             Superfamily = c("Str", "Str", "Str", "Str", "Str", "Str", "Str", 
                                             "Sra", "Sra", "Sra", "Sra", 
                                             "Srg", "Srg", "Srg", "Srg", "Srg", "Srg",
                                             "Solo", "Solo", "Solo", "Solo", "Solo"))

df.p <- left_join(df.p, chemo_families, by = "Family")

# change genus_species to G. species
species <- df.p$Full %>%
  str_replace("(^.)(.*_)", "\\1. ") %>%
  str_replace("(^.)", toupper)

df.p$Full <- species

plot <- ggplot(df.p, aes(Full, Family)) +
  geom_tile(aes(fill = Count_Norm), height = 1, width = 1) +
  facet_grid(Superfamily ~ Clade, scales = "free", space = "free") +
  labs(x = "Species", y = "Chemoreceptor Family") +
  # theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.25, "line"),
        panel.spacing.y = unit(0.25, "line"),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        strip.text = element_text(face = "bold")) +
  scale_fill_gradient(name = "Normalized Receptor Count",low = "white", high = "palegreen4") +
  # scale_fill_viridis_c(name = "Normalized Gene Count", option = "plasma") +
  NULL
plot

save_plot("CHEMO_HEATMAP.pdf", plot, base_width = 10, base_height = 6)


# # change facet labels
# gtab <- ggplot_gtable(ggplot_build(plot))
# stript <- which(grepl('strip-t', gtab$layout$name))
# fillt <- c("#a361c7", "#58a865", "#c65c8a", "gray", "#9b9c3b", "#648ace", "#c98443", "palegreen", "#cb4f42")
# k <- 1
# for (i in stript) {
#   j <- which(grepl('rect', gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
#   gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fillt[k]
#   k <- k+1
# }
# stripr <- which(grepl('strip-r', gtab$layout$name))
# fills <- c("gray", "#708cc9", "#c8615d", "#5fa271")
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
#   gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid::grid.draw(gtab)

# ggsave("CHEMO_HEATMAP.pdf", gtab, width = 15, height = 15)

# bar plot
df.b <- group_by(df.p, Species) %>%
  summarise(Sum = sum(Count)) %>%
  left_join(., species_info)

bar <- ggplot(df.b, aes(x = Species, y = Sum)) +
  geom_bar(aes(fill = Clade), stat = "identity", width = 1) + 
  scale_fill_manual(limits = c("I", "IIIa", "IIIb", "IVa", "IVb", "Va", "Vb", "Vd", "Ve"), values = c("#a361c7", "#58a865", "#c65c8a", "gray", "#9b9c3b", "#648ace", "#c98443", "palegreen", "#cb4f42")) +
  facet_grid(. ~ Clade, scales = "free", space = "free") +
  labs(y = "Receptor Count") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.25, "line"),
        panel.spacing.y = unit(0.25, "line"),
        axis.title.y = element_text(face = "bold"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  NULL
bar

save_plot("CHEMO_BARPLOT.pdf", bar, base_width = 10, base_height = 4)

 #########################################################################################################
######################                                                             ######################
######################                  Supplementary Figure 2                     ######################
######################                                                             ######################
#########################################################################################################

plot2 <- ggplot(df.p, aes(Full, Family)) +
  geom_tile(aes(fill = Count)) +
  facet_grid(Superfamily ~ Clade, scales = "free", space = "free_x") +
  labs(x = "Species", y = "Chemoreceptor Family") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        strip.text = element_text(face = "bold")) +
  scale_fill_viridis(name = "Gene Count", option = "plasma")
# plot

gtab2 <- ggplot_gtable(ggplot_build(plot2))
stript2 <- which(grepl('strip-t', gtab2$layout$name))
fillt2 <- c("#a361c7", "#58a865", "#c65c8a", "gray", "#9b9c3b", "#648ace", "#c98443", "palegreen", "#cb4f42")
k <- 1
for (i in stript2) {
  j <- which(grepl('rect', gtab2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gtab2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fillt2[k]
  k <- k+1
}
stripr2 <- which(grepl('strip-r', gtab$layout$name))
fills2 <- c("gray", "#708cc9", "#c8615d", "#5fa271")
k <- 1
for (i in stripr2) {
  j <- which(grepl('rect', gtab2$grobs[[i]]$grobs[[1]]$childrenOrder))
  gtab2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills2[k]
  k <- k+1
}
grid::grid.draw(gtab2)

ggsave("CHEMO_HEATMAP_raw.pdf", gtab2, width = 15, height = 15)

#########################################################################################################
######################                                                             ######################
######################                                                             ######################
######################                                                             ######################
#########################################################################################################

df3 <- group_by(df.p, Superfamily, Full, Species, Clade) %>%
  summarise(SF_Count = sum(Count)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(SF_Count_Norm = SF_Count / sum(SF_Count))

plot3 <- ggplot(df3, aes(Full, Superfamily)) +
  geom_tile(aes(fill = SF_Count_Norm), alpha = 0.5) +
  facet_grid(Superfamily ~ Clade, scales = "free", space = "free_x") +
  labs(x = "Species", y = "Chemoreceptor Superfamily") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_fill_viridis(name = "Normalized Gene Count", option = "plasma")
# plot3

gtab <- ggplot_gtable(ggplot_build(plot3))
stript <- which(grepl('strip-t', gtab$layout$name))
fillt <- c("#a361c7", "#58a865", "#c65c8a", "gray", "#9b9c3b", "#648ace", "#c98443", "palegreen", "#cb4f42")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
  gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fillt[k]
  k <- k+1
}
stripr <- which(grepl('strip-r', gtab$layout$name))
fills <- c("gray", "#708cc9", "#c8615d", "#5fa271")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
  gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(gtab)

ggsave("CHEMO_HEATMAP_SF.pdf", gtab, width = 15, height = 5)

#########################################################################################################
######################                                                             ######################
######################                                                             ######################
######################                                                             ######################
#########################################################################################################

plot4 <- ggplot(df3, aes(Full, Superfamily)) +
  geom_tile(aes(fill = SF_Count)) +
  facet_grid(Superfamily ~ Clade, scales = "free", space = "free_x") +
  labs(x = "Species", y = "Chemoreceptor Superfamily") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_fill_viridis(name = "Gene Count", option = "plasma")
# plot4

gtab <- ggplot_gtable(ggplot_build(plot4))
stript <- which(grepl('strip-t', gtab$layout$name))
fillt <- c("#a361c7", "#58a865", "#c65c8a", "gray", "#9b9c3b", "#648ace", "#c98443", "palegreen", "#cb4f42")
k <- 1
for (i in stript) {
  j <- which(grepl('rect', gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
  gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fillt[k]
  k <- k+1
}
stripr <- which(grepl('strip-r', gtab$layout$name))
fills <- c("gray", "#708cc9", "#c8615d", "#5fa271")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
  gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(gtab)

ggsave("CHEMO_HEATMAP_SF_raw.pdf", gtab, width = 15, height = 5)

###########################################
df4 <- left_join(df, species_info) %>%
  left_join(., chemo_families) %>%
  select(-Clade, -Species) %>%
  rename(Species = Full) %>%
  filter(!grepl("pep", Family)) %>%
  select(Gene_ID, Species, Family, Superfamily)

write.csv(df4, file="family_assignment.csv", row.names = FALSE, col.names = TRUE)
write.csv(df4, file="../../synteny/family_assignment.csv", row.names = FALSE, col.names = TRUE)
write.csv(df4, file="~/Box Sync/ZamanianLab/Manuscripts/2017-ChemoR/Scripts/family_assignment.csv", row.names = FALSE, col.names = TRUE)


