library(tidyverse)
library(viridis)

setwd("~/Box Sync/GHdata/50HGI/NChR/phylo/ML/")

tree <- read_csv("family_clades.csv")

# load in reference file matching species <-> clade
reference <- read.csv("~/Box Sync/ZamanianLab/Data/Phylogenetics/ChemoR/clade_species_phylum.csv", header = FALSE) %>%
  rename(Species = V1, Clade = V2, Phylum = V3)
species_info <- read.csv("~/Box Sync/GitHub/50HGI/auxillary/species_info.csv", header = TRUE) %>%
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

plot <- ggplot(df.p, aes(Full, Family)) +
  geom_tile(aes(fill = Count_Norm)) +
  facet_grid(Superfamily ~ Clade, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        strip.text = element_text(face = "bold")) +
  scale_fill_viridis(name = "Normalized Gene Count", option = "plasma")
# plot

gtab <- ggplot_gtable(ggplot_build(plot))
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

ggsave("CHEMO_HEATMAP.pdf", gtab, width = 15, height = 15)



plot2 <- ggplot(df.p, aes(Full, Family)) +
  geom_tile(aes(fill = Count)) +
  facet_grid(Superfamily ~ Clade, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        strip.text = element_text(face = "bold")) +
  scale_fill_viridis(name = "Raw Gene Count", option = "plasma")
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



###########################################
df3 <- left_join(df, species_info) %>%
  left_join(., chemo_families) %>%
  select(-Clade, -Species) %>%
  rename(Species = Full) %>%
  filter(!grepl("pep", Family)) %>%
  select(Gene_ID, Species, Family, Superfamily)

write.csv(df3, file="family_assignment.csv", row.names = FALSE, col.names = TRUE)
write.csv(df3, file="../../synteny/family_assignment.csv", row.names = FALSE, col.names = TRUE)
write.csv(df3, file="~/Box Sync/ZamanianLab/Manuscripts/2017-ChemoR/Scripts/family_assignment.csv", row.names = FALSE, col.names = TRUE)


