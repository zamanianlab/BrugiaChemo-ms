library(tidyverse)

setwd("~/Box Sync/GHdata/50HGI/NChR/phylo/ML/")

tree <- read_csv("family_clades.csv")

# load in reference file matching species <-> clade
reference <- read.csv("~/Box Sync/ZamanianLab/Data/Phylogenetics/ChemoR/clade_species_phylum.csv", header = FALSE) %>%
  rename(Species = V1, Clade = V2, Phylum = V3)

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
  ungroup

chemo_families <- data.frame(Family = c("srh", "str", "sri", "srd", "srj", "srm", "srn", 
                                        "sre", "sra", "srab", "srb",
                                        "srx", "srt", "srg", "sru", "srv", "srxa", 
                                        "srw", "srz", "srbc", "srsx", "srr"),
                             Superfamily = c("Str", "Str", "Str", "Str", "Str", "Str", "Str", 
                                             "Sra", "Sra", "Sra", "Sra", 
                                             "Srg", "Srg", "Srg", "Srg", "Srg", "Srg",
                                             "Solo", "Solo", "Solo", "Solo", "Solo"))

df.p <- left_join(df.p, chemo_families, by = "Family")

plot <- ggplot(df.p, aes(Species, Family)) +
  geom_tile(aes(fill = Count_Norm)) +
  facet_grid(Superfamily ~ Clade, scales = "free") + 
  theme(axis.text.x = element_text(size = rel(0.7), angle = 60, hjust = 1))
plot

ggsave("CHEMO_HEATMAP.pdf", plot, width = 15, height = 15)

write.csv(df, file="family_assignment.csv", row.names = FALSE, col.names = TRUE)



