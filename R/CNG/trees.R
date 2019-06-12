library(ape)
library(phytools)
library(phangorn)
library(seqinr)
library(tidyverse)
library(ggtree)
library(geiger)
library(conflicted)
library(here)
library(cowplot)

conflict_prefer("scale_fill_viridis_c", "ggplot2")
conflict_prefer("ggsave", "cowplot")
conflict_prefer("filter", "dplyr")

# Reading and rooting -----------------------------------------------------

# read in tree
bayes_file <- here("data", "S3_Data.nxs")
bayes <- ape::read.nexus(bayes_file)
bayes <- bayes[[1]]

# basic tree with node labels
unrooted <- ggtree(bayes, branch.length = "none", layout = "circular") +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.3) +
  geom_tiplab2() +
  NULL
# save_plot(here("plots", "unrooted.pdf"), unrooted, base_height = 40)

# reroot
bayes <- phytools::reroot(bayes, 119)
rooted <- ggtree(bayes, size = 1.3, layout = "circular", branch.length = "none") + 
  geom_label2(aes(subset = !isTip, label = label, fill = as.numeric(label))) +
  geom_tiplab2(size = 12) +
  geom_text2(aes(subset = !isTip, label = node), size = 16, hjust = -.3) +
  # scale_fill_viridis_c(limits = c(0, 100)) +
  NULL
# save_plot(here("plots", "rooted.pdf"), rooted, base_height = 40)

# Annotating --------------------------------------------------------------

# load in reference file matching species <-> clade
species_info <- read_delim("../../auxillary/species_info.csv",
                           col_names = TRUE, 
                           delim = ",") %>%
  select(-BioProject)

# change genus_species to G. species
species <- species_info$Full %>%
  str_replace("(^.)(.*_)", "\\1. ") %>%
  str_replace("(^.)", toupper)
species_info$Full <- species

# create a new tip.label column that uses Protein_IDs for C. elegans (where applicable) and Transcript_IDs for everything else
data <- as_tibble(rooted$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(Species = case_when(Species %in% c("acant", "acani") ~ Species,
                             !Species %in% c("acant", "acani") ~ str_extract(Species, "^.{4}"))) %>%
  left_join(., species_info) 

# select nodes to display bootstrap
bootstrap_families <- c(
  tax2 = 168,
  tax4 = 142,
  cng1 = 125,
  cng3 = 195,
  iiicng2 = 204,
  ivcng2 = 217,
  vcng2 = 114
)

bootstrap_families <- tibble::enframe(bootstrap_families) %>%
  rename(bootstrap_family = name, node = value) %>%
  select(node, bootstrap_family)

data <- left_join(data, bootstrap_families)

# Tree plotting -----------------------------------------------------------

t <- ggtree(bayes, size = 1.3, layout = "circular", branch.length = "none")

# attach family node/label information to the tree, for later labeling of collapsed nodes
t <- t %<+% data

# family highlighting
t_ann <- t +
  geom_cladelabel(node = 168, label = "TAX-2", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 142, label = "TAX-4", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 125, label = "CNG-1", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 195, label = "CNG-3", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 204, label = "Clade III CNG-2/CHE-6", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 217, label = "Clade IV CNG-2/CHE-6", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 114, label = "Clade V CNG-2/CHE-6", barsize = 10, color = "black", offset = 3, fontsize = 16) +
  NULL

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
# t_ann$data <- t_ann$data %>%
#   mutate(original.label = label) %>%
#   mutate(label = "")

t_final <- t_ann +
  geom_tippoint(aes(color = Clade, angle = angle), size = 12) +
  geom_label2(aes(subset = !is.na(bootstrap_family), label = round(as.numeric(label) * 100, digits = 0), fill = round(as.numeric(label) * 100)), size = 16) +
  labs(fill = "Posterior Probability") +
  scale_fill_viridis_c(limits = c(0, 100)) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  theme(
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  NULL
# t_final
# ggsave(here("plots", "Fig3B_raw.pdf"), t_final, width = 40, height = 40, useDingbats = FALSE)

# saveRDS(t_final, here("data", "final.tree"))

t_label <- t %<+% data +
  geom_tippoint(aes(color = Clade), size = 12) +
  geom_tiplab2(hjust = -0.06) +
  geom_label2(aes(subset = !is.na(bootstrap_family), label = round(as.numeric(label) * 100, digits = 0), fill = round(as.numeric(label) * 100)), size = 16) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  scale_fill_viridis_c(limits = c(0, 100)) +
  NULL

t_label$data <- t_label$data %>%
  mutate(label = str_glue("{Transcript_ID} ({Full})")) 
# t_label

ggsave(here("plots", "Fig3B_label.pdf"), t_label, width = 40, height = 40, useDingbats = FALSE)


# Subtrees ----------------------------------------------------------------

# TAX-4 -------------------------------------------------------------------

tax4 <- extract.clade(bayes, node = 142)

tax4.tree <- ggtree(tax4, layout = "rectangular", branch.length = "none")

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
tax4.data <- as_tibble(tax4.tree$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(Species = case_when(Species %in% c("acant", "acani") ~ Species,
                             !Species %in% c("acant", "acani") ~ str_extract(Species, "^.{4}"))) %>%
  left_join(., species_info) 

# add point colors for each clade
tax4.tree <- tax4.tree %<+% tax4.data

tax4.tree$data <- tax4.tree$data %>%
  mutate(label = str_glue("{Transcript_ID} ({Full})"))

tax4.final <- tax4.tree +
  geom_tippoint(aes(color = Clade), size = 6) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  geom_tiplab(fontface = "bold") +
  xlim(NA, 24) +
  NULL
# tax4.final

# saveRDS(tax4.final, here("data", "tax4.tree"))

# ggsave(here("plots", "Fig3D_raw.pdf"), tax4.final, width = 10, height = 10, useDingbats = FALSE)


# TAX-2 -------------------------------------------------------------------

tax2 <- extract.clade(bayes, node = 168)

tax2.tree <- ggtree(tax2, layout = "rectangular", branch.length = "none")

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
tax2.data <- as_tibble(tax2.tree$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(Species = case_when(Species %in% c("acant", "acani") ~ Species,
                             !Species %in% c("acant", "acani") ~ str_extract(Species, "^.{4}"))) %>%
  left_join(., species_info) 

# add point colors for each clade
tax2.tree <- tax2.tree %<+% tax2.data

tax2.tree$data <- tax2.tree$data %>%
  mutate(label = str_glue("{Transcript_ID} ({Full})"))

tax2.final <- tax2.tree +
  geom_tippoint(aes(color = Clade), size = 6) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  geom_tiplab(fontface = "bold") +
  xlim(NA, 24) +
  NULL
# tax2.final

# saveRDS(tax2.final, here("data", "tax2.tree"))

# ggsave(here("plots", "Fig3F_raw.pdf"), tax2.final, width = 10, height = 10, useDingbats = FALSE)

