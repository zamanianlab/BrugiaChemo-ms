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

#read in tree
bayes_file <- here("data", "S2_Data.nxs")
bayes <- ape::read.nexus(bayes_file)
bayes <- bayes[[1]]
bayes_drop <- drop.tip(bayes, c("celeg-B0212.5",
                                "celeg-C05C12.3",
                                "celeg-C29E6.2a",
                                "celeg-F28H7.10a",
                                "celeg-F54D1.5",
                                "celeg-K01A11.4a",
                                "celeg-M05B5.6",
                                "celeg-R06B10.4a",
                                "celeg-R13A5.1a",
                                "celeg-T01H8.5a",
                                "celeg-T09A12.3a",
                                "celeg-T10B10.7",
                                "celeg-Y40C5A.2",
                                "celeg-Y71A12B.4",
                                "celeg-Y73F8A.1",
                                "celeg-ZC21.2a",
                                "celeg-ZK512.3"))


# basic tree with node labels
unrooted <- ggtree(bayes_drop, branch.length = "none", layout = "circular") +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.3) +
  geom_tiplab2() +
  NULL
# save_plot(here("plots", "unrooted.pdf"), unrooted, base_height = 40)

# reroot with CUP-5 clade as outgroup
bayes <- phytools::reroot(bayes_drop, 400)
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

data <- as_tibble(rooted$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(Species = case_when(Species %in% c("acant", "acani") ~ Species,
                             !Species %in% c("acant", "acani") ~ str_extract(Species, "^.{4}"))) %>%
  left_join(., species_info) 

# select nodes to display bootstrap
bootstrap_families <- c(
  cup5 = 708, 
  pkd2 = 694,
  trpa1 = 661, 
  trpa2 = 668, 
  osm9 = 635, 
  ocr4 = 631, 
  ocr3 = 619,
  ocr2 = 583,
  ocr1 = 587,
  focr1 = 608, 
  focr2 = 599,
  focr = 593,
  ced11 = 447, 
  gtl2 = 420, 
  fgon2 = 404, 
  gon2 = 399,
  gtl1 = 386, 
  trp4 = 544, 
  spe41 = 510, 
  trp1 = 495, 
  trp2 = 472
)

bootstrap_families <- tibble::enframe(bootstrap_families) %>%
  rename(bootstrap_family = name, node = value) %>%
  select(node, bootstrap_family)

data <- left_join(data, bootstrap_families)

# Tree plotting -----------------------------------------------------------

t <- ggtree(bayes, layout = "circular", branch.length="none") #+ geom_tiplab2(size = 2) + geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)

# attach family node/label information to the tree, for later labeling of collapsed nodes
t <- t %<+% data

# superfamily highlighting
t_ann <- t +
  geom_cladelabel(node = 708, label = "TRPML", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  geom_cladelabel(node = 694, label = "TRPP", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  geom_cladelabel(node = 659, label = "TRPA", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  geom_cladelabel(node = 574, label = "TRPV", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  geom_cladelabel(node = 376, label = "TRPM", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  geom_cladelabel(node = 544, label = "TRPN", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  geom_cladelabel(node = 469, label = "TRPC", barsize = 10, color = "black", offset = 5, fontsize = 18) +
  NULL
  
# family highlighting
t_ann <- t_ann + 
  geom_cladelabel(node = 708, label = "CUP-5", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 694, label = "PKD-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 661, label = "TRPA-1", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 668, label = "TRPA-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 635, label = "OSM-9", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 631, label = "OCR-4", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 619, label = "OCR-3", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 583, label = "OCR-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 587, label = "OCR-1", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 593, label = "Clade III OCR-1/2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 447, label = "CED-11", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 420, label = "GTL-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 404, label = "Clade III GON-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 399, label = "GON-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 386, label = "GTL-1", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 544, label = "TRP-4", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 510, label = "TRP-3/SPE-41", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 495, label = "TRP-1", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  geom_cladelabel(node = 472, label = "TRP-2", barsize = 10, color = "gray40", offset = 1, fontsize = 16) +
  NULL

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
# t_ann$data <- t_ann$data %>%
#   mutate(label = str_glue("{Gene_ID} ({Full})")) %>%
#   mutate(label = ifelse(Clade == "IIIb",final_label, Species)) # Species + Gene_ID + Clade for clade III, species for all other clades

t_final <- t_ann +
  geom_tippoint(aes(color = Clade, angle = angle), size = 8) +
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
# ggsave(here("plots", "Fig3A_raw.pdf"), t_final, width = 40, height = 40, useDingbats = FALSE)

# saveRDS(t_final, here("data", "final.tree"))

t_label <- t %<+% data +
  geom_tippoint(aes(color = Clade), size = 6) +
  geom_tiplab2(hjust = -0.06) +
  geom_label2(aes(subset = !is.na(bootstrap_family), label = round(as.numeric(label) * 100, digits = 0), fill = round(as.numeric(label) * 100)), size = 16) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  scale_fill_viridis_c(limits = c(0, 100)) +
  NULL

t_label$data <- t_label$data %>%
  mutate(label = str_glue("{Transcript_ID} ({Full})")) 
# t_label

# ggsave(here("plots", "Fig3A_label.pdf"), t_label, width = 40, height = 40, useDingbats = FALSE)

# Subtrees ----------------------------------------------------------------

# OSM-9 -------------------------------------------------------------------

osm9 <- extract.clade(bayes, node = 635)

osm9.tree <- ggtree(osm9, layout = "rectangular", branch.length = "none") 

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
osm9.data <- as_tibble(osm9.tree$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(Species = case_when(Species %in% c("acant", "acani") ~ Species,
                             !Species %in% c("acant", "acani") ~ str_extract(Species, "^.{4}"))) %>%
  left_join(., species_info) 

# add point colors for each clade
osm9.tree <- osm9.tree %<+% osm9.data

osm9.tree$data <- osm9.tree$data %>%
  mutate(label = str_glue("{Transcript_ID} ({Full})")) 

osm9.final <- osm9.tree +
  geom_tippoint(aes(color = Clade), size = 6) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  geom_tiplab(fontface = "bold") +
  xlim(NA, 30) +
  NULL
# osm9.final

# saveRDS(osm9.final, here("data", "osm9.tree"))

# ggsave(here("plots", "Fig3C_raw.pdf"), osm9.final, width = 10, height = 10, useDingbats = FALSE)


# OCR-1/2 -----------------------------------------------------------------

ocr <- extract.clade(bayes, node = 581)

ocr.tree <- ggtree(ocr, layout = "rectangular", branch.length = "none") 

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
ocr.data <- as_tibble(ocr.tree$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(Species = case_when(Species %in% c("acant", "acani") ~ Species,
                             !Species %in% c("acant", "acani") ~ str_extract(Species, "^.{4}"))) %>%
  left_join(., species_info) 

# add point colors for each clade
ocr.tree <- ocr.tree %<+% ocr.data

ocr.tree$data <- ocr.tree$data %>%
  mutate(label = str_glue("{Transcript_ID} ({Full})")) 

ocr.final <- ocr.tree +
  geom_tippoint(aes(color = Clade), size = 6) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  geom_tiplab(fontface = "bold") +
  xlim(NA, 30) +
  NULL
# ocr.final

# saveRDS(ocr.final, here("data", "ocr.tree"))

# ggsave(here("plots", "Fig3E_raw.pdf"), ocr.final, width = 10, height = 10, useDingbats = FALSE)

