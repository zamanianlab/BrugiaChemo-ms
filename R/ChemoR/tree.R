library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(seqinr)
library(ggtree)
library(tidyverse)
library(magrittr)
library(cowplot)
library(here)


# Reading and rooting -----------------------------------------------------

# parse the iqtree output
iqtree_file <- here("data", "S1_Data.nwk")
iqtree <- read.newick(iqtree_file)

# reroot on srw
unrooted <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none") + 
  geom_label2(aes(subset = !isTip, label = label, fill = label)) + geom_tiplab2(size = 2) + 
  geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
  NULL
# save_plot("ggtree/unrooted.pdf", unrooted, base_height = 40)

iqtree <- phytools::reroot(iqtree, 2882)
rooted <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none") + 
  geom_label2(aes(subset = !isTip, label = label, fill = label)) + geom_tiplab2(size = 2) + 
  geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
  NULL
# save_plot(here("plots", "rooted.pdf"), rooted, base_height = 40)

# Annotating --------------------------------------------------------------

# load in reference file matching species <-> clade
reference <- read_delim(here("../../auxillary", "species_info.csv"),
                           col_names = TRUE, 
                           delim = ",") %>%
  select(-BioProject)

# load in reference file matching C. elegans Transcript_ID to Protein_ID
cele_reference <- read_csv(here("../../auxillary/ChemoR/celegans_chemor.csv"), col_names = TRUE) %>%
  rename(Transcript_ID = "Sequence_Name") %>%
  select(-Sequence)

# create a new tip.label column that uses Protein_IDs for C. elegans (where applicable) and Transcript_IDs for everything else
data <- as_tibble(rooted$data) %>%
  mutate(Species = case_when(str_detect(label, "-") == TRUE ~ label,
                             str_detect(label, "-") == FALSE ~ as.character(NA))) %>%
  separate(Species, into = c("Species", "Transcript_ID"), sep = "-", remove = TRUE) %>%
  left_join(., reference) %>%
  left_join(., cele_reference) %>%
  mutate(tip.label = case_when(isTip == TRUE & Species != "cele" ~ str_c(Species, Transcript_ID, Clade, sep = "-"),
                               isTip == TRUE & Species == "cele" & !is.na(Protein_ID) ~ str_c(Species, Protein_ID, Clade, sep = "-"),
                               isTip == TRUE & Species == "cele" & is.na(Protein_ID) ~ str_c(Species, Transcript_ID, Clade, sep = "-"),
                               isTip == FALSE ~ "NA"))

rooted <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none") %<+% data +
  geom_tiplab2(aes(label = tip.label)) +
  geom_label2(aes(subset = !isTip, label = label, fill = label)) + geom_tiplab2(size = 2) + 
  geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
  NULL
# save_plot(here("plots", "rooted.pdf"), rooted, base_height = 40)

# select nodes to display bootstrap
data %<>% mutate(bootstrap = case_when(str_detect(label, "-") == TRUE ~ as.numeric(NA),
                                       str_detect(label, "-") == FALSE ~ as.numeric(label)))

bootstrap_families <- c(
  srw = 4846,
  srw = 4961,
  srxa = 4628,
  srbc = 4687,
  srxa_srbc = 4627,
  srsx = 4469,
  sro = 4453,
  srz = 4412,
  srv = 4367,
  srv = 4313,
  sru = 4172,
  srg = 4196,
  srv_sru_srg = 4150, 
  srr = 4137,
  srt = 3968,
  srx = 3767,
  srb = 2858,
  sre = 2765,
  sra = 2626,
  srab = 2674,
  srn = 3252,
  sri = 3076,
  srh = 2888,
  srd = 3635,
  # srm = 3601,
  srj = 3484,
  sr4 = 3557,
  sr3 = 3577,
  str = 3288,
  Str = 2882,
  Sra = 2623,
  Sra = 2624,
  Sra_Str = 2622,
  Srg = 3765,
  deep = 2621
)

bootstrap_families <- tibble::enframe(bootstrap_families) %>%
  rename(bootstrap_family = name, node = value) %>%
  select(node, bootstrap_family)

data <- left_join(data, bootstrap_families)

# choose C. elegans nodes to collapse
collapse <- c(
  srw = 4986,
  srxa = 4648,
  srbc = 4710,
  srbc = 4736,
  srsx = 4557,
  srsx = 4494,
  srz = 4412,
  srv = 4389,
  sru = 4172,
  srg = 4203,
  srg = 4239,
  srg = 4263,
  srr = 4137,
  srt = 4001,
  srt = 4010,
  srt = 3982,
  srx = 3795,
  srx = 3776,
  srx = 3874,
  srx = 3887,
  srb = 2872,
  sre = 2789,
  sra = 2643,
  srab = 2683,
  sri = 3140,
  srh = 2888,
  srd = 3639,
  srd = 3686,
  srj = 3487,
  str = 3289
)

# convert to data frame
collapse <- tibble::enframe(collapse) %>%
  rename(collapse_family = name, node = value)

# choose P. pacificus nodes to collapse
ppa_collapse <- c(
  srw = 5078,
  srsx = 4512,
  srv = 4294,
  srg = 4254,
  srt = 3296,
  srx = 3955,
  srx = 3833,
  sre = 2825,
  sra = 2666,
  sri = 3078,
  srd = 3715,
  srj = 3270,
  srj = 3519,
  str = 3452
)

# convert to data frame
ppa_collapse <- tibble::enframe(ppa_collapse) %>%
  rename(collapse_family = name, node = value)

# choose P. redivivus nodes to collapse
pre_collapse <- c(
  sru = 4157,
  srg = 4231,
  sri = 3191
)

# convert to data frame
pre_collapse <- tibble::enframe(pre_collapse) %>%
  rename(collapse_family = name, node = value)

# choose srw nodes to collapse
srw_collapse <- c(
  srw = 4852,
  srw = 5147,
  srw = 5143,
  srw = 5155,
  srw = 5208,
  srw = 5228
)

# convert to data frame
srw_collapse <- tibble::enframe(srw_collapse) %>%
  rename(collapse_family = name, node = value)

# bind data frames
collapse <- bind_rows(collapse, ppa_collapse, pre_collapse, srw_collapse)

# add to all data
data <- left_join(data, collapse)

# Tree plotting -----------------------------------------------------------

t <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none")

# attach family node/label information to the tree, for later labeling of collapsed nodes
t <- t %<+% data

# collapse node with loop
for (i in 1:length(collapse$node)) {
  t <- ggtree::collapse(t, collapse$node[i])
}

# superfamily hilighting
t_ann <- t +
  geom_cladelabel(node = 2882, label = "Str", barsize = 10, color = "black", offset = 9, fontsize = 18) + # Str
  
  geom_cladelabel(node = 2623, label = "Sra", barsize = 10, color = "black", offset = 9, fontsize = 18) + # Sra

  geom_cladelabel(node = 3765, label = "Srg", barsize = 10, color = "black", offset = 9, fontsize = 18) + # Srg
  
  geom_cladelabel(node = 4453, label = "sro", barsize = 10, color = "gray20", offset = 6, fontsize = 16) + # sro
  geom_cladelabel(node = 4628, label = "srxa", barsize = 10, color = "gray20", offset = 6, fontsize = 16) + # srxa
  geom_cladelabel(node = 4469, label = "srsx", barsize = 10, color = "gray20", offset = 6, fontsize = 16) + # srsx
  geom_cladelabel(node = 4687, label = "srbc", barsize = 10, color = "gray20", offset = 6, fontsize = 16) + # srbc
  geom_cladelabel(node = 4846, label = "srw", barsize = 10, color = "gray20", offset = 6, fontsize = 16) + # srw
  # geom_cladelabel(node = 4412, label = "srz", barsize = 10, color = "black", offset = 4) + # srz
  NULL

# family highlighting
t_ann <- t_ann +
  geom_cladelabel(node = 4196, label = "srg", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 4367, label = "srv", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 4313, label = "srv", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 4156, label = "sru", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 4196, label = "srg", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  # geom_cladelabel(node = 4137, label = "srr", barsize = 10, color = "gray40", offset = 2) +
  geom_cladelabel(node = 3968, label = "srt", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 3767, label = "srx", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 2858, label = "srb", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 2765, label = "sre", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 2626, label = "sra", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 2674, label = "srab", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 3076, label = "sri", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 3229, label = "srh", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  # geom_cladelabel(node = 3191, label = "srh", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  # geom_cladelabel(node = 2888, label = "srh", barsize = 10, color = "gray40", offset = 2) +
  geom_cladelabel(node = 3592, label = "srm", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 3635, label = "srd", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 3484, label = "srj", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  geom_cladelabel(node = 3288, label = "str", barsize = 10, color = "gray40", offset = 3, fontsize = 16) +
  NULL

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
t_ann$data <- t_ann$data %>%
  mutate(original.label = label) %>%
  mutate(label = collapse_family) %>%
  mutate(isTip = case_when(node %in% collapse$node ~ TRUE,
                           !is.na(Species) ~ TRUE))

# final tree
t_final <- t_ann +
  geom_tiplab2(aes(subset = !is.na(collapse_family), label = label), align = TRUE, size = 12, linesize = 1) +
  geom_tippoint(aes(color = Clade, angle= angle), size = 5) +
  geom_point2(aes(subset = !is.na(bootstrap_family)), size = 4, fill = "black") +
  geom_point2(aes(subset = !is.na(collapse_family)), size = 6, shape = 21, fill = "#ACA4E2") +
  geom_label2(aes(subset = !is.na(bootstrap_family), label = bootstrap, fill = bootstrap), size = 16) +
  labs(fill = "Bootstrap Support") +
  scale_fill_viridis_c(limits = c(0, 100)) +
  scale_color_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  theme(
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  NULL
t_final
# save_plot(here("plots", "Fig1A_raw.pdf"), t_final, base_height = 40)

# saveRDS(t_final, here("data", "final.tree"))


# Export data for heatmap generation --------------------------------------

family_nodes <- c(
  srw = 4961,
  pep = 4852,
  pep = 5143,
  pep = 5147,
  pep = 5155,
  pep = 5208,
  pep = 5228,
  srxa = 4628,
  srbc = 4687,
  srsx = 4469,
  sro = 4453,
  sro = 4451,
  srz = 4412,
  srv = 4369,
  srv = 1035,
  srv = 1081,
  srv = 1082,
  srv = 1083,
  srv = 4367,
  srv = 4358,
  srv = 4354,
  srv = 4313,
  srv = 4292,
  sru = 4156,
  srg = 4196,
  srr = 4137,
  srt = 3968,
  srx = 3767,
  srb = 2858,
  sre = 2765,
  sra = 2626,
  srab = 2674,
  srn = 3252,
  srh = 3251,
  srh = 3229,
  srh = 3191,
  srh = 2888,
  sri = 3076,
  srd = 3635,
  srd = 2007,
  four = 3631, 
  four = 3630,
  four = 3629,
  four = 3624,
  four = 3622,
  four = 3621,
  four = 1857,
  four = 3610,
  four = 1855,
  srm = 3592,
  three_four = 3577,
  three_four = 3557,
  srj = 3270,
  srj = 3544,
  srj = 3519,
  srj = 3484,
  str = 3288
)

species <- unique(t_final$data$Species)

des <- list()
n <- 1
for (x in family_nodes) {
  y <- geiger::tips(iqtree, x)
  des[n] <- list(y)
  n <- n + 1
}
names(des) <- names(family_nodes)
# get lengths of lists in list-des
len <- sapply(des, length)
n <- max(len)
len <- n - len
# write csv
csv <- mapply(function(x, y) c(x, rep(NA, y)), des, len)
csv.m <- pivot_longer(as.data.frame(csv), cols = everything(), names_to = "Family", values_to = "Transcript_ID") %>%
  select(-.copy) %>%
  separate(Transcript_ID, into = c("Species", "Transcript_ID"), sep = "-", extra = "merge") %>%
  filter(!is.na(Transcript_ID))

write.csv(csv.m, file = here("data", "tree_clades.csv"), row.names = FALSE, col.names = TRUE)
