library(ape)
library(phytools)
library(phangorn)
library(seqinr)
library(ggtree)
library(tidyverse)
library(geiger)

setwd("~/Box/ZamanianLab/Data/Genomics/Phylogenetics/ChemoR/tree")

# parse RAxML output
raxml_file <- "RAxML_bipartitionsBranchLabels.final_ML"
raxml <- read.raxml(raxml_file)

# reroot based on the output of the ggtree command
# ggtree(raxml, size = 1.3,  layout = "circular", branch.length="none") + geom_label(aes(label=bootstrap, fill=bootstrap)) +  geom_tiplab2(size = 2) + geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)
raxml@phylo <-  phytools::reroot(raxml@phylo, 2269)

# load in reference file matching species <-> clade
reference <- read.csv("clade_species_phylum.csv", header = FALSE) %>%
  rename(Species = V1, Clade = V2, Phylum = V3)

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
# retain Transcript_ID for C. elegans genes
tip_labels <- as.data.frame(raxml@phylo$tip.label) %>%
  rename(tip_label = "raxml@phylo$tip.label") %>%
  tidyr::separate(tip_label, into = c("Species", "Gene_ID"), sep = "-", extra = "merge", remove = FALSE) %>%
  mutate(Transcript_ID = ifelse(Species %in% reference$Species, "", Species)) %>%
  mutate(Species = ifelse(Species %in% reference$Species, Species, "celeg")) %>% # Species + Gene_ID + Clade for clade III, species for all other clades 
  rename(label = tip_label) %>%
  left_join(., reference, by="Species") %>%
  dplyr::select(-Phylum)

# basic tree with node labels
unrooted <- ggtree(raxml, branch.length = "none", layout = "circular") + 
  geom_text2(aes(subset =! isTip, label = node), hjust = -.3) + 
  geom_tiplab2() +
  NULL

# ggsave("unrooted.pdf", unrooted, width = 40, height = 40)

# plot tree
t <- ggtree(raxml, size = 1.3,  layout = "circular", branch.length="none") #+ 
  # geom_tiplab2(size = 2) +
  # geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)

node_labels <- t +
  geom_tiplab2(size = 2) +
  geom_text2(aes(subset =! isTip, label = node), size = 2, hjust=-.3) +
  NULL

# ggsave("node_labels.pdf", node_labels, width = 40, height = 40)

bootstrap <- t +
  geom_tiplab2(size = 2) +
  geom_label2(na.rm = TRUE, aes(label = bootstrap, fill = bootstrap)) +
  NULL

# ggsave("bootstrap.pdf", bootstrap, width = 40, height = 40)

# attach annotation df to ggtree object
data <- left_join(t$data, tip_labels, by = "label") %>%
  dplyr::select(-parent, -x, -y)

# choose celeg nodes to collapse
families <- c(srw = 3798,
              srsx = 3758,
              srbc = 2206,
              sre = 2161,
              srb = 2051,
              sra = 2007,
              srab = 2092,
              srxa = 2511,
              srx = 2417,
              srx = 2364,
              srt = 3630,
              srz = 2549,
              srv = 2769,
              sru = 2644,
              srg = 2691,
              sri = 3116,
              srh = 2854,
              srd = 3565,
              srj = 3194,
              str = 3255
)

# convert to data frame
celeg_clusters <- tibble::enframe(families) %>%
  rename(celeg_family = name, node = value)

# add to all data
data <- left_join(data, celeg_clusters)
data$celeg_family[is.na(data$celeg_family)] <- ""

# select nodes to display bootstrap
bs <- raxml@data

bootstrap_families <- c(srw = 3796,
                        srsx = 3692,
                        srbc = 2199,
                        sre = 2138,
                        srb = 2046,
                        sra = 2006,
                        srab = 2092,
                        srxa = 2511,
                        srx = 2294,
                        # srx = 2364,
                        srt = 3606,
                        srz = 2549,
                        srv = 2763,
                        sru = 2644,
                        srg = 2691,
                        sri = 3114,
                        srh = 2817,
                        srd = 3520,
                        srj = 3192,
                        str = 3255
)

bootstrap_families <- tibble::enframe(bootstrap_families) %>%
  rename(bootstrap_family = name, node = value) %>%
  dplyr::select(node = node, bootstrap_family = bootstrap_family)

bs2 <- left_join(bootstrap_families, bs) %>%
  rename(BS = bootstrap)
bs2$BS <- as.numeric(bs2$BS)

data <- left_join(data, bs2)

# attach family node/label information to the tree, for later labeling of collapsed nodes
t_ann <- t %<+% data 

# collapse node with loop
for(i in 1:length(celeg_clusters$node)){
  t_ann <- ggtree::collapse(t_ann, celeg_clusters$node[i])
}

# superfamily hilighting
t_ann <- t_ann + 
  geom_cladelabel(node=2814, label = "", barsize = 5, color = "black", offset = 2) + #Str
  geom_cladelabel(node=2003, label = "", barsize = 5, color = "black", offset = 2) + #Sra
  
  geom_cladelabel(node=2642, label = "", barsize = 5, color = "black", offset = 2) + #Srg
  geom_cladelabel(node=3606, label = "", barsize = 5, color = "black", offset = 2) + #Srg
  geom_cladelabel(node=2293, label = "", barsize = 5, color = "black", offset = 2) + #Srg
  geom_cladelabel(node=3274, label = "", barsize = 5, color = "black", offset = 2) + #Srg
  
  geom_cladelabel(node=3692, label = "", barsize = 5, color = "black", offset = 2) + #srsx
  geom_cladelabel(node=2199, label = "", barsize = 5, color = "black", offset = 2) + #srbc
  geom_cladelabel(node=3784, label = "", barsize = 5, color = "black", offset = 2) + #srw
  geom_cladelabel(node=2549, label = "", barsize = 5, color = "black", offset = 2) + #srz
  NULL

# filarid clusters
# t_ann <- t_ann + 
#   geom_cladelabel(node = 3779, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 3756, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 3707, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 3727, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2202, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2080, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2119, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2535, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2496, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2506, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2310, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2353, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 2312, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   # geom_cladelabel(node = 2317, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   # geom_cladelabel(node = 2323, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   geom_cladelabel(node = 3624, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   # geom_cladelabel(node = 2845, label = "", barsize = 5, color = "#c65c8a", offset = 1) +
#   NULL

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
t_ann$data <- t_ann$data %>%
  mutate(original_label = label) %>%
  mutate(label = celeg_family) 

# final tree
t_final <- t_ann +
  geom_tiplab2(size = 6, align = TRUE, linetype = 5, hjust = -0.25) +
  geom_tippoint(aes(color = Clade), size = 6) +
  geom_point2(aes(subset = celeg_family != ""), size = 15, shape = 21, fill = "yellow2") +
  geom_label2(aes(na.rm = TRUE, label = BS, fill = BS), size = 15) +
  labs(fill = "Nodal Support") +
  scale_fill_viridis_c(limits = c(0, 100)) + 
  scale_color_manual(values = c(I = "#a361c7", IIIa = "#58a865", IIIb = "#c65c8a", IVb = "#9b9c3b", Va = "#648ace", Vb = "#c98443", Ve = "#cb4f42")) +
  theme(legend.position = "right") +
  NULL
# t_final

ggsave("CHEMO_ML_TREE.pdf", t_final, width = 40, height = 40, useDingbats = FALSE)

# add node labels
t_final_nodes <- t_final + geom_text2(aes(subset=!isTip, label = node), hjust = -0.3)
ggsave("CHEMO_ML_TREE_nodes.pdf", t_final_nodes, width = 40, height = 40)


############################################
## Get taxon labels for a particular node ##
############################################

filarid_nodes <- c(srsx1 = 3779,
          srsx2 = 3756,
          srsx3 = 3707,
          srsx4 = 3727,
          srbc = 2202,
          srab1 = 2080,
          srab2 = 2119,
          srxa = 2535,
          srx1 = 2496,
          srx2 = 2506,
          srx3 = 2310,
          srx4 = 2353,
          srx5 = 2317,
          srx6 = 2323,
          srt = 3624,
          srh = 2845)

# create list of lists that include gene names grouped by clade
des <- list()
n <- 1
for(x in filarid_nodes) {
  y <- geiger::tips(raxml@phylo, x)
  des[n] <- list(y)
  n <- n + 1
}
names(des) <- names(filarid_nodes)
# get lengths of lists in list-des
len <- sapply(des, length)
n <- max(len)
len <- n - len
# write csv
csv <- mapply(function(x,y) c( x , rep( NA , y ) ) , des , len )
csv.m <- gather(as.data.frame(csv), key = "Clade", value = "Gene_ID", na.rm = TRUE)
write.csv(csv.m, file="filarid_clades.csv", row.names = FALSE, col.names = TRUE)

family_nodes <- c(srw = 3797,
                  pep1 = 3991,
                  pep2 = 3792,
                  pep3 = 3785,
                  srsx = 3692,
                  srbc = 2199,
                  sre = 2137,
                  srb = 2046,
                  sra = 2006,
                  srab = 2063,
                  srxa = 2510,
                  srx = 2294,
                  srt = 3606,
                  srz = 2549,
                  srv = 2763,
                  sru = 2644,
                  srg = 2691,
                  sri = 3114,
                  srh = 2817,
                  srd = 3520,
                  srj = 3192,
                  str = 3255
)

des <- list()
n <- 1
for(x in family_nodes) {
  y <- geiger::tips(raxml@phylo, x)
  des[n] <- list(y)
  n <- n + 1
}
names(des) <- names(family_nodes)
# get lengths of lists in list-des
len <- sapply(des, length)
n <- max(len)
len <- n - len
# write csv
csv <- mapply(function(x,y) c( x , rep( NA , y ) ) , des , len )
csv.m <- gather(as.data.frame(csv), key = "Clade", value = "Gene_ID", na.rm = TRUE)
write.csv(csv.m, file="family_clades.csv", row.names = FALSE, col.names = TRUE)

















