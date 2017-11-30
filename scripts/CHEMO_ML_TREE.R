library(ape)
library(phytools)
library(phangorn)
library(seqinr)
library(ggtree)
library(tidyverse)
library(geiger)

setwd("~/Box Sync/GHdata/50HGI/NChR/phylo/ML/")

# parse RAxML output
raxml_file <- "RAxML_bipartitionsBranchLabels.final_ML"
raxml <- read.raxml(raxml_file)

# reroot based on the output of the ggtree command
# ggtree(raxml, size = 1.3,  layout = "circular", branch.length="none") + geom_label(aes(label=bootstrap, fill=bootstrap)) +  geom_tiplab2(size = 2) + geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)
raxml <-  ggtree::reroot(raxml, 2269)


# load in reference file matching species <-> clade
reference <- read.csv("~/Box Sync/ZamanianLab/Data/Phylogenetics/ChemoR/clade_species_phylum.csv", header = FALSE) %>%
  rename(Species = V1, Clade = V2, Phylum = V3)

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
# retain Transcript_ID for C. elegans genes
tip_labels <- as.data.frame(raxml@phylo$tip.label) %>%
  rename(tip_label = "raxml@phylo$tip.label")
tip_labels <- tidyr::separate(tip_labels, tip_label, into = c("Species", "Gene_ID"), sep = "-", remove = TRUE, extra = "merge")
tip_labels <- tip_labels %>%
  mutate(Transcript_ID = ifelse(Species %in% reference$Species, "", Species)) %>%
  mutate(Species = ifelse(Species %in% reference$Species, Species, "celeg")) # Species + Gene_ID + Clade for clade III, species for all other clades
tip_labels <- left_join(tip_labels, reference, by="Species")
tip_labels <- tip_labels %>%
  mutate(final_label = paste0(tip_labels$Species, "-", tip_labels$Transcript_ID, "-", tip_labels$Gene_ID, "-", tip_labels$Clade)) %>%
  select(-Phylum)
tip_labels <- data.frame(lapply(tip_labels, function(x) {
  gsub("--", "-", x)
}))

# replace the original tip labels
raxml@phylo$tip.label <- as.character(tip_labels$final_label)

# plot tree
t <- ggtree(raxml, size = 1.3,  layout = "circular", branch.length="none") #+ 
  # geom_tiplab2(size = 2) +
  # geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)

# attach annotation df to ggtree object
df <- tip_labels %>%
  mutate(label = final_label)
df <- left_join(t$data, df, by="label")
df <- select(df, -parent, -x, -y)

# add point colors for each clade
t_ann <- t %<+% df + 
  scale_color_manual(values = c(I = "#a361c7", IIIa = "#58a865", IIIb = "#c65c8a", IVb = "#9b9c3b", Va = "#648ace", Vb = "#c98443", Ve = "#cb4f42"))

# choose nodes to collapse, family names are below
nodes <- c(3798,
           3758,
           2206,
           2161,
           2051,
           2007,
           2092,
           2511,
           2417,
           2364,
           3630,
           2549,
           2769,
           2644,
           2691,
           3116,
           2854,
           3565,
           3194,
           3255
)









labels <- c("srw",
            "srsx",
            "srbc",
            "sre",
            "srb",
            "sra",
            "srab",
            "srxa",
            "srx",
            "srx",
            "srt",
            "srz",
            "srv",
            "sru",
            "srg",
            "sri",
            "srh",
            "srd",
            "srj",
            "str"
)

# collapse node with loop
for(i in 1:length(nodes)){
  t_ann <- ggtree::collapse(t_ann, nodes[i])
}

# create data frame from family nodes and labels
c_df <- data.frame(node = nodes, family = labels)

# attach family node/label information to the tree, for later labeling of collapsed nodes
t_ann <- t_ann %<+% c_df 

# superfamily hilighting
t_ann <- t_ann + 
  geom_cladelabel(node=2814, label = "", barsize = 5, color = "#5fa271", offset = 2) + #Str

  geom_cladelabel(node=2003, label = "", barsize = 5, color = "#708cc9", offset = 2) + #Sra

  geom_cladelabel(node=2642, label = "", barsize = 5, color = "#c8615d", offset = 2) + #Srg
  geom_cladelabel(node=3606, label = "", barsize = 5, color = "#c8615d", offset = 2) + #Srg
  geom_cladelabel(node=2293, label = "", barsize = 5, color = "#c8615d", offset = 2) + #Srg
  geom_cladelabel(node=3274, label = "", barsize = 5, color = "#c8615d", offset = 2) + #Srg

  geom_cladelabel(node=3692, label = "", barsize = 5, color = "#88a83d", offset = 2) + #srsx
  geom_cladelabel(node=2199, label = "", barsize = 5, color = "#9b63ca", offset = 2) + #srbc
  geom_cladelabel(node=3784, label = "", barsize = 5, color = "#c55d93", offset = 2) + #srw
  geom_cladelabel(node=2549, label = "", barsize = 5, color = "#c77332", offset = 2) #srz

# filarid clusters
t_ann <- t_ann + 
  geom_cladelabel(node = 3779, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3756, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3707, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3727, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2202, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2080, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2119, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2535, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2496, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2506, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2310, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2353, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2317, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2323, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3624, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2845, label = "", offset = 1, barsize = 5, color = "#c65c8a") 

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
t_ann$data <- t_ann$data %>%
  # mutate(final_label = label) %>%
  # mutate(label = ifelse(Clade == "IIIb", as.character(final_label), as.character(Species))) # Species + Gene_ID + Clade for clade III, species for all other clades
  mutate(label = ifelse(Clade == "IIIb", as.character(Gene_ID), "")) # Species + Gene_ID + Clade for clade III, blank for all other clades

# selected nodes to display bootstrap
bs <- raxml@bootstrap
nodes <- c(3784, #srw
        3692, #srsx
        2199, #srbc
        2002, #Sra
        2137, #sre
        2049, #srb
        2006, #sra
        2063, #srab
        2510, #srxa
        2482,
        2496,
        2502,
        2293, #Srg
        2295, #srx
        3606, #srt
        2549, #srz
        2763, #srv
        2644, #sru
        2691, #srg
        2814, #Str
        3116, #sri
        2817, #srh
        3522, #srd
        3192, #srj
        3255) #str

bs2 <- bs[bs$node %in% nodes,] %>%
  rename(BS = bootstrap)
bs2$BS <- as.numeric(bs2$BS)
t_bs <- t_ann %<+% bs2

# final tree
t_final <- t_bs +
  geom_tiplab2(aes(label = family), size = 6, align = TRUE, linetype = 5, hjust = -0.25) +
  # geom_tiplab2(size = 4, align = TRUE, hjust = -1) +
  theme(legend.position="right", legend.title = element_text(size = 40), legend.text = element_text(size = 35), legend.key.size = unit(2, "cm")) +
  geom_tippoint(aes(color = Clade), size = 5) +
  geom_point2(aes(subset=(is.na(family) == FALSE)), size = 8, shape = 21, fill = "yellow2") +
  geom_label2(aes(na.rm=TRUE, label = BS, fill = BS)) +
  # geom_label(aes(label=bootstrap, fill=bootstrap)) +
  scale_fill_continuous(low='darkgreen', high='red')
# t_final

ggsave("CHEMO_ML_TREE_label.pdf", t_final, width = 40, height = 40)

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

















