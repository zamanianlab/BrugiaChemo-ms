library(ape)
library(phytools)
library(phangorn)
library(seqinr)
library(ggtree)
library(tidyverse)
library(geiger)

setwd("~/Box Sync/GHdata/50HGI/NChR/phylo/")
  
# read in alignment, calculate distance matrix, infer neighbor-joining tree
aln <- read.alignment("final.trim.filter.aln", format = "fasta")
dat <- as.phyDat(aln, type = "AA")
# mt <- modelTest(dat, model = "all", FREQ = TRUE, multicore = TRUE, mc.cores = 4)
aa_dist <- dist.ml(dat, model = "VT") # can change model, VT was chosen by RAxML for the C. elegans tree. can also test models with command above
NJ  <- NJ(aa_dist)

# reroot based on the output of the ggtree command
ggtree(NJ, size = 1.3,  layout = "circular", branch.length="none") +  geom_tiplab2(size = 2) + geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)
NJ <-  phytools::reroot(NJ, 2716)


library("polysat")
# 100 NJ trees
for(i in 1:100){
  
}

# load in reference file matching species <-> clade
reference <- read.csv("~/Box Sync/ZamanianLab/Data/Phylogenetics/ChemoR/clade_species_phylum.csv", header = FALSE) %>%
  rename(Species = V1, Clade = V2, Phylum = V3)

# add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
tip_labels <- as.data.frame(NJ[]$tip.label) %>%
  rename(tip_label = "NJ[]$tip.label")
tip_labels <- tidyr::separate(tip_labels, tip_label, into = c("Species", "Gene_ID"), sep = "-", remove = TRUE, extra = "merge")
tip_labels <- left_join(tip_labels,reference,by="Species")
tip_labels <- tip_labels %>%
  mutate(final_label = paste0(tip_labels$Species,"-",tip_labels$Gene_ID,"-",tip_labels$Clade)) %>%
  select(-Phylum)

# replace the original tip labels
NJ[]$tip.label <- tip_labels$final_label

# plot tree
t <- ggtree(NJ, size = 1.3,  layout = "circular", branch.length="none") + 
  geom_tiplab2(size = 2) +
  geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)

# attach annotation df to ggtree object
df <- tip_labels %>%
  mutate(label = final_label)
df <- left_join(t$data, df, by="label")
df <- select(df, -parent, -branch.length, -x, -y)

# add point colors for each clade
t_ann <- t %<+% df + 
  scale_color_manual(values = c(I = "#a361c7", IIIa = "#58a865", IIIb = "#c65c8a", IVb = "#9b9c3b", Va = "#648ace", Vb = "#c98443", Ve = "#cb4f42"))

# choose nodes to collapse, family names are below
nodes <- c(3905,
           2766,
           2479,
           2422,
           2330,
           2063,
           3831,
           3773,
           3716,
           3647,
           3628,
           3452,
           3594,
           3518,
           3329,
           3235,
           2861,
           3192,
           3086,
           2962
           )

labels <- c("srw",
            "sri",
            "srh",
            "srd",
            "srj",
            "str",
            "srt",
            "srv",
            "sru",
            "srg",
            "srb",
            "sre",
            "sra",
            "srab",
            "srz",
            "srbc",
            "srsx",
            "srxa",
            "srx",
            "srx"
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
  geom_hilight(node=2057, fill = "yellow" , alpha = 0.5) + #Str
  geom_hilight(node=3423, fill = "blue" , alpha = 0.5) + #Sra
  
  geom_hilight(node=3644, fill = "red" , alpha = 0.5) + #Srg
  geom_hilight(node=3647, fill = "red" , alpha = 0.5) + #Srg
  geom_hilight(node=3527, fill = "red" , alpha = 0.5) + #Srg
  geom_hilight(node=3274, fill = "red" , alpha = 0.5) + #Srg
  geom_hilight(node=3212, fill = "red" , alpha = 0.5) + #Srg
  geom_hilight(node=2893, fill = "green" , alpha = 0.5) + #srsx
  geom_hilight(node=3119, fill = "green" , alpha = 0.5) + #srbc
  geom_hilight(node=3784, fill = "green" , alpha = 0.5) + #srw
  geom_hilight(node=2795, fill = "green" , alpha = 0.5) #srz

# filarid clusters
t_ann <- t_ann + geom_cladelabel(node = 3460, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3430, label = "", offset = 1, barsize = 5, color = "#c65c8a") + 
  geom_cladelabel(node = 3424, label = "", offset = 1, barsize = 5, color = "#c65c8a") + 
  geom_cladelabel(node = 3420, label = "", offset = 1, barsize = 5, color = "#c65c8a") + 
  geom_cladelabel(node = 3346, label = "", offset = 1, barsize = 5, color = "#c65c8a") + 
  geom_cladelabel(node = 3358, label = "", offset = 1, barsize = 5, color = "#c65c8a") + 
  geom_cladelabel(node = 2927, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2965, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 2952, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3208, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3052, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3056, label = "", offset = 1, barsize = 5, color = "#c65c8a") +
  geom_cladelabel(node = 3011, label = "", offset = 1, barsize = 5, color = "#c65c8a") 

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
t_ann$data <- t_ann$data %>%
  # mutate(final_label = label) %>%
  mutate(label = ifelse(Clade == "IIIb",final_label, Species)) # Species + Gene_ID + Clade for clade III, species for all other clades

# final tree
t_final <- t_ann +
  geom_tiplab2(aes(label = family), size = 6, align = TRUE, linetype = 5, hjust = -0.25) +
  # geom_tiplab2(size = 4, align = TRUE, linetype = 5, hjust = -.5) +
  theme(legend.position="right", legend.title = element_text(size = 40), legend.text = element_text(size = 35), legend.key.size = unit(2, "cm")) +
  geom_tippoint(aes(color = Clade), size = 5) +
  geom_point2(aes(subset=(is.na(family) == FALSE)), size = 8, shape = 21, fill = "yellow2")
t_final

ggsave("CHEMO_NJ_TREE_label.pdf", t_final, width = 40, height = 40)

# add node labels
t_final_nodes <- t_final + geom_text2(aes(subset=!isTip, label = node), hjust = -0.3)
ggsave("CHEMO_NJ_TREE_nodes.pdf", t_final_nodes, width = 40, height = 40)


############################################
## Get taxon labels for a particular node ##
############################################



node <- c(S1 = 3713,
          S2 = 3460,
          S3 = 3430,
          S4 = 3424,
          S5 = 3420,
          S6 = 3341,
          L1 = 3346,
          L2 = 3358,
          S7 = 2911,
          L3 = 2927,
          S8 = 2965,
          S9 = 2952,
          S10 = 3208,
          S11 = 3050,
          L4 = 3011,
          S12 = 2285)

# create list of lists that include gene names grouped by clade
des <- list()
n <- 1
for(x in node) {
  y <- geiger::tips(NJ, x)
  des[n] <- list(y)
  n <- n + 1
}
names(des) <- names(node)
# get lengths of lists in list-des
len <- sapply(des, length)
n <- max(len)
len <- n - len
# write csv
csv <- mapply(function(x,y) c( x , rep( NA , y ) ) , des , len )
csv.m <- gather(as.data.frame(csv), key = "Clade", value = "Gene_ID", na.rm = TRUE)
write.csv(csv.m, file="clades.csv", row.names = FALSE, col.names = TRUE)




















