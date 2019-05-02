library(ape)
library(phytools)
library(phangorn)
library(seqinr)
library(tidyverse)
library(ggtree)
library(geiger)
library(conflicted)

setwd("~/Box Sync/ZamanianLab/Data/Genomics/Phylogenetics/CNG/")

#read in tree
bayes_file <- "CNG.con.tre"
bayes <- ape::read.nexus(bayes_file)
bayes <- bayes[[1]]

# basic tree with node labels
unrooted <- ggtree(bayes, branch.length = "none", layout = "circular") + 
  geom_text2(aes(subset =! isTip, label = node), hjust = -.3) + 
  geom_tiplab2() +
  NULL

# ggsave("unrooted.pdf", unrooted, width = 40, height = 40)

# reroot
bayes <- phytools::reroot(bayes, 119)

#load in reference file matching species <-> clade
reference <- read.csv("clade_species_phylum.csv", header = FALSE)
colnames(reference) <- c("Species", "Clade", "Phylum")

#add species name to tip label, options to include: Species."-".Gene_ID."-".Clade."-".Phylum
tip_labels <- as.data.frame(bayes[]$tip.label) %>%
  rename(tip_label = "bayes[]$tip.label")
tip_labels <- tidyr::separate(tip_labels, tip_label, into = c("Species", "Gene_ID"), sep = "-", remove = TRUE, extra = "merge")
tip_labels <- left_join(tip_labels, reference, by="Species")
tip_labels <- tip_labels %>%
  mutate(final_label = paste0(tip_labels$Species,"-",tip_labels$Gene_ID,"-",tip_labels$Clade)) %>%
  select(-Phylum)

#replace the original tip labels
bayes[]$tip.label <- tip_labels$final_label


#Create a df containing internal node numbers and probs, as well as terminal node numbers and NA (for labels)
node_prob <- data.frame(node = seq(length(tip_labels$final_label) + 1 , length(tip_labels$final_label) + bayes$Nnode), label = bayes$node.label)
node_prob$label[node_prob$label == "Root"] <- NA
node <- seq(1, length(tip_labels$final_label))
label <- c(rep(NA, length(tip_labels$final_label)))
prob2 <- data.frame(node, label)
prob <- rbind(prob2, node_prob) %>%
  rename(prob = label)

#############################
#create entire ggtree object
#############################

t <- ggtree(bayes, layout = "circular", branch.length="none") #+ geom_tiplab2(size = 2) + geom_text2(aes(subset=!isTip, label=node), size = 2, hjust=-.3)

#attach annotation df to ggtree object
df <- tip_labels %>%
  mutate(label = final_label)
df <- left_join(t$data, df, by="label")
df <- select(df, -parent, -branch.length, -x, -y)

# add point colors for each clade
t_ann <- t %<+% df + 
  scale_color_manual(values = c(I = "#a361c7", IIIa = "#58a865", IIIb = "#c65c8a", IVa = "black", IVb = "#9b9c3b", Va = "#648ace", Vb = "#c98443", Vd = "yellow", Ve = "#cb4f42"))

t_ann$data <- t_ann$data %>%
  #mutate(final_label = label) #%>%
  mutate(label = ifelse(Clade == "IIIb", final_label, " "))

#add support values to nodes (for manual selection of nodes)
tbs <- t %<+% prob + 
  geom_tiplab2() +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  geom_label2(aes(subset=!is.na(prob), label = prob, fill = prob), size = 2) +
  NULL

# ggsave("bootstrap.pdf", tbs, width = 40, height = 40)

t_nodes <- t +
  geom_tiplab2() +
  geom_text2(aes(subset =! isTip, label = node), hjust = -.3) + 
  NULL

# ggsave("node_labels.pdf", t_nodes, width = 40, height = 40)

#selected nodes based on support values
c_df <- c(tax2 = 168,
          tax4 = 142,
          cng1 = 125,
          cng3 = 195,
          iiicng2 = 204,
          ivcng2 = 217,
          vcng2 = 114)
c_df <- as.data.frame(c_df)

prob3 <- prob[prob$node %in% c_df$c_df,]

# superfamily hilighting
# t_ann <- t_ann +
#   geom_cladelabel(node=708, label = "TRPML", barsize = 5, color = "#84ad3e", offset = 4, fontsize = 16, offset.text = 4) +
#   geom_cladelabel(node=694, label = "TRPP", barsize = 5, color = "#9d63cb", offset = 4, fontsize = 16, offset.text = 4) +
#   geom_cladelabel(node=659, label = "TRPA", barsize = 5, color = "#4aac8b", offset = 4, fontsize = 16, offset.text = 4) +
#   geom_cladelabel(node=574, label = "TRPV", barsize = 5, color = "#c95a80", offset = 4, fontsize = 16, offset.text = 4) +
#   geom_cladelabel(node=376, label = "TRPM", barsize = 5, color = "#8f8c4c", offset = 4, fontsize = 16, offset.text = 4) +
#   geom_cladelabel(node=544, label = "TRPN", barsize = 5, color = "#7787c5", offset = 4, fontsize = 16, offset.text = 4) +
#   geom_cladelabel(node=469, label = "TRPC", barsize = 5, color = "#ce6c38", offset = 4, fontsize = 16, offset.text = 4) +
#   NULL
  
# family highlighting
t_ann <- t_ann + 
  geom_cladelabel(node = 168, label = "TAX-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 142, label = "TAX-4", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 125, label = "CNG-1", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 195, label = "CNG-3", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 204, label = "Clade III CNG-2/CHE-6", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 217, label = "Clade IV CNG-2/CHE-6", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 114, label = "Clade V CNG-2/CHE-6", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  NULL

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
t_ann$data <- t_ann$data %>%
  # mutate(final_label = label) %>%
  mutate(label = ifelse(Clade == "IIIb", Species, "")) # Species + Gene_ID + Clade for clade III, species for all other clades

t_final <- t_ann %<+% prob3 + 
  geom_label2(aes(na.rm = TRUE, label = round(as.numeric(prob) * 100), fill = as.numeric(prob) * 100), size = 8) + 
  geom_tiplab2(hjust = -1) +
  scale_fill_viridis_c(limits = c(0, 100)) + 
  geom_tippoint(aes(color=Clade), size = 6) +
  labs(fill = "Posterior Probability") +
  theme(legend.position="right", legend.title = element_text(size = 40), legend.text = element_text(size = 35), legend.key.size = unit(2, "cm")) +
  NULL
t_final


ggsave("CNG_MB_TREE.pdf", t_final, width = 40, height = 40)



