library(ape)
library(phytools)
library(phangorn)
library(seqinr)
library(tidyverse)
library(ggtree)
library(geiger)
library(conflicted)

setwd("~/Box Sync/ZamanianLab/Data/Genomics/Phylogenetics/TRP/")

#read in tree
bayes_file <- "TRP.con.tre"
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
  geom_text2(aes(subset =! isTip, label = node), hjust = -.3) + 
  geom_tiplab2() +
  NULL

# ggsave("unrooted.pdf", unrooted, width = 40, height = 40)

# reroot with CUP-5 clade as outgroup
bayes <- phytools::reroot(bayes_drop, 400)

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
c_df <- c(cup5 = 708, 
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
          trp2 = 472)
c_df <- as.data.frame(c_df)

prob3 <- prob[prob$node %in% c_df$c_df,]

# superfamily hilighting
t_ann <- t_ann +
  geom_cladelabel(node=708, label = "TRPML", barsize = 5, color = "#84ad3e", offset = 4, fontsize = 16, offset.text = 4) +
  geom_cladelabel(node=694, label = "TRPP", barsize = 5, color = "#9d63cb", offset = 4, fontsize = 16, offset.text = 4) +
  geom_cladelabel(node=659, label = "TRPA", barsize = 5, color = "#4aac8b", offset = 4, fontsize = 16, offset.text = 4) +
  geom_cladelabel(node=574, label = "TRPV", barsize = 5, color = "#c95a80", offset = 4, fontsize = 16, offset.text = 4) +
  geom_cladelabel(node=376, label = "TRPM", barsize = 5, color = "#8f8c4c", offset = 4, fontsize = 16, offset.text = 4) +
  geom_cladelabel(node=544, label = "TRPN", barsize = 5, color = "#7787c5", offset = 4, fontsize = 16, offset.text = 4) +
  geom_cladelabel(node=469, label = "TRPC", barsize = 5, color = "#ce6c38", offset = 4, fontsize = 16, offset.text = 4) +
  NULL
  
# family highlighting
t_ann <- t_ann + 
  geom_cladelabel(node = 708, label = "CUP-5", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 694, label = "PKD-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 661, label = "TRPA-1", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 668, label = "TRPA-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 635, label = "OSM-9", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 631, label = "OCR-4", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 619, label = "OCR-3", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 583, label = "OCR-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 587, label = "OCR-1", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 593, label = "Clade III OCR-1/2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 447, label = "CED-11", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 420, label = "GTL-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 404, label = "Clade III GON-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 399, label = "GON-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 386, label = "GTL-1", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 544, label = "TRP-4", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 510, label = "TRP-3/SPE-41", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 495, label = "TRP-1", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  geom_cladelabel(node = 472, label = "TRP-2", barsize = 5, color = "black", offset = 1, fontsize = 10, offset.text = 2) +
  NULL

# geom_tiplab2 uses column label by default, so copy whatever data you want to the label column
t_ann$data <- t_ann$data %>%
  # mutate(final_label = label) %>%
  mutate(label = ifelse(Clade == "IIIb",final_label, Species)) # Species + Gene_ID + Clade for clade III, species for all other clades

t_final <- t_ann %<+% prob3 + 
  geom_label2(aes(na.rm = TRUE, label = round(as.numeric(prob) * 100), fill = as.numeric(prob) * 100), size = 8) + 
  scale_fill_viridis_c(limits = c(0, 100)) + 
  geom_tippoint(aes(color=Clade), size = 6) +
  labs(fill = "Posterior Probability") +
  theme(legend.position="right", legend.title = element_text(size = 40), legend.text = element_text(size = 35), legend.key.size = unit(2, "cm")) +
  NULL
# t_final


ggsave("TRP_MB_TREE.pdf", t_final, width = 40, height = 40)



