library(tidyverse)
library(cowplot)
library(ggdendro)
library(reshape2)
library(scales)
library(gridExtra)
library(viridis)
library(here)

here()

########################################################
####### Load in HT data and merge with UGA data, filter for TPM
########################################################

#### LOAD HT EXPRESSION FILE
load(here("data", "RNAsamples_HT.Rda"))
RNAsamples.HT <- RNAsamples.melt

#### LOAD UGA EXPRESSION FILE
load(here("data", "RNAsamples_UGA.Rda"))
RNAsamples.UGA <- RNAsamples.melt

#### FILTER for TPM (gene-specific measurements)
RNAsamples.HT <- RNAsamples.HT %>%
  filter(Metric == "TPM_g")
RNAsamples.UGA <- RNAsamples.UGA %>%
  filter(Metric == "TPM_g") %>%
  filter(Condition == "CON", Stage != "MF") %>%
  group_by(Transcript_ID,Stage) %>%
  summarize(Whole_Expression = mean(Expression)) #whole-body expression

RNAsamples.HT.UGA <- inner_join(RNAsamples.HT,RNAsamples.UGA)

#Made new dataframe of just whole body and re-merge (so that UGA data gets own line)
temp <- RNAsamples.HT.UGA
temp <- temp %>%
  mutate(Expression = Whole_Expression) %>%
  select(-Whole_Expression) %>%
  mutate(Location = "WB") %>%
  mutate(SID = ifelse(grepl("Bm-F-Head", SID),"Bm-F-WB", SID)) %>%
  mutate(SID = ifelse(grepl("Bm-M-Head", SID),"Bm-M-WB", SID)) %>%
  mutate(SID = ifelse(grepl("Bm-M-Tail", SID),"Bm-M-WB", SID)) %>%
  group_by(Transcript_ID, SID) %>%
  distinct(.keep_all=True) %>% ungroup()
temp2 <- RNAsamples.HT.UGA %>% select(-Whole_Expression)
RNAsamples.HT.UGA <- rbind(temp,temp2)

#load in new Bm ChemoR list (includes)
ChemoR_list = read.csv("~/Box/ZamanianLab/Data/Genomics/ChemoR/phylo/family_assignment.csv",header = TRUE) 
ChemoR_list <- ChemoR_list %>%
  filter(Species == "brugia_malayi") %>%
  select(-Species)
gene_list <- as.character(ChemoR_list$Transcript_ID)
genecount<-as.integer(length(gene_list))

##########################
# PLOT male vs female head
##########################

RNAsamples.MF <- RNAsamples.HT.UGA %>%
  filter(Location == "H") %>%
  filter(Transcript_ID %in% gene_list) %>%
  group_by(Transcript_ID, Stage, Location) %>%
  distinct(.keep_all=TRUE) %>% ungroup() %>%
  spread(Stage, Expression)
RNAsamples.MF.M <- RNAsamples.MF %>% select(-AF) %>% filter(!is.na(AM))
RNAsamples.MF.F <- RNAsamples.MF %>% select(-AM) %>% filter(!is.na(AF)) %>% select(Transcript_ID, AF)
RNAsamples.MF <- left_join(RNAsamples.MF.M,RNAsamples.MF.F)

#join chemoR family data
RNAsamples.MF <- left_join(RNAsamples.MF,ChemoR_list)
RNAsamples.MF$Family <- factor(RNAsamples.MF$Family, levels = c("srw","srsx","srbc","srab","srt","srx","srxa","srh"))

RNAsamples.MF.2 <- RNAsamples.MF %>%
  mutate(FC = (AF+.1)/(AM+.1)) %>%
  group_by(Gene_ID) %>%
  distinct(Gene_ID,.keep_all=TRUE) %>%
  mutate(Max = max(AF,AM)) %>%
  filter((AF+AM) > 1) #total TPM greater than one


MFplot <- ggplot(RNAsamples.MF.2)+
  aes(x = Superfamily, y = log2(FC))+ 
  geom_jitter(aes(size=Max,colour = Superfamily), alpha = 0.75, width=0.05) +
  theme_minimal() +
  #geom_hline(aes(yintercept=5.6), linetype="solid", colour = "black", size = 0.5) +
  #geom_hline(aes(yintercept=-3.5), linetype="solid", colour = "black", size = 0.5) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour = "grey37", size = 0.5) +
  scale_y_continuous(breaks=c(-3, -2, -1, 0, 1, 2, 3, 4, 5)) +
  theme(axis.title.x= element_blank(),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12),
        axis.ticks.x= element_blank()) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_line(size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank()) +
  scale_size_continuous(name="Expression (TPM)") +
  scale_colour_discrete(name="Superfamily") + # set matching colors later
  geom_text(mapping=aes(x='', y=1, label=c("♀")), colour = "black", fontface = "bold", size=9) +
  geom_text(mapping=aes(x='', y=2, label=c("head")), colour = "black", fontface = "bold", size=5) +
  geom_text(mapping=aes(x='', y=2.5, label=c("↑")), colour = "black", fontface = "bold", size=6) +
  geom_text(mapping=aes(x='', y=-1, label=c("♂")), colour = "black", fontface = "bold", size=9) +
  geom_text(mapping=aes(x='', y=-2, label=c("head")), colour = "black", fontface = "bold", size=5) +
  geom_text(mapping=aes(x='', y=-2.5, label=c("↓")), colour = "black", fontface = "bold", size=6) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  theme(legend.position="right")
MFplot

#Identify transcripts enriched in male head (FC < 1)
Mlist <- RNAsamples.MF.2 #%>%
  #filter(FC < 1 )
Mlist <- Mlist$Transcript_ID

##########################
# PLOT male head vs tail 
##########################

RNAsamples.HT <- RNAsamples.HT.UGA %>%
  filter(Stage == "AM") %>% filter(Location != "WB") %>%
  #filter(Transcript_ID %in% gene_list) %>%
  filter(Transcript_ID %in% Mlist) %>% #just look at Male enriched
  group_by(Transcript_ID, Stage, Location) %>%
  distinct(.keep_all=TRUE) %>% ungroup() %>%
  spread(Location, Expression)

RNAsamples.HT.H <- RNAsamples.HT %>% select(-T) %>% filter(!is.na(H))
RNAsamples.HT.T <- RNAsamples.HT %>% select(-H) %>% filter(!is.na(T)) %>% select(Transcript_ID, T)
RNAsamples.HT <- left_join(RNAsamples.HT.H,RNAsamples.HT.T)

#join chemoR family data
RNAsamples.HT <- left_join(RNAsamples.HT,ChemoR_list)
RNAsamples.HT$Family <- factor(RNAsamples.HT$Family, levels = c("srw","srsx","srbc","srab","srt","srx","srxa","srh"))

RNAsamples.HT <- RNAsamples.HT %>%
  mutate(FC = (H+0.1)/(T+0.1)) %>%
  group_by(Gene_ID) %>%
  distinct(Gene_ID,.keep_all=TRUE) %>%
  mutate(Max = max(H,T)) #%>%
  #filter((H+T) > 1) #total TPM less than one

HTplot <- ggplot(RNAsamples.HT)+
  aes(x = log2(FC), y =  '♂')+ 
  geom_jitter(aes(size=Max,colour = Superfamily), alpha = 0.75, width=0.1) +
  theme_minimal() +
  #geom_hline(aes(yintercept=5.6), linetype="solid", colour = "black", size = 0.5) +
  #geom_hline(aes(yintercept=-3.5), linetype="solid", colour = "black", size = 0.5) +
  geom_vline(aes(xintercept=0), linetype="dashed", colour = "grey37", size = 0.5) +
  scale_x_continuous(breaks=c(-3, -2, -1, 0, 1, 2, 3)) +
  theme(axis.title.x= element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=30),
        axis.title.y = element_text(face="bold", size=0),
        axis.ticks.x= element_blank()) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank()) +
  scale_size_continuous(name="Expresion (TPM)") +
  scale_colour_discrete(name="Superfamily") + # set matching colors later
  geom_text(mapping=aes(y='', x=0.8, label=c("head")), colour = "black", fontface = "bold", size=5) +
  geom_text(mapping=aes(y='', x=1.3, label=c("→")), colour = "black", fontface = "bold", size=6) +
  geom_text(mapping=aes(y='', x=0, label=c("")), colour = "black", fontface = "bold", size=9) +
  geom_text(mapping=aes(y='', x=-0.8, label=c("tail")), colour = "black", fontface = "bold", size=5) +
  geom_text(mapping=aes(y='', x=-1.3, label=c("←")), colour = "black", fontface = "bold", size=6) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  theme(legend.position="none")
HTplot



##############################
###### Supplemental figure ##########
##############################


HTseq_supp <- plot_grid(MFplot + theme(plot.margin = unit(c(0.8,0.8, 0.8, 0.8), "cm")),
                   HTplot + theme(plot.margin = unit(c(0.8,0.8, 0.8, 0.8), "cm")), 
                   nrow=2, rel_heights = c(2.5,1), labels = c('A','B'))
HTseq_supp

ggsave(here("plots", "S2_Figure.tiff"), HTseq_supp, width = 12, height =7, units = "in")




##############################
###### MAIN HEATMAP / CLUSTERING FIGURE ##########
##############################

#Filter RNA-Seq data for specific gene set (optional)
RNAsamples.cluster <- RNAsamples.HT.UGA %>%
  filter(Location != "WB") %>%
  filter(Transcript_ID %in% gene_list) %>%
  group_by(Transcript_ID, Stage, Location) %>%  
  distinct(.keep_all=TRUE) %>% ungroup() 
  
#Select relevant columns, filter by min expression, spread data
RNAsamples.cluster <- RNAsamples.cluster %>%
  select(Gene_ID, SID, Expression, Transcript_ID) %>%
  group_by(Gene_ID,SID) %>%
  distinct(Gene_ID,SID,Expression) %>%
  spread(SID,Expression) %>%
  ungroup() %>%
  dplyr::filter(`Bm-M-Head` > 1 | `Bm-M-Tail` > 1 | `Bm-F-Head` > 1) %>%
  rename(Bm_M_Head = "Bm-M-Head", Bm_M_Tail ="Bm-M-Tail",Bm_F_Head = "Bm-F-Head") 

#Set rownames to Gene_ID (will soon be deprecated)
rownames(RNAsamples.cluster) <- RNAsamples.cluster$Gene_ID
RNAsamples.cluster <- select(RNAsamples.cluster, -Gene_ID) 

#Convert to matrix 
y <- data.matrix(RNAsamples.cluster, rownames.force = NA)
ind <- apply(y, 1, var) == 0  #not really necessary after filtering for >10 reads
y <- y[!ind,]
#Optional: add pseudocount and log2-transform
#y <- log2(y+1)
#optional: scale normalize
y <- scale(t(y))
y<- t(y)

#cluster 
hc <- hclust( dist(y, method = "euclidean"), method = "ward.D" )
# plot(hc)

# plotting customized
op = par(bg = "#ffffff")
# plot(hc, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", 
#      col.axis = "#F38630", lwd = 2.5, lty = 1, sub = "", hang = -1, axes = FALSE, cex = 0.5, ann=FALSE)
# add axis
axis(side = 2, at = seq(0, 400, 100), col = "#F38630", labels = FALSE, 
     lwd = 2)
# add text in margin
mtext(seq(0, 400, 100), side = 2, at = seq(0, 400, 100), line = 1, 
      col = "#A38630", las = 2)

# plot dendogram
#dend <- ggdendrogram(hc, rotate = FALSE, size = 2)
dend <- ggdendrogram(hc, rotate = FALSE, size = 2) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) 
# dend

# saveRDS(dend, "~/Box/ZamanianLab/Manuscripts/2019-ChemoR/v1/Figs/Fig2/dend.plot")

# plot heatmap
ord <- hc$order
pd <- as.data.frame(y)  
pd$id <- seq_along(pd$'Bm_F_Head') #change input here
pd$id2 <- rownames(pd)
pd <- left_join(pd, select(RNAsamples.MF, Gene_ID, Superfamily, Family), by = c("id2" = "Gene_ID"))
pd <- unique(pd)
pd <- melt(pd, id.vars = c("id", "id2", "Superfamily", "Family"))
pd$id <- factor(pd$id, levels = c(ord))
pd$id2 <- factor(pd$id2, levels = c(pd$id2[ord]))

pd <- mutate(pd, row = "1")

family.map <- ggplot(pd, aes(x = id2, y = row)) +
  geom_tile(aes(fill = Superfamily)) +
  scale_fill_manual(limits = c("Sra", "Srg", "Str", "Solo"), values = c("#FB6A4A", "#6BAED6", "#74C476", "#FD8D3C")) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 9)) +
  NULL
family.map

heatmap <- ggplot(pd, aes(id2, variable)) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis() +
  scale_y_discrete(limits = c("Bm_F_Head", "Bm_M_Tail", "Bm_M_Head"), labels = c("Female Head", "Male Tail", "Male Head")) +
  labs(y = "Adult Tissue", x = "", fill = "Relative Expression") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(panel.border = element_blank(),
        # panel.grid = element_blank(),
        axis.title.x  = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(),
        axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 9),
        legend.text = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 9)) +
  NULL
heatmap

heatmap_plot <- grid.arrange(dend + theme(plot.margin=unit(c(1,2.3,-2,0.8), "cm")), heatmap + 
                             theme(axis.text.y = element_blank()) + theme(axis.ticks.y = element_blank()) + 
                             theme(plot.margin=unit(c(1,1,1,1), "cm")), ncol=1, heights=c(0.3,1))


final.heatmap <- plot_grid(dend + theme(plot.margin = margin(b = -35)), 
                           family.map,
                           heatmap + theme(axis.title.x = element_blank(),
                                           axis.text.x = element_text(size = 6, face = "plain", hjust = 1, vjust = 1.25),
                                           plot.margin = margin(t = -20)),
                           nrow = 3, align = "v", axis = "rl", 
                           scale = c(1.05, 1, 1), rel_heights = c(0.5, 0.25, 2))
# final.heatmap

save_plot(here("plots", "Fig2B_raw.pdf"), final.heatmap, base_width = 8, base_height = 6)

