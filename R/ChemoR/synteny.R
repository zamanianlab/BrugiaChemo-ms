library(tidyverse)
library(conflicted)
library(here)
library(cowplot)

here()

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# read in karyotypes, which were generated from GFF3 files
bm.karyotype <- read_delim(here("data", "karyotype.bmalayi.txt"), delim = " ", col_names = c("Number", "Name", "Length", "Species"))
ce.karyotype <- read_delim(here("data", "karyotype.celegans.txt"), delim = " ", col_names = c("Number", "Name", "Length", "Species"))
karyotype <- bind_rows(bm.karyotype, ce.karyotype) %>%
  mutate(Number = factor(Number, levels = c("1", "2", "3", "4", "5", "X", "007", "023"))) %>%
  mutate(Species = factor(Species)) %>%
  mutate(Species.Order = case_when(
    Species == "caenorhabditis_elegans" ~ 2,
    Species == "brugia_malayi" ~ 1,
  ))

# read in GTFs
bm.gtf <- read_csv(here("data", "brugia_malayi_prot_table.txt"), col_names = FALSE)
ce.gtf <- read_csv(here("data", "caenorhabditis_elegans_prot_table.txt"), col_names = FALSE)

chemor.gtf <- bind_rows(bm.gtf, ce.gtf) %>%
  select(-X4, Chromosome = X1, Start = X2, Stop = X3, Gene_ID = X5, Transcript_ID = X6) %>%
  mutate(Transcript_ID = str_remove(Transcript_ID, "\\.[0-9]$"
  ))

bm.chemor <- read_csv(here("data", "brugia_malayi_ids.txt"), col_names = c("Transcript_ID"))
ce.chemor <- read_csv(here("data", "caenorhabditis_elegans_ids.txt"), col_names = c("Transcript_ID"))
all.chemor <- bind_rows(bm.chemor, ce.chemor) %>%
  left_join(., read_csv(here("data", "family_assignment.csv"), col_names = TRUE)) %>%
  filter(Species != 0) %>%
  mutate(Species = factor(Species)) %>%
  left_join(chemor.gtf)

# Highlights --------------------------------------------------------------
# read in chemoR data and get start/stop coordinates
all.chemor <- mutate(all.chemor, Length = Stop - Start) %>%
  arrange(Gene_ID, desc(Length)) %>%
  distinct(Gene_ID, .keep_all = TRUE) %>% # remove isoforms
  mutate(Number = case_when(
    Chromosome == "I" ~ "1",
    Chromosome == "II" ~ "2",
    Chromosome == "III" ~ "3",
    Chromosome == "IV" ~ "4",
    Chromosome == "V" ~ "5",
    Chromosome == "X" ~ "X",
    Chromosome == "Bm_v4_Chr1_scaffold_001" ~ "1",
    Chromosome == "Bm_v4_Chr2_contig_001" ~ "2",
    Chromosome == "Bm_v4_Chr3_scaffold_001" ~ "3",
    Chromosome == "Bm_v4_Chr4_scaffold_001" ~ "4",
    Chromosome == "Bm_v4_ChrX_scaffold_001" ~ "X",
    Chromosome == "Bm_007" ~ "007",
    Chromosome == "Bm_023" ~ "023"
  )) %>%
  mutate(Species.Order = case_when(
    Species == "caenorhabditis_elegans" ~ 2,
    Species == "brugia_malayi" ~ 1,
  ))

# first row --> chromosomes with highlights

chromosome_plot <- ggplot() +
  geom_rect(data = filter(karyotype, !Number %in% c("007", "023")),
            aes(xmin = 0, 
                xmax = Length / 1000000,
                ymin = Species.Order + 0.1,
                ymax = Species.Order + 0.35),
            fill = "white",
            color = "black") +
  geom_rect(data = filter(all.chemor, !Number %in% c("007", "023")),
            aes(xmin = (Start - 50000) / 1000000,
                xmax = (Stop + 50000) / 1000000,
                ymin = Species.Order + 0.1,
                ymax = Species.Order + 0.35),
            fill = "black") +
  facet_grid(cols = vars(Number), scales = "free") +
  NULL
# chromosome_plot

# second row --> superfamily colors

# get a list of chromosomes and calculate chunks so that each gene will have equidistant area and the heatmap will span the entire chromosome/ideogram
chromosomes <- all.chemor %>%
  filter(!is.na(Superfamily)) %>%
  group_by(Chromosome) %>%
  mutate(n = n()) %>%
  distinct(Chromosome, .keep_all = TRUE) %>%
  select(Species, Chromosome, Number, n) %>%
  left_join(., karyotype, by = c("Species", "Number")) %>%
  rename(Chromosome.Length = Length) %>%
  filter(Chromosome.Length != 0) %>%
  mutate(Chunk = ifelse(Chromosome.Length / n > 3000000, 2000000, Chromosome.Length / n)) %>%
  select(Chromosome, Chunk, Chromosome.Length)

all.chemor <- all.chemor %>%
  filter(!is.na(Superfamily)) %>%
  arrange(Chromosome, Start) %>%
  left_join(., chromosomes, by = "Chromosome") %>%
  mutate(Chunk = floor(Chunk)) %>%
  group_by(Chromosome) %>%
  mutate(Chunk.Stop = accumulate(.x = Chunk, `+`)) %>%
  mutate(Chunk.Start = Chunk.Stop - Chunk + 1) %>%
  ungroup()

# sf.colors <- tibble(Superfamily = c("Sra", "Srg", "Str", "Solo"),
#                     Color = c("red", "blue", "green", "orange"))  

superfamily_plot <- chromosome_plot +
  geom_rect(data = filter(all.chemor, Species == "brugia_malayi", !Number %in% c("007", "023")),
            aes(xmin = Chunk.Start / 1000000, 
                xmax = Chunk.Stop / 1000000,
                ymin = Species.Order - 0.025, 
                ymax = Species.Order - 0.225, 
                fill = Superfamily), color = "black", size = 0.5) +
  geom_segment(data = filter(all.chemor, Species == "brugia_malayi", !Number %in% c("007", "023")),
               aes(x = (Start + Stop) / 2000000,
                   xend = (Chunk.Start + Chunk.Stop) / 2000000,
                   y = Species.Order + 0.1,
                   yend = Species.Order - 0.025),
               size = 0.35) +
  geom_rect(data = filter(all.chemor, Species == "caenorhabditis_elegans", !is.na(Superfamily)),
            aes(xmin = Chunk.Start / 1000000, 
                xmax = Chunk.Stop / 1000000, 
                ymin = Species.Order - 0.025, 
                ymax = Species.Order - 0.225, 
                fill = Superfamily,
                color = Superfamily)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  # scale_x_continuous(labels = scales::scientific, breaks = c(5000000, 10000000, 15000000, 20000000, 25000000)) +
  NULL
# superfamily_plot

# following rows --> staged RNA-seq
rnaseq <- read_csv(here("data", "RNAseq.data.csv"), col_names = TRUE)

# add GTF metadata to RNAseq data
chemor.rnaseq <- filter(all.chemor, Species == "brugia_malayi") %>% 
  left_join(., rnaseq)

# merge with RNAseq data and get rid of extraneous columns
heatmaps <- chemor.rnaseq %>%
  select(Chromosome, Number, Chunk.Start, Chunk.Stop, Mean_Expression, Sample_Name) %>%
  mutate(Mean_Expression = log2(Mean_Expression + 1)) %>%
  mutate(Order = case_when(
    Sample_Name == "E" ~ 0.75,
    Sample_Name == "MF_60DPI" ~ 0.65,
    Sample_Name == "MFL1_V18HPI" ~ 0.55,
    Sample_Name == "L2_V4DPI" ~ 0.45,
    Sample_Name == "L3_V8DPI" ~ 0.35,
    Sample_Name == "L3_1DPI" ~ 0.25,
    Sample_Name == "L4M_20DPIM" ~ 0.15,
    Sample_Name == "L4F_20DPIF" ~ 0.05,
    Sample_Name == "AM" ~ -0.05,
    Sample_Name == "AF" ~ -0.15,
  )) %>%
  filter(!is.na(Mean_Expression))

text_labels <- distinct(heatmaps, Sample_Name, Order) %>%
  mutate(Number = 5, x = 10000000, Order = Order - 0.05) %>%
  mutate(Label = case_when(
    Sample_Name == "E" ~ "Embryo",
    Sample_Name == "MF_60DPI" ~ "Microfilaria, 60 DPI",
    Sample_Name == "MFL1_V18HPI" ~ "mf/L1, 18 HPI",
    Sample_Name == "L2_V4DPI" ~ "L2, 4 DPI",
    Sample_Name == "L3_V8DPI" ~ "L3, 8 DPI",
    Sample_Name == "L3_1DPI" ~ "L3, 1 DPI",
    Sample_Name == "L4M_20DPIM" ~ "L4 Male, 20 DPI",
    Sample_Name == "L4F_20DPIF" ~ "L4 Female, 20 DPI",
    Sample_Name == "AM" ~ "Adult Male",
    Sample_Name == "AF" ~ "Adult Female",
  )) 

rnaseq_plot <- 
  superfamily_plot +
  ggnewscale::new_scale_fill() +
  geom_rect(data = filter(heatmaps, !Number %in% c("007", "023")), 
            aes(xmin = Chunk.Start / 1000000,
                xmax = Chunk.Stop / 1000000,
                ymin = Order,
                ymax = Order - 0.1, 
                fill = Mean_Expression)) +
  geom_text(data = text_labels, 
            aes(label = Label, x = x / 1000000, y = Order),
            fontface = "bold",
            size = 4) +
  scale_fill_gradient(low = "white", high = "black", name = "Log2(TPM)") 
# rnaseq_plot

# following rows --> head/tail RNA-seq
rnaseq <- read_csv(here("data", "RNAseq.data.csv"), col_names = TRUE)

# add GTF metadata to RNAseq data
chemor.rnaseq <- filter(all.chemor, Species == "brugia_malayi") %>% 
  left_join(., rnaseq)

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
  group_by(Transcript_ID, Stage) %>%
  summarize(Whole_Expression = mean(Expression)) # whole-body expression

RNAsamples.HT.UGA <- inner_join(RNAsamples.HT, RNAsamples.UGA)

# Made new dataframe of just whole body and re-merge (so that UGA data gets own line)
temp <- RNAsamples.HT.UGA %>%
  mutate(Expression = Whole_Expression) %>%
  select(-Whole_Expression) %>%
  mutate(Location = "WB") %>%
  mutate(SID = ifelse(grepl("Bm-F-Head", SID), "Bm-F-WB", SID)) %>%
  mutate(SID = ifelse(grepl("Bm-M-Head", SID), "Bm-M-WB", SID)) %>%
  mutate(SID = ifelse(grepl("Bm-M-Tail", SID), "Bm-M-WB", SID)) %>%
  group_by(Transcript_ID, SID) %>%
  distinct(.keep_all = True) %>%
  ungroup()
temp2 <- RNAsamples.HT.UGA %>% select(-Whole_Expression)
RNAsamples.HT.UGA <- rbind(temp, temp2) %>%
  select(Transcript_ID, Gene_ID, SID, Metric, Expression, Sample_Name = Stage, Location) %>%
  filter(Metric == "TPM_g")

ht.rnaseq <- filter(RNAsamples.HT.UGA, Gene_ID %in% all.chemor$Gene_ID, Location != "WB") %>%
  left_join(., distinct(select(chemor.rnaseq, Number, Gene_ID, Chunk.Start, Chunk.Stop))) %>%
  mutate(Mean_Expression = log2(Expression + 1)) %>%
  mutate(Order = case_when(
    SID == "Bm-F-Head" ~ -0.3,
    SID == "Bm-M-Head" ~ -0.4,
    SID == "Bm-M-Tail" ~ -0.5
  )) %>%
  filter(!is.na(Mean_Expression))

ht_labels <- distinct(ht.rnaseq, SID, Order) %>%
  mutate(Number = 5, x = 10000000, Order = Order - 0.05) %>%
  mutate(Label = case_when(
    SID == "Bm-F-Head" ~ "Female Head",
    SID == "Bm-M-Head" ~ "Male Head",
    SID == "Bm-M-Tail" ~ "Male Tale"
  )) 

ht_plot <- 
  rnaseq_plot +
  ggnewscale::new_scale_fill() +
  geom_rect(data = filter(ht.rnaseq, !Number %in% c("007", "023")),
            aes(xmin = Chunk.Start / 1000000,
                xmax = Chunk.Stop / 1000000,
                ymin = Order,
                ymax = Order - 0.1, 
                fill = Mean_Expression)) +
  geom_text(data = ht_labels, 
            aes(label = Label, x = x / 1000000, y = Order),
            fontface = "bold",
            size = 4) +
  scale_fill_viridis_c(option = "plasma", name = "Log2(TPM)") 
# ht_plot

# add rectangles and labels
rectangles <- select(all.chemor, Number, Species, Length, Chunk.Stop) %>%
  filter(Species == "brugia_malayi", !Number %in% c("007", "023")) %>%
  group_by(Number) %>%
  summarize(xmax = max(Chunk.Stop)) %>%
  mutate(xmin = 0)
  
# host labels
host_labels <- tibble(Label = c("Human Stage", "Vector Stage", "Human Stage"),
                      y = c(0.65, 0.4, 0.095),
                      Number = c(3, 3, 3),
                      x = c(9800000, 9800000, 9800000))

species_labels <- tibble(Label = c("C. elegans", "B. malayi"),
                         y = c(2.05, 1.05),
                         Number = c(1, 1),
                         x = c(-1500000, -1500000))

final_plot <- ht_plot +
  geom_hline(yintercept = 0.5455, color = "grey", size = 0.5, linetype = 2) +
  geom_hline(yintercept = 0.2545, color = "grey", size = 0.5, linetype = 2) +
  geom_rect(data = rectangles,
            aes(xmin = xmin / 1000000, xmax = xmax / 1000000),
            ymin = -0.25, ymax = 0.75, color = "black", fill = NA, size = 0.5) +
  geom_rect(data = rectangles,
            aes(xmin = xmin / 1000000, xmax = xmax / 1000000),
            ymin = -0.3, ymax = -0.6, color = "black", fill = NA, size = 0.5) +
  # geom_rect(data = rectangles, 
  #           aes(xmin = xmin / 1000000, xmax = xmax / 1000000),
  #           ymin = -0.25, ymax = 0.34, color = "steelblue", size = 1, fill = NA, linetype = 2) +
  geom_text(data = host_labels, 
            aes(x = x / 1000000, y = y, label = Label),
            fontface = "bold",
            size = 3.9) +
  geom_text(data = species_labels,
            aes(x = x / 1000000, y = y, label = Label),
            fontface = "bold.italic",
            size = 3.9,
            angle = 90) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  labs(x = "Chromosome (Mb)", y = "") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, vjust = 113),
    axis.text.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 12, vjust = 105),
    strip.text.x = element_text(face = "bold", size = 20),
    panel.grid = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)) +
  NULL
final_plot

save_plot(here("plots", "Fig2.pdf"), final_plot, base_height = 6, base_width = 12)




