library(tidyverse)
library(conflicted)
library(magrittr)
library(here)

here()

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

# read in karyotypes, which were generated from GFF3 files
bm.karyotype <- read_delim(here("data", "karyotype", "karyotype.bmalayi.txt"), delim = " ", col_names = FALSE)
ce.karyotype <- read_delim(here("data", "karyotype", "karyotype.celegans.txt"), delim = " ", col_names = FALSE)
karyotype <- bind_rows(bm.karyotype, ce.karyotype) %>%
  select(Circos.Chrom = X3, Length = X6) %>%
  mutate(Chromosome = c("Bm_v4_Chr1_scaffold_001", "Bm_v4_Chr2_contig_001", "Bm_v4_Chr3_scaffold_001", "Bm_v4_Chr4_scaffold_001", "Bm_v4_ChrX_scaffold_001", "Bm_007",
                        "I", "II", "III", "IV", "V", "X"))

# trim bm.karyotype for later
bm.karyotype <- filter(karyotype, str_detect(Circos.Chrom, "bm"))

# read in GTFs
bm.gtf <- read_delim("~/Box/ZamanianLab/Data/Genomics/ChemoR/gtf/brugia_malayi_prot_table.txt", delim = ",", col_names = FALSE)
ce.gtf <- read_delim("~/Box/ZamanianLab/Data/Genomics/ChemoR/gtf/caenorhabditis_elegans_prot_table.txt", delim = ",", col_names = FALSE)

gtf <- bind_rows(bm.gtf, ce.gtf) %>%
  select(-X4, Chromosome = X1, Start = X2, Stop = X3, Gene_ID = X5, Transcript_ID = X6) %>%
  mutate(Transcript_ID = str_remove(Transcript_ID, "\\.[0-9]$"
  ))

bm.chemoR <- read_csv("~/Box/ZamanianLab/Data/Genomics/ChemoR/list/brugia_malayi_ids.txt", col_names = c("Transcript_ID"))
ce.chemoR <- read_csv("~/Box/ZamanianLab/Data/Genomics/ChemoR/list/caenorhabditis_elegans_ids.txt", col_names = c("Transcript_ID"))
all.chemoR <- bind_rows(bm.chemoR, ce.chemoR) %>%
  left_join(., read_delim("../R/ChemoR/data/family_assignment.csv", delim = ",", col_names = TRUE)) %>%
  filter(Species != 0)

# Highlights --------------------------------------------------------------
# read in chemoR data and get start/stop coordinates
all.chemoR <- left_join(all.chemoR, gtf) %>%
  left_join(., karyotype) %>%
  mutate(Length = Stop - Start) %>%
  arrange(Gene_ID, desc(Length)) %>%
  distinct(Gene_ID, .keep_all = TRUE) # remove isoforms

highlights <- select(all.chemoR, Circos.Chrom, Start, Stop) %>%
  arrange(Circos.Chrom, Start)

# get a list of chromosomes and calculate chunks so that each gene will have equidistant area and the heatmap will span the entire chromosome/ideogram
chromosomes <- all.chemoR %>%
  group_by(Chromosome) %>%
  mutate(n = n()) %>%
  distinct(Chromosome, .keep_all = TRUE) %>%
  select(Chromosome, n) %>%
  left_join(., karyotype) %>%
  filter(Length != 0) %>%
  mutate(Chunk = ifelse(Length / n > 3000000, 2000000, Length / n)) %>%
  select(Chromosome, Chunk, Length)

sf.colors <- tibble(Superfamily = c("Sra", "Srg", "Str", "Solo"),
                    Color = c("red", "blue", "green", "orange"))  

expand <- all.chemoR %>%
  select(Chromosome, Circos.Chrom, Start, Stop, Superfamily) %>%
  arrange(Chromosome, Start) %>%
  left_join(., chromosomes) %>%
  mutate(Chunk = floor(Chunk)) %>%
  group_by(Chromosome) %>%
  mutate(Chunk.Stop = accumulate(.x = Chunk, `+`)) %>%
  mutate(Chunk.Start = Chunk.Stop - Chunk + 1) %>%
  ungroup()
  
highlights.expand <- expand %>%
  select(Circos.Chrom, Chunk.Start, Chunk.Stop, Superfamily) %>%
  left_join(., sf.colors) %>%
  select(-Superfamily) %>%
  mutate(Fill.Color = str_glue("fill_color={Color},stroke_color={Color}")) %>%
  select(-Color)

# Heatmaps ----------------------------------------------------------------

# read in RNAseq data
rnaseq <- read_delim(here("data", "RNAseq.data.csv"), delim = ",", col_names = TRUE)

# add GTF metadata to RNAseq daata
chemoR.RNAseq <- filter(all.chemoR, Species == "brugia_malayi") %>% 
  left_join(bm.karyotype) %>%
  left_join(., rnaseq)

# create start/stop values based on chunk size for each chromosome
circos.chemoR <- filter(all.chemoR, Species == "brugia_malayi") %>%
  select(Chromosome, Start, Stop) %>%
  arrange(Chromosome, Start) %>%
  left_join(., chromosomes) %>%
  mutate(Chunk = floor(Chunk)) %>%
  group_by(Chromosome) %>%
  mutate(Heatmap.Stop = accumulate(.x = Chunk, `+`)) %>%
  mutate(Heatmap.Start = Heatmap.Stop - Chunk + 1)

# merge with RNAseq data and get rid of extraneous columns
heatmaps <- select(chemoR.RNAseq, Circos.Chrom, Chromosome, Start, Stop, Mean_Expression, Sample_Name) %>%
  group_by(Circos.Chrom) %>%
  left_join(., circos.chemoR) %>%
  arrange(Start) %>%
  select(Circos.Chrom, Heatmap.Start, Heatmap.Stop, Mean_Expression, Sample_Name) %>%
  mutate(Mean_Expression = log2(Mean_Expression + 1))

legend <- plot(cowplot::get_legend(ggplot(heatmaps) + geom_tile(aes(x = Sample_Name, y = Heatmap.Start, height = Heatmap.Stop - Heatmap.Start, fill = Mean_Expression)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = "Purples")) +
  labs(fill = "Log2(Transcripts Per Million)") +
  theme(legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10))
  )
  )

# Connectors --------------------------------------------------------------
# make connectors from the ideogram highlights to the first heatmap
connectors <- ungroup(expand) %>%
  mutate(Highlight.Connector = round((Start + Stop) / 2)) %>%
  mutate(Chunk.Connector = round((Chunk.Start + Chunk.Stop) / 2)) %>%
  left_join(., karyotype) %>%
  ungroup() %>%
  select(Circos.Chrom, Chunk.Connector, Highlight.Connector) %>%
  arrange(Circos.Chrom, Chunk.Connector)
  
stages <- unique(chemoR.RNAseq$Sample_Name)

# Writing ------------------------------------------------------------

# write out so that each stage has a separate file
write_data <- function(stage) {
  
  write_delim(highlights, here("data", "highlights", "all_highlights.txt"), delim = " ", col_names = FALSE)
  write_delim(highlights.expand, here("data", "highlights", "highlights_superfamily.txt"), delim = " ", col_names = FALSE)
  temp <- filter(heatmaps, Sample_Name == stage) %>%
    select(-Sample_Name)
  write_delim(temp, here("data", "heatmap", stage), delim = " ", col_names = FALSE)
  write_delim(connectors, here("data", "connector", "connectors.txt"), delim = " ", col_names = FALSE)
  
}

temp <- apply(tibble(stages), 1, write_data)
