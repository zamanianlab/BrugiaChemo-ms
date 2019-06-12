library(tidyverse)
library(cowplot)
library(ggrepel)
library(chroma)
library(here)
library(conflicted)

conflict_prefer("scale_color_viridis_d", "ggplot2")
conflicted::conflict_prefer("filter", "dplyr")

# Tidy data -----------------------------------------------------------
# get clades as defined by the ML tree
tree <- read_delim(here("data", "tree_clades.csv"), delim = ",") %>%
  distinct(.keep_all = TRUE)

# load in reference file matching species <-> clade
species_info <- read_delim(here("../../auxillary", "species_info.csv"),
                           col_names = TRUE, 
                           delim = ",") %>%
  select(-BioProject)

# get blast results, which blasted predicted ChemoRs from species not included in the tree against
# the ChemoRs included in the tree
blast <- read_delim(here("data", "ChemoR.blastout"),
                    col_names = c("Species", "Transcript_ID", "Best_Hit", "qseqid", "sseqid", "pident", "ppos", "length", "mismatch", "evalue", "bitscore"),
                    delim = "\t") %>%
  distinct(Transcript_ID, .keep_all = TRUE) %>%
  select(Species, Transcript_ID, Best_Hit) %>%
  filter(!Transcript_ID %in% tree$Transcript_ID) %>%
  left_join(., tree, by = c("Best_Hit" = "Transcript_ID")) %>%
  select(-Best_Hit, -Species.y) %>%
  rename(Species = Species.x) %>%
  select(Family, Species, Transcript_ID)

# add blasted ChemoRs to ChemoRs from tree
df <- bind_rows(tree, blast) %>%
  left_join(., species_info, by = "Species") %>%
  distinct(Transcript_ID, .keep_all = TRUE)

# remove isoforms 
dedup <- mutate(df, Temp = case_when(
  Species == "acani" ~ Transcript_ID,
  Species == "acant" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "acey" ~ str_remove(Transcript_ID, "\\.t[0-9]$"),
  Species == "aduo" ~ Transcript_ID,
  Species == "alum" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "asuu" ~ str_remove(Transcript_ID, "\\_t[0-9]*$"),
  Species == "bmal" ~ str_remove(Transcript_ID, "[a-z]$"),
  Species == "bpah" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "btim" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "cbri" ~ str_remove(Transcript_ID, "[a-z]$"),
  Species == "cele" ~ str_remove(Transcript_ID, "[a-z]$"),
  Species == "dimm" ~ Transcript_ID,
  Species == "dmed" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "dviv" ~ Transcript_ID,
  Species == "hcon" ~ str_remove(Transcript_ID, "\\.[0-9*]$"),
  Species == "hpol" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "lloa" ~ Transcript_ID,
  Species == "lsig" ~ Transcript_ID,
  Species == "mhap" ~ Transcript_ID,
  Species == "name" ~ Transcript_ID,
  Species == "nbra" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "ooch" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "ovol" ~ Transcript_ID,
  Species == "ppac" ~ Transcript_ID,
  Species == "pred" ~ str_remove(Transcript_ID, "\\.t[0-9]$"),
  Species == "puni" ~ str_remove(Transcript_ID, "\\_t[0-9]*$"),
  Species == "rcul" ~ Transcript_ID,
  Species == "rkr3" ~ str_remove(Transcript_ID, "\\.[0-9]*$"),
  Species == "scar" ~ str_remove(Transcript_ID, "\\.t[0-9]$"),
  Species == "smur" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "srat" ~ str_remove(Transcript_ID, "[a-z]$"),
  Species == "sste" ~ str_remove(Transcript_ID, "\\.[0-9]*$"),
  Species == "sven" ~ str_remove(Transcript_ID, "\\.[0-9]*$"),
  Species == "tcal" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$"),
  Species == "tcan" ~ str_remove(Transcript_ID, "\\.[0-9]*$"),
  Species == "tmur" ~ Transcript_ID,
  Species == "tspi" ~ Transcript_ID,
  Species == "tsui" ~ Transcript_ID,
  Species == "wban" ~ str_remove(Transcript_ID, "-mRNA-[0-9]$")
)) %>%
  arrange(Species, Transcript_ID) %>%
  distinct(Temp, .keep_all = TRUE) %>%
  select(-Temp)

# count number of ChemoRs in each family for each species
# remove peptide families
dedup.p <- dedup %>%
  group_by(Family, Species, Clade) %>%
  summarise(n()) %>%
  rename(Count = "n()") %>%
  ungroup() %>%
  filter(!grepl("pep", Family)) %>%
  group_by(Species) %>%
  mutate(Count_Norm = Count / sum(Count)) %>%
  ungroup %>%
  left_join(., species_info)

# create a data.frame with superfamily info
chemo_families <- data.frame(Family = c("srh", "str", "sri", "srd", "srj", "srm", "srn", "four", "three_four",
                                        "sre", "sra", "srab", "srb",
                                        "srx", "srt", "srg", "sru", "srv", "srxa",
                                        "srw", "srz", "srbc", "srsx", "srr", "sro"),
                             Superfamily = c("Str", "Str", "Str", "Str", "Str", "Str", "Str", "Str", "Str",
                                             "Sra", "Sra", "Sra", "Sra",
                                             "Srg", "Srg", "Srg", "Srg", "Srg", "Srg",
                                             "Solo", "Solo", "Solo", "Solo", "Solo", "Solo"))

# add superfamily info
dedup.p <- left_join(dedup.p, chemo_families, by = "Family") %>%
  mutate(Family = factor(Family, levels = c("srh", "str", "sri", "srd", "srj", "srm", "srn", "four", "three_four",
                                            "sre", "sra", "srab", "srb",
                                            "srx", "srt", "srg", "sru", "srv", "srxa",
                                            "srw", "srz", "srbc", "srsx", "srr", "sro"),
                         labels = c("srh", "str", "sri", "srd", "srj", "srm", "srn", "IV-Specific", "III/IV-Specific",
                                    "sre", "sra", "srab", "srb",
                                    "srx", "srt", "srg", "sru", "srv", "srxa",
                                    "srw", "srz", "srbc", "srsx", "srr", "sro")))

# change genus_species to G. species
species <- dedup.p$Full %>%
  str_replace("(^.)(.*_)", "\\1. ") %>%
  str_replace("(^.)", toupper)
dedup.p$Full <- species

# Numbers/Stats -----------------------------------------------------------

life.history <- tibble(Type = c("Free-Living", "Free-Living", "Free-Living", "Free-Living", "Free-Living",
                                "Facultative", "Facultative", "Facultative",
                                "Skin-Penetrating", "Skin-Penetrating", "Skin-Penetrating", "Skin-Penetrating", "Skin-Penetrating", "Skin-Penetrating", "Skin-Penetrating", "Skin-Penetrating",
                                "Ingested Larvae", "Ingested Larvae", "Ingested Larvae",
                                "Ingested Egg", "Ingested Egg", "Ingested Egg", "Ingested Egg", "Ingested Egg", "Ingested Egg", "Ingested Egg",
                                "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted", "Vector Transmitted",
                                "Host-Contained",
                                "Plant-Parasitic"),
                       Species = c("cbri", "cele", "ppac", "rkr3", "pred",
                                   "srat", "sste", "sven",
                                   "rcul", "name", "acani", "nbra", "dmed", "scar", "acey", "aduo",
                                   "dviv", "hcon", "hpol",
                                   "tmur", "tcan", "asuu", "smur", "tsui", "alum", "puni",
                                   "bmal", "bpah", "btim", "dimm", "lloa", "lsig", "ooch", "ovol", "tcal", "wban", "acant",
                                   "tspi",
                                   "mhap"))
type <- left_join(dedup.p, life.history) %>%
  filter(Species != "rculi") %>% # doesn't fit well in any of the groups
  group_by(Species, Type, Full) %>%
  summarise(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(Type = factor(Type, levels = c("Free-Living", "Plant-Parasitic", "Facultative", "Skin-Penetrating", "Ingested Larvae", "Ingested Egg", "Vector Transmitted", "Host-Contained")))

rect <- as.data.frame(unique(type$Type)) %>%
  mutate(fill = seq(1, 7, 1)) %>%
  rename(Type = 1)

gradient <- matrix(viridis_colors(8), nrow = 50, ncol = length(viridis_colors(8)), byrow = TRUE)

pos <- position_jitter(width = 0.25, seed = 1)

# Plotting ----------------------------------------------------------------

# heatmap -----------------------------------------------------------------

heatmap <- ggplot(dedup.p, aes(Full, Family)) +
  geom_tile(aes(fill = Count_Norm), height = 1, width = 1) +
  facet_grid(Superfamily ~ Clade, scales = "free", space = "free") +
  scale_fill_gradient(name = paste0("Fraction Total", "\n", "(per species)"), low = "white", high = "palegreen4") +
  labs(x = "", y = "Chemoreceptor Family") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.25, "line"),
        panel.spacing.y = unit(0.25, "line"),
        axis.line = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(face = "bold", angle = 90, size = 12),
        axis.text.x = element_text(size = 9.5, face = "bold.italic", angle = 60, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 8),
        strip.text = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white", size = 0.5),
        legend.text = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 9)) +
  NULL
# heatmap

# saveRDS(heatmap, file = here("data", "heatmap.rds"))

# save_plot(here("plots", "Fig1B_heatmap.pdf"), heatmap, base_width = 10, base_height = 6)

# bar plot -----------------------------------------------------------------

# sum all ChemoRs by species
df.b <- group_by(dedup.p, Species) %>%
  summarise(Sum = sum(Count)) %>%
  left_join(., species_info) %>%
  mutate(All = "")

bar <- ggplot(df.b, aes(x = Species, y = Sum)) +
  geom_bar(aes(fill = Clade), stat = "identity", width = 1) + 
  scale_fill_manual(values = c(I = "#ABB065", IIIa = "#E495A5", IIIb = "#E495A5", IIIc = "#E495A5", IVa = "#39BEB1", IVb = "#39BEB1", Va = "#ACA4E2", Vb = "#ACA4E2", Vc = "#ACA4E2", Vo = "#ACA4E2")) +
  facet_grid(All ~ Clade, scales = "free", space = "free") +
  labs(y = "Receptor Count") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.25, "line"),
        panel.spacing.y = unit(0.25, "line"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y  = element_text(face = "bold", angle = 90, size = 12),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 12)) +
  NULL
# bar

# saveRDS(bar, file = here("data", "bar.rds"))

# save_plot(here("plots", "Fig1B_bar.pdf"), bar, base_width = 10, base_height = 3)

# type plot ---------------------------------------------------------------

rho_test <- cor.test(as.numeric(type$Type), log2(type$Total), data = type, method = "spearman")
estimate <- rho_test$estimate
p <- rho_test$p.value

type.plot <- ggplot(type, aes(x = Type)) +
  geom_abline(slope = estimate, intercept = log2(1500), linetype = 18, alpha = 0.75, size = 0.5) +
  geom_label_repel(aes(label = Full, y = log2(Total)), fill = "white", position = pos, 
                   size = 3, label.size = 0.2, segment.size = 0.2, fontface = "italic") + # two labels so that only the box fill is transparent
  geom_point(aes(y = log2(Total), color = Type), size = 4, position = pos) +
  # geom_label_repel(aes(label = Full, fill = Type, y = log2(Total)), alpha = 0.4, box.padding = 0.5, seed = 1234, size = 3, fontface = "italic") +
  annotation_raster(gradient, 0.5, 8.5, 1.25, 1.5) +
  annotate("text", x = 7, y = 11, label = expression(paste(rho, " = ", -0.813)), size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = 7, y = 10.3, label = expression(paste("p", " = ", 3.24e-10)), size = 4, fontface = "bold", hjust = 0) +
  annotate("text", 2.25, 2, label = "Decreasing Environmental Exposure", size = 4, fontface = "bold") +
  annotate("segment", x = 4, xend = 4.5, y = 2, yend = 2, colour = "black", size = 0.5, arrow = arrow(angle = 25, type = "closed")) +
  ylim(1.75, 11.25) +
  scale_color_viridis_d() +
  labs(x = "", y = expression(paste("log"[2], "(Total Chemoreceptors)"))) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(panel.spacing.x = unit(0.25, "line"),
        panel.spacing.y = unit(0.25, "line"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(color="black", size = 0.25),
        axis.line.y = element_line(color="black", size = 0.25),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(face = "bold", angle = 90, size = 12),
        axis.text.x = element_text(size = 9.5, face = "bold", angle = 35, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 9),
        strip.text = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  NULL
# type.plot

# saveRDS(type.plot, here("data", "type.plot"))

save_plot(here("plots", "Fig1C_raw.pdf"), type.plot, base_width = 14, base_height = 10)

# Export family assignments provided from the tree and BLAST --------------

df4 <- left_join(df, species_info) %>%
  left_join(., chemo_families) %>%
  select(-Clade, -Species) %>%
  rename(Species = Full) %>%
  filter(!grepl("pep", Family)) %>%
  select(Transcript_ID, Species, Family, Superfamily) %>%
  separate(Transcript_ID, into = c("Transcript_ID", "TEMP"), sep = "-", remove = FALSE, extra = "drop") %>%
  select(-TEMP)

# write.csv(df4, file = here("data", "family_assignment.csv"), row.names = FALSE, col.names = TRUE)
