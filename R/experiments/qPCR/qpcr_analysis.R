library(tidyverse)
library(cowplot)
library(conflicted)
library(magrittr)
library(lubridate)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)
library(Hmisc)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("here", "here")

# Import and tidy data ----------------------------------------------------

results <- read_delim(here("data", "results.txt"), col_names = TRUE, delim = "\t") %>%
  filter(!is.na(`Sample Name`)) %>%
  mutate(Replicate = c(rep(c(rep.int(1, 6), rep(2, 6), rep(3, 6)), 3)))

amplification_data <- read_delim(here("data", "amplification_data.txt"), col_names = TRUE, delim = "\t", col_types = "cicdd") %>%
  filter(!is.na(`Target Name`)) %>%
  mutate(log.dRn = log10(dRn)) %>%
  pivot_longer(cols = Rn:log.dRn, names_to = "Measurement", values_to = "Values")

melt_region_temperature_data <- read_delim(here("data", "melt_region_temperature_data.txt"), col_names = TRUE, delim = "\t") %>%
  pivot_longer(cols = contains("Reading"), names_to = "Reading", values_to = "Temperature") %>%
  mutate(Reading = str_replace(Reading, "Reading ", "")) %>%
  mutate(Reading = as.integer(Reading))

melt_region_derivative_data <- read_delim(here("data", "melt_region_derivative_data.txt"), col_names = TRUE, delim = "\t") %>%
  pivot_longer(cols = contains("Reading"), names_to = "Reading", values_to = "Derivative") %>%
  mutate(Reading = str_replace(Reading, "Reading ", "")) %>%
  mutate(Reading = as.integer(Reading)) 

melt_region_normalized_data <- read_delim(here("data", "melt_region_normalized_data.txt"), col_names = TRUE, delim = "\t") %>%
  pivot_longer(cols = contains("Reading"), names_to = "Reading", values_to = "Value") %>%
  mutate(Reading = str_replace(Reading, "Reading ", "")) %>%
  mutate(Reading = as.integer(Reading)) %>%
  mutate(Well = as.integer(Well)) %>%
  mutate(Value = as.numeric(Value))

melt_region_normalized_data <- left_join(melt_region_normalized_data, melt_region_temperature_data) %>%
  left_join(., melt_region_derivative_data) %>%
  select(-Reading, -`Reporter Dye`, -`Well Location`) %>%
  pivot_longer(cols = c(Value, Derivative), names_to = "Measurement", values_to = "Value")
  
# Plotting and analysis ---------------------------------------------------

melt_curves <- ggplot(melt_region_normalized_data, aes(x = Temperature, y = Value)) +
  geom_line(aes(color = as.factor(Well))) +
  facet_grid(rows = vars(Target), cols = vars(Measurement), scales = "free", space = "free") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none") +
  NULL
melt_curves

ct <- select(results, `Sample Name`, `Target Name`, `Ct`, `Ct Threshold`) %>%
  mutate(Measurement = "dRn") %>%
  distinct()

amplification_plot <- ggplot(filter(amplification_data, Measurement != "Rn"), aes(x = Cycle, y = Values)) +
  geom_line(aes(color = as.factor(Well))) +
  geom_hline(data = ct, aes(yintercept = `Ct Threshold`)) +
  facet_wrap(vars(`Target Name`, Measurement), nrow = 3, scales = "free_y") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "none") +
  NULL
amplification_plot

# Ct ----------------------------------------------------------------------

ct_df <- select(results, Replicate, `Sample Name`, `Target Name`, Ct) %>%
  mutate(Ct = as.numeric(Ct)) %>%
  group_by(Replicate, `Sample Name`, `Target Name`) %>%
  summarise(Tech.Mean.Ct = mean(Ct))

ct_plot <- ggplot(ct_df, aes(x = `Sample Name`, y = Tech.Mean.Ct)) +
  geom_point(aes(color = as.factor(Replicate)), size = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  scale_x_discrete(limits = c("Bp L3 FRESH", "Bp L3 37C 4 HR", "Bp L3 RT 4 HR"), labels = c("Freshly Extracted", paste0("4 Hr. 37", "\u00b0", "C"), "4 Hr. Room Temp")) +
  labs(x = "", y = expression(C[T]), title = " ") +
  facet_grid(cols = vars(`Target Name`)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(angle = 90, vjust = 0.5, size = 0),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none") +
  NULL
ct_plot

# ddCt --------------------------------------------------------------------

summary <- select(results, Replicate, `Sample Name`, `Target Name`, Ct) %>%
  filter(`Target Name` != "TAX-4") %>%
  mutate(Ct = as.numeric(Ct)) %>%
  group_by(Replicate, `Sample Name`, `Target Name`) %>%
  summarize(Mean.Ct = mean(Ct)) %>%
  pivot_wider(names_from = `Target Name`, values_from = Mean.Ct)

ddCt_df <- summary %>%
  mutate(dCt = `OSM-9` - GAPDH) %>%
  select(-GAPDH, -`OSM-9`) %>%
  pivot_wider(names_from = `Sample Name`, values_from = dCt) %>%
  mutate(RTvFresh = `Bp L3 RT 4 HR` - `Bp L3 FRESH`, TSvFresh = `Bp L3 37C 4 HR` - `Bp L3 FRESH`, RTvTS = `Bp L3 RT 4 HR` - `Bp L3 37C 4 HR`) %>%
  select(RTvFresh, TSvFresh, RTvTS) %>%
  pivot_longer(cols = everything(), names_to = "Comparison", values_to = "ddCt")

ddCt_summary <- group_by(ddCt_df, Comparison) %>%
  summarise(Mean = mean(ddCt), SEM = sd(ddCt) / sqrt(length(ddCt)))

ddCt <- ggplot(ddCt_df) +
  # geom_point() +
  geom_col(data = ddCt_summary, aes(x = Comparison, y = Mean), fill = "black") +
  geom_errorbar(data = ddCt_summary, aes(x = Comparison, ymax = Mean + SEM, ymin = 0.1), alpha = 0.5, width = 0.25) +
  scale_x_discrete(limits = c("TSvFresh", "RTvFresh"), labels = c(paste0("4 HPE, 37", "\u00b0", "C \n vs. 0 HPE"), paste0("4 HPE. 21", "\u00b0", "C \n vs. 0 HPE"))) +
  labs(x = "", y = expression(Delta*Delta*C[T]), title = " ") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(angle = 90, vjust = 0.5, size = 0),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none") +
  NULL
ddCt

save_plot(here("plots", "Fig4C_raw.pdf"), ddCt, base_width = 2)
