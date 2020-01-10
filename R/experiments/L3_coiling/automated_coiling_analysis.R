library(tidyverse)
library(cowplot)
library(conflicted)
library(ggbeeswarm)
library(ggpubr)
library(Hmisc)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("here", "here")

# Import data -------------------------------------------------------------

# get all experiments and create data frame of experiment metadata
experiments <- tibble(File = list.files(pattern = ".*hr.csv$", recursive = TRUE)) %>%
  separate(File, c("Data", "Replicate", "Hour"), sep = "/", remove = FALSE, extra = "drop", fill = "warn") %>%
  mutate(Exp.File = File) %>%
  select(-File)

plate.designs <- tibble(File = list.files(pattern = "plate_design.csv$", recursive = TRUE)) %>%
  separate(File, c("Data", "Replicate"), sep = "/", remove = FALSE, extra = "drop", fill = "warn") %>%
  mutate(Plate.File = File) %>%
  select(-File)

experiments <- left_join(experiments, plate.designs) %>%
  select(-Data)

# Tidy data ---------------------------------------------------------------

tidy_data <- function(experiment) {
  
  plate.design <- read_csv(experiment[4], col_names = c("row", "0001", "0002", "0003", "0004", "0005", "0006", "0007", "0008", "0009", "0010", "0011", "0012"))

  plate.m <- plate.design %>%
    slice(2:9) %>%
    gather(col, Treatment, 2:13) %>%
    mutate(Well = paste0(row, col)) %>%
    select(-row, -col) %>%
    separate(Treatment, c("Drug", "Dose", "Rep"), sep = "_", remove = FALSE) %>%
    filter(Drug != "NA") %>%
    mutate(Replicate = experiment[1])

  data <- read.csv(experiment[3], header = TRUE, sep = ",") %>%
    mutate(Hour = experiment[2])%>%
    mutate(Replicate = experiment[1])

  data <- left_join(data, plate.m)
}

tidy.data <- apply(experiments, 1, tidy_data) %>%
  bind_rows() %>%
  select(Well, Hour, Treatment, Drug, Dose, Replicate, Tech.Rep = Rep, everything(), Motility = Normalized.Motility) %>%
  filter(Hour %in% c("24hr", "48hr"))

# summary
summary <- group_by(tidy.data, Replicate, Hour, Drug, Dose) %>%
  summarise(Mean = mean(Motility))

total.summary <- group_by(tidy.data, Replicate, Hour, Drug, Dose) %>%
  summarise(Mean = mean(Total.Motility))

# control summary
control.summary <- ungroup(summary) %>%
  filter(Drug == "CON") %>%
  select(-Drug, -Dose, Control.Mean = Mean)

control.total.summary <- ungroup(total.summary) %>%
  dplyr::filter(Drug == "CON") %>%
  select(-Drug, -Dose, Control.Mean = Mean)

# add control to data for normalization and normalize to control
tidy.data <- left_join(tidy.data, control.summary) %>%
  mutate(Normalized.Motility = Motility / Control.Mean) %>%
  pivot_longer(cols = Total.Motility:Normalized.Motility, names_to = "Measure", values_to = "Value")

# prepare for plotting
lookup <- data.frame(
  Dose = c("NA", "10-3", "10-4", "10-5", "10-6", "10-7"),
  DEC = c(0.01, 1000, 100, 10, 1, 0.1)
)

tidy.data <- left_join(tidy.data, lookup) %>%
  select(Well, Hour, -Treatment, Drug, Dose, Replicate, DEC, Measure, Value)

summary_stats <- filter(tidy.data, Measure == "Total.Motility") %>%
  group_by(DEC, Hour) %>%
  summarize(Mean = mean(Value), SE = sd(Value)/sqrt(n()))

### 24 hr stats

aov.24 <- aov(data = filter(tidy.data, Measure == "Total.Motility", Hour == "24hr"),
              formula = Value ~ as.factor(DEC))
summary(aov.24)
# one-sided t-test with Holm's adjustment
pairwise.t.test(pluck(filter(tidy.data, Measure == "Total.Motility", Hour == "24hr"), "Value"),
                pluck(filter(tidy.data, Measure == "Total.Motility", Hour == "24hr"), "DEC"),
                alternative = "greater",
                p.adjust.method = "none")

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

significance <- tribble(~Hour, ~DEC, ~Value, ~label,
                        "24hr", "0.1", 8e7, "ns",
                        "24hr", "1", 8e7, "ns",
                        "24hr", "10", 8e7, "*",
                        "24hr", "100", 8e7, "**",
                        "24hr", "1000", 8e7, "**")

### 48 hr stats

aov.48 <- aov(data = filter(tidy.data, Measure == "Total.Motility", Hour == "48hr"),
              formula = Value ~ as.factor(DEC))
summary(aov.48)
# one-sided t-test with Holm's adjustment
pairwise.t.test(pluck(filter(tidy.data, Measure == "Total.Motility", Hour == "48hr"), "Value"),
                pluck(filter(tidy.data, Measure == "Total.Motility", Hour == "48hr"), "DEC"),
                alternative = "greater",
                p.adjust.method = "none")

significance <- bind_rows(significance, tribble(~Hour, ~DEC, ~Value, ~label,
                                                "48hr", "0.1", 8e7, "ns",
                                                "48hr", "1", 8e7, "*",
                                                "48hr", "10", 8e7, "**",
                                                "48hr", "100", 8e7, "**",
                                                "48hr", "1000", 8e7, "***"))

final.plot <- ggplot(filter(tidy.data, Measure == "Total.Motility"), aes(x = as.factor(DEC), y = Value / 1e7)) + # y scale is arbitrary, divide units by 1e7
  geom_pointrange(data = summary_stats, 
                  aes(x = as.factor(DEC), y = Mean / 1e7, ymin = (Mean - SE) / 1e7, ymax = (Mean + SE) / 1e7),
                  color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = significance, aes(label = label)) +
  scale_x_discrete(limits = c("0.01", "0.1", "1", "10", "100", "1000"),
                   labels = c("Control", "0.1", "1", "10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 8)) +
  facet_grid(rows = vars(Hour)) +
  labs(x = "NAM (µM)", y = "Mean Motility Units") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, vjust = .5),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.y = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    axis.ticks.y = element_line(size = 0.5),
    legend.position = "none"
  ) +
  NULL
final.plot

save_plot(here("plots", "Fig5D.pdf"), final.plot, base_width = 3, base_height = 6)

