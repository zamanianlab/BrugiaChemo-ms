library(tidyverse)
library(cowplot)
library(conflicted)
library(ggbeeswarm)
library(ggpubr)
library(Hmisc)
library(here)

conflict_prefer("filter", "dplyr")


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
    dplyr::filter(Drug != "NA") %>%
    mutate(Replicate = experiment[1])

  data <- read.csv(experiment[3], header = TRUE, sep = ",") %>%
    mutate(Hour = experiment[2])%>%
    mutate(Replicate = experiment[1])

  data <- left_join(data, plate.m)
}

tidy.data <- apply(experiments, 1, tidy_data) %>%
  bind_rows() %>%
  select(Well, Hour, Treatment, Drug, Dose, Replicate, Tech.Rep = Rep, everything(), Motility = Normalized.Motility)

# remove outliers
motility.quantile <- quantile(tidy.data$Motility, probs = c(.25, .75))
trimmed.data <- filter(tidy.data, Motility > motility.quantile[1] - 1.5 * IQR(tidy.data$Motility), Motility < motility.quantile[2] + 1.5 * IQR(tidy.data$Motility))

total.quantile <- quantile(tidy.data$Total.Motility, probs = c(.25, .75))
trimmed.data <- filter(trimmed.data, Total.Motility > total.quantile[1] - 1.5 * IQR(tidy.data$Total.Motility), Total.Motility < total.quantile[2] + 1.5 * IQR(tidy.data$Total.Motility))

norm_quantile <- quantile(tidy.data$Normalization.Factor, probs = c(.25, .75))
trimmed.data <- filter(trimmed.data, Normalization.Factor > norm_quantile[1] - 1.5 * IQR(tidy.data$Normalization.Factor), Normalization.Factor < norm_quantile[2] + 1.5 * IQR(tidy.data$Normalization.Factor))

# summary
summary <- group_by(trimmed.data, Replicate, Hour, Drug, Dose) %>%
  summarise(Mean = mean(Motility))

total.summary <- group_by(trimmed.data, Replicate, Hour, Drug, Dose) %>%
  summarise(Mean = mean(Total.Motility))

# control summary
control.summary <- ungroup(summary) %>%
  dplyr::filter(Drug == "CON") %>%
  select(-Drug, -Dose, Control.Mean = Mean)

control.total.summary <- ungroup(total.summary) %>%
  dplyr::filter(Drug == "CON") %>%
  select(-Drug, -Dose, Control.Mean = Mean)

# add control to data for normalization and normalize to control
trimmed.data <- left_join(trimmed.data, control.summary) %>%
  mutate(Normalized.Motility = Motility / Control.Mean) %>%
  pivot_longer(cols = Total.Motility:Normalized.Motility, names_to = "Measure", values_to = "Value")

# prepare for plotting
lookup <- data.frame(
  Dose = c("NA", "10-3", "10-4", "10-5", "10-6", "10-7"),
  DEC = c(0.01, 1000, 100, 10, 1, 0.1)
)

trimmed.data <- left_join(trimmed.data, lookup) %>%
  # mutate(DEC = as.factor(DEC)) %>%
  select(Well, Hour, -Treatment, Drug, Dose, Replicate, DEC, Measure, Value)

comparisons <- list(c("0.01", "0.1"), c("0.01", "1"), c("0.01", "10"), c("0.01", "100"), c("0.01", "1000"))

all.plot <- ggplot(filter(trimmed.data, Hour %in% c("24hr", "48hr"), Measure == "Normalized.Motility"), aes(x = as.character(DEC), y = Value)) +
  # geom_beeswarm(aes(color = Replicate)) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.75, size = 0.5) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif") +
  scale_x_discrete(limits = c("0.01", "0.1", "1", "10", "100", "1000"),
                   labels = c("Control", "0.1", "1", "10", "100", "1000")) +
  # scale_y_continuous(limits = c(-1, 3.5), breaks = seq(0, 3, 1)) +
  facet_grid(rows = vars(Hour)) +
  labs(x = "NAM (µM)", y = expression(bold(paste("log"[2], "(Mean Movement Units) + 1")))) +
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
    legend.position = "none"
  ) +
  NULL
all.plot

aov <- aov(Value ~ as.factor(DEC), data = filter(trimmed.data, Hour == "48hr", Measure == "Total.Motility"))
TukeyHSD(aov)


save_plot(here("plots", "Fig5D.pdf"), all.plot, base_width = 3, base_height = 6)

