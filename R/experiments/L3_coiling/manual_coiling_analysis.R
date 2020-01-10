library(tidyverse)
library(Hmisc)
library(cowplot)
library(here)
library(conflicted)

conflict_prefer("here", "here")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# get all experiments and create data frame of experiment metadata
experiments <- tibble(File = list.files(pattern = "key.csv$", recursive = TRUE)) %>%
  separate(File, c("Data", "Rep", "Hours"), sep = "/", remove = FALSE, extra = "drop", fill = "warn") %>%
  select(-Data, Data.File = File)

plates <- tibble(File = list.files(pattern = "plate_design.csv$", recursive = TRUE)) %>%
  separate(File, c("Data", "Rep"), sep = "/", remove = FALSE, extra = "drop", fill = "warn") %>%
  select(-Data, Plate.File = File)

experiments <- left_join(experiments, plates)

# Date Tidying ------------------------------------------------------------

tidy_data <- function(experiments) {
  
  plate <- read_csv(experiments[4], col_names = c("row", "0001", "0002", "0003", "0004", "0005", "0006", "0007", "0008", "0009", "0010", "0011", "0012")) %>%
    slice(2:9) %>%
    gather(col, Treatment, 2:13) %>% 
    mutate(Well = paste0(row, col)) %>%
    select(-row, -col) %>%
    separate(Treatment ,c("Drug", "Dose", "Rep"), sep = "_", remove = FALSE) %>%
    dplyr::filter(Drug != "NA")
  
  lookup <- data.frame(Dose = c("NA", "10-3", "10-4", "10-5", "10-6", "10-7"),
                        DEC = c(0.01, 1000, 100, 10, 1, 0.1))
  
  data <- read_csv(experiments[1], col_names = TRUE) %>%
    select(-Code, -"MD5-lite") %>%
    mutate(Hours = experiments[3]) %>%
    gather(Scorer, Score, -Well, -Hours) %>%
    left_join(., plate) %>%
    select(-Rep) %>%
    left_join(., lookup) %>%
    mutate(Bio.Rep = experiments[2])
  
}

tidy.data <- apply(experiments, 1, tidy_data) %>%
  bind_rows() %>%
  filter(Hours == "48hr" | Scorer != "Tran") # remove scores from scorer who used incorrect scale

# summarize
summary <- group_by(tidy.data, Bio.Rep, Drug, Dose, Hours) %>%
  summarise(Mean = mean(Score))

# control summary
control.summary <- ungroup(summary) %>%
  dplyr::filter(Drug == "CON") %>%
  select(-Drug, -Dose, Control.Score = Mean)

# add control to data for normalization and normalize to control
tidy.data <- left_join(tidy.data, control.summary) %>%
  mutate(Normalized.Score = Score / Control.Score) %>%
  mutate(Dose = as.factor(Dose))

# Plotting ----------------------------------------------------------------

summary_stats <- tidy.data %>%
  group_by(DEC, Hours) %>%
  summarize(Mean = mean(Score), SE = sd(Score)/sqrt(n()))

### 24 hr stats

aov.24 <- aov(data = filter(tidy.data, Hours == "24hr"),
              formula = Score ~ as.factor(DEC))
summary(aov.24)
# one-sided t-test with Holm's adjustment
pairwise.t.test(pluck(filter(tidy.data, Hours == "24hr"), "Score"),
                pluck(filter(tidy.data, Hours == "24hr"), "DEC"),
                alternative = "greater",
                p.adjust.method = "none")

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

significance <- tribble(~Hours, ~DEC, ~Score, ~label,
                        "24hr", "0.1", 5, "***",
                        "24hr", "1", 5, "****",
                        "24hr", "10", 5, "****",
                        "24hr", "100", 5, "****",
                        "24hr", "1000", 5, "****")

### 48 hr stats

aov.48 <- aov(data = filter(tidy.data, Hours == "48hr"),
              formula = Score ~ as.factor(DEC))
summary(aov.48)
# one-sided t-test with Holm's adjustment
pairwise.t.test(pluck(filter(tidy.data, Hours == "48hr"), "Score"),
                pluck(filter(tidy.data, Hours == "48hr"), "DEC"),
                alternative = "greater",
                p.adjust.method = "none")

significance <- bind_rows(significance, tribble(~Hours, ~DEC, ~Score, ~label,
                                                "48hr", "0.1", 5, "ns",
                                                "48hr", "1", 5, "**",
                                                "48hr", "10", 5, "****",
                                                "48hr", "100", 5, "****",
                                                "48hr", "1000", 5, "****"))



final.plot <- ggplot(tidy.data, aes(x = as.character(DEC), y = Score, group = DEC)) + 
  geom_pointrange(data = summary_stats, 
                  aes(x = as.factor(DEC), y = Mean, ymin = Mean - SE, ymax = Mean + SE),
                  color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = significance, aes(label = label)) +
  scale_x_discrete(limits = c("0.01", "0.1", "1", "10", "100", "1000"), labels = c("Control", "0.1", "1", "10", "100", "1000")) +
  # scale_x_log10(
  #   limits = c(0.0035, 3003),
  #   breaks = c(0.005, 0.01, 0.1, 1, 10, 100, 1000),
  #   labels = c("", "Control", "0.1", "1", "10", "100", "1000")
  # ) +
  scale_y_continuous(limits = c(0, 5.1), breaks = c(0, 1, 2, 3, 4, 5)) +
  facet_grid(rows = vars(Hours)) +
  labs(x = "NAM (µM)", y = "Coiling Score") +
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
final.plot

save_plot(here("plots", "Fig5C.pdf"), final.plot, base_width = 3, base_height = 6)
