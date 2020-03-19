library(tidyverse)
library(conflicted)
library(ggbeeswarm)
library(cowplot)
library(lubridate)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("get_legend", "cowplot")


# Import and tidy ---------------------------------------------------------

import_tidy <- function(experiments) {

  # print(experiments[1])
  data <- read_csv(experiments[1], col_names = TRUE) %>%
    mutate(Date = experiments[2])
}

tidy_data <- readRDS(here("data", "ce.nosetouch.tidy.data")) %>%
  filter(Strain != "ZAM23") %>%
  mutate(Date = as.character(Date)) %>%
  mutate(Reversal_Rate = Reversals / Touches)

# Plotting ----------------------------------------------------------------

pairwise.t.test(pluck(
  filter(
    tidy_data,
    Date %in% c("2019-12-18", "2019-12-19", "2020-01-09")
  ),
  "Reversal_Rate"
),
pluck(
  filter(
    tidy_data,
    Date %in% c("2019-12-18", "2019-12-19", "2020-01-09")
  ),
  "Strain"
),
alternative = "g",
p.adjust.method = "none"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

nose.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 1.1, "***",
  "ZAM22", 1.1, "****",
  "ZAM24", 1.1, "***"
)

reversal_rate <- ggplot(
  filter(
    tidy_data,
    Date %in% c("2019-12-18", "2019-12-19", "2020-01-09"),
    Strain %in% c("N2", "CX10", "ZAM22", "ZAM24")
  ),
  aes(x = Strain, y = Reversal_Rate)
) +
  geom_violin(alpha = 0.6, color = "black", size = 0.7, fill = "grey90") +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = nose.significance, aes(x = Strain, y = Value, label = label)) +
  scale_x_discrete(
    limits = c("N2", "CX10", "ZAM24", "ZAM22"),
    labels = c("N2", "osm-9(ky10)", "osm-9(ky10);\ncel-osm-9", "osm-9(ky10);\nbma-osm-9")
  ) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold.italic", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none"
  ) +
  labs(x = "Strain", y = "Reversal Rate", title = "Nose Touch Reversal (ASH)") +
  NULL
reversal_rate

save_plot(here("plots", "Fig8B_raw.pdf"), reversal_rate, base_height = 6, base_width = 8)

pairwise.t.test(pluck(
  filter(
    tidy_data,
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    Date %in% c("2019-11-12", "2019-11-13", "2020-01-14", "2020-01-15")
  ),
  "Reversal_Rate"
),
pluck(
  filter(
    tidy_data,
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    Date %in% c("2019-11-12", "2019-11-13", "2020-01-14", "2020-01-15")
  ),
  "Strain"
),
alternative = "g",
p.adjust.method = "none"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

supp.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 1.1, "****",
  "ZAM18", 1.1, "ns",
  "ZAM13", 1.1, "ns",
  "ZAM17", 1.1, "ns"
)

supp_plot <- ggplot(
  filter(
    tidy_data,
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    Date %in% c("2019-11-12", "2019-11-13", "2020-01-14", "2020-01-15")
  ),
  aes(x = Strain, y = Reversal_Rate)
) +
  geom_violin(alpha = 0.6, color = "black", size = 0.7, fill = "grey90") +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = supp.significance, aes(x = Strain, y = Value, label = label)) +
  scale_x_discrete(
    limits = c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    labels = c(
      "N2",
      "osm-9(ky10)",
      "osm-9(ky10);\ncel-osm-9::\nunc-54 3' UTR",
      "osm-9(ky10);\nbma-osm-9::\nunc-54 3' UTR",
      "osm-9(ky10);\nbma-osm-9::\nunc-54 3' UTR +\nbma-ocr-1/2a::\nunc-54 3' UTR"
    )
  ) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold.italic", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none"
  ) +
  labs(x = "Strain", y = "Reversal Rate", title = "Nose Touch Reversal (ASH)") +
  NULL
supp_plot

save_plot(here("plots", "S11C_Figure_raw.pdf"), supp_plot, base_height = 6, base_width = 8)
