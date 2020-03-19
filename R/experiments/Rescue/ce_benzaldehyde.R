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

tidy_data <- readRDS(here("data", "ce.benzaldehyde.tidy.data")) %>%
  filter(Strain != "ZAM23") %>%
  mutate(Date = as.character(Date))

# Plotting ----------------------------------------------------------------

pairwise.t.test(pluck(filter(tidy_data, Notes == "Unseeded plate", Date %in% c("2019-12-17", "2019-12-18")), "Time_To_Reversal"),
  pluck(filter(tidy_data, Notes == "Unseeded plate", Date %in% c("2019-12-17", "2019-12-18")), "Strain"),
  alternative = "less",
  p.adjust.method = "none"
)

benz.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 22, "****",
  "ZAM22", 22, "****",
  "ZAM24", 22, "****"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

time_reversal_plot <- ggplot(
  filter(
    tidy_data,
    Notes == "Unseeded plate",
    Strain %in% c("N2", "CX10", "ZAM22", "ZAM24"),
    Date %in% c("2019-12-17", "2019-12-18")
  ),
  aes(x = Strain, y = Time_To_Reversal)
) +
  geom_violin(alpha = 0.6, color = "black", size = 0.7, fill = "grey90") +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = benz.significance, aes(x = Strain, y = Value, label = label)) +
  scale_x_discrete(
    limits = c("N2", "CX10", "ZAM24", "ZAM22"),
    labels = c("N2", "osm-9(ky10)", "osm-9(ky10);\ncel-osm-9", "osm-9(ky10);\nbma-osm-9")
  ) +
  scale_y_continuous(limits = c(0, 22)) +
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
  labs(x = "Strain", y = "Time to Reversal", title = "Benzaldehyde Avoidance (ASH)") +
  NULL
time_reversal_plot

save_plot(here("plots", "Fig8A_raw.pdf"), time_reversal_plot, base_height = 6, base_width = 8)

pairwise.t.test(pluck(
  filter(
    tidy_data, Notes == "Unseeded plate",
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    # Person != "Tran",
    Date %in% c("2020-01-13", "2020-01-14", "2019-11-14", "2019-11-15", "2020-01-14")
  ),
  "Time_To_Reversal"
),
pluck(
  filter(
    tidy_data, Notes == "Unseeded plate",
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    # Person != "Tran",
    Date %in% c("2020-01-13", "2020-01-14", "2019-11-14", "2019-11-15", "2020-01-14")
  ),
  "Strain"
),
alternative = "less",
p.adjust.method = "none"
)

supp.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 22, "*",
  "ZAM18", 22, "ns",
  "ZAM13", 22, "ns",
  "ZAM17", 22, "ns"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

supp_plot <- ggplot(
  filter(
    tidy_data,
    Notes == "Unseeded plate",
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    # Person != "Tran",
    Date %in% c("2020-01-13", "2020-01-14", "2019-11-14", "2019-11-15", "2020-01-14")
  ),
  aes(x = Strain, y = Time_To_Reversal)
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
  scale_y_continuous(limits = c(0, 22)) +
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
  labs(x = "Strain", y = "Time to Reversal", title = "Benzaldehyde Avoidance (ASH)") +
  NULL
supp_plot

save_plot(here("plots", "S11B_Figure_raw.pdf"), supp_plot, base_height = 6, base_width = 8)
