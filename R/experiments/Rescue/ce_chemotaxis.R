library(tidyverse)
library(conflicted)
library(ggbeeswarm)
library(cowplot)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("get_legend", "cowplot")

# Import and tidy ---------------------------------------------------------

tidy_data <- readRDS(here("data", "ce.chemotaxis.tidy.data")) %>%
  bind_rows()

# # Plotting ----------------------------------------------------------------

pairwise.t.test(
  pluck(filter(
    tidy_data, Cue == "Diacetyl",
    Region == "CI",
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    Notes == "10 cm CTX",
    Date %in% c("20181217", "20181218", "20181221", "20190112", "20190115", "20190118", "20191002", "20191014", "20191015", "20191102", "20191106", "20200128")
  ), "Value"),
  pluck(filter(
    tidy_data, Cue == "Diacetyl",
    Region == "CI",
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    Notes == "10 cm CTX",
    Date %in% c("20181217", "20181218", "20181221", "20190112", "20190115", "20190118", "20191002", "20191014", "20191015", "20191102", "20191106", "20200128")
  ), "Strain")
)

s11.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 1, "****",
  "ZAM13", 1, "ns",
  "ZAM17", 1, "ns",
  "ZAM18", 1, "ns"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

s11.plot <- ggplot(
  filter(
    tidy_data,
    Cue == "Diacetyl",
    Region == "CI",
    Strain %in% c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17"),
    Notes == "10 cm CTX",
    Date %in% c("20181217", "20181218", "20181221", "20190112", "20190115", "20190118", "20191002", "20191014", "20191015", "20191102", "20191106", "20200128")
  ),
  aes(x = Strain, y = Value)
) +
  geom_violin(alpha = 0.6, color = "black", size = 0.7, fill = "grey90") +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = s11.significance, aes(x = Strain, y = Value, label = label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
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
  scale_y_continuous(limits = c(-0.6, 1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Strain", y = "Chemotaxis Index", title = "Diacetyl Chemotaxis (AWA)") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold.italic", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(angle = 90, vjust = 0.5, size = 0),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none"
  ) +
  NULL
s11.plot

save_plot(here("plots", "S11A_Figure.pdf"), s11.plot, base_height = 6, base_width = 8)


pairwise.t.test(
  pluck(filter(
    tidy_data, Cue == "Diacetyl",
    Region == "CI",
    Strain %in% c("N2", "CX10", "ZAM22", "ZAM24"),
    Notes == "10 cm CTX",
    Date %in% c("20200121", "20200122", "20200125")
  ), "Value"),
  pluck(filter(
    tidy_data, Cue == "Diacetyl",
    Region == "CI",
    Strain %in% c("N2", "CX10", "ZAM22", "ZAM24"),
    Notes == "10 cm CTX",
    Date %in% c("20200121", "20200122", "20200125")
  ), "Strain")
)

s12.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 1, "**",
  "ZAM22", 1, "ns",
  "ZAM24", 1, "ns"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

s12.plot <- ggplot(
  filter(
    tidy_data,
    Cue == "Diacetyl",
    Region == "CI",
    Strain %in% c("N2", "CX10", "ZAM22", "ZAM24"),
    Notes == "10 cm CTX",
    Date %in% c("20200121", "20200122", "20200125")
  ),
  aes(x = Strain, y = Value)
) +
  geom_violin(alpha = 0.6, color = "black", size = 0.7, fill = "grey90") +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = s12.significance, aes(x = Strain, y = Value, label = label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
  scale_x_discrete(
    limits = c("N2", "CX10", "ZAM22", "ZAM24"),
    labels = c(
      "N2",
      "osm-9(ky10)",
      "osm-9(ky10);\ncel-osm-9::\nosm-9 3' UTR",
      "osm-9(ky10);\nbma-osm-9::\nosm-9 3' UTR"
    )
  ) +
  scale_y_continuous(limits = c(-0.6, 1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Strain", y = "Chemotaxis Index", title = "Diacetyl Chemotaxis (AWA)") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold.italic", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(angle = 90, vjust = 0.5, size = 0),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none"
  ) +
  NULL
s12.plot

save_plot(here("plots", "S12_Figure.pdf"), s12.plot, base_height = 6, base_width = 6)

pairwise.t.test(pluck(filter(tidy_data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), "Value"),
  pluck(filter(tidy_data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), "Strain"),
  p.adjust.method = "none"
)

tax4.significance <- tribble(
  ~Strain, ~Value, ~label,
  "N2", 1.1, "****",
  "ZAM14", 1.1, "*",
  "ZAM21", 1.1, "****"
)

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

tax4.plot <- ggplot(filter(tidy_data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), aes(x = Strain, y = Value)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.7, fill = "grey90") +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  geom_text(data = tax4.significance, aes(x = Strain, y = Value, label = label)) +
  scale_x_discrete(
    limits = c("N2", "PR678", "ZAM21", "ZAM14"),
    labels = c("N2", "tax-4(p678)", "tax-4(p678);\ncel-tax-4", "tax-4(p678);\nbma-tax-4")
  ) +
  # facet_grid(. ~ Cue, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Strain", y = "Chemotaxis Index", title = "Isoamyl Alcohol Chemotaxis (AWC)") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold.italic", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(angle = 90, vjust = 0.5, size = 0),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "none"
  ) +
  NULL
tax4.plot

save_plot(here("plots", "Fig9.pdf"), tax4.plot, base_height = 6, base_width = 8)
