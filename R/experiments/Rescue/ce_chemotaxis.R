library(tidyverse)
library(conflicted)
library(ggbeeswarm)
library(cowplot)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("get_legend", "cowplot")

# Import and tidy ---------------------------------------------------------

tidy_data <- readRDS(here("data", "ce.chemotaxis.tidy.data")) %>%
  bind_rows() %>%
  mutate(Date = factor(Date))

# Plotting ----------------------------------------------------------------

genes <- tribble(~Strain, ~Element, ~Start, ~Stop,
                 "N2", NA, 0, 0,
                 "CX10", NA, 0, 0,
                 "ZAM18", "cel-osm-9p", 3.5, 4.5,
                 "ZAM18", "cel-osm-9", 4.5, 5.5,
                 "ZAM18", "cel-unc-54", 5.5, 6.5,
                 "ZAM13", "cel-osm-9p", 3.5, 4.5,
                 "ZAM13", "bma-osm-9", 4.5, 5.5,
                 "ZAM13", "cel-unc-54", 5.5, 6.5,
                 "ZAM17", "cel-osm-9p", 3.5, 4.5,
                 "ZAM17", "bma-osm-9", 4.5, 5.5,
                 "ZAM17", "cel-unc-54", 5.5, 6.5,
                 "ZAM17", "cel-osm-9p", 0, 1,
                 "ZAM17", "bma-ocr-1/2a", 1, 2,
                 "ZAM17", "cel-unc-54", 2, 3,
                 "ZAM24", "cel-osm-9p", 3.5, 4.5,
                 "ZAM24", "cel-osm-9", 4.5, 5.5,
                 "ZAM24", "cel-osm-9", 5.5, 6.5,
                 "ZAM22", "cel-osm-9p", 3.5, 4.5,
                 "ZAM22", "bma-osm-9", 4.5, 5.5,
                 "ZAM22", "cel-osm-9", 5.5, 6.5
)

model <- ggplot(filter(genes, !is.na(Element))) +
  geom_tile(aes(x = Start, y = Strain, fill = Element), width = 1, height = 0.3, color = "black", size = 0.3, alpha = 0.75) +
  annotate("text", x = 2.75, y = 3, label = "+") +
  scale_y_discrete(limits = rev(c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17", "ZAM24", "ZAM22"))) +
  scale_fill_manual(limits = c("cel-osm-9p", "cel-osm-9", "bma-osm-9", "bma-ocr-1/2a", "cel-unc-54"), 
                    values = c("grey", "steelblue", "red", "darkgreen", "purple")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank()) +
  NULL
# model

pairwise.t.test(pluck(filter(tidy_data, Cue == "Diacetyl", Region == "CI", Notes == "10 cm CTX"), "Value"),
                pluck(filter(tidy_data, Cue == "Diacetyl", Region == "CI", Notes == "10 cm CTX"), "Strain"))

osm9.significance <- tribble(~Strain, ~Value, ~label,
                        "N2", 1, "****",
                        "ZAM13", 1, "ns",
                        "ZAM17", 1, "ns",
                        "ZAM18", 1, "ns",
                        "ZAM22", 1, "ns",
                        "ZAM24", 1, "ns")

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

osm9.plot <- ggplot(filter(tidy_data, Cue == "Diacetyl", Region == "CI", Notes == "10 cm CTX"), aes(x = Strain, y = Value)) +
  geom_boxplot() +
  geom_beeswarm(size = 1.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  geom_text(data = osm9.significance, aes(x = Strain, y = Value, label = label)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
  scale_x_discrete(limits = rev(c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17", "ZAM24", "ZAM22")), labels = rev(c("N2", "osm-9(ky10)", "", "",  "", "", ""))) +
  facet_grid(. ~ Cue, scales = "free_x") +
  scale_y_continuous(limits = c(-0.6, 1), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.x = element_text(face = "plain", size = 10),
        axis.text.y = element_text(face = "italic", size = 10),
        axis.title.x  = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.75, colour = "black"),
        legend.position = "none") +
  labs(x = "Strain", y = "Chemotaxis Index") +
  coord_flip() +
  NULL
osm9.plot

legend <- get_legend(model)

osm <- plot_grid(model + theme(legend.position = "none", plot.margin = margin(r = -70, l = 10)), osm9.plot, axis = "tb", align = "h", rel_widths = c(0.15, 1))
final.osm <- plot_grid(osm, legend, nrow = 2, rel_heights = c(1, 0.1))

save_plot(here("plots", "Fig8A_raw.pdf"), final.osm, base_height = 6, base_width = 8)

pairwise.t.test(pluck(filter(tidy_data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), "Value"),
                pluck(filter(tidy_data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), "Strain"),
                p.adjust.method = "none")

tax4.significance <- tribble(~Strain, ~Value, ~label,
                        "N2", 1.1, "****",
                        "ZAM14", 1.1, "*",
                        "ZAM21", 1.1, "****")

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

tax4.plot <- ggplot(filter(tidy_data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), aes(x = Strain, y = Value)) +
  geom_boxplot() +
  geom_beeswarm(size = 2.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  geom_text(data = tax4.significance, aes(x = Strain, y = Value, label = label)) +
  scale_x_discrete(limits = rev(c("N2", "PR678", "ZAM21", "ZAM14")), labels = rev(c("N2", "tax-4(p678)", "tax-4(p678); tax-4p::ce-tax-4", "tax-4p::bm-tax-4"))) +
  facet_grid(. ~ Cue, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.x = element_text(face = "plain", size = 10),
        axis.text.y = element_text(face = "italic", size = 10),
        axis.title.x  = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.75, colour = "black"),
        legend.position = "none") +
  labs(x = "Strain", y = "Chemotaxis Index") +
  coord_flip() +
  NULL
tax4.plot

save_plot(here("plots", "Fig9.pdf"), tax4.plot, base_height = 6, base_width = 8)

