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
  mutate(Date = ymd(Date)) %>%
  mutate(Reversal_Rate = Reversals / Touches)

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

pairwise.t.test(pluck(tidy_data, "Reversal_Rate"),
                pluck(tidy_data, "Strain"),
                alternative = "g")

# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

significance <- tribble(~Strain, ~Value, ~label,
                        "N2", 1.1, "****",
                        "ZAM13", 1.1, "ns",
                        "ZAM17", 1.1, "ns",
                        "ZAM18", 1.1, "ns",
                        "ZAM22", 1.1, "***",
                        "ZAM24", 1.1, "*")

reversal_rate <- ggplot(filter(tidy_data, Date != "2019-12-18"), aes(x = Strain, y = Reversal_Rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 1.5, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  geom_text(data = significance, aes(x = Strain, y = Value, label = label)) +
  scale_x_discrete(
    limits = rev(c("N2", "CX10", "ZAM18", "ZAM13", "ZAM17", "ZAM24", "ZAM22")),
    labels = rev(c("N2", "osm-9(ky10)", "", "",  "", "", ""))
  ) +
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
  labs(x = "Strain", y = "Reversal Rate") +
  coord_flip() +
  NULL
# reversal_rate

legend <- get_legend(model)

osm <- plot_grid(model + theme(legend.position = "none", plot.margin = margin(r = -70, l = 10)), reversal_rate, axis = "tb", align = "h", rel_widths = c(0.15, 1))
final.osm <- plot_grid(osm, legend, nrow = 2, rel_heights = c(1, 0.1))

save_plot(here("plots", "Fig8C_raw.pdf"), final.osm, base_height = 6, base_width = 8)
