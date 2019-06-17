library(tidyverse)
library(conflicted)
library(magrittr)
library(ggbeeswarm)
library(Hmisc)
library(cowplot)
library(lubridate)
library(here)
library(ggpubr)

conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
here()


# Date Tidying ------------------------------------------------------------

tidy.data <- read_csv(here("data", "data.csv"), col_names = TRUE) %>%
  mutate(CI = (Cue_N - Control_N) / Total) %>%
  select(Date, Strain, Cue, Control, Outside, Notes, everything()) %>%
  pivot_longer(Cue_N:CI, names_to = "Region", values_to = "Value") %>%
  mutate(Date = factor(Date))

# Plotting ----------------------------------------------------------------

osm9.comparisons <- list(c("N2", "CX10"), c("N2", "ZAM13"), c("CX10", "ZAM13"))

osm9.plot <- ggplot(dplyr::filter(tidy.data, Cue == "Diacetyl", Region == "CI", Notes == "10 cm CTX"), aes(x = Strain, y = Value)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(size = 3.25, alpha = 0.6, groupOnX = TRUE, cex = 3.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
  scale_x_discrete(limits = c("N2", "CX10", "ZAM13"), labels = c("N2", "osm-9(ky10)", "osm-9(ky10);\n osm-9p::bm-osm-9")) +
  facet_grid(. ~ Cue, scales = "free_x") +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1.5) +
  stat_compare_means(comparisons = osm9.comparisons, method = "t.test", label = "p.signif", label.y = c(1.05, 1.2, 1.35)) +
  scale_y_continuous(limits = c(-0.6, 1.4), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.x = element_text(face = "bold.italic", size = 10, angle = 45, vjust = .5),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.x  = element_text(angle = 90, vjust = 0.5, size = 0), 
        axis.title.y  = element_text(face = "bold", angle=90, size=12),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.75, colour = "black"),
        legend.position = "none") +
  labs(x = "Strain", y = "Chemotaxis Index") +
  NULL
osm9.plot

save_plot(here("plots", "Fig8A_raw.pdf"), osm9.plot, base_width = 4, base_height = 5)

tax4.comparisons <- list(c("N2", "PR678"), c("N2", "ZAM14"), c("PR678", "ZAM14"))

tax4.plot <- ggplot(dplyr::filter(tidy.data, Cue == "Isoamyl alcohol", Region == "CI", Notes == "10 cm CTX"), aes(x = Strain, y = Value)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(size = 3.25, alpha = 0.6, groupOnX = TRUE, cex = 3.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
  scale_x_discrete(limits = c("N2", "PR678", "ZAM14"), labels = c("N2", "tax-4(p678)", "tax-4(p678);\n tax-4p::bm-tax-4")) +
  scale_color_viridis_d() +
  facet_grid(. ~ Cue, scales = "free_x") +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1.5) +
  stat_compare_means(comparisons = tax4.comparisons, method = "t.test", label = "p.signif", label.y = c(1.05, 1.2, 1.35)) +
  scale_y_continuous(limits = c(-0.6, 1.4), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.x = element_text(face = "bold.italic", size = 10, angle = 45, vjust = .5),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.x  = element_text(angle = 90, vjust = 0.5, size = 0), 
        axis.title.y  = element_text(face = "bold", angle=90, size=12),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.75, colour = "black")) +
  labs(x = "Strain", y = "Chemotaxis Index", color = "Experiment Date") +
  NULL
tax4.plot

save_plot(here("plots", "Fig8B_raw.pdf"), tax4.plot, base_width = 4, base_height = 5)