library(tidyverse)
library(conflicted)
library(magrittr)
library(ggbeeswarm)
library(Hmisc)
library(cowplot)
library(here)
library(ggpubr)

here()

# Import and tidy data Tidying --------------------------------------------
tidy.data <- read_csv(here("data", "data.csv")) %>%
  mutate(CI = (T_n - C_n)/(T_n + C_n + M_n + O_n)) %>%
  mutate(CI2 = (T_n - C_n)/(T_n + C_n)) %>%
  mutate(CI3 = (T_n - C_n)/(T_n + C_n + O_n)) %>%
  mutate(X = (M_n + O_n)/(T_n + C_n + M_n + O_n)) %>%
  mutate(Mf = (M_n)/(T_n + C_n + M_n + O_n)) %>%
  select(-Notes) 

# change NaN > 0 for all CIs
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
tidy.data[is.nan(tidy.data)] <- 0

# Plotting ----------------------------------------------------------------

comparisons <- list(c("lacZ", "tax-4"), c("lacZ", "osm-9"))

dot.plot <- ggplot(dplyr::filter(tidy.data, T_n + C_n > 2), aes(x = dsRNA, y = CI3)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", alpha = 0.8, dotsize = 0.7) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  # stat_summary(fun.y = mean, geom = "point", color = "red", shape = 18, size = 5) +
  # geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  # geom_beeswarm(aes(color = as.factor(Date)), size = 3.25, alpha = 0.6, groupOnX = TRUE, cex = 1.5) +
  geom_hline(aes(yintercept = 0), linetype="dashed", colour = "grey37", size = 0.5) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  scale_x_discrete(limits = c("lacZ", "osm-9", "tax-4"), labels = c("lacZ(RNAi)", "osm-9(RNAi)", "tax-4(RNAi)")) +
  scale_y_continuous(limits = c(-1.2, 1.5), breaks = c(-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1)) +
  # facet_grid(. ~ Injection.Day) +
  ylab("Fraction leaving center") + xlab("") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold.italic"),
        axis.text.x = element_text(face = "bold", size = 10, angle = 45, vjust = .5),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.x  = element_text(angle = 90, vjust = 0.5, size = 0), 
        axis.title.y  = element_text(face = "bold", angle = 90, size = 12),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.75, colour = "black"),
        legend.position = "none") +
  labs(color = "Trial", x = "Target, Injection Date", y = "Chemotaxis Index (CI)") +
  NULL
dot.plot

save_plot(here("plots", "Fig7B_raw.pdf"), dot.plot, base_width = 5)

move.plot <- ggplot(dplyr::filter(tidy.data, T_n + C_n > 2), aes(x = dsRNA, y = 1 - Mf)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", alpha = 0.8, dotsize = 0.7) +
  # stat_summary(fun.y = mean, geom = "point", color = "red", shape = 18, size = 5) +
  # geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  # geom_beeswarm(aes(color = as.factor(Date)), size = 3.25, alpha = 0.6, groupOnX = TRUE, cex = 1.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  # scale_x_discrete(limits = c("lacZ-9", "tax-4-9", "tax-4-10", "osm-9-10"), labels = c("lacZ 9 DPI", "tax-4 9 DPI", "tax-4 10 DPI", "osm-9 10 DPI")) +
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
  scale_colour_manual(values = c("royalblue2","royalblue3","steelblue4", "steelblue3", "steelblue2","steelblue1")) +
  ylab("Fraction leaving center") + xlab("") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold.italic"),
        axis.text.x = element_text(face = "bold", size = 10, angle = 45, vjust = .5),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.x  = element_text(angle = 90, vjust = 0.5, size = 0), 
        axis.title.y  = element_text(face = "bold", angle = 90, size = 12),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.75, colour = "black"),
        legend.position = "none") +
  labs(color = "Trial", x = "Target, Injection Date", y = "Fraction Leaving Center") +
  NULL
move.plot

save_plot(here("plots", "Fig7C_raw.pdf"), move.plot, base_width = 5)
