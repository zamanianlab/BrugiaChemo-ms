library(tidyverse)
library(conflicted)
library(ggbeeswarm)
library(Hmisc)
library(cowplot)
library(here)
library(ggpubr)

conflict_prefer("here", "here")
conflict_prefer("filter", "dplyr")

here()

# Import Data -------------------------------------------------------------

tidy.data <- readRDS(here("data", "NAM.chemotaxis.tidy.data")) %>%
  mutate(CI = (T_n - C_n) / (T_n + C_n + M_n + O_n)) %>%
  mutate(CI2 = (T_n - C_n) / (T_n + C_n)) %>%
  mutate(CI3 = (T_n - C_n) / (T_n + C_n + O_n)) %>%
  mutate(X = (M_n + O_n) / (T_n + C_n + M_n + O_n)) %>%
  mutate(Mf = (M_n) / (T_n + C_n + M_n + O_n)) 

# Plot data ----------------------------------------------------------------

comparisons <- list(c("Nicotinamide", "Untreated"))

# movement of fresh worms away from center
move.plot <- ggplot(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi"), aes(x = Group, y = 1 - Mf)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(color = "royalblue2", size = 1.25, alpha = 0.6, groupOnX = TRUE, cex = 2) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif", label.x = 1.5, label.y = 1.1) +
  scale_x_discrete(limits = c("Untreated", "Nicotinamide"), labels = c("Untreated", "NAM")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
  labs(title = "0 DPE", x = "", y = "Fraction leaving center") + 
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
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
move.plot

# save_plot(here("plots", "Fig4B_raw.pdf"), move.plot, base_width = 3)

FBS_mean <- mean(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5 & Group == "Untreated" & Treatment == "FBS")$CI3)
NAM_FBS_mean <- mean(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5 & Group == "Nicotinamide" & Treatment == "FBS")$CI3)
NaCl_mean <- mean(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5 & Group == "Untreated" & Treatment == "NaCl")$CI3)
NAM_NaCl_mean <- mean(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5 & Group == "Nicotinamide" & Treatment == "NaCl")$CI3)
MB_mean <- mean(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5 & Group == "Untreated" & Treatment == "3-methyl-1-butanol")$CI3)
NAM_MB_mean <- mean(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5 & Group == "Nicotinamide" & Treatment == "3-methyl-1-butanol")$CI3)
means <- tibble(
  Treatment = factor(x = rep(c("FBS", "NaCl", "3-methyl-1-butanol"), 2), levels = c("FBS", "NaCl", "3-methyl-1-butanol")),
  CI3 = c(FBS_mean, NaCl_mean, MB_mean, NAM_FBS_mean, NAM_NaCl_mean, NAM_MB_mean),
  Group = c(rep("Untreated", 3), rep("Nicotinamide", 3))
)

# chemotaxis index of fresh worms
fresh.ci.plot <- ggplot(dplyr::filter(tidy.data, Age == "Fresh" & Species == "Bpahangi" & T_n + C_n + O_n >= 5), aes(x = Group, y = CI3)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(color = "royalblue2", size = 2, alpha = 0.6, groupOnX = TRUE, cex = 3.5) +
  geom_text(data = means, mapping = aes(x = Group, y = 1.15, label = paste0("CI= ", round(CI3, digits = 2))), size = 3.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif", label.x = 1.5, label.y = 1.3) +
  scale_x_discrete(limits = c("Untreated", "Nicotinamide"), labels = c("Untreated", "NAM")) +
  scale_y_continuous(limits = c(-1, 1.4), breaks = c(-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1)) +
  labs(title = "0 DPE", x = "", y = "Chemotaxis Index (CI)") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
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
  facet_grid(. ~ Treatment) +
  NULL
fresh.ci.plot

save_plot(here("plots", "Fig4A_raw.pdf"), fresh.ci.plot, base_width = 8)

# movement of old worms away from center
old.move.plot <- ggplot(dplyr::filter(tidy.data, Age == "Old" & Species == "Bpahangi"), aes(x = Group, y = 1 - Mf)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(color = "darkgreen", size = 2, alpha = 0.6, groupOnX = TRUE, cex = 2.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
  labs(title = "1 DPE", x = "", y = "Fraction leaving center") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
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
old.move.plot

# save_plot(here("plots", "Fig4D_right.pdf"), old.move.plot, base_width = 2)

# chemotaxis index of old worms
FBS_mean <- mean(filter(tidy.data, Age == "Old" & Species == "Bpahangi" & Group == "Untreated" & Treatment == "FBS")$CI3, na.rm = TRUE)

old.ci.plot <- ggplot(dplyr::filter(tidy.data, Age == "Old" & Species == "Bpahangi" & Treatment == "FBS"), aes(x = Group, y = CI3)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(color = "darkgreen", size = 2, alpha = 0.6, groupOnX = TRUE, cex = 2.5) +
  annotate("text", x = 1, y = 1.15, label = paste0("Mean CI: ", round(FBS_mean, digits = 2)), size = 3.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey37", size = 0.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  scale_y_continuous(limits = c(-1, 1.2), breaks = c(-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1)) +
  labs(title = "1 DPE", x = "", y = "Chemotaxis Index (CI)") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
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
old.ci.plot

save_plot(here("plots", "Fig4D_left.pdf"), old.ci.plot, base_width = 2)
