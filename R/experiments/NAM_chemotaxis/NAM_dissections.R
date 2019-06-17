library(tidyverse)
library(conflicted)
library(magrittr)
library(ggbeeswarm)
library(Hmisc)
library(cowplot)
library(here)
library(ggpubr)
library(magick)
library(ggridges)

conflict_prefer("filter", "dplyr")
conflict_prefer("get_legend", "cowplot")

here()

# Import and tidy data ----------------------------------------------------

feeding.dr <- read_csv(here("data", "feeding_dose_response.csv")) %>%
  gather(Rep, Proportion.Fed, -NAM_mM) %>%
  rename(Dose = NAM_mM)

feeding.choice <- read_csv(here("data", "feeding_5_25.csv")) %>%
  rename(Proportion.Fed = Proportion_fed) %>%
  mutate(Group = factor(Group, levels = c("Control", "5mM_NAM", "25mM_NAM"), labels = c("None", "5 mM", "25 mM")))

feeding.size <- read_csv(here("data", "feeding_size.csv"))

tidy.data <- read_csv(here("data", "dissection_data.csv"), col_names = TRUE, ) %>%
  pivot_longer(cols = Head:Abdomen, names_to = "Tissue", values_to = "Count") %>%
  dplyr::filter(!is.na(Count)) %>%
  mutate(Tissue = fct_relevel(Tissue, "Abdomen", "Thorax", "Head")) %>%
  mutate(Date = factor(Date))

total <- group_by(tidy.data, Date, Treatment, Concentration, Number) %>%
  summarise(Total = sum(Count))

total.summary <- group_by(total, Treatment) %>%
  summarise(Mean = mean(Total), SEM = sd(Total) / sqrt(length(Total)))

tidy.data %<>% left_join(., total) %>%
  mutate(Proportion = Count / Total)

# remove outliers
quantile <- quantile(dplyr::filter(tidy.data, Concentration == "5 mM")$Total)
trimmed.data <- dplyr::filter(tidy.data, Total < quantile[4] + 2.5 * IQR(tidy.data$Total))

summary <- group_by(trimmed.data, Tissue, Treatment, Concentration) %>%
  summarise(Mean = mean(Count), SEM = sd(Count) / sqrt(length(Count)))

total.data <- select(trimmed.data, -Tissue, -Proportion, -Count) %>%
  unique()

# Plot data ----------------------------------------------------------------

# feeding dose reponse
dr.plot <- ggplot(dplyr::filter(feeding.dr, Dose != 0), aes(x = Dose, y = Proportion.Fed)) +
  stat_summary(fun.y = "mean", color = "black", geom = "point", shape = 18, alpha = 0.5, size = 5) +
  stat_summary(fun.y = "mean", color = "black", geom = "line", alpha = 0.5) +
  geom_vline(xintercept = 5, linetype = 2, color = "firebrick1", size = 0.5, alpha = 0.7) +
  geom_vline(xintercept = 25, linetype = 2, color = "firebrick4", size = 0.5, alpha = 0.7) +
  geom_point() + 
  scale_x_log10(breaks = c(10^-4, 10^-2, 1, 5, 25, 100), labels = c(0.0001, 0.01, 1, 5, 25, 100)) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", vjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)) +
  labs(x = "NAM (mM)", y = "Proportion Fed") +
  NULL
dr.plot

save_plot(here("plots", "Fig6A_raw.pdf"), dr.plot, base_width = 4)

# feeding of selected concentrations
comparisons <- list(c("None", "5 mM"), c("None", "25 mM"), c("5 mM", "25 mM"))

feeding.plot <- ggplot(feeding.choice, aes(x = Group, y = Proportion.Fed, color = Group)) +
  geom_boxplot() +
  geom_point(alpha = 0.7) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  scale_y_continuous(breaks = c(70, 80, 90, 100)) +
  scale_x_discrete(limits = c("None", "5 mM", "25 mM"), labels = c("0", "5", "25")) +
  scale_color_manual(values = c("black", "firebrick1", "firebrick4")) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", vjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)) +
  labs(x = "NAM (mM)", y = "Proportion Fed") +
  NULL
feeding.plot

save_plot(here("plots", "Fig6B_raw.pdf"), feeding.plot, base_width = 4)

# size of distended abdomens after eating
size.plot <- ggplot(feeding.size, aes(x = Width, y = Length, color = Group)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = c("firebrick4", "firebrick1", "black", "grey")) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", vjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)) +
  labs(x = "Abdomen Width (mm)", y = "Abdomen Lenth (mm)") +
  NULL
size.plot

save_plot(here("plots", "Fig6C_raw.pdf"), size.plot, base_width = 4)

# number of worms per mosquito
total.lm <- lm(Total ~ Concentration, data = total.data)
total.aov <- aov(total.lm)
summary(total.aov)
tukey <- TukeyHSD(total.aov)
tukey.df <- as.data.frame(tukey$Concentration) %>%
  slice(2:3) %>%
  mutate(Concentration = c("25 mM", "5 mM")) %>%
  mutate(p = `p adj`)

density.plot <- ggplot(total.data, aes(x = Total, y = Concentration)) +
  geom_density_ridges(aes(fill = Concentration, color = Concentration,
                          point_shape = Concentration, point_fill = Concentration, point_color = Concentration),
                      jittered_points = TRUE,
                      quantile_lines = TRUE,
                      alpha = 0.5, point_alpha = 0.7, point_size = 2, point_color = "black", scale = 0.95, size = 0.1) +
  geom_text(data = tukey.df, aes(y = Concentration), x = 25, label = "p < 0.005", vjust = -5) +
  scale_color_manual(values = c("firebrick4", "firebrick1", "black")) +
  scale_fill_manual(values = c("firebrick4", "firebrick1", "black")) +
  # stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  ggridges::scale_discrete_manual(aesthetics = "point_shape", values = c(21, 21, 21)) +
  ggridges::scale_discrete_manual(aesthetics = "point_fill", values = c("firebrick4", "firebrick1", "black")) +
  scale_x_continuous(expand = c(0.01, 0)) + 
  scale_y_discrete(expand = c(0.01, 0), breaks = c("None", "5 mM", "25 mM")) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", vjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)) +
  labs(x = "Total L3s per Mosquito", y = "Treatment") +
  NULL
density.plot

save_plot(here("plots", "Fig6D_raw.pdf"), density.plot, base_width = 8)

# mean worms per mosquito
bar.plot <- ggplot(ungroup(summary), aes(x = Concentration, y = Mean, fill = Tissue)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("Head" = "#CB891D", "Thorax" = "#4C7E14", "Abdomen" = "#407CC3")) +
  # scale_fill_viridis_d() +
  scale_x_discrete(limits = c("None", "5 mM", "25 mM"), labels = c("0", "5", "25")) +
  # scale_y_continuous(limits = c(0, 16), breaks = c(0, 5, 10, 15)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", vjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)) +
  labs(x = "NAM (mM)", y = "Mean L3s per Mosquito") +
  # facet_grid(. ~ Concentration) +
  NULL
bar.plot

save_plot(here("plots", "Fig6E_raw.pdf"), bar.plot, base_width = 4)

# proportion of worms in each tissue
trimmed.data %<>% mutate(Tissue = fct_relevel(Tissue, "Head", "Thorax", "Abdomen"))

dot.plot <- ggplot(trimmed.data, aes(x = Concentration, y = Proportion, color = Tissue)) +
  geom_violin(alpha = 0.6, color = "black", size = 0.75, fill = "grey90") +
  geom_beeswarm(size = 1, alpha = 0.6, groupOnX = TRUE, cex = 1.5) +
  scale_color_manual(values = c("Head" = "#CB891D", "Thorax" = "#4C7E14", "Abdomen" = "#407CC3")) +
  facet_grid(. ~ Tissue) +
  # stat_summary(fun.y = "mean", geom = "point", color = "red", shape = 18, alpha = 0.5, size = 5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1), color = "red", shape = 18, alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = list(c("None", "5 mM"), c("None", "25 mM")), method = "t.test", label = "p.signif", label.y = c(1.1, 1.2, 1.3)) +
  scale_x_discrete(limits = c("None", "5 mM", "25 mM"), labels = c("0", "5", "25")) +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.5, 1)) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", vjust = 0.5, size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 11),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.line = element_line(size = 0.75, colour = "black")
  ) +
  labs(x = "NAM (mM)", y = "Proportion of L3s in Tissue") +
  NULL
dot.plot

save_plot(here("plots", "Fig6F_raw.pdf"), dot.plot, base_width = 14)

# Supplemental ------------------------------------------------------------

# scatter plot with regression
facetRegression <- function(dat, xvar, yvar, group) {
  fml <- paste(yvar, "~", xvar)
  
  group <- rlang::sym(group)
  wrap_fml <- rlang::new_formula(rhs = group, lhs = NULL)
  dot <- rlang::quo(-!!group)
  
  dat %>%
    nest(!!dot) %>%
    mutate(
      model = map(data, ~ lm(fml, data = .x)),
      adj.r.squared = map_dbl(model, ~ signif(summary(.x)$adj.r.squared, 5)),
      intercept = map_dbl(model, ~ signif(.x$coef[[1]], 5)),
      slope = map_dbl(model, ~ signif(.x$coef[[2]], 5)),
      pvalue = map_dbl(model, ~ signif(summary(.x)$coef[2, 4], 5))
    ) %>%
    select(-data, -model) %>%
    left_join(dat)
}

temp <- mutate(trimmed.data, Proportion = Count / Total) %>%
  filter(Tissue == "Thorax")

lm.df <- facetRegression(temp, "Total", "Proportion", "Concentration") %>%
  select(Concentration:pvalue) %>%
  unique() %>%
  mutate(Concentration = factor(Concentration, levels = c("None", "5 mM", "25 mM")))

scatter.data <- dplyr::filter(trimmed.data, Tissue == "Thorax") %>%
  mutate(Concentration = factor(Concentration, levels = c("None", "5 mM", "25 mM")))

scatter.plot <- ggplot(scatter.data, aes(x = Total, y = Count / Total, color = Concentration, group = Concentration)) +
  geom_jitter(size = 2, height = 0, width = 1, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = lm.df, aes(label = paste0("R^2 = ", adj.r.squared)), x = 8.5, y = 1) +
  scale_color_manual(values = c("black", "firebrick1", "firebrick4")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.5, 1)) +
  facet_grid(Concentration ~ .) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", angle = 90, size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    strip.text.y = element_text(face = "bold", size = 12),
    axis.ticks = element_line(size = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.75, colour = "black"),
    legend.position = "None"
  ) +
  labs(y = "Proportion L3s in Thorax", x = "Total L3s per Mosquito") +
  NULL
scatter.plot

save_plot(here("plots", "S4_Figure.pdf"), scatter.plot, base_width = 5)
