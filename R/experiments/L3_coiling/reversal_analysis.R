library(tidyverse)
library(here)
library(conflicted)
library(cowplot)

conflicted::conflict_prefer("filter", "dplyr")

tidy_data <- readRDS(here("data", "coiling_reversal.data")) %>%
  rename(Treatment = Well) %>%
  separate(Treatment, c("Stage", "Temperature", "Step", "Serum", "Replicate"), "_") %>%
  filter(Serum == "FBS") %>%
  select(Date, Temperature, Step, Replicate, Normalized.Motility)
  
start_motility <- filter(tidy_data, Step == "A") %>%
  rename(A.Value = Normalized.Motility) %>%
  select(-Step, -Temperature)

tidy_data <- left_join(tidy_data, start_motility) %>%
  mutate(Control.Normalized = Normalized.Motility / A.Value) %>%
  filter(Step != "A")

final.plot <- ggplot(tidy_data, aes(x = Step, y = Control.Normalized * 100)) +
  geom_point() +
  geom_line(aes(group = interaction(Replicate, Date))) +
  scale_x_discrete(labels = c("Room Temp.", "37C°")) +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  theme(axis.text.x = element_text(face="bold", size=10, angle = 45, vjust = .5),
        axis.text.y = element_text(face="bold", size=9),
        axis.title.x  = element_text(angle=90, vjust=0.5, size=0), 
        axis.title.y  = element_text(face="bold", angle=90, size=12),
        strip.text.x = element_text(face="bold", size=12),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.75, colour = "black")) +
  labs(y = "Percent Change") +
  NULL
final.plot

save_plot(here("plots", "S4_Figure.pdf"), final.plot, base_height = 6, base_width = 5)
