
library(tidyverse)

# Load predicted ages
UC_preds <- readRDS("output/UC_clock_predictions.rds")
WH_preds <- readRDS("output/WH_combined_ages.rds")

# Function to print MAE and correlation for either clock
print_accuracy <- function(df, age, pred) {
  clock.mae <- ie2misc::mae(df[[pred]], df[[age]])
  clock.corr <- cor.test(df[[pred]], df[[age]])$estimate
  cat('MAE: ', clock.mae, ' ; ', 'correlation: ', clock.corr)
}

# Assess accuracy of UC clocks
print_accuracy(df = UC_preds, pred = 'DNAmAgeClock2', age = 'Age')
print_accuracy(df = UC_preds, pred = 'DNAmAgeClock3', age = 'Age')

# Assess accuracy of western Hudson Bay clock
print_accuracy(df = WH_preds, pred = 'AgePredict', age = 'Age')

# Plot UC clock 3
ggplot(UC_preds, aes(x = Age, y = DNAmAgeClock3)) +
  geom_abline(intercept = 0, slope = 1, colour = "#E8E8E8", alpha = 0.5) +
  geom_point(aes(colour = Spec), size = 2) +
  scale_colour_manual(values = c("#B23A3A", "#6E7F8D")) +
  geom_smooth(method = 'lm', se = F, colour = "#E8E8E8", linewidth = 2, lineend = "round") +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), "cm"),
        panel.background = element_rect(fill = "#32404A"),
        plot.background = element_rect(fill = "#32404A", colour = NA),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "cm"),
        axis.text = element_text(colour = '#E8E8E8', size = 18),
        axis.title.x = element_text(colour = '#E8E8E8', size = 18, vjust = -5),
        axis.title.y = element_text(colour = '#E8E8E8', size = 18, vjust = 5),
        strip.background = element_rect(fill = "#32404A"),
        strip.text = element_text(colour = '#E8E8E8', size = 18),
        legend.position = "None") +
  labs(x = "Chronological age", y = "Epigenetic age")

# Save the plot
ggsave("figures/presentation/UC_clock3.svg", device = "svg", width = 15, height = 13, units = "cm")

# Plot the WH clock (axes same as UC clock)
ggplot(WH_preds, aes(x = Age, y = AgePredict)) +
  geom_abline(intercept = 0, slope = 1, colour = "#E8E8E8", alpha = 0.5) +
  geom_point(colour = "#E8E8E8", size = 2) +
  geom_smooth(method = 'lm', se = F, colour = "#E8E8E8", linewidth = 2, lineend = "round") +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), "cm"),
        panel.background = element_rect(fill = "#32404A"),
        plot.background = element_rect(fill = "#32404A", colour = NA),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "cm"),
        axis.text = element_text(colour = '#E8E8E8', size = 18),
        axis.title.x = element_text(colour = '#E8E8E8', size = 18, vjust = -5),
        axis.title.y = element_text(colour = '#E8E8E8', size = 18, vjust = 5),
        strip.background = element_rect(fill = "#32404A"),
        strip.text = element_text(colour = '#E8E8E8', size = 18),
        legend.position = "None") +
  ylim(0, 50) +
  labs(x = "Chronological age", y = "Epigenetic age")

# Save the plot
ggsave("figures/presentation/WH_clock.svg", device = "svg", width = 15, height = 13, units = "cm")
