
# 1 Load limmas ====

limma_bl_a <- readRDS('output/limma_blood_age.rds')
limma_sk_a <- readRDS('output/limma_skin_age.rds')
limma_t_s <- readRDS('output/limma_tissue_sex.rds')

# 2 Volcano plots ====

# Volcano plots

# EWAS for blood Cpgs significantly associated with age 
bl_age_ewas <- ggplot(limma_bl_a, aes(x = age_coeff, y = -log10(age_pval_adj), col = sig)) + 
  geom_point() + 
  scale_colour_brewer(palette = 'Greys', direction = 1) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 18),
        axis.text.x = element_text(colour = 'black', size = 18),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),
        legend.position = 'none') +
  xlab('Coefficient') + ylab('-log10(Adjusted p-value)')

# EWAS for skin Cpgs significantly associated with age 
sk_age_ewas <- ggplot(limma_sk_a, aes(x = age_coeff, y = -log10(age_pval_adj), col = sig)) + 
  geom_point() + 
  scale_colour_brewer(palette = 'Greys', direction = 1) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 18),
        axis.text.x = element_text(colour = 'black', size = 18),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),
        legend.position = 'none') +
  xlab('Coefficient') + ylab('-log10(Adjusted p-value)')

# EWAS for Cpgs not significantly associated with sex
sex_ewas <- ggplot(limma_t_s, aes(x = sexM_coeff, y = -log10(sexM_pval_adj), col = sig_sex)) +
  geom_point() +
  scale_colour_brewer(palette = 'Greys', direction = -1) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 18),
        axis.text.x = element_text(colour = 'black', size = 18),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5'),
        legend.position = 'none') +
  xlab('Coefficient') + ylab('-log10(Adjusted p-value)')

# Put panels together
panel_ewas <- plot_grid(bl_age_ewas, sk_age_ewas, sex_ewas, ncol = 3,
                        labels = c('A', 'B', 'C'), label_size = 18, align = 'hv')

# X and Y labels for ewas panel plot
Ylab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = '-log10(Adjusted p-value)', size = 7,angle = 90) + theme_void()
Xlab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Coefficient', size = 7, hjust = 0.1) + theme_void()

# Add axis labels
panel_ewas_with_y <- plot_grid(Ylab, panel_ewas, rel_widths = c(0.1, 1))
panel_ewas_with_x_y <- plot_grid(panel_ewas_with_y, Xlab, rel_heights = c(1, 0.05), ncol = 1)

# 3 Save plot ====

ggsave('ewas_panel.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/supplementary/', dpi = 300, width = 30, height = 12, units = 'cm', bg = 'white')

