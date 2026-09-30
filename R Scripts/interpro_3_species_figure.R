library(tidyverse)

# Build the master 3-species dataset
ficus_data <- tibble(
  Metric = rep(c("Annotated (≥1 Match)", "InterPro Entry (IPR)", "GO Mapped", "Pathways Mapped"), 3),
  Species = c(
    rep("F. craterostoma", 4),
    rep("F. benjamina", 4),
    rep("F. salicifolia", 4)
  ),
  Percentage = c(
    80.50, 61.62, 49.16, 50.03,  # F. craterostoma
    76.35, 59.43, 46.01, 47.61,  # F. benjamina
    70.28, 48.91, 38.42, 38.90   # F. salicifolia
  ),
  Count = c(
    30065, 23016, 18361, 18685,  # F. craterostoma
    27167, 21145, 16372, 16940,  # F. benjamina
    62514, 43506, 34174, 34603   # F. salicifolia
  )
)

# Order species and metrics for the plot
ficus_data <- ficus_data %>%
  mutate(
    Species = factor(Species, levels = c("F. craterostoma", "F. benjamina", "F. salicifolia")),
    Metric = factor(Metric, levels = c("Annotated (≥1 Match)", "InterPro Entry (IPR)", "GO Mapped", "Pathways Mapped"))
  )

# Plot percentage comparison across species
plot_3sp <- ggplot(ficus_data, aes(x = Metric, y = Percentage, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.2) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = position_dodge(width = 0.8),
    vjust = -0.5, size = 3.2, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "F. craterostoma" = "#2b5c8f",
    "F. benjamina"    = "#408080",
    "F. salicifolia"  = "#d95f02"
  )) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Comparative InterProScan Functional Annotation Across Ficus Species",
    x = NULL,
    y = "Percentage of Total Proteome (%)"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  )

print(plot_3sp)
ggsave("ficus_three_species_annotation.png", plot = plot_3sp, width = 9, height = 5.5, dpi = 300)