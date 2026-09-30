library(tidyverse)

stats <- read.delim(
  "Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv",
  check.names = FALSE
)

names(stats)[1] <- "Metric"

stats_long <- stats %>%
  pivot_longer(
    cols = -Metric,               # Pivot everything EXCEPT the Metric column
    names_to = "Species",         # The column names (species) go here
    values_to = "Value"           # The actual numbers go here
  ) %>%
  mutate(Value = as.numeric(Value))

adaptation_metrics <- c(
  "Number of species-specific orthogroups",
  "Number of genes in species-specific orthogroups"
)

plot_data <- stats_long %>% 
  filter(Metric %in% adaptation_metrics)

# 3. PLOT: Create a clean, high-contrast bar chart with value labels
ggplot(plot_data, aes(x = Species, y = Value, fill = Metric)) +
  # We set a fixed dodge width (e.g., 0.9) so the bars and text align perfectly
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", lwd = 0.4) +
  
  # --- THE VALUE LABELS FIX ---
  geom_text(
    aes(label = scales::comma(Value)), # Formats numbers nicely (e.g., 10,000 instead of 10000)
    position = position_dodge(width = 0.9), # Matches the bar layout exactly
    vjust = -0.5,                          # Pushes the text slightly ABOVE the top of the bar
    size = 4,                              # Font size of the labels
    fontface = "bold"
  ) +
  
  # Using the clean Magenta and Blue from your trusted colorblind palette
  scale_fill_manual(
    values = c("#D81B60", "#1E88E5"),
    labels = c("Genes in Specific Groups", "Unique Orthogroups")
  ) + 
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Lineage-Specific Lineage Innovation in Ficus",
    subtitle = "Comparison of unique gene families and total specific genes",
    x = "Species",
    y = "Count",
    fill = "Evolutionary Metric"
  ) +
  
  # Expand the y-axis slightly so the labels at the very top don't get cut off
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical"
  )
