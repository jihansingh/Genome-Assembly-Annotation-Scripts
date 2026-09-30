# Install required packages if you haven't already:
# install.packages(c("tidyverse", "readr", "ggplot2"))

library(tidyverse)

# ----------------------------------------------------------------------
# 1. READ & PREPARE THE DATA
# ----------------------------------------------------------------------
# Replace the file paths below with your actual file locations
csv_files <- c(
  "Ficus craterostoma"     = "ficus_craterostoma/double_trouble/duplication_mode_summary.csv",
  "Ficus salicifolia" = "ficus_salicifolia/double_trouble/duplication_mode_summary.csv",
  "Ficus benjamina"    = "ficus_benjamina/double_trouble/duplication_mode_summary.csv"
)

# Read all 3 files and combine them automatically into one tidy data frame
df <- map_df(csv_files, ~ read_csv(.x, show_col_types = FALSE), .id = "Species")

# Inspect the column names to make sure they match expectations
# DoubleTrouble TSV files typically use column names like "dup_mode" and "count"
print(head(df))



library(tidyverse)
library(scales) # For comma formatting

# Define label mappings for DoubleTrouble duplication modes
mode_labels <- c(
  "DD"  = "Dispersed Duplication",
  "PD"  = "Proximal Duplication",
  "SD"  = "Segmental Duplication",
  "TD"  = "Tandem Duplication"
)

# Generate Plot
plot_dup_facet <- ggplot(df, aes(x = Species, y = Freq, fill = Species)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = comma(Freq)), 
    vjust = -0.5, 
    size = 3, 
    fontface = "bold"
  ) +
  facet_wrap(
    ~ Var1, 
    scales = "free_y", 
    labeller = labeller(Var1 = mode_labels)
  ) +
  scale_y_continuous(
    name = "Number of Gene Pairs",
    labels = comma, 
    expand = expansion(mult = c(0, 0.2))
  ) +
  scale_fill_manual(values = c(
    "Ficus craterostoma" = "#2b5c8f", 
    "Ficus salicifolia"  = "#d95f02", 
    "Ficus benjamina"    = "#2e8b57"
  )) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text = element_text(face = "italic")
  )

ggsave(
  filename = "Ficus_duplication_modes_faceted.png",
  plot = plot_dup_facet,
  width = 8.5,
  height = 6,
  dpi = 600,
  bg = "white"
)
ggsave(
  filename = "Ficus_duplication_modes_faceted.png",
  plot = plot_dup_facet,
  width = 8.5,
  height = 6,
  dpi = 600,
  bg = "white"
)
