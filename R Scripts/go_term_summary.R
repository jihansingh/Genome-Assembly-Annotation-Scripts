# ==============================================================================
# Pipeline: 3-Species Comparative GO Term Mapping & Visualization
# ==============================================================================

library(tidyverse)
library(GO.db)
library(AnnotationDbi)

# ------------------------------------------------------------------------------
# 1. Function to Parse InterProScan TSV and Calculate GO Category Rates
# ------------------------------------------------------------------------------
process_interpro <- function(tsv_path, species_name, total_genes) {
  message(paste("Processing data for:", species_name, "..."))
  
  tsv_data <- read_tsv(
    tsv_path, 
    col_names = FALSE, 
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  # 1. Filter rows with GO terms in Column 14
  go_rows <- tsv_data %>%
    filter(!is.na(X14) & X14 != "" & X14 != "-") %>%
    dplyr::select(Protein_ID = X1, GO_Terms = X14) %>%
    distinct()
  
  # 2. Separate pipe-delimited terms AND extract pure GO IDs (GO: followed by 7 digits)
  unnested_go <- go_rows %>%
    separate_rows(GO_Terms, sep = "\\|") %>%
    mutate(GO_Terms = str_extract(GO_Terms, "GO:\\d{7}")) %>% # <--- CLEAN EXTRACT HERE
    filter(!is.na(GO_Terms)) %>%
    distinct()
  
  # 3. Map GO IDs to Ontology Categories using AnnotationDbi::Ontology()
  unique_gos <- unique(unnested_go$GO_Terms)
  go_ontologies <- AnnotationDbi::Ontology(unique_gos)
  
  go_map <- data.frame(
    GO_Terms = names(go_ontologies),
    ONTOLOGY = as.character(go_ontologies),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(ONTOLOGY))
  
  # 4. Join back to proteins and summarize
  protein_ontology <- unnested_go %>%
    inner_join(go_map, by = "GO_Terms") %>%
    dplyr::select(Protein_ID, ONTOLOGY) %>%
    distinct()
  
  ontology_names <- c(
    "BP" = "Biological Process",
    "MF" = "Molecular Function",
    "CC" = "Cellular Component"
  )
  
  summary_stats <- protein_ontology %>%
    group_by(ONTOLOGY) %>%
    summarise(Unique_Proteins = n_distinct(Protein_ID), .groups = "drop") %>%
    mutate(
      Species = species_name,
      Category = ontology_names[ONTOLOGY],
      Total_Genes = total_genes,
      Percentage = (Unique_Proteins / total_genes) * 100
    ) %>%
    filter(!is.na(Category))
  
  return(summary_stats)
}
# ------------------------------------------------------------------------------
# 2. Process Data for All Three Species
# ------------------------------------------------------------------------------
# UPDATE THESE FILE PATHS TO YOUR ACTUAL TSV FILES!
file_paths <- list(
  "F. craterostoma" = "ficus_craterostoma/interproscan.tsv",
  "F. benjamina"    = "ficus_benjamina/interproscan.tsv",
  "F. salicifolia"   = "ficus_salicifolia/interproscan.tsv"
)

# Total predicted gene counts for percentage calculation
gene_counts <- list(
  "F. craterostoma" = 37349,
  "F. benjamina"    = 35303,
  "F. salicifolia"   = 88952
)

# Run loop over all 3 species and combine results into one dataframe
combined_go_data <- map2_dfr(
  file_paths, 
  names(file_paths), 
  ~ process_interpro(tsv_path = .x, species_name = .y, total_genes = gene_counts[[.y]])
)

# Refine Factor Levels for Plotting Order
combined_go_data <- combined_go_data %>%
  mutate(
    Species = factor(Species, levels = c("F. craterostoma", "F. benjamina", "F. salicifolia")),
    Category = factor(Category, levels = c("Biological Process", "Molecular Function", "Cellular Component"))
  )

# Print Summary Table to Console
print("Calculated Summary Table:")
print(combined_go_data %>% select(Species, Category, Unique_Proteins, Total_Genes, Percentage))

# ------------------------------------------------------------------------------
# 3. Generate Comparative Plot
# ------------------------------------------------------------------------------
ficus_colors <- c(
  "F. craterostoma" = "#2B5C8F",  # Dark Blue
  "F. benjamina"    = "#3A7D7C",  # Teal
  "F. salicifolia"   = "#D95F02"   # Orange
)

p <- ggplot(combined_go_data, aes(x = Category, y = Percentage, fill = Species)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = position_dodge(width = 0.8),
    vjust = -0.6,
    fontface = "bold",
    size = 3.6
  ) +
  scale_fill_manual(values = ficus_colors) +
  scale_y_continuous(
    limits = c(0, max(combined_go_data$Percentage, na.rm = TRUE) * 1.15),
    expand = c(0, 0)
  ) +
  labs(
    title = "Comparative GO Term Category Distribution Across Ficus Species",
    x = NULL,
    y = "Percentage of Total Proteome (%)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 15)),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 11, color = "black", margin = margin(t = 5)),
    axis.text.y  = element_text(size = 10, color = "#444444"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.4),
    legend.position = "top",
    legend.text = element_text(size = 10, face = "italic"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Render plot
print(p)
# Save as a standard CSV file
write_csv(combined_go_data, "Ficus_GO_Category_Summary.csv")
ggsave(
  filename = "Ficus_GO_Category_Distribution.png",
  plot = p,                  # Points to your ggplot object 'p'
  width = 8,                 # Width in inches
  height = 5.5,              # Height in inches
  dpi = 300                  # 300 DPI for high-resolution print standard
)
