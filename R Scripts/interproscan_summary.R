# ==============================================================================
# InterProScan Summary Statistics and Data Visualization Script
# ==============================================================================

# Install required packages if not already installed
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 1. User Settings & File Inputs
# ------------------------------------------------------------------------------
# Set paths and total predicted gene/protein count
tsv_file <- "ficus_craterostoma/interproscan.tsv"  # Path to your InterProScan TSV output
total_proteins <- 37349                  # Replace with your actual predicted protein count

# Standard InterProScan TSV column names
col_names <- c(
  "Protein_ID", "MD5", "Length", "Analysis", "Signature_Acc",
  "Signature_Desc", "Start", "Stop", "Score", "Status", "Date",
  "InterPro_Acc", "InterPro_Desc", "GO_Terms", "Pathways"
)

# ------------------------------------------------------------------------------
# 2. Load and Clean Data
# ------------------------------------------------------------------------------
# Load TSV (handling missing values in unmapped columns)
interpro_df <- read_tsv(
  file = tsv_file,
  col_names = col_names,
  show_col_types = FALSE,
  na = c("-", "", "NA")
)

# ------------------------------------------------------------------------------
# 3. Compute High-Level Metrics & Table
# ------------------------------------------------------------------------------
# Unique annotated proteins overall
annotated_count <- n_distinct(interpro_df$Protein_ID)
unannotated_count <- total_proteins - annotated_count

# Proteins with assigned InterPro IDs (IPRxxxxx)
ipr_count <- interpro_df %>%
  filter(!is.na(InterPro_Acc)) %>%
  summarise(n = n_distinct(Protein_ID)) %>%
  pull(n)

# Proteins with GO Terms
go_count <- interpro_df %>%
  filter(!is.na(GO_Terms)) %>%
  summarise(n = n_distinct(Protein_ID)) %>%
  pull(n)

# Proteins with Pathways (KEGG, MetaCyc, Reactome)
pathway_count <- interpro_df %>%
  filter(!is.na(Pathways)) %>%
  summarise(n = n_distinct(Protein_ID)) %>%
  pull(n)

# Signal Peptides and Transmembrane Helices
signalp_count <- interpro_df %>%
  filter(str_detect(Analysis, regex("SignalP", ignore_case = TRUE))) %>%
  summarise(n = n_distinct(Protein_ID)) %>%
  pull(n)

tmhmm_count <- interpro_df %>%
  filter(str_detect(Analysis, regex("TMHMM|Phobius", ignore_case = TRUE))) %>%
  summarise(n = n_distinct(Protein_ID)) %>%
  pull(n)

# Summary table output to console
cat("\n=== INTERPROSCAN ANNOTATION SUMMARY ===\n")
cat(sprintf("Total Proteins Analyzed:   %d\n", total_proteins))
cat(sprintf("Annotated (≥1 Match):      %d (%.2f%%)\n", annotated_count, (annotated_count / total_proteins) * 100))
cat(sprintf("Unannotated:               %d (%.2f%%)\n", unannotated_count, (unannotated_count / total_proteins) * 100))
cat(sprintf("InterPro Entry Assigned:   %d (%.2f%%)\n", ipr_count, (ipr_count / total_proteins) * 100))
cat(sprintf("GO Terms Mapped:           %d (%.2f%%)\n", go_count, (go_count / total_proteins) * 100))
cat(sprintf("Pathways Mapped:           %d (%.2f%%)\n", pathway_count, (pathway_count / total_proteins) * 100))
cat(sprintf("Signal Peptides (SignalP): %d (%.2f%%)\n", signalp_count, (signalp_count / total_proteins) * 100))
cat(sprintf("Transmembrane (TMHMM):    %d (%.2f%%)\n", tmhmm_count, (tmhmm_count / total_proteins) * 100))
cat("=========================================\n\n")

# ------------------------------------------------------------------------------
# Export Summary Metrics to a TSV File
# ------------------------------------------------------------------------------

# 1. Combine all key metrics into a clean table
export_table <- tibble(
  Metric = c(
    "Total Predicted Proteins",
    "Annotated (≥1 Match)",
    "Unannotated Fraction",
    "InterPro Entry Assigned (IPR)",
    "Gene Ontology (GO) Mapped",
    "Pathways Mapped (KEGG/Reactome/MetaCyc)",
    "Predicted Signal Peptides (SignalP)",
    "Predicted Transmembrane Helices (TMHMM/Phobius)"
  ),
  Protein_Count = c(
    total_proteins,
    annotated_count,
    unannotated_count,
    ipr_count,
    go_count,
    pathway_count,
    signalp_count,
    tmhmm_count
  )
) %>%
  mutate(
    Percentage = round((Protein_Count / total_proteins) * 100, 2)
  )

# 2. Write to TSV file
write_tsv(export_table, "interpro_summary_metrics.tsv")

cat("\nSummary metrics successfully saved to 'interpro_summary_metrics.tsv'!\n")

# ------------------------------------------------------------------------------
# 4. Figure 1: Member Database Breakdown Plot
# ------------------------------------------------------------------------------
# Count unique proteins per member database
db_counts <- interpro_df %>%
  group_by(Analysis) %>%
  summarise(Unique_Proteins = n_distinct(Protein_ID)) %>%
  mutate(Percentage = (Unique_Proteins / total_proteins) * 100) %>%
  arrange(desc(Unique_Proteins))

# Plot
plot_db <- ggplot(db_counts, aes(x = reorder(Analysis, Unique_Proteins), y = Unique_Proteins)) +
  geom_col(fill = "#2b5c8f", width = 0.7) +
  geom_text(
    aes(label = sprintf("%s (%.1f%%)", comma(Unique_Proteins), Percentage)),
    hjust = -0.1, size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.25))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Proteins Annotated Across InterPro Member Databases",
    x = "Member Database",
    y = "Number of Annotated Proteins"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

print(plot_db)
ggsave("interpro_db_breakdown.png", plot = plot_db, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Figure 2: Functional Annotation Summary
# ------------------------------------------------------------------------------
summary_df <- tibble(
  Category = c(
    "Overall Annotated", "InterPro Entry (IPR)", "GO Terms", 
    "Pathways", "Transmembrane", "Signal Peptide"
  ),
  Count = c(
    annotated_count, ipr_count, go_count, 
    pathway_count, tmhmm_count, signalp_count
  )
) %>%
  mutate(
    Percentage = (Count / total_proteins) * 100,
    Category = factor(Category, levels = rev(Category))
  )

plot_summary <- ggplot(summary_df, aes(x = Category, y = Count)) +
  geom_col(fill = "#408080", width = 0.6) +
  geom_text(
    aes(label = sprintf("%s (%.1f%%)", comma(Count), Percentage)),
    hjust = -0.1, size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.25))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Functional Annotation Overview",
    x = NULL,
    y = "Number of Proteins"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

print(plot_summary)
ggsave("interpro_annotation_summary.png", plot = plot_summary, width = 8, height = 5, dpi = 300)