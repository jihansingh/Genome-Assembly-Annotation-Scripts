# ==============================================================================
# LINEAGE-SPECIFIC DUPLICATION GO ENRICHMENT PIPELINE
# ==============================================================================
# Purpose: Identify and plot GO terms for lineage-specific duplicated orthologs
# Target Species: Ficus salicifolia

# 1. DEPENDENCIES & INSTALLATION (Uncomment if needed)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("GO.db")

library(tidyverse)
library(clusterProfiler)
library(GO.db)

# ==============================================================================
# PART 1: REGENERATE THE ORTHOGROUP MAPPING & Wrangle Data
# ==============================================================================

# Read the raw wide Orthogroups file from OrthoFinder
raw_orthogroups <- read.delim("Orthofinder/Orthogroups/Orthogroups.tsv", sep = "\t")

# Reshape the wide matrix into a tidy, two-column mapping table
og_mapping <- raw_orthogroups %>%
  pivot_longer(
    cols = -Orthogroup, 
    names_to = "Species", 
    values_to = "Gene_ID"
  ) %>%
  filter(Gene_ID != "", !is.na(Gene_ID)) %>%
  separate_rows(Gene_ID, sep = ", ")

# Standardize species names to match your formatting style perfectly
og_mapping <- og_mapping %>%
  mutate(Species = case_when(
    Species == "ficus_salicifolia" ~ "ficus_salicifolia",
    # Add other species translations here if their column names look different:
    # str_detect(Species, "benjamina") ~ "ficus_benjamina",
    TRUE ~ Species
  ))

# ==============================================================================
# PART 2: LOAD INTERPROSCAN AND MERGE DATA
# ==============================================================================

# Define the standard 15 columns outputted by InterProScan
ipr_cols <- c(
  "Gene_ID", "Sequence_MD5", "Sequence_Length", "Analysis", 
  "Signature_Accession", "Signature_Description", "Start", "Stop", 
  "Score", "Status", "Date", "InterPro_Accession", 
  "InterPro_Description", "GO_annotations", "Pathways"
)

# Read raw InterProScan TSV file
ipr_raw <- read.delim(
  "ficus_salicifolia/interProScan.tsv", 
  header = FALSE, 
  col.names = ipr_cols, 
  na.strings = c("", "NA", "-")
)

# Create the final unified cross-referenced dataset
combined_data <- inner_join(og_mapping, ipr_raw, by = "Gene_ID", relationship = "many-to-many")

# ==============================================================================
# PART 3: BUILD CUSTOM GO TERM DEFINITIONS (SALICIFOLIA)
# ==============================================================================

go_term2gene <- combined_data %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::select(Gene_ID, GO_annotations) %>% 
  drop_na(GO_annotations) %>%
  distinct() %>%
  separate_rows(GO_annotations, sep = "\\|") %>%
  dplyr::select(GO_annotations, Gene_ID) %>%
  distinct()

# Extract a master list of all unique GO IDs present in your dataset
unique_go_ids <- unique(go_term2gene$GO_annotations)

# Query GO.db database to retrieve formal text labels and structural categories
go_term2name <- AnnotationDbi::select(
  GO.db, 
  keys = unique_go_ids, 
  columns = c("TERM", "ONTOLOGY"), 
  keytype = "GOID"
) %>%
  dplyr::select(GOID, TERM)

# ==============================================================================
# PART 4: STRATIFY BY VENN CATEGORY & DUPLICATION STATUS
# ==============================================================================

# Step 4.1: Rebuild orthogroup classes based on species presence matrix
orthogroup_classes <- og_mapping %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::summarize(
    has_benj = "ficus_benjamina" %in% Species,
    has_crat = "ficus_craterostoma" %in% Species,
    has_sali = "ficus_salicifolia" %in% Species
  ) %>%
  dplyr::mutate(Venn_Category = case_when(
    has_benj & has_crat & has_sali ~ "Core (Shared by all 3)",
    has_sali & !has_benj & !has_crat ~ "Salicifolia Unique",
    has_crat & (has_benj | has_sali) ~ "Accessory (Shared with 1 other)",
    TRUE ~ "Other"
  ))

# Step 4.2: Combine stratification information
stratified_data <- combined_data %>%
  dplyr::filter(Species == "ficus_salicifolia") %>%
  dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")

# Step 4.3: Calculate copy counts inside Ficus salicifolia to pinpoint duplications
sali_copy_counts <- og_mapping %>%
  dplyr::filter(Species == "ficus_salicifolia") %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::tally(name = "Internal_Copies")

sali_duplication_biology <- stratified_data %>%
  dplyr::left_join(sali_copy_counts, by = "Orthogroup") %>%
  dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))

# ==============================================================================
# PART 5: EXTRACT TARGET GENES (FIXED FOR LINEAGE-SPECIFIC DUPLICATES ONLY)
# ==============================================================================

# The universal background genes for F. salicifolia
sali_universe <- combined_data %>%
  dplyr::filter(Species == "ficus_salicifolia") %>%
  dplyr::pull(Gene_ID) %>%
  unique()

# TARGET SPECIFICATION: Only duplicated genes belonging exclusively to the F. salicifolia lineage
sali_target_genes <- sali_duplication_biology %>%
  dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Salicifolia Unique") %>%
  dplyr::pull(Gene_ID) %>%
  unique()

# ==============================================================================
# PART 6: HYPERGEOMETRIC ENRICHMENT TEST
# ==============================================================================

sali_enrichment_results <- enricher(
  gene = sali_target_genes,
  universe = sali_universe,
  TERM2GENE = go_term2gene,
  TERM2NAME = go_term2name,
  pvalueCutoff = 1,          # Keeps thresholds wide open to parse all available metrics
  qvalueCutoff = 1,          # Disables FDR truncation step
  pAdjustMethod = "none"     # Ranks primarily using raw p-value depth
)

sali_enrichment_table <- as.data.frame(sali_enrichment_results)
write.csv(sali_enrichment_table, "Salicifolia_Lineage_Specific_Duplicated_GO_Enrichment.csv", row.names = FALSE)

# ==============================================================================
# PART 7: GENERATE & SAVE THE REVISED DOTPLOT
# ==============================================================================

dotplot(
  sali_enrichment_results, 
  showCategory = 15, 
  font.size = 10,
  title = "GO Enrichment of Lineage-Specific Duplications in Ficus salicifolia"
) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(color = "black"),
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10, unit = "pt")
  )

ggsave("Figure_Salicifolia_Lineage_Specific_GO_Dotplot.pdf", width = 8.5, height = 6.5)
