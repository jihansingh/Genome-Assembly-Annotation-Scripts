# ==============================================================================
# LINEAGE-SPECIFIC DUPLICATION GO ENRICHMENT PIPELINE (FINAL CRATEROSTOMA RUN)
# ==============================================================================
# Target Species: Ficus craterostoma
# Features: Dynamic text-scrubbing for appended source tags in GO annotations

library(tidyverse)
library(clusterProfiler)
library(GO.db)

# ==============================================================================
# PART 1: REGENERATE THE ORTHOGROUP MAPPING 
# ==============================================================================
raw_orthogroups <- read.delim("Orthofinder/Orthogroups/Orthogroups.tsv", sep = "\t", stringsAsFactors = FALSE)

og_mapping <- raw_orthogroups %>%
  pivot_longer(
    cols = -Orthogroup, 
    names_to = "Species", 
    values_to = "Gene_ID"
  ) %>%
  filter(Gene_ID != "", !is.na(Gene_ID)) %>%
  separate_rows(Gene_ID, sep = ", ")

og_mapping <- og_mapping %>%
  mutate(Species = case_when(
    Species == "ficus_craterostoma" ~ "ficus_craterostoma",
    TRUE ~ Species
  ))

# ==============================================================================
# PART 2: LOAD INTERPROSCAN AND MERGE DATA
# ==============================================================================
ipr_cols <- c(
  "Gene_ID", "Sequence_MD5", "Sequence_Length", "Analysis", 
  "Signature_Accession", "Signature_Description", "Start", "Stop", 
  "Score", "Status", "Date", "InterPro_Accession", 
  "InterPro_Description", "GO_annotations", "Pathways"
)

ipr_raw <- read.delim(
  "ficus_craterostoma/interProScan.tsv", 
  header = FALSE, 
  col.names = ipr_cols, 
  na.strings = c("", "NA", "-"),
  stringsAsFactors = FALSE
)

combined_data <- inner_join(og_mapping, ipr_raw, by = "Gene_ID", relationship = "many-to-many")

# ==============================================================================
# PART 3: BUILD CUSTOM GO TERM DEFINITIONS (WITH TEXT CLEANING)
# ==============================================================================

# 1. Isolate rows that contain active data
go_annotated_rows <- combined_data %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::filter(!is.na(GO_annotations) & GO_annotations != "" & GO_annotations != "-")

# 2. Process, melt, and strip trailing parentheses text e.g., "GO:0005515(InterPro)"
go_term2gene <- go_annotated_rows %>%
  dplyr::select(Gene_ID, GO_annotations) %>% 
  distinct() %>%
  separate_rows(GO_annotations, sep = "\\|") %>%
  mutate(GO_annotations = trimws(GO_annotations)) %>%
  # Strip anything within parentheses and the parentheses themselves
  mutate(GO_annotations = gsub("\\(.*\\)", "", GO_annotations)) %>%
  mutate(GO_annotations = trimws(GO_annotations)) %>%
  # Keep only standard, clean GO identifiers
  dplyr::filter(grepl("^GO:[0-9]+$", GO_annotations)) %>%
  dplyr::select(GO_annotations, Gene_ID) %>%
  distinct()

# Extract and verify unique structural database keys
unique_go_ids <- unique(go_term2gene$GO_annotations)
message(paste("Successfully processed and loaded", length(unique_go_ids), "clean unique GO terms."))

# 3. Query GO.db database to fetch descriptive text definitions
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
orthogroup_classes <- og_mapping %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::summarize(
    has_benj = "ficus_benjamina" %in% Species,
    has_crat = "ficus_craterostoma" %in% Species,
    has_sali = "ficus_salicifolia" %in% Species
  ) %>%
  dplyr::mutate(Venn_Category = case_when(
    has_benj & has_crat & has_sali ~ "Core (Shared by all 3)",
    has_crat & !has_benj & !has_sali ~ "Craterostoma Unique",
    has_sali & (has_benj | has_crat) ~ "Accessory (Shared with 1 other)",
    TRUE ~ "Other"
  ))

stratified_data <- combined_data %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")

crat_copy_counts <- og_mapping %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::tally(name = "Internal_Copies")

crat_duplication_biology <- stratified_data %>%
  dplyr::left_join(crat_copy_counts, by = "Orthogroup") %>%
  dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))

# ==============================================================================
# PART 5: EXTRACT TARGET GENES (LINEAGE-SPECIFIC DUPLICATES ONLY)
# ==============================================================================
crat_universe <- combined_data %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::pull(Gene_ID) %>%
  unique()

# Targeted filtering profiles
crat_target_genes <- crat_duplication_biology %>%
  dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Craterostoma Unique") %>%
  dplyr::pull(Gene_ID) %>%
  unique()

message(paste("Analyzing", length(crat_target_genes), "lineage-specific duplicated target genes."))

# ==============================================================================
# PART 6: HYPERGEOMETRIC ENRICHMENT TEST
# ==============================================================================
crat_enrichment_results <- enricher(
  gene = crat_target_genes,
  universe = crat_universe,
  TERM2GENE = go_term2gene,
  TERM2NAME = go_term2name,
  pvalueCutoff = 1,      
  qvalueCutoff = 1,      
  pAdjustMethod = "none" 
)

crat_enrichment_table <- as.data.frame(crat_enrichment_results)
write.csv(crat_enrichment_table, "Craterostoma_Lineage_Specific_Duplicated_GO_Enrichment.csv", row.names = FALSE)

# ==============================================================================
# PART 7: GENERATE & SAVE THE REVISED DOTPLOT
# ==============================================================================
dotplot_fig <- dotplot(
  crat_enrichment_results, 
  showCategory = 15, 
  font.size = 10,
  title = "GO Enrichment of Lineage-Specific Duplications in Ficus craterostoma"
) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(color = "black"),
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10, unit = "pt")
  )

print(dotplot_fig)
ggsave("Figure_Craterostoma_Lineage_Specific_GO_Dotplot.pdf", plot = dotplot_fig, width = 8.5, height = 6.5)

# ==============================================================================
# LINEAGE-SPECIFIC DUPLICATION GO ENRICHMENT PIPELINE
# ==============================================================================
# Purpose: Identify and plot GO terms for lineage-specific duplicated orthologs
# Target Species: Ficus benjamina

library(tidyverse)
library(clusterProfiler)
library(GO.db)

# ==============================================================================
# PART 1: REGENERATE THE ORTHOGROUP MAPPING 
# ==============================================================================
raw_orthogroups <- read.delim("Orthofinder/Orthogroups/Orthogroups.tsv", sep = "\t", stringsAsFactors = FALSE)

og_mapping <- raw_orthogroups %>%
  pivot_longer(
    cols = -Orthogroup, 
    names_to = "Species", 
    values_to = "Gene_ID"
  ) %>%
  filter(Gene_ID != "", !is.na(Gene_ID)) %>%
  separate_rows(Gene_ID, sep = ", ")

og_mapping <- og_mapping %>%
  mutate(Species = case_when(
    Species == "ficus_benjamina" ~ "ficus_benjamina",
    TRUE ~ Species
  ))

# ==============================================================================
# PART 2: LOAD INTERPROSCAN AND MERGE DATA
# ==============================================================================
ipr_cols <- c(
  "Gene_ID", "Sequence_MD5", "Sequence_Length", "Analysis", 
  "Signature_Accession", "Signature_Description", "Start", "Stop", 
  "Score", "Status", "Date", "InterPro_Accession", 
  "InterPro_Description", "GO_annotations", "Pathways"
)

# Read raw InterProScan TSV file for Ficus benjamina
ipr_raw <- read.delim(
  "ficus_benjamina/interProScan.tsv", 
  header = FALSE, 
  col.names = ipr_cols, 
  na.strings = c("", "NA", "-"),
  stringsAsFactors = FALSE
)

combined_data <- inner_join(og_mapping, ipr_raw, by = "Gene_ID", relationship = "many-to-many")

# ==============================================================================
# PART 3: BUILD CUSTOM GO TERM DEFINITIONS (WITH TEXT CLEANING)
# ==============================================================================

# 1. Isolate rows that contain active data
go_annotated_rows <- combined_data %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::filter(!is.na(GO_annotations) & GO_annotations != "" & GO_annotations != "-")

# 2. Process, melt, and strip trailing parentheses text e.g., "GO:0005515(InterPro)"
go_term2gene <- go_annotated_rows %>%
  dplyr::select(Gene_ID, GO_annotations) %>% 
  distinct() %>%
  separate_rows(GO_annotations, sep = "\\|") %>%
  mutate(GO_annotations = trimws(GO_annotations)) %>%
  # Strip source tags in parentheses
  mutate(GO_annotations = gsub("\\(.*\\)", "", GO_annotations)) %>%
  mutate(GO_annotations = trimws(GO_annotations)) %>%
  # Filter for structural GO identification keys
  dplyr::filter(grepl("^GO:[0-9]+$", GO_annotations)) %>%
  dplyr::select(GO_annotations, Gene_ID) %>%
  distinct()

# Extract and verify unique structural database keys
unique_go_ids <- unique(go_term2gene$GO_annotations)
message(paste("Successfully processed and loaded", length(unique_go_ids), "clean unique GO terms."))

# 3. Query GO.db database to fetch descriptive text definitions
go_term2name <- AnnotationDbi::select(
  GO.db, 
  keys = unique_go_ids, 
  columns = c("TERM", "ONTOLOGY"), 
  keytype = "GOID"
) %>%
  dplyr::select(GOID, TERM)

# === THE FIX: Manually inject the missing proline biosynthesis term description ===
if ("GO:0006561" %in% go_term2name$GOID) {
  go_term2name <- go_term2name %>%
    mutate(TERM = if_else(GOID == "GO:0006561" & (is.na(TERM) | TERM == ""), 
                          "proline biosynthetic process", TERM))
} else {
  go_term2name <- rbind(go_term2name, data.frame(GOID = "GO:0006561", TERM = "proline biosynthetic process"))
}

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
    has_benj & !has_crat & !has_sali ~ "Benjamina Unique",
    has_crat & (has_benj | has_sali) ~ "Accessory (Shared with 1 other)",
    TRUE ~ "Other"
  ))

# Step 4.2: Combine stratification information
stratified_data <- combined_data %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")

# Step 4.3: Calculate copy counts inside Ficus benjamina to pinpoint duplications
benj_copy_counts <- og_mapping %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::tally(name = "Internal_Copies")

benj_duplication_biology <- stratified_data %>%
  dplyr::left_join(benj_copy_counts, by = "Orthogroup") %>%
  dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))

# ==============================================================================
# PART 5: EXTRACT TARGET GENES (LINEAGE-SPECIFIC DUPLICATES ONLY)
# ==============================================================================

# Universal genomic background framework (all mapped genes)
benj_universe <- combined_data %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::pull(Gene_ID) %>%
  unique()

# Targeted selection: Only duplicated genes belonging exclusively to the F. benjamina lineage
benj_target_genes <- benj_duplication_biology %>%
  dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Benjamina Unique") %>%
  dplyr::pull(Gene_ID) %>%
  unique()

message(paste("Analyzing", length(benj_target_genes), "lineage-specific duplicated target genes."))

# ==============================================================================
# PART 6: HYPERGEOMETRIC ENRICHMENT TEST
# ==============================================================================
benj_enrichment_results <- enricher(
  gene = benj_target_genes,
  universe = benj_universe,
  TERM2GENE = go_term2gene,
  TERM2NAME = go_term2name,
  pvalueCutoff = 1,      
  qvalueCutoff = 1,      
  pAdjustMethod = "none" 
)

benj_enrichment_table <- as.data.frame(benj_enrichment_results)
write.csv(benj_enrichment_table, "Benjamina_Lineage_Specific_Duplicated_GO_Enrichment.csv", row.names = FALSE)

# ==============================================================================
# PART 7: GENERATE & SAVE THE REVISED DOTPLOT
# ==============================================================================
dotplot_fig <- dotplot(
  benj_enrichment_results, 
  showCategory = 15, 
  font.size = 10,
  title = "GO Enrichment of Lineage-Specific Duplications in Ficus benjamina"
) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(color = "black"),
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10, unit = "pt")
  )

print(dotplot_fig)
ggsave("Figure_Benjamina_Lineage_Specific_GO_Dotplot.png", plot = dotplot_fig, width = 8.5, height = 6.5)

