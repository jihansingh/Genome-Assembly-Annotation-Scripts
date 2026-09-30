library(tidyverse)
library(clusterProfiler)
library(GO.db)

# ==========================================
# 1. PARSE & PROCESS RAW ORTHOFINDER MAP
# ==========================================
# Read and shape raw wide Orthogroups matrix into a tidy long format
raw_orthogroups <- read.delim("Orthofinder/Orthogroups/Orthogroups.tsv", sep = "\t")

og_mapping <- raw_orthogroups %>%
  pivot_longer(
    cols = -Orthogroup, 
    names_to = "Species", 
    values_to = "Gene_ID"
  ) %>%
  filter(Gene_ID != "", !is.na(Gene_ID)) %>%
  separate_rows(Gene_ID, sep = ", ")

# Standardize species formatting style
og_mapping <- og_mapping %>%
  mutate(Species = case_when(
    Species == "ficus_benjamina" ~ "ficus_benjamina",
    TRUE ~ Species
  ))

# ==========================================
# 2. LOAD INTERPROSCAN AND MERGE DATA
# ==========================================
ipr_cols <- c(
  "Gene_ID", "Sequence_MD5", "Sequence_Length", "Analysis", 
  "Signature_Accession", "Signature_Description", "Start", "Stop", 
  "Score", "Status", "Date", "InterPro_Accession", 
  "InterPro_Description", "GO_annotations", "Pathways"
)

ipr_raw <- read.delim(
  "ficus_benjamina/interProScan.tsv", 
  header = FALSE, 
  col.names = ipr_cols, 
  na.strings = c("", "NA", "-")
)

# Core combination join using many-to-many argument
combined_data <- inner_join(og_mapping, ipr_raw, by = "Gene_ID", relationship = "many-to-many")

# ==========================================
# 3. BUILD CUSTOM GO TERM DEFINITIONS (SALICIFOLIA)
# ==========================================
# Parse out pipe-separated GO strings specifically for F. salicifolia
library(dplyr)
library(tidyr)

go_term2gene <- combined_data %>%
  # 1. Filter for F. salicifia (per your comment)
  filter(Species == "ficus_benjamina") %>%
  
  # 2. Use explicit dplyr:: namespace to avoid the collision
  dplyr::select(Gene_ID, GO_annotations) %>% 
  drop_na(GO_annotations) %>%
  distinct() %>%
  
  # 3. Split the pipe-separated terms
  separate_rows(GO_annotations, sep = "\\|") %>%
  
  # 4. Use explicit dplyr:: namespace again
  dplyr::select(GO_annotations, Gene_ID) %>%
  distinct()

# Fetch formal descriptive names using GO.db
unique_go_ids <- unique(go_term2gene$GO_annotations)

# Ensure unique_go_ids is a clean, flat character vector with no NAs
clean_go_ids <- as.character(na.omit(unique_go_ids))

go_term2name <- AnnotationDbi::select(
  GO.db, 
  keys = clean_go_ids, 
  columns = c("TERM", "ONTOLOGY"), 
  keytype = "GOID"
) %>%
  dplyr::select(GOID, TERM)

# ==========================================
# 4. DEFINE REFERENCE UNIVERSE & TARGET LIST
# ==========================================
# Background universe: All annotated genes in F. salicifolia
benj_universe <- combined_data %>%
  filter(Species == "ficus_benjamina") %>%
  pull(Gene_ID) %>%
  unique()

# Calculate copy counts per orthogroup inside salicifolia to identify duplications
benj_copy_counts <- og_mapping %>%
  filter(Species == "ficus_benjamina") %>%
  group_by(Orthogroup) %>%
  tally(name = "Internal_Copies")

# Target query list: Genes belonging to orthogroups with expansions (> 1 copy)
benj_target_genes <- og_mapping %>%
  filter(Species == "ficus_benjamina") %>%
  left_join(benj_copy_counts, by = "Orthogroup") %>%
  filter(Internal_Copies > 1) %>%
  pull(Gene_ID) %>%
  unique()

# ==========================================
# 5. RUN HYPERGEOMETRIC ENRICHMENT
# ==========================================
benj_enrichment_results <- enricher(
  gene = benj_target_genes,
  universe = benj_universe,
  TERM2GENE = go_term2gene,
  TERM2NAME = go_term2name,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Export structured table of outcomes
benj_enrichment_table <- as.data.frame(benj_enrichment_results)
write.csv(benj_enrichment_table, "Benjamina_Duplicated_GO_Enrichment.csv", row.names = FALSE)
benj_enrichment_results@result <- benj_enrichment_results@result %>%
  mutate(Description = ifelse(ID == "GO:0006561", "proline biosynthetic process", Description))
# ==========================================
# 6. GENERATE THE MANUSCRIPT DOTPLOT
# ==========================================
dotplot(
  benj_enrichment_results, 
  showCategory = 15, 
  font.size = 10,
  title = "Enriched Gene Ontology Terms in Ficus benjamina Duplications"
) +
theme_bw(base_size = 11) +
theme(
  plot.title = element_text(face = "bold", size = 11),
  axis.text.y = element_text(color = "black"),
  # Added padding to right margin so legend text never gets cut off
  plot.margin = margin(t = 10, r = 25, b = 10, l = 10, unit = "pt")
)

# Save high-resolution publication vector layout
ggsave("Figure_Benjamina_GO_Enrichment_Dotplot.pdf", width = 8.5, height = 6.5)
