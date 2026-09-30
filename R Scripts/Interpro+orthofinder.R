library(tidyverse)

# ==============================================================================
# PART A: LOAD ORTHOFINDER GENE MAPPING
# ==============================================================================
# Read wide orthogroup configuration matrix and pivot to clean layout
og_mapping <- read.delim("Orthofinder/Orthogroups/Orthogroups.tsv", sep = "\t", check.names = FALSE) %>%
  pivot_longer(
    cols = -Orthogroup, 
    names_to = "Species", 
    values_to = "Gene_ID"
  ) %>%
  # Separate multiple genes in a single cell if an orthogroup has duplicates
  separate_rows(Gene_ID, sep = ", ") %>% 
  filter(!is.na(Gene_ID) & Gene_ID != "")

# Standardize species name attributes to guarantee filter matches
og_mapping <- og_mapping %>%
  mutate(Species = case_when(
    Species == "ficus_salicifolia" ~ "ficus_salicifolia",
    TRUE ~ Species
  ))

# ==============================================================================
# PART B: LOAD INTERPROSCAN FUNCTIONAL ANNOTATIONS
# ==============================================================================
ipr_raw <- read.delim("ficus_salicifolia/interproscan.tsv", header = FALSE, sep = "\t")

ipr_clean <- ipr_raw %>%
  dplyr::select(
    Gene_ID = V1,          # Column 1: Match Key
    Database = V4,         # Column 4: Source Database
    Domain_ID = V5,        # Column 5: Unique Identifier Key
    Domain_Name = V6       # Column 6: Functional Structural Name Text
  ) %>%
  dplyr::distinct()

# ==============================================================================
# PART C: INTEGRATE GENOMIC DATASTREAMS
# ==============================================================================
combined_data <- inner_join(og_mapping, ipr_clean, by = "Gene_ID", relationship = "many-to-many")

# Generate global total baseline profile for validation work
sal_global_profile <- combined_data %>%
  filter(Species == "ficus_salicifolia") %>%
  group_by(Domain_ID, Domain_Name) %>%
  tally(name = "Total_Occurrences") %>%
  arrange(desc(Total_Occurrences))

print("Global functional database profile loaded:")
head(sal_global_profile, 10)

# ==============================================================================
# PART D: COMPUTE VENN CLASSIFICATIONS & DUPLICATIONS
# ==============================================================================
# 1. Isolate strict lineage specificity rules
# ==============================================================================
# PART D: COMPUTE VENN CLASSIFICATIONS & DUPLICATIONS
# ==============================================================================
# 1. Isolate strict lineage specificity rules
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

# 2. Join environmental stratification back to base structural maps (FIXED WITH NAMESPACE)
stratified_data <- combined_data %>%
  dplyr::filter(Species == "ficus_salicifolia") %>%
  dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")

# 3. Calculate duplication status specifically inside F. salicifolia genomes
sali_copy_counts <- og_mapping %>%
  dplyr::filter(Species == "ficus_salicifolia") %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::tally(name = "Internal_Copies")

sali_duplication_biology <- stratified_data %>%
  dplyr::left_join(sali_copy_counts, by = "Orthogroup") %>%
  dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))

# ==============================================================================
# PART E: ISOLATE LINEAGE-SPECIFIC DUPLICATIONS (THE TARGET PHENOTYPE)
# ==============================================================================
# TARGET CONDITIONAL ACTION: Lineage Specific AND Mechanistically Duplicated
sali_duplicated_domains <- sali_duplication_biology %>%
  filter(Duplication_Status == "Duplicated" & Venn_Category == "Salicifolia Unique") %>%
  group_by(Domain_Name) %>%
  tally(name = "Count") %>%
  arrange(desc(Count))

print("Lineage-specific duplicated domains isolated:")
head(sali_duplicated_domains, 20)

# Extract core background groups for comparative table metrics
sali_unique_functions <- stratified_data %>%
  filter(Venn_Category == "Salicifolia Unique") %>%
  group_by(Domain_Name) %>%
  tally(name = "Gene_Count") %>%
  arrange(desc(Gene_Count))

sal_core_functions <- stratified_data %>%
  filter(Venn_Category == "Core (Shared by all 3)") %>%
  group_by(Domain_Name) %>%
  tally(name = "Gene_Count") %>%
  arrange(desc(Gene_Count))

# ==============================================================================
# PART F: OUTPUT CSV SUMMARY TRACKS
# ==============================================================================
write.csv(sal_global_profile, "Table_S1_Salicifolia_Global_Domain_Profile.csv", row.names = FALSE)
write.csv(sali_unique_functions, "Table_S2_Salicifolia_Lineage_Specific_Domains.csv", row.names = FALSE)
write.csv(sal_core_functions, "Table_S3_Salicifolia_Shared_Core_Domains.csv", row.names = FALSE)
write.csv(sali_duplicated_domains, "Table_S4_Salicifolia_Lineage_Specific_Duplicated_Domains.csv", row.names = FALSE)

# ==============================================================================
# PART G: GENERATE HIGH-RESOLUTION SUMMARY PLOT
# ==============================================================================
# Select the top enriched functional categories for plotting
plot_data_raw_flat <- sali_duplicated_domains %>%
  filter(!is.na(Domain_Name) & Domain_Name != "") %>%
  head(20)

ggplot(plot_data_raw_flat, aes(x = reorder(Domain_Name, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#E64A19", color = "black", lwd = 0.3, width = 0.75) +
  
  # Inject explicit gene sequence value integers to bar ends
  geom_text(aes(label = scales::comma(Count)), stat = "identity", hjust = -0.15, size = 2.5, color = "black") +
  
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw(base_size = 10) +
  labs(
    title = "Lineage-Specific Duplicated Domain Landscape",
    subtitle = "Ficus salicifolia Unique Genomic Alterations",
    x = NULL,
    y = "Total Duplicated Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, face = "italic", color = "grey30"),
    axis.text.y = element_text(color = "black", size = 7.5), 
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_Salicifolia_Lineage_Specific_Duplicated_Domains.pdf", width = 10, height = 5.5)
