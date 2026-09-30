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
# PART E: ISOLATE LINEAGE-SPECIFIC DUPLICATIONS & TEST ENRICHMENT
# ==============================================================================
# 1. Target conditional action: Lineage Specific AND Mechanistically Duplicated
sali_duplicated_domains <- sali_duplication_biology %>%
  dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Salicifolia Unique") %>%
  dplyr::group_by(Domain_ID, Domain_Name) %>%
  dplyr::tally(name = "Dup_Count") %>%
  dplyr::arrange(desc(Dup_Count))

# 2. Extract gene totals for Fisher's Exact 2x2 Contingency Table
total_sali_genes <- length(unique(stratified_data$Gene_ID))
total_dup_genes  <- length(unique(sali_duplication_biology %>% 
                                    dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Salicifolia Unique") %>% 
                                    dplyr::pull(Gene_ID)))

# 3. Perform Fisher's Exact Test against genome background
sali_enrichment_results <- sal_global_profile %>%
  dplyr::left_join(sali_duplicated_domains, by = c("Domain_ID", "Domain_Name")) %>%
  dplyr::mutate(
    Dup_Count = dplyr::coalesce(Dup_Count, 0L),
    BG_Count  = Total_Occurrences
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p_value = fisher.test(matrix(c(
      Dup_Count, 
      total_dup_genes - Dup_Count, 
      BG_Count - Dup_Count, 
      total_sali_genes - (BG_Count - Dup_Count)
    ), nrow = 2), alternative = "greater")$p.value
  ) %>%
  dplyr::ungroup() %>%
  # Correct for multiple testing across all domain tests
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::filter(Dup_Count > 0) %>%
  dplyr::arrange(FDR, p_value)

print("Statistically enriched lineage-specific duplicated domains:")
head(sali_enrichment_results, 20)

# Extract core background groups for comparative table metrics
sali_unique_functions <- stratified_data %>%
  dplyr::filter(Venn_Category == "Salicifolia Unique") %>%
  dplyr::group_by(Domain_Name) %>%
  dplyr::tally(name = "Gene_Count") %>%
  dplyr::arrange(desc(Gene_Count))

sal_core_functions <- stratified_data %>%
  dplyr::filter(Venn_Category == "Core (Shared by all 3)") %>%
  dplyr::group_by(Domain_Name) %>%
  dplyr::tally(name = "Gene_Count") %>%
  dplyr::arrange(desc(Gene_Count))

# ==============================================================================
# PART F: OUTPUT CSV SUMMARY TRACKS
# ==============================================================================
write.csv(sal_global_profile, "Table_S1_Salicifolia_Global_Domain_Profile.csv", row.names = FALSE)
write.csv(sali_unique_functions, "Table_S2_Salicifolia_Lineage_Specific_Domains.csv", row.names = FALSE)
write.csv(sal_core_functions, "Table_S3_Salicifolia_Shared_Core_Domains.csv", row.names = FALSE)
write.csv(sali_enrichment_results, "Table_S4_Salicifolia_Enriched_Duplicated_Domains.csv", row.names = FALSE)

# ==============================================================================
# PART G: GENERATE HIGH-RESOLUTION SUMMARY PLOT
# ==============================================================================
# Select the top enriched functional categories (FDR < 0.05) for plotting
plot_data_raw_flat <- sali_enrichment_results %>%
  dplyr::filter(!is.na(Domain_Name) & Domain_Name != "" & Domain_Name != "-" & FDR < 0.05) %>%
  dplyr::slice_head(n = 20) %>%
  # Construct a combined label: "Count (p = X.XXe-XX) ***"
  dplyr::mutate(
    FDR_Stars = case_when(
      FDR < 0.001 ~ " ***",
      FDR < 0.01  ~ " **",
      FDR < 0.05  ~ " *",
      TRUE        ~ ""
    ),
    Label_Text = paste0(
      scales::comma(Dup_Count), 
      " (p = ", formatC(p_value, format = "e", digits = 2), ")", 
      FDR_Stars
    )
  )

ggplot(plot_data_raw_flat, aes(x = reorder(Domain_Name, Dup_Count), y = Dup_Count)) +
  geom_bar(stat = "identity", fill = "#E64A19", color = "black", lwd = 0.3, width = 0.75) +
  
  # Inject count + formatted p-value + FDR stars to bar ends
  geom_text(aes(label = Label_Text), stat = "identity", hjust = -0.05, size = 2.4, color = "black") +
  
  coord_flip() +
  # Expanded multiplier (0.45) provides ample room for the longer text labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  theme_bw(base_size = 10) +
  labs(
    title = "Lineage-Specific Duplicated Domain Landscape",
    subtitle = "Ficus salicifolia Enriched Unique Genomic Alterations (* FDR < 0.05, ** < 0.01, *** < 0.001)",
    x = NULL,
    y = "Total Duplicated Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, face = "italic", color = "grey30"),
    axis.text.y = element_text(color = "black", size = 7.5), 
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_Salicifolia_Lineage_Specific_Duplicated_Domains.png", width = 11, height = 5.5)

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
    Species == "ficus_craterostoma" ~ "ficus_craterostoma",
    TRUE ~ Species
  ))

# ==============================================================================
# PART B: LOAD INTERPROSCAN FUNCTIONAL ANNOTATIONS
# ==============================================================================
ipr_raw <- read.delim("ficus_craterostoma/interproscan.tsv", header = FALSE, sep = "\t")

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
combined_data <- dplyr::inner_join(og_mapping, ipr_clean, by = "Gene_ID", relationship = "many-to-many")

# Generate global total baseline profile for validation work
crat_global_profile <- combined_data %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::group_by(Domain_ID, Domain_Name) %>%
  dplyr::tally(name = "Total_Occurrences") %>%
  dplyr::arrange(desc(Total_Occurrences))

print("Global functional database profile loaded:")
head(crat_global_profile, 10)

# ==============================================================================
# PART D: COMPUTE VENN CLASSIFICATIONS & DUPLICATIONS
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

# 2. Join environmental stratification back to base structural maps
stratified_data <- combined_data %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")

# 3. Calculate duplication status specifically inside F. craterostoma genomes
crat_copy_counts <- og_mapping %>%
  dplyr::filter(Species == "ficus_craterostoma") %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::tally(name = "Internal_Copies")

crat_duplication_biology <- stratified_data %>%
  dplyr::left_join(crat_copy_counts, by = "Orthogroup") %>%
  dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))

# ==============================================================================
# PART E: ISOLATE LINEAGE-SPECIFIC DUPLICATIONS & TEST ENRICHMENT
# ==============================================================================
# 1. Target conditional action: Lineage Specific AND Mechanistically Duplicated
crat_duplicated_domains <- crat_duplication_biology %>%
  dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Craterostoma Unique") %>%
  dplyr::group_by(Domain_ID, Domain_Name) %>%
  dplyr::tally(name = "Dup_Count") %>%
  dplyr::arrange(desc(Dup_Count))

# 2. Extract gene totals for Fisher's Exact 2x2 Contingency Table
total_crat_genes <- length(unique(stratified_data$Gene_ID))
total_dup_genes  <- length(unique(crat_duplication_biology %>% 
                                    dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Craterostoma Unique") %>% 
                                    dplyr::pull(Gene_ID)))

# 3. Perform Fisher's Exact Test against genome background
crat_enrichment_results <- crat_global_profile %>%
  dplyr::left_join(crat_duplicated_domains, by = c("Domain_ID", "Domain_Name")) %>%
  dplyr::mutate(
    Dup_Count = dplyr::coalesce(Dup_Count, 0L),
    BG_Count  = Total_Occurrences
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p_value = fisher.test(matrix(c(
      Dup_Count, 
      total_dup_genes - Dup_Count, 
      BG_Count - Dup_Count, 
      total_crat_genes - (BG_Count - Dup_Count)
    ), nrow = 2), alternative = "greater")$p.value
  ) %>%
  dplyr::ungroup() %>%
  # Correct for multiple testing across all domain tests
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::filter(Dup_Count > 0) %>%
  dplyr::arrange(FDR, p_value)

print("Statistically enriched lineage-specific duplicated domains:")
head(crat_enrichment_results, 20)

# Extract core background groups for comparative table metrics
crat_unique_functions <- stratified_data %>%
  dplyr::filter(Venn_Category == "Craterostoma Unique") %>%
  dplyr::group_by(Domain_Name) %>%
  dplyr::tally(name = "Gene_Count") %>%
  dplyr::arrange(desc(Gene_Count))

crat_core_functions <- stratified_data %>%
  dplyr::filter(Venn_Category == "Core (Shared by all 3)") %>%
  dplyr::group_by(Domain_Name) %>%
  dplyr::tally(name = "Gene_Count") %>%
  dplyr::arrange(desc(Gene_Count))

# ==============================================================================
# PART F: OUTPUT CSV SUMMARY TRACKS
# ==============================================================================
write.csv(crat_global_profile, "Table_S1_Craterostoma_Global_Domain_Profile.csv", row.names = FALSE)
write.csv(crat_unique_functions, "Table_S2_Craterostoma_Lineage_Specific_Domains.csv", row.names = FALSE)
write.csv(crat_core_functions, "Table_S3_Craterostoma_Shared_Core_Domains.csv", row.names = FALSE)
write.csv(crat_enrichment_results, "Table_S4_Craterostoma_Enriched_Duplicated_Domains.csv", row.names = FALSE)
# ==============================================================================
# PART G: GENERATE HIGH-RESOLUTION SUMMARY PLOT
# ==============================================================================
# Select the top enriched functional categories (FDR < 0.05) for plotting
plot_data_raw_flat <- crat_enrichment_results %>%
  dplyr::filter(!is.na(Domain_Name) & Domain_Name != "" & Domain_Name != "-" & FDR < 0.05) %>%
  dplyr::slice_head(n = 20) %>%
  # Construct a combined label: "Count (p = X.XXe-XX) ***"
  dplyr::mutate(
    FDR_Stars = case_when(
      FDR < 0.001 ~ " ***",
      FDR < 0.01  ~ " **",
      FDR < 0.05  ~ " *",
      TRUE        ~ ""
    ),
    Label_Text = paste0(
      scales::comma(Dup_Count), 
      " (p = ", formatC(p_value, format = "e", digits = 2), ")", 
      FDR_Stars
    )
  )

ggplot(plot_data_raw_flat, aes(x = reorder(Domain_Name, Dup_Count), y = Dup_Count)) +
  geom_bar(stat = "identity", fill = "#1E88E5", color = "black", lwd = 0.3, width = 0.75) +
  
  # Inject count + formatted p-value + FDR stars to bar ends
  geom_text(aes(label = Label_Text), stat = "identity", hjust = -0.05, size = 2.4, color = "black") +
  
  coord_flip() +
  # Expanded multiplier (0.45) provides ample room for the longer text labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  theme_bw(base_size = 10) +
  labs(
    title = "Lineage-Specific Duplicated Domain Landscape",
    subtitle = "Ficus craterostoma Enriched Unique Genomic Alterations (* FDR < 0.05, ** < 0.01, *** < 0.001)",
    x = NULL,
    y = "Total Duplicated Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, face = "italic", color = "grey30"),
    axis.text.y = element_text(color = "black", size = 7.5), 
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_Craterostoma_Lineage_Specific_Duplicated_Domains.png", width = 11, height = 5.5)

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
    Species == "ficus_benjamina" ~ "ficus_benjamina",
    TRUE ~ Species
  ))

# ==============================================================================
# PART B: LOAD INTERPROSCAN FUNCTIONAL ANNOTATIONS
# ==============================================================================
ipr_raw <- read.delim("ficus_benjamina/interproscan.tsv", header = FALSE, sep = "\t")

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
combined_data <- dplyr::inner_join(og_mapping, ipr_clean, by = "Gene_ID", relationship = "many-to-many")

# Generate global total baseline profile for validation work
benj_global_profile <- combined_data %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::group_by(Domain_ID, Domain_Name) %>%
  dplyr::tally(name = "Total_Occurrences") %>%
  dplyr::arrange(desc(Total_Occurrences))

print("Global functional database profile loaded:")
head(benj_global_profile, 10)

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
    has_benj & !has_crat & !has_sali ~ "Benjamina Unique",
    has_crat & (has_benj | has_sali) ~ "Accessory (Shared with 1 other)",
    TRUE ~ "Other"
  ))

# 2. Join environmental stratification back to base structural maps
stratified_data <- combined_data %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")

# 3. Calculate duplication status specifically inside F. benjamina genomes
benj_copy_counts <- og_mapping %>%
  dplyr::filter(Species == "ficus_benjamina") %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::tally(name = "Internal_Copies")

benj_duplication_biology <- stratified_data %>%
  dplyr::left_join(benj_copy_counts, by = "Orthogroup") %>%
  dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))

# ==============================================================================
# PART E: ISOLATE LINEAGE-SPECIFIC DUPLICATIONS & TEST ENRICHMENT
# ==============================================================================
# 1. Target conditional action: Lineage Specific AND Mechanistically Duplicated
benj_duplicated_domains <- benj_duplication_biology %>%
  dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Benjamina Unique") %>%
  dplyr::group_by(Domain_ID, Domain_Name) %>%
  dplyr::tally(name = "Dup_Count") %>%
  dplyr::arrange(desc(Dup_Count))

# 2. Extract gene totals for Fisher's Exact 2x2 Contingency Table
total_benj_genes <- length(unique(stratified_data$Gene_ID))
total_dup_genes  <- length(unique(benj_duplication_biology %>% 
                                    dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == "Benjamina Unique") %>% 
                                    dplyr::pull(Gene_ID)))

# 3. Perform Fisher's Exact Test against genome background
benj_enrichment_results <- benj_global_profile %>%
  dplyr::left_join(benj_duplicated_domains, by = c("Domain_ID", "Domain_Name")) %>%
  dplyr::mutate(
    Dup_Count = dplyr::coalesce(Dup_Count, 0L),
    BG_Count  = Total_Occurrences
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p_value = fisher.test(matrix(c(
      Dup_Count, 
      total_dup_genes - Dup_Count, 
      BG_Count - Dup_Count, 
      total_benj_genes - (BG_Count - Dup_Count)
    ), nrow = 2), alternative = "greater")$p.value
  ) %>%
  dplyr::ungroup() %>%
  # Correct for multiple testing across all domain tests
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::filter(Dup_Count > 0) %>%
  dplyr::arrange(FDR, p_value)

print("Statistically enriched lineage-specific duplicated domains:")
head(benj_enrichment_results, 20)

# Extract core background groups for comparative table metrics
benj_unique_functions <- stratified_data %>%
  dplyr::filter(Venn_Category == "Benjamina Unique") %>%
  dplyr::group_by(Domain_Name) %>%
  dplyr::tally(name = "Gene_Count") %>%
  dplyr::arrange(desc(Gene_Count))

benj_core_functions <- stratified_data %>%
  dplyr::filter(Venn_Category == "Core (Shared by all 3)") %>%
  dplyr::group_by(Domain_Name) %>%
  dplyr::tally(name = "Gene_Count") %>%
  dplyr::arrange(desc(Gene_Count))

# ==============================================================================
# PART F: OUTPUT CSV SUMMARY TRACKS
# ==============================================================================
write.csv(benj_global_profile, "Table_S1_Benjamina_Global_Domain_Profile.csv", row.names = FALSE)
write.csv(benj_unique_functions, "Table_S2_Benjamina_Lineage_Specific_Domains.csv", row.names = FALSE)
write.csv(benj_core_functions, "Table_S3_Benjamina_Shared_Core_Domains.csv", row.names = FALSE)
write.csv(benj_enrichment_results, "Table_S4_Benjamina_Enriched_Duplicated_Domains.csv", row.names = FALSE)

# ==============================================================================
# PART G: GENERATE HIGH-RESOLUTION SUMMARY PLOT
# ==============================================================================
# Select the top enriched functional categories (FDR < 0.05) for plotting
plot_data_raw_flat <- benj_enrichment_results %>%
  dplyr::filter(!is.na(Domain_Name) & Domain_Name != "" & Domain_Name != "-" & FDR < 0.05) %>%
  dplyr::slice_head(n = 20) %>%
  # Construct a combined label: "Count (p = X.XXe-XX) ***"
  dplyr::mutate(
    FDR_Stars = case_when(
      FDR < 0.001 ~ " ***",
      FDR < 0.01  ~ " **",
      FDR < 0.05  ~ " *",
      TRUE        ~ ""
    ),
    Label_Text = paste0(
      scales::comma(Dup_Count), 
      " (p = ", formatC(p_value, format = "e", digits = 2), ")", 
      FDR_Stars
    )
  )

ggplot(plot_data_raw_flat, aes(x = reorder(Domain_Name, Dup_Count), y = Dup_Count)) +
  geom_bar(stat = "identity", fill = "#4CAF50", color = "black", lwd = 0.3, width = 0.75) +
  
  # Inject count + formatted p-value + FDR stars to bar ends
  geom_text(aes(label = Label_Text), stat = "identity", hjust = -0.05, size = 2.4, color = "black") +
  
  coord_flip() +
  # Expanded multiplier (0.45) provides ample room for the longer text labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  theme_bw(base_size = 10) +
  labs(
    title = "Lineage-Specific Duplicated Domain Landscape",
    subtitle = "Ficus benjamina Enriched Unique Genomic Alterations (* FDR < 0.05, ** < 0.01, *** < 0.001)",
    x = NULL,
    y = "Total Duplicated Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, face = "italic", color = "grey30"),
    axis.text.y = element_text(color = "black", size = 7.5), 
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_Benjamina_Lineage_Specific_Duplicated_Domains.png", width = 11, height = 5.5)


sali_df <- sali_enrichment_results %>% mutate(Species = "Ficus salicifolia")
crat_df <- crat_enrichment_results %>% mutate(Species = "Ficus craterostoma")
benj_df <- benj_enrichment_results %>% mutate(Species = "Ficus benjamina")

library(ggplot2)
library(dplyr)
library(tidytext) # Contains reorder_within for clean per-facet ordering

# ==============================================================================
# VERTICALLY FACETED SUMMARY PLOT (1 COLUMN LAYOUT)
# ==============================================================================

# 1. Combine data from all three species
all_species_enrichment <- bind_rows(sali_df, crat_df, benj_df)

# 2. Filter top 15 enriched domains per species & format labels
plot_data_faceted <- all_species_enrichment %>%
  filter(!is.na(Domain_Name) & Domain_Name != "" & Domain_Name != "-" & FDR < 0.05) %>%
  group_by(Species) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  mutate(
    FDR_Stars = case_when(
      FDR < 0.001 ~ " ***",
      FDR < 0.01  ~ " **",
      FDR < 0.05  ~ " *",
      TRUE        ~ ""
    ),
    Label_Text = paste0(
      scales::comma(Dup_Count), 
      " (p = ", formatC(p_value, format = "e", digits = 1), ")", 
      FDR_Stars
    ),
    # Create a unique species-domain key for clean per-facet ordering
    Domain_Unique = paste0(Domain_Name, "___", Species)
  )

# 3. Create Vertical Facet Plot
ggplot(plot_data_faceted, aes(x = reorder(Domain_Unique, Dup_Count), y = Dup_Count, fill = Species)) +
  geom_bar(stat = "identity", color = "black", lwd = 0.3, width = 0.75, show.legend = FALSE) +
  
  # Text labels anchored cleanly past the bar end
  geom_text(aes(label = Label_Text), stat = "identity", hjust = -0.05, size = 2.7, color = "black") +
  
  coord_flip() +
  
  # Custom colors matching your individual plots
  scale_fill_manual(values = c(
    "Ficus salicifolia"   = "#E64A19",
    "Ficus craterostoma" = "#1E88E5",
    "Ficus benjamina"    = "#4CAF50"
  )) +
  
  # Strips the internal species suffix so the Y-axis labels stay crisp
  scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
  
  # Stack panels into 1 column
  facet_wrap(~ Species, scales = "free", ncol = 1) +
  
  # Expands horizontal room so labels won't clip off the right edge
  scale_y_continuous(expand = expansion(mult = c(0, 0.40))) +
  
  theme_bw(base_size = 11) +
  labs(
    title = "Comparative Lineage-Specific Duplicated Domain Landscape",
    subtitle = "Enriched Genomic Alterations Across Ficus Lineages (* FDR < 0.05, ** < 0.01, *** < 0.001)",
    x = NULL,
    y = "Total Duplicated Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, face = "italic", color = "grey30"),
    strip.text = element_text(face = "bold.italic", size = 10.5), # Facet Header formatting
    strip.background = element_rect(fill = "grey92"),
    axis.text.y = element_text(color = "black", size = 8.5), 
    axis.text.x = element_text(color = "black", size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines") # Space between stacked facets
  )

# Save as a tall manuscript figure
ggsave("Figure_Combined_Ficus_Lineage_Specific_Duplications_Vertical.png", width = 10, height = 14)
