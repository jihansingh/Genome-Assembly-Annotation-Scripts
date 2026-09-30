library(tidyverse)
library(ggplot2)
library(tidytext)
library(GO.db) # Provides programmatic mapping from GO:XXXXXXX -> Term Name

# ==============================================================================
# PART A: LOAD ORTHOFINDER GENE MAPPING
# ==============================================================================
og_mapping <- read.delim("Orthofinder/Orthogroups/Orthogroups.tsv", sep = "\t", check.names = FALSE) %>%
  pivot_longer(
    cols = -Orthogroup, 
    names_to = "Species", 
    values_to = "Gene_ID"
  ) %>%
  separate_rows(Gene_ID, sep = ", ") %>% 
  filter(!is.na(Gene_ID) & Gene_ID != "")

# Compute genome-wide Venn classifications once for all species
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
    has_crat & !has_benj & !has_sali ~ "Craterostoma Unique",
    has_benj & !has_crat & !has_sali ~ "Benjamina Unique",
    TRUE ~ "Accessory (Shared with 1 other)"
  ))

# Map of GO Accessions to Terms globally using GO.db
go_term_map <- AnnotationDbi::select(GO.db, keys = keys(GO.db), columns = "TERM") %>%
  dplyr::rename(GO_ID = GOID, GO_Name = TERM)

# ==============================================================================
# HELPER FUNCTION: RUN GO ENRICHMENT PIPELINE PER SPECIES
# ==============================================================================
run_go_pipeline <- function(sp_name, venn_target, tsv_path, plot_color) {
  
  ipr_raw <- read.delim(tsv_path, header = FALSE, sep = "\t", fill = TRUE)
  
  # Target GO terms column (standard TSV col 14)
  go_col <- if (ncol(ipr_raw) >= 14) 14 else ncol(ipr_raw)
  
  ipr_clean <- ipr_raw %>%
    dplyr::select(Gene_ID = V1, GO_Raw = dplyr::all_of(go_col)) %>%
    dplyr::filter(!is.na(GO_Raw) & GO_Raw != "" & GO_Raw != "-") %>%
    # Expand multi-mapped GO terms
    tidyr::separate_rows(GO_Raw, sep = "\\|") %>%
    tidyr::separate_rows(GO_Raw, sep = ",") %>%
    # Isolate exact accession key: strips out "(InterPro)", "(PANTHER)", etc.
    dplyr::mutate(GO_ID = stringr::str_extract(GO_Raw, "GO:\\d+")) %>%
    dplyr::filter(!is.na(GO_ID)) %>%
    dplyr::select(Gene_ID, GO_ID) %>%
    dplyr::distinct() %>%
    # Join human-readable GO names from GO.db map
    dplyr::left_join(go_term_map, by = "GO_ID") %>%
    dplyr::mutate(
      GO_Name = dplyr::coalesce(GO_Name, GO_ID) # Fallback to GO_ID if name absent
    )
  
  # PART C: INTEGRATE DATASTREAMS
  combined_data <- dplyr::inner_join(og_mapping, ipr_clean, by = "Gene_ID", relationship = "many-to-many")
  
  global_profile <- combined_data %>%
    dplyr::filter(Species == sp_name) %>%
    dplyr::group_by(GO_ID, GO_Name) %>%
    dplyr::tally(name = "Total_Occurrences") %>%
    dplyr::arrange(desc(Total_Occurrences))
  
  # PART D: STRATIFICATION & DUPLICATION COUNTS
  stratified_data <- combined_data %>%
    dplyr::filter(Species == sp_name) %>%
    dplyr::left_join(dplyr::select(orthogroup_classes, Orthogroup, Venn_Category), by = "Orthogroup")
  
  copy_counts <- og_mapping %>%
    dplyr::filter(Species == sp_name) %>%
    dplyr::group_by(Orthogroup) %>%
    dplyr::tally(name = "Internal_Copies")
  
  duplication_biology <- stratified_data %>%
    dplyr::left_join(copy_counts, by = "Orthogroup") %>%
    dplyr::mutate(Duplication_Status = if_else(Internal_Copies > 1, "Duplicated", "Strict Single-Copy"))
  
  # PART E: FISHER'S EXACT ENRICHMENT TEST
  duplicated_terms <- duplication_biology %>%
    dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == venn_target) %>%
    dplyr::group_by(GO_ID, GO_Name) %>%
    dplyr::tally(name = "Dup_Count") %>%
    dplyr::arrange(desc(Dup_Count))
  
  total_sp_genes  <- length(unique(stratified_data$Gene_ID))
  total_dup_genes <- length(unique(duplication_biology %>% 
                                     dplyr::filter(Duplication_Status == "Duplicated" & Venn_Category == venn_target) %>% 
                                     dplyr::pull(Gene_ID)))
  
  enrichment_results <- global_profile %>%
    dplyr::left_join(duplicated_terms, by = c("GO_ID", "GO_Name")) %>%
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
        total_sp_genes - (BG_Count - Dup_Count)
      ), nrow = 2), alternative = "greater")$p.value
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::filter(Dup_Count > 0) %>%
    dplyr::arrange(FDR, p_value)
  
  # PART F: OUTPUT CSV SUMMARY TRACKS
  clean_sp <- gsub("ficus_", "", sp_name)
  write.csv(enrichment_results, paste0("Table_S4_", tools::toTitleCase(clean_sp), "_Enriched_Duplicated_GO_Terms.csv"), row.names = FALSE)
  
  # PART G: PLOT INDIVIDUAL SPECIES
  plot_data <- enrichment_results %>%
    dplyr::filter(FDR < 0.05) %>%
    dplyr::slice_head(n = 20) %>%
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
  
  if(nrow(plot_data) > 0) {
    p <- ggplot(plot_data, aes(x = reorder(GO_Name, Dup_Count), y = Dup_Count)) +
      geom_bar(stat = "identity", fill = plot_color, color = "black", lwd = 0.3, width = 0.75) +
      geom_text(aes(label = Label_Text), stat = "identity", hjust = -0.05, size = 2.4, color = "black") +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
      theme_bw(base_size = 10) +
      labs(
        title = "Lineage-Specific Duplicated GO Landscape",
        subtitle = paste0("Ficus ", clean_sp, " Enriched Unique GO Terms (* FDR < 0.05, ** < 0.01, *** < 0.001)"),
        x = NULL,
        y = "Total Duplicated Gene Count"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, face = "italic", color = "grey30"),
        axis.text.y = element_text(color = "black", size = 7.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    ggsave(paste0("Figure_", tools::toTitleCase(clean_sp), "_Lineage_Specific_Duplicated_GO.png"), p, width = 11, height = 5.5)
  }
  
  return(enrichment_results)
}

# ==============================================================================
# RUN PIPELINE FOR ALL SPECIES
# ==============================================================================
sali_enrichment_results <- run_go_pipeline("ficus_salicifolia", "Salicifolia Unique", "ficus_salicifolia/interproscan.tsv", "#E64A19")
crat_enrichment_results <- run_go_pipeline("ficus_craterostoma", "Craterostoma Unique", "ficus_craterostoma/interproscan.tsv", "#1E88E5")
benj_enrichment_results <- run_go_pipeline("ficus_benjamina", "Benjamina Unique", "ficus_benjamina/interproscan.tsv", "#4CAF50")

# ==============================================================================
# COMBINED VERTICAL FACETED PLOT WITH ACTUAL GO TERM NAMES
# ==============================================================================
sali_df <- sali_enrichment_results %>% mutate(Species = "Ficus salicifolia")
crat_df <- crat_enrichment_results %>% mutate(Species = "Ficus craterostoma")
benj_df <- benj_enrichment_results %>% mutate(Species = "Ficus benjamina")

all_species_go <- bind_rows(sali_df, crat_df, benj_df)

plot_data_faceted <- all_species_go %>%
  filter(!is.na(GO_Name) & GO_Name != "" & GO_Name != "-" & FDR < 0.05) %>%
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
    GO_Unique = paste0(GO_Name, "___", Species)
  )

ggplot(plot_data_faceted, aes(x = reorder(GO_Unique, Dup_Count), y = Dup_Count, fill = Species)) +
  geom_bar(stat = "identity", color = "black", lwd = 0.3, width = 0.75, show.legend = FALSE) +
  geom_text(aes(label = Label_Text), stat = "identity", hjust = -0.05, size = 2.7, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c(
    "Ficus salicifolia"   = "#E64A19",
    "Ficus craterostoma" = "#1E88E5",
    "Ficus benjamina"    = "#4CAF50"
  )) +
  scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
  facet_wrap(~ Species, scales = "free", ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.40))) +
  theme_bw(base_size = 11) +
  labs(
    title = "Comparative Lineage-Specific Duplicated GO Landscape",
    subtitle = "Enriched Gene Ontology Terms Across Ficus Lineages (* FDR < 0.05, ** < 0.01, *** < 0.001)",
    x = NULL,
    y = "Total Duplicated Gene Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, face = "italic", color = "grey30"),
    strip.text = element_text(face = "bold.italic", size = 10.5),
    strip.background = element_rect(fill = "grey92"),
    axis.text.y = element_text(color = "black", size = 8.5), 
    axis.text.x = element_text(color = "black", size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

ggsave("Figure_Combined_Ficus_Lineage_Specific_GO_Vertical.png", width = 12, height = 14)
