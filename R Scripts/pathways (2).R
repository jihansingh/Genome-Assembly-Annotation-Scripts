library(data.table)
library(tidyverse)

# ==============================================================================
# 1. SETUP FILE PATHS & SPECIES LIST
# ==============================================================================
# Update these paths to point to your actual InterProScan TSV files
species_files <- list(
  "Ficus_craterostoma" = "ficus_craterostoma/interproscan.tsv",
  "Ficus_salicifolia"   = "ficus_salicifolia/interproscan.tsv",
  "Ficus_benjamina"     = "ficus_benjamina/interproscan.tsv"
)

# Lists to store output tables per species
reactome_by_species <- list()
cyc_by_species      <- list()
summary_list        <- list()

# ==============================================================================
# 2. PROCESS EACH SPECIES
# ==============================================================================
for (sp in names(species_files)) {
  
  file_path <- species_files[[sp]]
  cat("\n--------------------------------------------------\n")
  cat("Processing:", sp, "\n")
  cat("Loading file:", file_path, "\n")
  
  # Read InterProScan TSV fast with data.table (only reading Gene ID and Pathway col 15)
  # X1 = Sequence/Gene ID, X15 = Pathway Annotations
  tsv <- fread(
    file_path, 
    sep = "\t", 
    header = FALSE, 
    select = c(1, 15), 
    col.names = c("GeneID", "Pathways"),
    quote = ""
  )
  
  # Step A: Filter out empty/unannotated pathway rows
  tsv_pathways <- tsv[Pathways != "-" & !is.na(Pathways) & Pathways != ""]
  
  # Step B: Fast split of pipe-delimited strings
  pathway_tokens <- tsv_pathways[, .(Pathway = unlist(strsplit(Pathways, "\\|"))), by = .(GeneID)]
  pathway_tokens <- unique(pathway_tokens)
  
  # Step C: Filter Plant Reactome (R-ATH = Arabidopsis plant reference pathways)
  reactome <- pathway_tokens[like(Pathway, "R-ATH")]
  reactome[, `:=`(
    PathwayID = sub("^Reactome:", "", Pathway),
    Species   = sp
  )]
  
  # Step D: Filter PlantCyc / MetaCyc pathways
  cyc <- pathway_tokens[like(Pathway, "MetaCyc")]
  cyc[, `:=`(
    PathwayID = sub("^MetaCyc:", "", Pathway),
    Species   = sp
  )]
  
  # Save to species lists
  reactome_by_species[[sp]] <- reactome
  cyc_by_species[[sp]]      <- cyc
  
  # Log summary stats
  summary_list[[sp]] <- data.table(
    Species              = sp,
    Reactome_Pairs       = nrow(reactome),
    Reactome_Unique_Genes = uniqueN(reactome$GeneID),
    PlantCyc_Pairs       = nrow(cyc),
    PlantCyc_Unique_Genes = uniqueN(cyc$GeneID)
  )
  
  # Clear memory before next iteration
  rm(tsv, tsv_pathways, pathway_tokens)
  gc()
}

# ==============================================================================
# 3. COMBINE DATASETS & PRINT SUMMARY
# ==============================================================================
# Print overall summary breakdown table
summary_df <- rbindlist(summary_list)
cat("\n==================================================\n")
cat("          PLANT PATHWAY SUMMARY OVERVIEW          \n")
cat("==================================================\n")
print(summary_df)

# Combine all species into master data.tables (for multi-species comparative work)
all_reactome <- rbindlist(reactome_by_species)
all_cyc      <- rbindlist(cyc_by_species)

library(clusterProfiler)

# Prepare PlantCyc TERM2GENE for each species
t2g_craterostoma <- cyc_by_species[["Ficus_craterostoma"]][, .(PathwayID, GeneID)]
t2g_salicifolia  <- cyc_by_species[["Ficus_salicifolia"]][, .(PathwayID, GeneID)]
t2g_benjamina    <- cyc_by_species[["Ficus_benjamina"]][, .(PathwayID, GeneID)]

# Helper function to create clean Pathway Names from MetaCyc IDs
get_t2n <- function(t2g_table) {
  unique_terms <- unique(t2g_table$PathwayID)
  # Reformat MetaCyc IDs into readable labels
  data.frame(
    PathwayID = unique_terms,
    Name      = gsub("-", " ", unique_terms)
  )
}

t2n_craterostoma <- get_t2n(t2g_craterostoma)
t2n_benjamina <- get_t2n(t2g_benjamina)
t2n_salicifolia <- get_t2n(t2g_salicifolia)

orthogroups_tsv  <- "Orthofinder/Orthogroups/Orthogroups.tsv"
gene_counts_tsv  <- "Orthofinder/Orthogroups/Orthogroups.GeneCount.tsv"

ipr_fcrat        <- "ficus_craterostoma/interproscan.tsv"
ipr_fsal         <- "ficus_salicifolia/interproscan.tsv"
ipr_fbenj        <- "ficus_benjamina/interproscan.tsv"

outdir           <- "GO_enrichment_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Column names as they appear in Orthogroups.tsv — adjust to match your headers
col_fcrat        <- "ficus_craterostoma"
col_fsal         <- "ficus_salicifolia"
col_fbenj        <- "ficus_benjamina"

p_cutoff         <- 0.05   # after BH correction
min_gene_count   <- 1      # minimum genes annotated with a GO term to test it

# ── 1. Libraries ──────────────────────────────────────────────────────────────
library(tidyverse)

# ── 2. Load OrthoFinder tables and extract species-specific OG gene lists ─────
og_counts <- read_tsv(gene_counts_tsv, show_col_types = FALSE)
og_genes  <- read_tsv(orthogroups_tsv,  show_col_types = FALSE)

count_cols <- c(col_fcrat, col_fsal, col_fbenj)

og_sp <- og_counts %>%
  dplyr::select(Orthogroup, all_of(count_cols)) %>%
  mutate(
    sp_specific = case_when(
      .data[[col_fcrat]] > 0 & .data[[col_fsal]] == 0 & .data[[col_fbenj]] == 0 ~ "Fcraterostoma",
      .data[[col_fsal]]  > 0 & .data[[col_fcrat]] == 0 & .data[[col_fbenj]] == 0 ~ "Fsalicifolia",
      .data[[col_fbenj]] > 0 & .data[[col_fcrat]] == 0 & .data[[col_fsal]]  == 0 ~ "Fbenjamina",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sp_specific))

message("Species-specific orthogroup counts:")
print(table(og_sp$sp_specific))

# Pull gene IDs from species-specific OGs
extract_genes <- function(og_df, og_ids, sp_col) {
  og_df %>%
    filter(Orthogroup %in% og_ids) %>%
    pull(all_of(sp_col)) %>%
    na.omit() %>%
    str_split(", ") %>%
    unlist() %>%
    str_trim() %>%
    unique()
}

sp_genes <- list(
  Fcraterostoma = extract_genes(og_genes,
                                og_sp %>% filter(sp_specific == "Fcraterostoma") %>% pull(Orthogroup),
                                col_fcrat),
  Fsalicifolia  = extract_genes(og_genes,
                                og_sp %>% filter(sp_specific == "Fsalicifolia")  %>% pull(Orthogroup),
                                col_fsal),
  Fbenjamina    = extract_genes(og_genes,
                                og_sp %>% filter(sp_specific == "Fbenjamina")    %>% pull(Orthogroup),
                                col_fbenj)
)

message("\nGenes in species-specific OGs:")
print(sapply(sp_genes, length))

library(clusterProfiler)
library(data.table)

# --- 1. Align Species Names & Map to PlantCyc ---
# Note: sp_genes uses "Fcraterostoma", "Fsalicifolia", "Fbenjamina"
# Ensure target IDs match your t2g table gene ID format (protein_id vs gene_id)

target_fcrat <- sp_genes$Fcraterostoma
target_fsal  <- sp_genes$Fsalicifolia
target_fbenj <- sp_genes$Fbenjamina

# --- 2. Run PlantCyc Pathway Enrichment ---

# Ficus craterostoma
res_craterostoma <- enricher(
  gene          = target_fcrat,
  universe      = unique(t2g_craterostoma$GeneID), # Background universe from PlantCyc
  TERM2GENE     = t2g_craterostoma,
  TERM2NAME     = t2n_craterostoma,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# Ficus salicifolia
res_salicifolia <- enricher(
  gene          = target_fsal,
  universe      = unique(t2g_salicifolia$GeneID),
  TERM2GENE     = t2g_salicifolia,
  TERM2NAME     = t2n_salicifolia,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# Ficus benjamina
res_benjamina <- enricher(
  gene          = target_fbenj,
  universe      = unique(t2g_benjamina$GeneID),
  TERM2GENE     = t2g_benjamina,
  TERM2NAME     = t2n_benjamina,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# Load libraries
library(data.table)

# Ensure clean pathway descriptions are updated across all species
df_craterostoma[, Species := "Ficus craterostoma"]
df_salicifolia[, Species := "Ficus salicifolia"]
df_benjamina[, Species := "Ficus benjamina"]

# Combine all species results into a single comprehensive data.table
all_pathways_df <- rbindlist(list(df_craterostoma, df_salicifolia, df_benjamina), use.names = TRUE, fill = TRUE)

# Sort logically by species and p.adjust
setorder(all_pathways_df, Species, p.adjust)

# -----------------------------------------------------------------------------
# Option 1: Save as individual species TSV files
# -----------------------------------------------------------------------------
fwrite(df_craterostoma, file = "Ficus_craterostoma_PlantCyc_Enrichment.tsv", sep = "\t")
fwrite(df_salicifolia, file = "Ficus_salicifolia_PlantCyc_Enrichment.tsv", sep = "\t")
fwrite(df_benjamina, file = "Ficus_benjamina_PlantCyc_Enrichment.tsv", sep = "\t")

# -----------------------------------------------------------------------------
# Option 2: Save combined table with all 3 species in one TSV
# -----------------------------------------------------------------------------
fwrite(all_pathways_df, file = "All_Ficus_Species_PlantCyc_Enrichment.tsv", sep = "\t")

library(ggplot2)
library(data.table)
library(stringr)

# -----------------------------------------------------------------------------
# 1. Update Clean Pathway Names
# -----------------------------------------------------------------------------
pathway_names <- c(
  "PWY-8270" = "Cycloartenol biosynthesis",
  "PWY-101"  = "Photosynthesis light reactions",
  "PWY-7980" = "ATP synthesis / energy turnover",
  "PWY-5381" = "Pyridine nucleotide cycling (NAD salvage)",
  "PWY-7953" = "Gallic acid biosynthesis (tannins)",
  "PWY-6387" = "Monolignol biosynthesis (late steps)",
  "PWY-6386" = "Monolignol biosynthesis (early steps)",
  "PWY-6466" = "General phenylpropanoid biosynthesis",
  "PWY-7947" = "Proanthocyanidin (tannin) biosynthesis",
  "PWY-8187" = "Pheophorbide a oxygenase pathway",
  "PWY-3341" = "Chlorophyll a degradation",
  "PWY-5083" = "NAD/NADP biosynthesis (de novo)",
  "PWY-8148" = "Secondary metabolite glycosylation"
)

# Replace cryptic IDs with descriptive pathway names
plot_dt[, Clean_Name := ifelse(ID %in% names(pathway_names), pathway_names[ID], Description)]

# Wrap text cleanly for y-axis readability
plot_dt[, Clean_Name_wrapped := str_wrap(Clean_Name, width = 38)]

# Calculate -log10 p-adj
plot_dt[, log_padj := -log10(p.adjust)]

# Format p.adjust into clean scientific notation for display on bars
plot_dt[, p_label := sprintf("p.adj = %.2e", p.adjust)]

# Ensure species factor order matches previous figures
plot_dt[, Species := factor(Species, levels = c("F. craterostoma", "F. salicifolia", "F. benjamina"))]

# -----------------------------------------------------------------------------
# 2. Build Faceted Bar Chart with Direct p-value Labels
# -----------------------------------------------------------------------------
p_bar_clean <- ggplot(plot_dt, aes(x = reorder(Clean_Name_wrapped, log_padj), y = log_padj)) +
  geom_col(fill = "#2c3e50", color = "black", width = 0.65) +
  geom_text(aes(label = p_label), 
            hjust = -0.15, 
            size = 3.2, 
            color = "black", 
            fontface = "bold") +
  coord_flip() +
  facet_wrap(~ Species, scales = "free_y", ncol = 1) +
  # Expand upper axis boundary slightly to prevent text clipping
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(
    x = NULL,
    y = expression(-log[10] ~ "(" * p[adj] * ")"),
    title = "PlantCyc Pathway Enrichment Across Ficus Species",
    subtitle = "Top enriched metabolic pathways identified in species-specific orthogroups"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold.italic", size = 11),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.spacing = unit(1.0, "lines"),
    panel.grid.major.y = element_blank() # Clean horizontal grid lines
  )

print(p_bar_clean)

# -----------------------------------------------------------------------------
# 3. Save Plot
# -----------------------------------------------------------------------------
ggsave("Ficus_Species_PlantCyc_Barplot_With_Pvalues.pdf", plot = p_bar_clean, width = 8.5, height = 9, dpi = 300)
ggsave("Ficus_Species_PlantCyc_Barplot_With_Pvalues.png", plot = p_bar_clean, width = 8.5, height = 9, dpi = 300)
