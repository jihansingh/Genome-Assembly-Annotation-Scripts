# GO Enrichment for Species-Specific Orthogroups
# Source: GO terms from InterProScan TSV (col 14), no topGO
# Test:   Fisher's exact test per GO term, BH correction
# ─────────────────────────────────────────────────────────────────────────────
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
# ── 0. Paths ──────────────────────────────────────────────────────────────────
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

# ── 3. Parse InterProScan TSV → long gene-GO table ───────────────────────────
# Col 1  = protein/gene ID
# Col 14 = GO terms, pipe-separated (e.g. GO:0003674|GO:0008150), or "-"
parse_ipr_go <- function(ipr_path, sp_label) {
  message("  Parsing: ", ipr_path)
  raw <- read_tsv(
    ipr_path,
    col_names     = FALSE,
    show_col_types = FALSE,
    comment       = "#"
  )
  if (ncol(raw) < 14) stop(paste("Fewer than 14 columns in", ipr_path,
                                 "-- was --goterms passed to InterProScan?"))
  raw %>%
    select(gene_id = X1, go_raw = X14) %>%
    distinct() %>%
    filter(!is.na(go_raw), go_raw != "-") %>%
    mutate(
      gene_id = str_trim(gene_id),
      go_list = str_split(go_raw, "\\|")
    ) %>%
    select(gene_id, go_list) %>%
    unnest(go_list) %>%
    rename(go_term = go_list) %>%
    filter(str_detect(go_term, "^GO:\\d{7}$")) %>%
    distinct() %>%
    mutate(species = sp_label)
}

parse_ipr_go <- function(path, species_label) {
  message("  Parsing: ", path)
  
  df <- data.table::fread(
    path,
    sep        = "\t",
    header     = FALSE,
    fill       = TRUE,
    quote      = "",
    col.names  = c("protein_id", "md5", "length", "analysis",
                   "sig_acc", "sig_desc", "start", "stop",
                   "score", "status", "date", "ipr_acc",
                   "ipr_desc", "go_terms", "pathways"),
    colClasses = "character",
    na.strings = c("", "-", "NA")
  )
  
  df <- df[!is.na(df$go_terms), ]
  df <- tidyr::separate_rows(df, go_terms, sep = "\\|")
  df$go_terms <- sub("\\(.*\\)$", "", df$go_terms)
  df <- df[grepl("^GO:\\d{7}$", df$go_terms), ]
  
  # Derive gene_id by stripping -mRNA-N suffix
  df$gene_id <- sub("-mRNA-\\d+$", "", df$protein_id)
  
  df$species <- species_label
  df[, c("gene_id", "protein_id", "go_terms", "species")]
}

ipr_all <- bind_rows(
  parse_ipr_go(ipr_fcrat, "Fcraterostoma"),
  parse_ipr_go(ipr_fsal,  "Fsalicifolia"),
  parse_ipr_go(ipr_fbenj, "Fbenjamina")
)

# Check gene_id format matches sp_genes
head(ipr_all$gene_id)
head(sp_genes$Fcraterostoma)

# Should now be > 0
length(intersect(ipr_all$gene_id[ipr_all$species == "Fcraterostoma"],
                 sp_genes$Fcraterostoma))

message("Total unique gene-GO pairs across all species: ", nrow(ipr_all))

# ── 4. Fisher's exact test enrichment function ────────────────────────────────
# Background = all annotated genes for that species (from InterProScan)
# Foreground = genes in species-specific OGs that have any GO annotation
#
# For each GO term, 2x2 contingency table:
#
#                  | In term | Not in term |
# Foreground genes |    a    |      b      |
# Background genes |    c    |      d      |

run_fisher_enrichment <- function(sp_label, ipr_df, interesting_genes) {
  
  sp_ipr <- ipr_df %>% filter(species == sp_label)
  
  # Universe: all genes with at least one GO annotation for this species
  all_annotated <- unique(sp_ipr$gene_id)
  N <- length(all_annotated)
  
  # Foreground: interesting genes that are also annotated
  fg_genes <- intersect(interesting_genes, all_annotated)
  n_fg <- length(fg_genes)
  
  message("\n  [", sp_label, "] ",
          n_fg, " / ", length(interesting_genes),
          " species-specific OG genes have GO annotations",
          " (background N = ", N, ")")
  
  if (n_fg == 0) {
    message("  No annotated foreground genes — skipping.")
    return(NULL)
  }
  
  # Summarise GO term counts in foreground and background
  go_fg <- sp_ipr %>%
    filter(gene_id %in% fg_genes) %>%
    count(go_term, name = "fg_count")
  
  go_bg <- sp_ipr %>%
    count(go_term, name = "bg_count")
  
  go_tbl <- go_bg %>%
    left_join(go_fg, by = "go_term") %>%
    mutate(fg_count = replace_na(fg_count, 0)) %>%
    filter(bg_count >= min_gene_count)    # drop very rare terms
  
  message("  Testing ", nrow(go_tbl), " GO terms with >= ", min_gene_count, " annotated genes")
  
  # Run Fisher's exact test row-wise
  run_fisher <- function(fg_in, bg_in, n_fg, N) {
    bg_out <- N - bg_in
    fg_out <- n_fg - fg_in
    mat <- matrix(c(fg_in, fg_out, bg_in - fg_in, bg_out - fg_out),
                  nrow = 2,
                  dimnames = list(c("fg", "bg"), c("in_term", "not_in_term")))
    # Clamp negatives from rounding edge cases
    mat[mat < 0] <- 0
    fisher.test(mat, alternative = "greater")$p.value
  }
  
  go_tbl <- go_tbl %>%
    mutate(
      p_value    = mapply(run_fisher, fg_count, bg_count,
                          MoreArgs = list(n_fg = n_fg, N = N)),
      p_adjusted = p.adjust(p_value, method = "BH"),
      enrichment = (fg_count / n_fg) / (bg_count / N),
      species    = sp_label
    ) %>%
    filter(p_adjusted < p_cutoff) %>%
    arrange(p_adjusted)
  
  message("  Significant GO terms (BH < ", p_cutoff, "): ", nrow(go_tbl))
  return(go_tbl)
}
run_fisher_enrichment <- function(sp_label, ipr_df, interesting_genes) {
  
  sp_ipr <- ipr_df %>% filter(species == sp_label)
  
  all_annotated <- unique(sp_ipr$protein_id)
  N <- length(all_annotated)
  
  fg_genes <- intersect(interesting_genes, all_annotated)
  n_fg <- length(fg_genes)
  
  message("\n  [", sp_label, "] ",
          n_fg, " / ", length(interesting_genes),
          " species-specific OG genes have GO annotations",
          " (background N = ", N, ")")
  
  if (n_fg == 0) {
    message("  No annotated foreground genes — skipping.")
    return(NULL)
  }
  
  go_fg <- sp_ipr %>%
    filter(protein_id %in% fg_genes) %>%
    count(go_terms, name = "fg_count")
  
  go_bg <- sp_ipr %>%
    count(go_terms, name = "bg_count")
  
  go_tbl <- go_bg %>%
    left_join(go_fg, by = "go_terms") %>%
    mutate(fg_count = replace_na(fg_count, 0)) %>%
    filter(bg_count >= min_gene_count)
  
  message("  Testing ", nrow(go_tbl), " GO terms with >= ", min_gene_count, " annotated genes")
  
  run_fisher <- function(fg_in, bg_in, n_fg, N) {
    bg_out <- N - bg_in
    fg_out <- n_fg - fg_in
    mat <- matrix(c(fg_in, fg_out, bg_in - fg_in, bg_out - fg_out),
                  nrow = 2,
                  dimnames = list(c("fg", "bg"), c("in_term", "not_in_term")))
    mat[mat < 0] <- 0
    fisher.test(mat, alternative = "greater")$p.value
  }
  
  go_tbl <- go_tbl %>%
    mutate(
      p_value    = mapply(run_fisher, fg_count, bg_count,
                          MoreArgs = list(n_fg = n_fg, N = N)),
      p_adjusted = p.adjust(p_value, method = "BH"),
      enrichment = (fg_count / n_fg) / (bg_count / N),
      species    = sp_label
    ) %>%
    filter(p_adjusted < p_cutoff) %>%
    arrange(p_adjusted)
  
  message("  Significant GO terms (BH < ", p_cutoff, "): ", nrow(go_tbl))
  return(go_tbl)
}

results_list <- list(
  run_fisher_enrichment("Fcraterostoma", ipr_all, sp_genes$Fcraterostoma),
  run_fisher_enrichment("Fsalicifolia",  ipr_all, sp_genes$Fsalicifolia),
  run_fisher_enrichment("Fbenjamina",    ipr_all, sp_genes$Fbenjamina)
)
# ── 5. Run for all three species ──────────────────────────────────────────────
results_list <- list(
  run_fisher_enrichment("Fcraterostoma", ipr_all, sp_genes$Fcraterostoma),
  run_fisher_enrichment("Fsalicifolia",  ipr_all, sp_genes$Fsalicifolia),
  run_fisher_enrichment("Fbenjamina",    ipr_all, sp_genes$Fbenjamina)
)

results_all <- bind_rows(results_list)

if (nrow(results_all) == 0) {
  message("\nNo significant GO terms found at BH < ", p_cutoff, ".")
} else {
  message("\nTotal significant GO terms across all species: ", nrow(results_all))
}

# ── 6. Write results ──────────────────────────────────────────────────────────
out_path <- file.path(outdir, "species_specific_GO_enrichment_IPR.tsv")
write_tsv(results_all, out_path)
message("Results written to: ", out_path)

# Get GO term descriptions
library(GO.db)

go_descriptions <- AnnotationDbi::select(
  GO.db,
  keys    = unique(results_all$go_terms),
  columns = c("GOID", "TERM", "ONTOLOGY"),
  keytype = "GOID"
) %>% rename(go_terms = GOID, go_name = TERM)

# ── 7. Bar plot — top 15 GO terms per species (by enrichment score) ───────────
# ── 7. Bar plot — top 15 GO terms per species (by enrichment score) ───────────
library(ggplot2)
library(scales)

sp_colours <- c(
  Fcraterostoma = "#D4547A",
  Fsalicifolia  = "#6B8EC2",
  Fbenjamina    = "#E8C547"
)

# Helper function to nicely format tiny p-values
format_p_val <- function(p) {
  ifelse(p == 0, "p < 1e-300", paste0("p = ", formatC(p, format = "e", digits = 1)))
}

plot_df <- results_all %>%
  left_join(go_descriptions, by = "go_terms") %>%
  mutate(
    go_name = ifelse(is.na(go_name), go_terms, go_name),
    go_name = ifelse(go_terms == "GO:0007205",
                     "phospholipase C-activating GPCR signaling pathway",
                     go_name)
  ) %>%
  group_by(species) %>%
  slice_min(p_adjusted, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # Increased wrapping width from 35 to 55 to prevent label line collision
    label     = str_wrap(go_name, width = 55),
    label     = fct_reorder(label, enrichment),
    p_label   = format_p_val(p_adjusted)
  )

p <- ggplot(plot_df, aes(x = enrichment, y = label, fill = species)) +
  geom_col(show.legend = FALSE, width = 0.75) +
  geom_text(aes(label = p_label),
            hjust = -0.1, size = 2.6) +
  facet_wrap(~ species, scales = "free_y", ncol = 1,
             labeller = labeller(species = c(
               Fcraterostoma = "F. craterostoma",
               Fsalicifolia  = "F. salicifolia",
               Fbenjamina    = "F. benjamina"
             ))) +
  scale_fill_manual(values = sp_colours) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) + # Extra room on right for p-labels
  labs(
    x     = "Fold enrichment",
    y     = NULL,
    title = "GO enrichment — species-specific orthogroups"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text  = element_text(face = "italic"),
    axis.text.y = element_text(size = 8.5, lineheight = 0.8) # Adjusted lineheight to avoid overlap
  )

plot_path <- file.path(outdir, "species_specific_GO_enrichment_IPR.pdf")

# Increased height from 14 to 16 to give wrapped y-axis labels vertical room
ggsave(plot_path, p, width = 10, height = 16)
message("Plot saved to: ", plot_path)