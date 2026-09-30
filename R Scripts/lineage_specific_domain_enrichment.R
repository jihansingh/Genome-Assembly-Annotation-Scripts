# ─────────────────────────────────────────────────────────────────────────────
# Domain Enrichment for Species-Specific Orthogroups
# Source: InterPro domain annotations from InterProScan TSV
# Test:   Fisher's exact test per domain, BH correction
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

outdir           <- "domain_enrichment_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

col_fcrat        <- "ficus_craterostoma"
col_fsal         <- "ficus_salicifolia"
col_fbenj        <- "ficus_benjamina"

p_cutoff         <- 0.05
min_gene_count   <- 1

# ── 1. Libraries ──────────────────────────────────────────────────────────────
library(tidyverse)
library(ggtext)
library(data.table)
library(patchwork)

# ── 2. Load OrthoFinder tables and extract species-specific OG gene lists ─────
og_counts <- read_tsv(gene_counts_tsv, show_col_types = FALSE)
og_genes  <- read_tsv(orthogroups_tsv, show_col_types = FALSE)

count_cols <- c(col_fcrat, col_fsal, col_fbenj)

og_sp <- og_counts %>%
  dplyr::select(Orthogroup, all_of(count_cols)) %>%
  mutate(
    sp_specific = case_when(
      .data[[col_fcrat]] > 0 & .data[[col_fsal]]  == 0 & .data[[col_fbenj]] == 0 ~ "Fcraterostoma",
      .data[[col_fsal]]  > 0 & .data[[col_fcrat]] == 0 & .data[[col_fbenj]] == 0 ~ "Fsalicifolia",
      .data[[col_fbenj]] > 0 & .data[[col_fcrat]] == 0 & .data[[col_fsal]]  == 0 ~ "Fbenjamina",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sp_specific))

message("Species-specific orthogroup counts:")
print(table(og_sp$sp_specific))

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

# ── 3. Parse InterProScan TSV → long protein-domain table ────────────────────
parse_ipr_domains <- function(path, species_label) {
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
  
  # Keep only rows with a valid InterPro accession
  df <- df[!is.na(df$ipr_acc) & grepl("^IPR", df$ipr_acc), ]
  
  # One row per unique protein-domain combination
  df <- unique(df[, c("protein_id", "ipr_acc", "ipr_desc")])
  
  df$species <- species_label
  df
}

ipr_all <- bind_rows(
  parse_ipr_domains(ipr_fcrat, "Fcraterostoma"),
  parse_ipr_domains(ipr_fsal,  "Fsalicifolia"),
  parse_ipr_domains(ipr_fbenj, "Fbenjamina")
)

message("Total unique protein-domain pairs: ", nrow(ipr_all))
table(ipr_all$species)

# ── 4. Domain enrichment function ─────────────────────────────────────────────
run_domain_enrichment <- function(sp_label, ipr_df, interesting_genes) {
  
  sp_ipr <- ipr_df %>% filter(species == sp_label)
  
  all_annotated <- unique(sp_ipr$protein_id)
  N <- length(all_annotated)
  
  fg_genes <- intersect(interesting_genes, all_annotated)
  n_fg <- length(fg_genes)
  
  message("\n  [", sp_label, "] ",
          n_fg, " / ", length(interesting_genes),
          " species-specific OG genes have domain annotations",
          " (background N = ", N, ")")
  
  if (n_fg == 0) {
    message("  No annotated foreground genes — skipping.")
    return(NULL)
  }
  
  domain_fg <- sp_ipr %>%
    filter(protein_id %in% fg_genes) %>%
    distinct(protein_id, ipr_acc, ipr_desc) %>%
    count(ipr_acc, ipr_desc, name = "fg_count")
  
  domain_bg <- sp_ipr %>%
    distinct(protein_id, ipr_acc, ipr_desc) %>%
    count(ipr_acc, ipr_desc, name = "bg_count")
  
  domain_tbl <- domain_bg %>%
    left_join(domain_fg, by = c("ipr_acc", "ipr_desc")) %>%
    mutate(fg_count = replace_na(fg_count, 0)) %>%
    filter(bg_count >= min_gene_count)
  
  message("  Testing ", nrow(domain_tbl), " domains")
  
  run_fisher <- function(fg_in, bg_in, n_fg, N) {
    bg_out <- N - bg_in
    fg_out <- n_fg - fg_in
    mat <- matrix(c(fg_in, fg_out, bg_in - fg_in, bg_out - fg_out),
                  nrow = 2)
    mat[mat < 0] <- 0
    fisher.test(mat, alternative = "greater")$p.value
  }
  
  domain_tbl %>%
    mutate(
      p_value    = mapply(run_fisher, fg_count, bg_count,
                          MoreArgs = list(n_fg = n_fg, N = N)),
      p_adjusted = p.adjust(p_value, method = "BH"),
      enrichment = (fg_count / n_fg) / (bg_count / N),
      species    = sp_label
    ) %>%
    filter(p_adjusted < p_cutoff) %>%
    arrange(p_adjusted)
}

# ── 5. Run for all three species ───────────────────────────────────────────────
domain_results_list <- list(
  run_domain_enrichment("Fcraterostoma", ipr_all, sp_genes$Fcraterostoma),
  run_domain_enrichment("Fsalicifolia",  ipr_all, sp_genes$Fsalicifolia),
  run_domain_enrichment("Fbenjamina",    ipr_all, sp_genes$Fbenjamina)
)

domain_results_all <- bind_rows(domain_results_list) %>%
  mutate(p_adjusted = ifelse(p_adjusted == 0, .Machine$double.xmin, p_adjusted))

message("\nTotal significant domains across all species: ", nrow(domain_results_all))
print(table(domain_results_all$species))

# ── 6. Write results ───────────────────────────────────────────────────────────
write_tsv(domain_results_all,
          file.path(outdir, "species_specific_domain_enrichment_IPR.tsv"))
message("Results written to: ", file.path(outdir, "species_specific_domain_enrichment_IPR.tsv"))

# ── 7. Bar plots per species ───────────────────────────────────────────────────
sp_colours <- c(
  Fcraterostoma = "#D4547A",
  Fsalicifolia  = "#6B8EC2",
  Fbenjamina    = "#E8C547"
)

sp_labels <- c(
  Fcraterostoma = "*Ficus craterostoma*",
  Fsalicifolia  = "*Ficus salicifolia*",
  Fbenjamina    = "*Ficus benjamina*"
)

# Use str_wrap to preserve readable domain names without line collision
domain_plot_df <- domain_results_all %>%
  mutate(label = stringr::str_wrap(ipr_desc, width = 38)) %>%
  group_by(species) %>%
  slice_min(p_adjusted, n = 15, with_ties = FALSE) %>%
  ungroup()

plot_list <- list()

for (sp in c("Fcraterostoma", "Fsalicifolia", "Fbenjamina")) {
  
  sp_df <- domain_plot_df %>% 
    filter(species == sp) %>%
    mutate(label = fct_reorder(label, enrichment))
  
  if (nrow(sp_df) == 0) {
    p <- ggplot() + 
      theme_void() + 
      labs(title = sp_labels[sp]) +
      annotate("text", x = 0.5, y = 0.5, label = "No significantly enriched domains") +
      theme(plot.title = element_markdown(size = 11, face = "bold"))
  } else {
    
    # Calculate max x-value to give generous padding for p-value labels
    max_x <- max(sp_df$enrichment, na.rm = TRUE)
    
    p <- ggplot(sp_df, aes(x = enrichment, y = label, fill = species)) +
      geom_col(show.legend = FALSE, width = 0.65) + 
      
      # Add text with slight nudge and clipping turned off so it never overlaps the bar
      geom_text(
        aes(label = paste0("p=", signif(p_adjusted, 2))),
        hjust = -0.15, 
        size = 2.5
      ) +
      scale_fill_manual(values = sp_colours) +
      
      # Add 35% margin space on the right side of the x-axis for text
      scale_x_continuous(
        limits = c(0, max_x * 1.35),
        expand = c(0, 0)
      ) +
      coord_cartesian(clip = "off") + # Allows p-value labels to extend past axis edge if needed
      
      labs(x = "Fold enrichment", y = NULL, title = sp_labels[sp]) +
      theme_bw(base_size = 9) +
      theme(
        plot.title       = element_markdown(size = 11, face = "bold"),
        # Lineheight 0.8 keeps wrapped labels compact without colliding vertically
        axis.text.y      = element_text(size = 7, color = "black", lineheight = 0.8),
        axis.title.x     = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(t = 5, r = 25, b = 5, l = 5) # Added right margin
      )
  }
  
  plot_list[[sp]] <- p
}

# Assemble panels with vertical spacing
combined_figure <- plot_list$Fbenjamina / 
  plot_list$Fcraterostoma / 
  plot_list$Fsalicifolia +
  plot_annotation(
    title = "Domain Enrichment in Species-Specific Orthogroups",
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.tag   = element_text(size = 13, face = "bold")
    )
  )

# Save with height = 15 to give 15 rows in each panel room to breathe
ggsave(file.path(outdir, "Figure_Combined_Domain_Enrichment.pdf"),
       combined_figure, width = 8.5, height = 15)

ggsave(file.path(outdir, "Figure_Combined_Domain_Enrichment.png"),
       combined_figure, width = 8.5, height = 15, dpi = 600)

message("Combined figure updated and saved to: ", outdir)
