suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# ==============================================================================
# CONFIGURATION & CONSTANTS
# ==============================================================================
ORTHOFINDER_DUP_PATH <- "OrthoFinder/Gene_Duplication_Events/Duplications.tsv"

species_config <- list(
  craterostoma = list(
    ips_path  = "ficus_craterostoma/interproscan.tsv",
    node_name = "ficus_craterostoma"  # Must match 'Species Tree Node' in Duplications.tsv
  ),
  salicifolia = list(
    ips_path  = "ficus_salicifolia/interproscan.tsv",
    node_name = "ficus_salicifolia"
  ),
  benjamina = list(
    ips_path  = "ficus_benjamina/interproscan.tsv",
    node_name = "ficus_benjamina"
  )
)

FDR_THRESHOLD   <- 0.05
MIN_TERM_GENES  <- 2
MAX_TERM_GENES  <- 500
OUT_DIR         <- "domain_go_ora_results"
dir.create(OUT_DIR, showWarnings = FALSE)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

normalize_gene_id <- function(x) {
  x <- sub("^ficus_[a-z]+_", "", x)
  x <- sub("^F_[a-z]_", "", x)
  x <- sub("-mRNA-[0-9]+$", "", x)
  x
}

# Extracts InterPro Domains (Column 12 & 13) and GO Terms (Column 14) from InterProScan TSV
read_ips_annotations <- function(ips_path) {
  # InterProScan default columns:
  # Col 1: Gene/Protein ID
  # Col 12: InterPro Accession (e.g. IPR000001)
  # Col 13: InterPro Description
  # Col 14: GO Annotations (e.g. GO:0008152|GO:0003674)
  
  dt <- fread(ips_path, sep = "\t", header = FALSE, fill = TRUE,
              quote = "", na.strings = c("-", ""))
  
  # Ensure file has sufficient columns
  if (ncol(dt) < 14) {
    stop(sprintf("%s is missing required InterPro/GO columns (expected >=14 cols).", ips_path))
  }
  
  gene_ids <- normalize_gene_id(dt[[1]])
  
  # --- 1. Extract InterPro Domains ---
  ipr_acc  <- dt[[12]]
  ipr_desc <- dt[[13]]
  
  domain_dt <- data.table(
    gene_id     = gene_ids,
    annotation  = ifelse(!is.na(ipr_acc) & !is.na(ipr_desc), 
                         paste0(ipr_acc, ": ", ipr_desc), 
                         ipr_acc),
    type        = "Domain"
  )
  domain_dt <- domain_dt[!is.na(annotation) & annotation != ""]
  
  # --- 2. Extract GO Terms ---
  go_raw <- dt[[14]]
  go_dt  <- data.table(gene_id = gene_ids, go_str = go_raw)
  go_dt  <- go_dt[!is.na(go_str) & go_str != ""]
  
  go_long <- go_dt[, .(annotation = unlist(str_split(go_str, "\\|"))), by = gene_id]
  go_long[, annotation := str_trim(annotation)]
  go_long <- go_long[annotation != ""]
  # Keep only standard GO accessions
  go_long <- go_long[str_detect(annotation, "^GO:[0-9]+")]
  go_long[, type := "GO"]
  
  # Combine annotations and eliminate duplicates per gene
  combined <- rbind(domain_dt, go_long)
  unique(combined)
}

get_duplicated_genes <- function(dup_tsv_path, species_node_name, all_genes) {
  dt <- fread(dup_tsv_path, sep = "\t")
  
  avail_nodes <- unique(dt$`Species Tree Node`)
  matched_node <- species_node_name
  
  if (!species_node_name %in% avail_nodes) {
    matched <- avail_nodes[str_detect(avail_nodes, regex(species_node_name, ignore_case = TRUE))]
    if (length(matched) > 0) {
      message(sprintf("    [INFO] Auto-matching node name '%s' -> '%s'", species_node_name, matched[1]))
      matched_node <- matched[1]
    } else {
      warning(sprintf("    [WARNING] Node '%s' not found in Duplications.tsv!", species_node_name))
    }
  }
  
  message("    [INFO] Using node: '", matched_node, "'")
  
  node_dt <- dt[`Species Tree Node` == matched_node]
  
  g1 <- unlist(str_split(node_dt[["Genes 1"]], ",\\s*"))
  g2 <- unlist(str_split(node_dt[["Genes 2"]], ",\\s*"))
  
  raw_genes <- normalize_gene_id(c(g1, g2))
  raw_genes <- unique(raw_genes[raw_genes != "" & !is.na(raw_genes)])
  
  intersect(raw_genes, all_genes)
}

run_ora_enrichment <- function(annotation_long, duplicated_genes, all_genes,
                               min_size = MIN_TERM_GENES, max_size = MAX_TERM_GENES) {
  
  term_list <- split(annotation_long$gene_id, annotation_long$annotation)
  term_list <- lapply(term_list, function(g) intersect(unique(g), all_genes))
  
  sizes <- sapply(term_list, length)
  term_list <- term_list[sizes >= min_size & sizes <= max_size]
  
  n_dup_total <- length(duplicated_genes)
  n_bg_total  <- length(all_genes)
  
  res <- rbindlist(lapply(names(term_list), function(term) {
    term_genes    <- term_list[[term]]
    n_term        <- length(term_genes)
    n_dup_in_term <- length(intersect(term_genes, duplicated_genes))
    
    # Contingency Table:
    #                 in_term    not_in_term
    # duplicated         a            b
    # not duplicated     c            d
    a <- n_dup_in_term
    b <- n_dup_total - a
    c <- n_term - a
    d <- n_bg_total - n_term - b
    
    m <- matrix(c(a, b, c, d), nrow = 2)
    ft <- fisher.test(m, alternative = "greater")
    
    data.table(
      term           = term,
      size           = n_term,
      dup_in_term    = a,
      dup_total      = n_dup_total,
      pval           = ft$p.value,
      odds_ratio     = unname(ft$estimate)
    )
  }))
  
  if (nrow(res) > 0) {
    res[, padj := p.adjust(pval, method = "BH")]
    setorder(res, padj, pval)
  }
  return(res)
}

# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================

all_domain_results <- list()
all_go_results     <- list()

for (sp in names(species_config)) {
  message("Processing ", sp, " ...")
  cfg <- species_config[[sp]]
  
  # Load annotations
  annot_long <- read_ips_annotations(cfg$ips_path)
  
  # Retrieve complete gene background
  raw_ips_genes <- fread(cfg$ips_path, sep = "\t", header = FALSE, fill = TRUE, quote = "", select = 1)[[1]]
  all_gene_ids  <- unique(normalize_gene_id(raw_ips_genes))
  
  # Identify gene duplications for focal species
  duplicated_genes <- get_duplicated_genes(ORTHOFINDER_DUP_PATH, cfg$node_name, all_gene_ids)
  
  message("  ", sp, ": ", length(duplicated_genes), " duplicated genes out of ",
          length(all_gene_ids), " total background genes.")
  
  # --- 1. Domain ORA ---
  domain_long <- annot_long[type == "Domain"]
  res_domain  <- run_ora_enrichment(domain_long, duplicated_genes, all_gene_ids)
  res_domain[, species := sp]
  fwrite(res_domain, file.path(OUT_DIR, paste0(sp, "_ora_domain_enrichment.csv")))
  all_domain_results[[sp]] <- res_domain
  
  # --- 2. GO ORA ---
  go_long <- annot_long[type == "GO"]
  res_go  <- run_ora_enrichment(go_long, duplicated_genes, all_gene_ids)
  res_go[, species := sp]
  fwrite(res_go, file.path(OUT_DIR, paste0(sp, "_ora_go_enrichment.csv")))
  all_go_results[[sp]] <- res_go
}

# ==============================================================================
# CROSS-SPECIES COMPARISON & LINEAGE SPECIFICITY
# ==============================================================================

process_lineage_specificity <- function(results_list, category_name) {
  if (length(results_list) != length(species_config)) return(NULL)
  
  combined <- rbindlist(results_list, fill = TRUE)
  
  wide <- dcast(
    combined,
    term ~ species,
    value.var = c("pval", "padj", "odds_ratio", "size", "dup_in_term")
  )
  
  # Significance filter
  wide[, craterostoma_sig := !is.na(padj_craterostoma) & padj_craterostoma < FDR_THRESHOLD & odds_ratio_craterostoma > 1]
  wide[, salicifolia_sig  := !is.na(padj_salicifolia)  & padj_salicifolia < FDR_THRESHOLD  & odds_ratio_salicifolia > 1]
  wide[, benjamina_sig    := !is.na(padj_benjamina)    & padj_benjamina < FDR_THRESHOLD    & odds_ratio_benjamina > 1]
  
  # Lineage specificity checks
  wide[, craterostoma_lineage_specific := craterostoma_sig &
         (is.na(pval_salicifolia) | pval_salicifolia > 0.05 | odds_ratio_salicifolia <= 1) &
         (is.na(pval_benjamina) | pval_benjamina > 0.05 | odds_ratio_benjamina <= 1)]
  
  wide[, salicifolia_lineage_specific := salicifolia_sig &
         (is.na(pval_craterostoma) | pval_craterostoma > 0.05 | odds_ratio_craterostoma <= 1) &
         (is.na(pval_benjamina) | pval_benjamina > 0.05 | odds_ratio_benjamina <= 1)]
  
  wide[, benjamina_lineage_specific := benjamina_sig &
         (is.na(pval_craterostoma) | pval_craterostoma > 0.05 | odds_ratio_craterostoma <= 1) &
         (is.na(pval_salicifolia) | pval_salicifolia > 0.05 | odds_ratio_salicifolia <= 1)]
  
  out_file <- file.path(OUT_DIR, paste0("lineage_specific_ora_", category_name, ".csv"))
  fwrite(wide, out_file)
  
  message(sprintf("\n================ Lineage-Specific ORA Summary (%s) ================", category_name))
  message("F. craterostoma specific duplicated terms: ", sum(wide$craterostoma_lineage_specific, na.rm = TRUE))
  message("F. salicifolia specific duplicated terms:  ", sum(wide$salicifolia_lineage_specific, na.rm = TRUE))
  message("F. benjamina specific duplicated terms:    ", sum(wide$benjamina_lineage_specific, na.rm = TRUE))
  message("============================================================================\n")
}

process_lineage_specificity(all_domain_results, "domains")
process_lineage_specificity(all_go_results, "go_terms")

# Save craterostoma significant domain table as an example summary export
if (!is.null(all_domain_results$craterostoma)) {
  fwrite(
    all_domain_results$craterostoma[padj < FDR_THRESHOLD & odds_ratio > 1, 
                                    .(Domain = term, `Genes Duplicated` = dup_in_term, `Term Size` = size, `Odds Ratio` = round(odds_ratio, 2), `Adjusted p-value` = signif(padj, 3))], 
    file.path(OUT_DIR, "craterostoma_significant_domains_table.csv")
  )
}
