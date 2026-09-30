suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

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

normalize_gene_id <- function(x) {
  x <- sub("^ficus_[a-z]+_", "", x)
  x <- sub("^F_[a-z]_", "", x)
  x <- sub("-mRNA-[0-9]+$", "", x)
  x
}

PLANTCYC_ID_PATH    <- "plantcyc_pathway_ids.csv"
STRICT_PLANT_FILTER <- TRUE
PATHWAY_SOURCE   <- "MetaCyc"
FDR_THRESHOLD    <- 0.05
MIN_PATH_GENES   <- 2
MAX_PATH_GENES   <- 500
OUT_DIR <- "pathway_gsea_results"
dir.create(OUT_DIR, showWarnings = FALSE)

read_ips_pathways <- function(ips_path) {
  dt <- fread(ips_path, sep = "\t", header = FALSE, fill = TRUE,
              quote = "", na.strings = c("-", ""))
  
  if (ncol(dt) < 15) {
    stop(sprintf("%s missing pathways column (expected >=15 cols).", ips_path))
  }
  
  setnames(dt, c(1, 15), c("gene_id", "pathways"))
  dt <- dt[!is.na(pathways), .(gene_id, pathways)]
  
  long <- dt[, .(pathway = unlist(str_split(pathways, "\\|"))), by = gene_id]
  long[, pathway := str_trim(pathway)]
  long <- unique(long[pathway != ""])
  
  long[, source := fcase(
    str_detect(pathway, regex("MetaCyc", ignore_case = TRUE)), "MetaCyc",
    str_detect(pathway, regex("KEGG", ignore_case = TRUE)),    "KEGG",
    str_detect(pathway, regex("Reactome", ignore_case = TRUE)),"Reactome",
    default = "Other"
  )]
  long
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
  
  # Return only duplicated genes that are also present in the annotated background
  intersect(raw_genes, all_genes)
}

load_plantcyc_ids <- function(path) {
  ids <- fread(path)[[1]]
  unique(str_trim(as.character(ids)))
}

extract_pathway_id <- function(pathway_string) {
  res <- str_remove(pathway_string, "^MetaCyc:\\s*")
  str_trim(res)
}

run_ora_enrichment <- function(pathway_long, duplicated_genes, all_genes,
                               min_size = MIN_PATH_GENES, max_size = MAX_PATH_GENES) {
  pathway_list <- split(pathway_long$gene_id, pathway_long$pathway)
  pathway_list <- lapply(pathway_list, function(g) intersect(unique(g), all_genes))
  
  sizes <- sapply(pathway_list, length)
  pathway_list <- pathway_list[sizes >= min_size & sizes <= max_size]
  
  n_dup_total <- length(duplicated_genes)
  n_bg_total  <- length(all_genes)
  
  res <- rbindlist(lapply(names(pathway_list), function(p) {
    path_genes    <- pathway_list[[p]]
    n_path        <- length(path_genes)
    n_dup_in_path <- length(intersect(path_genes, duplicated_genes))
    
    # 2x2 contingency table:
    #                 in_pathway   not_in_pathway
    # duplicated         a               b
    # not duplicated      c               d
    a <- n_dup_in_path
    b <- n_dup_total - a
    c <- n_path - a
    d <- n_bg_total - n_path - b
    
    m <- matrix(c(a, b, c, d), nrow = 2)
    ft <- fisher.test(m, alternative = "greater")
    
    data.table(
      pathway       = p,
      size          = n_path,
      dup_in_path   = a,
      dup_total     = n_dup_total,
      pval          = ft$p.value,
      odds_ratio    = unname(ft$estimate)
    )
  }))
  
  res[, padj := p.adjust(pval, method = "BH")]
  setorder(res, padj, pval)
  res
}

all_results <- list()

for (sp in names(species_config)) {
  message("Processing ", sp, " ...")
  cfg <- species_config[[sp]]
  
  pathway_long <- read_ips_pathways(cfg$ips_path)
  pathway_long <- pathway_long[source == PATHWAY_SOURCE]
  
  if (nrow(pathway_long) == 0) {
    warning(sp, ": no ", PATHWAY_SOURCE, " pathway annotations found - skipping.")
    next
  }
  
  pathway_long[, pathway_id := extract_pathway_id(pathway)]
  pathway_long[, gene_id := normalize_gene_id(gene_id)]
  
  if (STRICT_PLANT_FILTER) {
    plantcyc_ids <- load_plantcyc_ids(PLANTCYC_ID_PATH)
    pathway_long <- pathway_long[pathway_id %in% plantcyc_ids]
  }
  
  raw_ips_genes <- fread(cfg$ips_path, sep = "\t", header = FALSE, fill = TRUE, quote = "", select = 1)[[1]]
  all_gene_ids  <- unique(normalize_gene_id(raw_ips_genes))
  
  duplicated_genes <- get_duplicated_genes(ORTHOFINDER_DUP_PATH, cfg$node_name, all_gene_ids)
  
  message("  ", sp, ": ", length(duplicated_genes), " duplicated genes identified out of ",
          length(all_gene_ids), " total background genes.")
  
  res <- run_ora_enrichment(pathway_long, duplicated_genes, all_gene_ids)
  res[, species := sp]
  
  out_path <- file.path(OUT_DIR, paste0(sp, "_ora_pathway_enrichment.csv"))
  fwrite(res, out_path)
  
  all_results[[sp]] <- res
}

if (length(all_results) == length(species_config)) {
  
  combined <- rbindlist(all_results, fill = TRUE)
  
  wide <- dcast(
    combined,
    pathway ~ species,
    value.var = c("pval", "padj", "odds_ratio", "size", "dup_in_path")
  )
  
  # Significant = FDR < threshold AND odds ratio > 1 (i.e. over-, not under-represented)
  wide[, craterostoma_sig := !is.na(padj_craterostoma) & padj_craterostoma < FDR_THRESHOLD & odds_ratio_craterostoma > 1]
  wide[, salicifolia_sig  := !is.na(padj_salicifolia)  & padj_salicifolia < FDR_THRESHOLD  & odds_ratio_salicifolia > 1]
  wide[, benjamina_sig    := !is.na(padj_benjamina)    & padj_benjamina < FDR_THRESHOLD    & odds_ratio_benjamina > 1]
  
  # Lineage specificity: significant in the focal species, and not even nominally
  # significant (raw pval > 0.05) or not enriched in the same direction elsewhere
  wide[, craterostoma_lineage_specific := craterostoma_sig &
         (is.na(pval_salicifolia) | pval_salicifolia > 0.05 | odds_ratio_salicifolia <= 1) &
         (is.na(pval_benjamina) | pval_benjamina > 0.05 | odds_ratio_benjamina <= 1)]
  
  wide[, salicifolia_lineage_specific := salicifolia_sig &
         (is.na(pval_craterostoma) | pval_craterostoma > 0.05 | odds_ratio_craterostoma <= 1) &
         (is.na(pval_benjamina) | pval_benjamina > 0.05 | odds_ratio_benjamina <= 1)]
  
  wide[, benjamina_lineage_specific := benjamina_sig &
         (is.na(pval_craterostoma) | pval_craterostoma > 0.05 | odds_ratio_craterostoma <= 1) &
         (is.na(pval_salicifolia) | pval_salicifolia > 0.05 | odds_ratio_salicifolia <= 1)]
  
  comparison_path <- file.path(OUT_DIR, "lineage_specific_ora_pathways.csv")
  fwrite(wide, comparison_path)
  
  message("\n================ Lineage-Specific ORA Summary ================")
  message("F. craterostoma specific duplicated pathways: ", sum(wide$craterostoma_lineage_specific, na.rm = TRUE))
  message("F. salicifolia specific duplicated pathways:  ", sum(wide$salicifolia_lineage_specific, na.rm = TRUE))
  message("F. benjamina specific duplicated pathways:    ", sum(wide$benjamina_lineage_specific, na.rm = TRUE))
  message("================================================================\n")
}
fwrite(all_results$craterostoma[padj < 0.05 & odds_ratio > 1, .(Pathway = pathway, `Genes Duplicated` = dup_in_path, `Pathway Size` = size, `Odds Ratio` = round(odds_ratio, 2), `Adjusted p-value` = signif(padj, 3))], "pathway_gsea_results/craterostoma_significant_pathways_table.csv")
