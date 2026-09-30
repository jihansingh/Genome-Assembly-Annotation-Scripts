library(tidyverse)
library(GO.db)

# paths
og_file  <- "C:/Users/jihan/Documents/Uni_Work/Masters/Orthofinder/Orthogroups/Orthogroups.tsv"
ips_benj <- "C:/Users/jihan/Documents/Uni_Work/Masters/ficus_benjamina/interproscan.tsv"
ips_crat <- "C:/Users/jihan/Documents/Uni_Work/Masters/ficus_craterostoma/interproscan.tsv"
ips_sali <- "C:/Users/jihan/Documents/Uni_Work/Masters/ficus_salicifolia/interproscan.tsv"
out_file <- "C:/Users/jihan/Documents/Uni_Work/Masters/go_orthogroup_summary.tsv"

# read one InterProScan TSV into gene_id, go_id
read_ips_go <- function(path) {
  read_tsv(path, col_names = FALSE, col_types = cols(.default = "c"), quote = "") %>%
    select(gene_id = X1, go_id = X14) %>%
    filter(!is.na(go_id), go_id != "-") %>%
    separate_rows(go_id, sep = "\\|") %>%
    mutate(go_id = str_remove(go_id, "\\(.*\\)$")) %>%
    distinct(gene_id, go_id)
}

go <- bind_rows(read_ips_go(ips_benj), read_ips_go(ips_crat), read_ips_go(ips_sali)) %>%
  distinct(gene_id, go_id)

# GO names
go_names <- AnnotationDbi::select(GO.db, keys = unique(go$go_id),
                                  columns = "TERM", keytype = "GOID") %>%
  rename(go_id = GOID, go_name = TERM)

# orthogroups: wide to long, one row per gene
og <- read_tsv(og_file, col_types = cols(.default = "c")) %>%
  pivot_longer(-Orthogroup, names_to = "species", values_to = "gene_id") %>%
  filter(!is.na(gene_id)) %>%
  separate_rows(gene_id, sep = ",\\s*")

# diagnostics: species labels and ID match rate
unique(og$species)
n_distinct(go$gene_id)
n_distinct(inner_join(og, go, by = "gene_id")$gene_id)

# genes per orthogroup within each species x GO term
per_og <- og %>%
  inner_join(go, by = "gene_id", relationship = "many-to-many") %>%
  count(species, go_id, Orthogroup, name = "genes_in_og")

summary_tbl <- per_og %>%
  group_by(species, go_id) %>%
  summarise(
    n_genes       = sum(genes_in_og),
    n_orthogroups = n(),
    largest_og    = Orthogroup[which.max(genes_in_og)],
    largest_og_n  = max(genes_in_og),
    top3_genes    = sum(sort(genes_in_og, decreasing = TRUE)[1:min(3, n())]),
    .groups = "drop"
  ) %>%
  mutate(
    genes_per_og  = n_genes / n_orthogroups,
    top3_fraction = top3_genes / n_genes
  ) %>%
  left_join(go_names, by = "go_id") %>%
  relocate(go_name, .after = go_id) %>%
  arrange(species, desc(n_genes))

write_tsv(summary_tbl, out_file)

# your terms of interest in F. salicifolia
terms <- c("protein binding", "proteolysis",
           "cysteine-type peptidase activity", "zinc ion binding")
summary_tbl %>%
  filter(str_detect(species, "salicifolia"), go_name %in% terms)