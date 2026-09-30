library(tidyverse)
library(pheatmap)

counts <- read.delim(
  "Orthogroups/Orthogroups.GeneCount.tsv"
)

mat <- counts %>%
  column_to_rownames("Orthogroup")

mat <- as.matrix(mat)

# most variable orthogroups
vars <- apply(
  mat,
  1,
  var
)

top <- names(
  sort(
    vars,
    decreasing = TRUE
  )[1:100]
)

pheatmap(
  mat[top, ],
  scale = "row",
  fontsize_row = 4,
  cluster_cols = TRUE,   
  cluster_rows = TRUE,   
  show_rownames = TRUE, 
  show_colnames = TRUE, 
  fontsize_row = 4,
  fontsize_col = 12,     
  angle_col = 45,        
  main = "Orthogroup Abundance Profiles Across Ficus Species",
  filename = "OrthogroupHeatmap_Final.pdf",
  width = 7,
  height = 9
)