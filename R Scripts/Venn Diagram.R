# ==============================================================================
# Orthogroup Venn Diagram Generation Script
# Finalized Layout: Maximized Circles & Styled Labels in Outer Whitespace
# ==============================================================================

library(VennDiagram)
library(grid)

# 1. Load Data
orth <- read.delim(
  "Orthogroups/Orthogroups.tsv",
  check.names = FALSE,
  colClasses = "character" 
)

# 2. Extract Data Vectors
craterostoma <- which(!is.na(orth$ficus_craterostoma) & trimws(orth$ficus_craterostoma) != "")
salicifolia  <- which(!is.na(orth$ficus_salicifolia)  & trimws(orth$ficus_salicifolia)  != "")
benjamina    <- which(!is.na(orth$ficus_benjamina)    & trimws(orth$ficus_benjamina)    != "")

# 3. Build Base Plot (Maximized circles, hidden default labels)
venn.plot <- venn.diagram(
  x = list(
    "F. craterostoma (2n)" = craterostoma,
    "F. salicifolia (4n)"  = salicifolia,
    "F. benjamina (2n)"    = benjamina
  ),
  category.names = c("", "", ""), 
  filename = NULL,
  fill = c("#D81B60", "#1E88E5", "#FFC107"), 
  alpha = c(0.5, 0.5, 0.5),
  col = c("#B0124C", "#1565C0", "#E6A100"),
  
  # Tiny margin ensures the circles expand to full size on the canvas
  margin = 0.02, 
  
  # Counts styling with a clean white backing glow
  label.col = "black",
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  label.box = TRUE,              
  label.box.fill = "white",      
  label.box.col = "transparent", 
  label.box.alpha = 0.6
)

# 4. Open PDF Canvas Device
pdf("Orthogroup_Venn_Perfect_Spacing.pdf", width = 8, height = 8)

# 5. Render the Venn Diagram Structure
grid.draw(venn.plot)

# 6. Overlay Italicized Species Labels at Absolute Boundaries
# Left-aligned and pushed out into the top-left corner whitespace
grid.text(
  label = expression(italic("F. craterostoma")~"(2n)"),
  x = unit(0.05, "npc"), y = unit(0.97, "npc"), 
  just = "left", 
  gp = gpar(cex = 1.2, fontfamily = "sans", fontface = "bold")
)

# Right-aligned and pushed out into the top-right corner whitespace
grid.text(
  label = expression(italic("F. salicifolia")~"(4n)"),
  x = unit(0.95, "npc"), y = unit(0.97, "npc"), 
  just = "right", 
  gp = gpar(cex = 1.2, fontfamily = "sans", fontface = "bold")
)

# Centered at the absolute bottom edge workspace boundary
grid.text(
  label = expression(italic("F. benjamina")~"(2n)"),
  x = unit(0.50, "npc"), y = unit(0.03, "npc"), 
  just = "center",
  gp = gpar(cex = 1.2, fontfamily = "sans", fontface = "bold")
)

# Close and finalize file save
dev.off()
