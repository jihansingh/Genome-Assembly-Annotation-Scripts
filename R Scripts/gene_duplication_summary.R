library(tidyverse)
library(ggtree)
library(treeio)

# 1. Read the Orthofinder duplication tree file
tree <- read.tree("Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt")

# Create a clean data frame to hold our extracted tip numbers
# This safely pulls the numbers out before we clean the tree labels
tip_data <- data.frame(
  label = tree$tip.label,
  # This regex isolates just the number at the end of the text string
  dup_count = as.numeric(str_extract(tree$tip.label, "\\d+$"))
)

# Clean up the tree's tip labels so they are JUST the species names
tree$tip.label <- tree$tip.label %>%
  str_remove("\\s*\\d+$") %>%   # Remove the numbers
  str_replace_all("_", " ")      # Replace underscores with clean spaces

# Re-assign clean names as the row identifier for matching
tip_data$clean_label <- tree$tip.label

# 2. Plot the polished tree
ggtree(tree, linewidth = 1.2) %<+% tip_data + # The "%<+%" operator links our data frame to the tree!
  
  # Clean, purely italicized species names
  geom_tiplab(aes(label = clean_label), fontface = "italic", size = 5, offset = 0.02) + 
  
  # Add the terminal branch duplication counts at the tips using your Gold/Yellow anchor color
  geom_tiplab(aes(label = scales::comma(dup_count)), color = "#1E88E5", 
              fontface = "bold", size = 5, offset = -0.05, hjust = 1) +
  
  # Clean up internal node names to display just the ancestral duplication number (e.g., 324 instead of N1_324)
  geom_nodelab(aes(label = str_remove(label, ".*_")), vjust = -0.6, hjust = 1.3, 
               color = "#D81B60", fontface = "bold", size = 5) +
  
  theme_tree2() + 
  labs(
    title = "Lineage-Specific Gene Duplication Events Map",
    subtitle = "Pink = Shared Ancestral Events | Blue = Terminal Species-Specific Expansions",
    x = "Evolutionary Distance"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4)))

# 3. Save the final publication tree
ggsave("SpeciesTree_Duplications_Final.pdf", width = 8, height = 5)
