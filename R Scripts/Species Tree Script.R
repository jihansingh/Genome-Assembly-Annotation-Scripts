library(ggtree)
library(ape)

tree <- read.tree("Species_Tree/SpeciesTree_rooted.txt") 

p <- ggtree(tree) + 
  geom_tiplab(size = 4) + 
  theme_tree2() +
  hexpand(0.3, direction = 1)


ggsave( "SpeciesTree.pdf", p, width = 8, height = 6 )
