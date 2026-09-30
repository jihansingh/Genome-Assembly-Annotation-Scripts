library(tidyverse)

counts <- read.delim(
  "Orthogroups/Orthogroups.GeneCount.tsv"
)

counts_long <- counts %>%
  pivot_longer(
    -Orthogroup,
    names_to = "Species",
    values_to = "Copies"
  )

ggplot(
  counts_long,
  aes(Copies)
) +
  geom_histogram(
    bins = 30
  ) +
  facet_wrap(
    ~ Species,
    scales = "free_y"
  ) +
  theme_bw()

ggsave(
  "CopyNumberDistribution.pdf",
  width = 10,
  height = 6
)
