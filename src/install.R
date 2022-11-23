#!/usr/bin/env Rscript
install.packages("BiocManager")

BiocManager::install(
  c("tidyverse", "tidybulk", "edgeR", "Rtsne", "sva", "DESeq2", "janitor", "ggrepel"),
  Ncpus = 8
)
