#!/usr/bin/env Rscript

library(tidyverse)
library(janitor)

dir.create("results/transcriptome/", showWarnings = FALSE, recursive = TRUE )

annotation_raw <-
  read_tsv(
    file = "data/reference/chrX.gtf",
    comment = "#",
    col_names = c(
      "seqid", "source", "type", "start", "end", "score", "strand",
      "phase", "attributes"
    )
  )


gene_annotation <-
  annotation_raw %>%
  select(attributes) %>%
  distinct() %>%
  mutate(
    attributes = attributes %>%
      str_remove_all("\\\"") %>%
      str_replace_all("  ", " ") %>%
      str_remove(";$") %>%
      str_replace_all("; ", ";"),
  ) %>%
  separate(attributes, into = c("gene_id", "transcript_id", "gene_name", "product"), sep = ";", remove = TRUE) %>%
  mutate(
    gene_id = gene_id %>% str_remove("^gene_id "),
    transcript_id = transcript_id %>% str_remove("^transcript_id "),
    gene_name = gene_name %>% str_remove("^gene_name "),
    product = product %>% str_remove("^product "),
  ) %>%
  select(-transcript_id) %>%
  rename(gene = gene_id) %>%
  distinct() %>%
  write_tsv("results/gene_annotation.tsv")
