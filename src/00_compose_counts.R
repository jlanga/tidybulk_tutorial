#!/usr/bin/env Rscript

library(tidyverse)

dir.create("results/")

files <- list.files(path = "data/quant", pattern = "tab", full.names = TRUE)

counts <-
  files %>%
  map(
    .f = function(x) {
      read_tsv(
        file = x,
        skip = 4,
        col_names = c("gene", "counts", "forward", "reverse"),
        show_col_types = FALSE
      ) %>%
      select(gene, counts) %>%
      mutate(
        sample = x %>%
          str_remove("data/quant/") %>%
          str_remove("ReadsPerGene") %>%
          str_remove(".out.tab")
      )
    }
  ) %>%
  bind_rows() %>%
  write_tsv("results/counts.tsv")
