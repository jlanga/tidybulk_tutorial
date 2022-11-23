#!/usr/bin/env Rscript

library(tidybulk)
library(tidyverse)

# 1 - Read tables ----

sample_metadata <- read_csv("data/reference/geuvadis_phenodata.csv")
sample_metadata <- sample_metadata %>% rename(sample = ids)
sample_metadata

counts <- read_tsv("results/counts.tsv")
counts

# Add sample info to counts
counts %>%
  left_join(sample_metadata)

# Convert to tidybulk
counts <-
  counts %>%
  left_join(sample_metadata) %>%
  tidybulk(.sample = sample, .transcript = gene, .abundance = counts)

counts


# 2 - Explore ----
## 2.1 PCA ----
pca <-
  counts %>%
  identify_abundant() %>%
  reduce_dimensions(method = "PCA")

ggplot(pca, aes(x = PC1, y = PC2, color = sex)) +
  geom_point()

ggplot(pca, aes(x = PC1, y = PC2, color = population)) +
  geom_point()

## 2.2 MDS ----
mds <-
  counts %>%
  identify_abundant() %>%
  reduce_dimensions(method = "MDS")

ggplot(mds, aes(x = Dim1, y = Dim2, color = sex)) +
  geom_point()

ggplot(mds, aes(x = Dim1, y = Dim2, color = population)) +
  geom_point()


## 02.3 t-SNE ----

# tsne <-
#   counts %>%
#   identify_abundant() %>%
#   reduce_dimensions(method = "tSNE")



# 3 - Normalize data ----

counts_normalized <-
  counts %>%
  identify_abundant(factor_of_interest =  sex, minimum_counts = 2, minimum_proportion = 0.7) %>%
  keep_abundant(factor_of_interest = sex) %>%
  scale_abundance() %>%
  adjust_abundance(~ sex + population)

data_normalization <-
  counts_normalized %>%
  pivot_longer(
    c(counts, counts_scaled, counts_scaled_adjusted),
    values_to = "count",
    names_to = "normalisation"
  )

plot_normalization <-
  data_normalization %>%
  ggplot(aes(count + 1, color = sample)) +
  geom_density(na.rm = TRUE) +
  scale_x_log10() +
  facet_grid(~normalisation)
plot_normalization


# 4 - Inspect normalized data ----

## 4.1 PCA ----

counts_normalized %>%
  reduce_dimensions(method = "PCA") %>%
  ggplot(aes(x = PC1, y = PC2, color = sex)) +
  geom_point()

## 4.2 MDS ----
counts_normalized %>%
  reduce_dimensions(method = "MDS") %>%
  ggplot(aes(x = Dim1, y = Dim2, color = sex)) +
  geom_point()


# 5 - Differential expression ----
differential_expression <-
  counts_normalized %>%
  test_differential_abundance(
    ~ 0 + sex + population,
    contrasts = "sexfemale - sexmale"
  ) %>%
  rename(
    logFC = `logFC___sexfemale - sexmale`,
    logCPM = `logCPM___sexfemale - sexmale`,
    Fvalue = `F___sexfemale - sexmale`,
    PValue = `PValue___sexfemale - sexmale`,
    FDR = `FDR___sexfemale - sexmale`
  ) %>%
  pivot_transcript()

differential_expression

# 5.1 MA Plot ----

data_de_plot <-
  differential_expression %>%
  mutate(
    significant = FDR < 0.05 & abs(logFC) > 1.2,
    direction = sign(logFC),
    label = if_else(significant, gene,  NA_character_),
    color = if_else(significant, "red", "black")
  )

ma_plot <-
  data_de_plot %>%
  ggplot(aes(x = logCPM, y = logFC, label = label)) +
  geom_point(aes(color = color, size = significant)) +
  ggrepel::geom_text_repel() +
  scale_size_discrete(range = c(0, 2))
ma_plot



volcano_plot <-
  data_de_plot %>%
  ggplot(aes(x = logFC, y = FDR, label = label)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  ggrepel::geom_text_repel() +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_y_continuous(trans = "log10_reverse")
volcano_plot



# Identify that gene

feature_metadata <- read_tsv("results/gene_annotation.tsv")
feature_metadata

data_de_plot %>%
  left_join(feature_metadata) %>%
  filter(significant)
