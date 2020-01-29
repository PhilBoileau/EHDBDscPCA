# t-SNE on AML patient 027

library(here)
library(tidyverse)
library(Rtsne)
library(SingleCellExperiment)

# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(87612)

# retain genes that contain 5 or more non-zero counts accross all cells
sce <- sce[(Matrix::rowSums(counts(sce) != 0) > 4), ]

# retain 1000 most variable genes
vars <- rowVars(as.matrix(log1p(counts(sce))))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
core <- sce[names(vars)[1:1000], ]

# split into target and background datasets
patient_027_sce <- core[, which(core$cell_prov %in%
                                  c("pre_trans_027", "post_trans_027"))]
patient_027 <- as.matrix(t(counts(patient_027_sce)))

# get tsne representation
bmmc_tsne_pca <- Rtsne(as.matrix(patient_027), perplexity = 30, max_iter = 1000,
                       pca = TRUE, theta = 0)
saveRDS(bmmc_tsne_pca, file = here("analyses/bmmc_data/data/bmmc_tsne_pca_027.rds"))
