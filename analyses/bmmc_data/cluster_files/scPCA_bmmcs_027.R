################################################################################
############################# scPCA on BMMCs Data ##############################
################################################################################

# AML027

################################################################################

library(here)
library(tidyverse)
library(scPCA)
library(SingleCellExperiment)
library(BiocParallel)

# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(7860973)

# retain genes that contain 5 or more non-zero counts accross all cells
sce <- sce[(Matrix::rowSums(counts(sce) != 0) > 4), ]

# retain 1000 most variable genes
vars <- rowVars(as.matrix(log1p(counts(sce))))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
core <- sce[names(vars)[1:1000], ]

# split into target and background datasets
background_sce <- core[, which(core$cell_prov %in% c("healthy_1", "healthy_2"))]
patient_027_sce <- core[, which(core$cell_prov %in%
                                c("pre_trans_027", "post_trans_027"))]

background <- as.matrix(t(counts(background_sce)))
patient_027 <- as.matrix(t(counts(patient_027_sce)))

# perform scpca
snowparam <- SnowParam(workers = 30, type = "SOCK")
BiocParallel::register(snowparam, default = TRUE)
bmmc_scpca <- scPCA(patient_027, background,
                    penalties = exp(seq(log(1e-9), log(1), length.out = 20)),
                    center = TRUE, scale = TRUE, n_centers = 2, parallel = TRUE,
                    max_iter = 1000)

saveRDS(bmmc_scpca, file = here("analyses/bmmc_data/data/bmmc_scpca_027.rds"))