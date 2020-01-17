################################################################################
############################# scPCA on BMMCs Data ##############################
################################################################################

# AML035

################################################################################

library(here)
library(tidyverse)
library(scPCA)
library(SingleCellExperiment)
library(BiocParallel)

# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(573322)

# retain genes that contain 5 or more non-zero counts accross all cells
sce <- sce[(Matrix::rowSums(counts(sce) != 0) > 4), ]

# retain 1000 most variable genes
vars <- rowVars(as.matrix(log1p(counts(sce))))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
core <- sce[names(vars)[1:1000], ]

# split into target and background datasets
background_sce <- core[, which(core$cell_prov %in% c("healthy_1", "healthy_2"))]
patient_035_sce <- core[, which(core$cell_prov %in%
                                  c("pre_trans_035", "post_trans_035"))]

background <- as.matrix(t(counts(background_sce)))
patient_035 <- as.matrix(t(counts(patient_035_sce)))

# perform scpca
snowparam <- SnowParam(workers = 8, type = "SOCK")
BiocParallel::register(snowparam, default = TRUE)
bmmc_scpca <- scPCA(patient_035, background, penalties = 0, center = TRUE,
                    scale = TRUE, n_centers = 2, parallel = TRUE, cv = 5,
                    max_iter = 1000)

saveRDS(bmmc_scpca, file = here("analyses/bmmc_data/data/bmmc_cpca_cv_035.rds"))