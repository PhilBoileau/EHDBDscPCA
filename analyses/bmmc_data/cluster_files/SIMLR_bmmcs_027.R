################################################################################
############################# SIMLR on BMMCs Data ##########################
################################################################################

# AML027

################################################################################

library(here)
library(tidyverse)
library(SIMLR)
library(SingleCellExperiment)

# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(235423)

# retain genes that contain 5 or more non-zero counts accross all cells
sce <- sce[(Matrix::rowSums(counts(sce) != 0) > 4), ]

# retain 1000 most variable genes
vars <- rowVars(as.matrix(log1p(counts(sce))))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
core <- sce[names(vars)[1:1000], ]

# consider only the taret dataset
patient_027_sce <- core[, which(core$cell_prov %in%
                                  c("pre_trans_027", "post_trans_027"))]
patient_027 <- as.matrix(counts(patient_027_sce))

# perform zinbwave
bmmc_simlr <- SIMLR(X = scale(patient_027), c = 2, no.dim = 2, k = 30)

saveRDS(bmmc_simlr,
        file = here("analyses/bmmc_data/data/bmmc_simlr_027.rds"))
