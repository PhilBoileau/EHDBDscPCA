################################################################################
############################# SIMLR on BMMCs Data ##########################
################################################################################

# AML035

################################################################################

library(here)
library(tidyverse)
library(SIMLR)
library(SingleCellExperiment)
library(parallel)

# detect the number of cores
num_cores <- detectCores()

# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(12431)

# retain genes that contain 5 or more non-zero counts accross all cells
sce <- sce[(Matrix::rowSums(counts(sce) != 0) > 4), ]

# retain 1000 most variable genes
vars <- rowVars(as.matrix(log1p(counts(sce))))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
core <- sce[names(vars)[1:1000], ]

# consider only the taret dataset
patient_035_sce <- core[, which(core$cell_prov %in%
                                  c("pre_trans_035", "post_trans_035"))]
patient_035 <- as.matrix(t(counts(patient_027_sce)))

# perform zinbwave
bmmc_simlr <- SIMLR(X = scale(patient_035), c = 2, no.dim = 2,
                    cores.ratio = num_cores)

saveRDS(bmmc_simlr,
        file = here("analyses/bmmc_data/data/bmmc_simlr_035.rds"))