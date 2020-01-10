################################################################################
############################# ZINB-WaVE on BMMCs Data ###########################
################################################################################

# AML035

################################################################################

library(here)
library(tidyverse)
library(zinbwave)
library(SingleCellExperiment)
library(BiocParallel)

# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(65928)

# retain genes that contain 5 or more non-zero counts accross all cells
sce <- sce[(Matrix::rowSums(counts(sce) != 0) > 4), ]

# retain 1000 most variable genes
vars <- rowVars(as.matrix(log1p(counts(sce))))
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
core <- sce[names(vars)[1:1000], ]

# consider only the background dataset
patient_035_sce <- core[, which(core$cell_prov %in%
                                c("pre_trans_035", "post_trans_035"))]
patient_035_se <- SummarizedExperiment(
  list(counts = as.matrix(counts(patient_035_sce)))
)

# perform zinbwave
multicore <- MulticoreParam(workers = 30)
BiocParallel::register(multicore, default = TRUE)
bmmc_zinbwave_sce <- zinbwave(
  patient_035_se,
  K = 2,
  verbose = TRUE
)

saveRDS(bmmc_zinbwave_sce,
        file = here("analyses/bmmc_data/data/zinbwave_035.rds"))