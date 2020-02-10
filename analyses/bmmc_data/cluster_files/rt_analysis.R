# scPCA with CV tuned hyperparameters

library(here)
library(tidyverse)
library(scPCA)
library(Rtsne)
library(umap)
library(SingleCellExperiment)
library(zinbwave)
library(SIMLR)
library(cluster)
library(microbenchmark)


# load the data
source(file = here("analyses/bmmc_data/helpers/load_count_data.R"))

# ensure reproducibility
set.seed(31231)

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

patient_035_se <- SummarizedExperiment(
  list(counts = as.matrix(counts(patient_035_sce)))
)

# start the microbenchmarking
snowparam <- SnowParam(workers = 4, type = "SOCK")
BiocParallel::register(snowparam, default = TRUE)

aml035_microbench <- microbenchmark(
  "PCA" = prcomp(patient_035, center = TRUE, scale. = TRUE),
  "cPCA" = scPCA(patient_035, background,
                 center = TRUE, scale = TRUE, penalties = 0, n_centers = 2,
                 max_iter = 1000),
  "scPCA (iterative)" = scPCA(patient_035, background,
                              penalties = exp(seq(log(1e-9), log(1), length.out = 20)),
                              center = TRUE, scale = TRUE, n_centers = 2,
                              parallel = TRUE, max_iter = 1000,
                              alg = "iterative"),
  "scPCA (var. proj.)" = scPCA(patient_035, background,
                               penalties = exp(seq(log(1e-9), log(1), length.out = 20)),
                               center = TRUE, scale = TRUE, n_centers = 2,
                               parallel = TRUE, max_iter = 1000,
                               alg = "var_proj"),
  "scPCA (rand. var. proj.)" = scPCA(patient_035, background,
                                     penalties = exp(seq(log(1e-9), log(1), length.out = 20)),
                                     center = TRUE, scale = TRUE, n_centers = 2,
                                     parallel = TRUE, max_iter = 1000,
                                     alg = "rand_var_proj"),
  "t-SNE (PCA)" = Rtsne(patient_035, perplexity = 30, max_iter = 1000,
                        pca = TRUE, theta = 0, num_threads = 4),
  "t-SNE" = Rtsne(patient_035, perplexity = 30, max_iter = 1000,
                  pca = FALSE, theta = 0, num_threads = 4),
  "UMAP" = umap(patient_035, n_neighors = 30, min_dist = 0.02),
  "ZINB-WaVE" = zinbwave(
    patient_035_se,
    K = 2
  ),
  "SIMLR" = SIMLR(X = scale(patient_035), c = 2, no.dim = 2, k = 30,
                  normalize = FALSE),
  times = 3,
  control = list(order = "block", warmup = 0)
)

# save the file for later
write_rds(aml035_microbench,
          path = here("analyses/bmmc_data/data/aml035_microbench.rds"))
