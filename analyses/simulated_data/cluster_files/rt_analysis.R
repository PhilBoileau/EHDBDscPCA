# scPCA with CV tuned hyperparameters

library(here)
library(tidyverse)
library(scPCA)
library(Rtsne)
library(umap)
library(SingleCellExperiment)
library(splatter)
library(scater)
library(zinbwave)
library(SIMLR)
library(cluster)
library(microbenchmark)


# generate the dataset
params <- newSplatParams(
  seed = 6757293,
  nGenes = 500,
  batchCells = c(150, 150),
  batch.facLoc = c(0.05, 0.05),
  batch.facScale = c(0.05, 0.05),
  group.prob = rep(1/3, 3),
  de.prob = c(0.1, 0.05, 0.1),
  de.downProb = c(0.1, 0.05, 0.1),
  de.facLoc = rep(0.2, 3),
  de.facScale = rep(0.2, 3)
)
sim_groups_sce <- splatSimulate(params, method = "groups")

# get the logcounts of the data
sim_groups_sce <- normalize(sim_groups_sce)

# remove all cells without variation in counts
sim_groups_sce <- sim_groups_sce[which(rowVars(counts(sim_groups_sce)) != 0), ]

# get the logcount matrices
tg <- t(logcounts(sim_groups_sce)[, which(sim_groups_sce$Group != "Group2")])
bg <- t(logcounts(sim_groups_sce)[, which(sim_groups_sce$Group == "Group2")])

set.seed(412321)

# start the microbenchmarking
snowparam <- SnowParam(workers = 4, type = "SOCK")
BiocParallel::register(snowparam, default = TRUE)

sim_microbench <- microbenchmark(
  "PCA" = prcomp(tg, center = TRUE, scale. = TRUE),
  "cPCA" = scPCA(tg, bg, center = TRUE, scale = TRUE,
                 penalties = 0, n_centers = 2, parallel = TRUE),
  "scPCA (iterative)" = scPCA(tg, bg, center = TRUE, scale = TRUE,
                              n_centers = 2, parallel = TRUE, alg = "iterative"),
  "scPCA (var. proj.)" = scPCA(tg, bg, center = TRUE, scale = TRUE,
                               n_centers = 2, parallel = TRUE, alg = "var_proj"),
  "scPCA (rand. var. proj.)" = scPCA(tg, bg, center = TRUE, scale = TRUE,
                                     n_centers = 2, parallel = TRUE,
                                     alg = "rand_var_proj"),
  "t-SNE (PCA)" = Rtsne(scale(tg), theta = 0, perplexity = 30, max_iter = 1000,
                        pca = TRUE, num_threads = 4),
  "t-SNE" = Rtsne(scale(tg), theta = 0, perplexity = 30, max_iter = 1000,
                  pca = FALSE, num_threads = 4),
  "UMAP" = umap(scale(tg), n_neighbors = 30, min_dist = 0.02),
  "ZINB-WaVE (no batch)" = zinbwave(
    sim_groups_sce[, which(sim_groups_sce$Group != "Group2")],
    K = 2
  ),
  "ZINB-WaVE (batch)" = zinbwave(
    sim_groups_sce[, which(sim_groups_sce$Group != "Group2")],
    K = 2,
    X = "~Batch"
  ),
  "SIMLR" = SIMLR(X = t(scale(tg)), c = 2, no.dim = 2, normalize = FALSE,
                  cores.ratio = 4),
  times = 5,
  control = list(order = "block")
)

# save the file for later
write_rds(sim_microbench,
          path = here("analyses/simulated_data/cluster_files/sim_microbench.rds"))
