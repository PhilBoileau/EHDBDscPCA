# scPCA with CV tuned hyperparameters

library(here)
library(tidyverse)
library(scPCA)
library(SingleCellExperiment)
library(splatter)
library(scater)


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

# get the batch and group labels
target_batch <- sim_groups_sce$Batch[which(sim_groups_sce$Group != "Group2")]
target_group <- sim_groups_sce$Group[which(sim_groups_sce$Group != "Group2")]
batch_mem <- if_else(target_batch == "Batch1", 1, 2)
group_mem <- if_else(target_group == "Group1", 1, 2)

set.seed(412321)

# perform scpca, tuning hyperparameters through 5-fold cv
sim_scpca_cv <- scPCA(tg, bg, center = TRUE, scale = TRUE, n_centers = 2)

# save the file for later
write_rds(sim_scpca_cv,
          path = here("analyses/simulations/sim_scpca_cv.rds"))
