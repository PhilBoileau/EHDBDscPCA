# File for running scPCA on dengue data
library(GEOquery)
library(genefilter)
library(scPCA)
library(BiocParallel)

set.seed(871234)

# load the data, already log2 transformed
ges <- getGEO("GSE51808")$GSE51808_series_matrix.txt.gz

# keep the 50% most variable genes
var_filt_ges <- varFilter(ges, var.cutoff = 1-500/nrow(exprs(ges)))

# extract the target and background datasets
control_label <- which(var_filt_ges$`status:ch1` == "control")
target <- t(exprs(var_filt_ges)[, -control_label])
background <- t(exprs(var_filt_ges)[, control_label])

# get the target data labels
dengue_class <- var_filt_ges$`status:ch1`[-control_label]

# repeat scPCA for n_centers = 2,3,4,5
snowparam <- SnowParam(workers = 32, type = "SOCK")
BiocParallel::register(snowparam, default = TRUE)

diff_centers_scpca <- lapply(2:5, function(x) {
  
  # fit scPCA for different number of centers
  scPCA(target, background, center = TRUE, n_centers = x, parallel = TRUE,
        max_iter = 1000)
  
})

# save the results
save(diff_centers_scpca, file = "diff_centers_scpca.RData")
