# File for running cpca on dengue data
library(GEOquery)
library(genefilter)
library(scPCA)

set.seed(12312)

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

# run cpca on the expression data
dengue_cpca_cv <- scPCA(target, background, center = TRUE,
                        penalties = 0, n_centers = 3, cv = 3)

# save the results
save(dengue_cpca_cv, file = "cpca_cv.RData")