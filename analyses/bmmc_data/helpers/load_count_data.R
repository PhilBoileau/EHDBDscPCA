# Load the .mtx data into .Rmd
pre_trans_027 <- Matrix::readMM(
  file = here("analyses/bmmc_data/data/AML027_pre/matrix.mtx")
)
pre_trans_035 <- Matrix::readMM(
  file = here("analyses/bmmc_data/data/AML035_pre/matrix.mtx")
)
post_trans_027 <- Matrix::readMM(
  file = here("analyses/bmmc_data/data/AML027_post/matrix.mtx")
)
post_trans_035 <- Matrix::readMM(
  file = here("analyses/bmmc_data/data/AML035_post/matrix.mtx")
)
healthy_1 <- Matrix::readMM(
  file = here("analyses/bmmc_data/data/healthy_1/matrix.mtx")
)
healthy_2 <- Matrix::readMM(
  file = here("analyses/bmmc_data/data/healthy_2/matrix.mtx")
)

# create a single cell experiment object #######################################

# combine the count matrices 
count_matrix <- cbind(
  pre_trans_027,
  pre_trans_035,
  post_trans_027,
  post_trans_035,
  healthy_1,
  healthy_2
)

# get the providence of each cell
prov_labels <- c(
  rep("pre_trans_027", ncol(pre_trans_027)),
  rep("pre_trans_035", ncol(pre_trans_035)),
  rep("post_trans_027", ncol(post_trans_027)),
  rep("post_trans_035", ncol(post_trans_035)),
  rep("healthy_1", ncol(healthy_1)),
  rep("healthy_2", ncol(healthy_2))
)
prov_df <- data.frame(cell_prov = prov_labels)

# get the gene names
gene_names_df <- read_tsv(
  file = here("analyses/bmmc_data/data/AML027_pre/genes.tsv"),
  col_names = FALSE
) %>%
  as.data.frame()

rownames(gene_names_df) <- gene_names_df[, 1]
gene_names_df[, 1] <- NULL

# create the sce object
sce <- SingleCellExperiment(
  assays = list(counts = count_matrix),
  colData = prov_df,
  rowData = gene_names_df
)
###############################################################################
# remove leftover objects
rm(pre_trans_027, pre_trans_035, post_trans_027, post_trans_035, healthy_1, healthy_2,
   count_matrix, prov_labels, prov_df, gene_names_df)
