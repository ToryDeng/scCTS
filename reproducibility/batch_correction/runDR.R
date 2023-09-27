library(scater)
library(tidyverse)
library(patchwork)
library(lisi)
library(BiocParallel)
library(sva)
library(limma)
library(scMerge)
library(batchelor)



pbmc_lupus <- readRDS("~/scCTS/reproducibility/batch_correction/pbmc_lupus_bc.rds")
pbmc_covid <- readRDS("~/scCTS/reproducibility/batch_correction/pbmc_covid_bc.rds")
save.path <- "~/scCTS/reproducibility/batch_correction"
BPPARAM <- MulticoreParam(workers = 20)

for (bc.method in c('origin', 'limma', 'fastMNN', 'combat', 'scMerge', 'scMerge2')){
  # create name of assay that stores corrected, PCA and UMAP data
  expr.assay.name <- ifelse(bc.method == 'origin', 'logcounts', paste(bc.method, 'logcounts', sep = '_'))
  pca.assay.name <- paste(bc.method, 'pca', sep='_')
  umap.assay.name <- paste(bc.method, 'umap', sep='_')

  cat('Start to compute principal components of', expr.assay.name, '\n')
  pbmc_lupus <- runPCA(pbmc_lupus, assay.type=expr.assay.name, name=pca.assay.name, subset_row=rownames(pbmc_lupus))
  pbmc_covid <- runPCA(pbmc_covid, assay.type=expr.assay.name, name=pca.assay.name, subset_row=rownames(pbmc_covid))
  cat('Principal components have been computed.\n')

  cat('Start to compute UMAP of', expr.assay.name, '\n')
  pbmc_lupus <- runUMAP(pbmc_lupus, dimred=pca.assay.name, name=umap.assay.name, BPPARAM=BPPARAM)
  pbmc_covid <- runUMAP(pbmc_covid, dimred=pca.assay.name, name=umap.assay.name, BPPARAM=BPPARAM)
  cat('UMAP have been computed.\n')
}

saveRDS(pbmc_lupus, file.path(save.path, "pbmc_lupus_bc_dr.rds"))
saveRDS(pbmc_covid, file.path(save.path, "pbmc_covid_bc_dr.rds"))
