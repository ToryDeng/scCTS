reticulate::use_python(required = T, python = '/usr/local/anaconda3/envs/bio/bin/')
suppressPackageStartupMessages(library(anndata))

load.path <- "~/scCTS/reproducibility/data_preprocessing/covid19_GSE158055_subset100k.h5ad"
save.path <- "~/scCTS/reproducibility/data_preprocessing"
ad <- read_h5ad(load.path)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(scMerge))
suppressPackageStartupMessages(library(scater))
source("~/scCTS/R/utils.R")




pbmc <- SingleCellExperiment(
  assays = list(counts=Matrix::t(ad$X)),
  colData = ad$obs
)

sce.subjects <- list()
for (sub in unique(colData(pbmc)[['sampleID']])){
  sce.subjects[[sub]] <- pbmc[,colData(pbmc)[['sampleID']] == sub]
}

# filter low expressed genes
for (sub in names(sce.subjects)){
  cat('Current subject: ', sub, '\n')
  sce <-sce.subjects[[sub]]
  sce.subjects[[sub]] <- sce[nexprs(sce, byrow=TRUE) > 10,]
}


# merge subjects
pbmc.covid19 <- sce_cbind(sce.subjects, cut_off_batch=0, cut_off_overall=0,
                          exprs=c('counts'), colData_names=names(colData(pbmc)))
colData(pbmc.covid19) <- colData(pbmc.covid19)[names(colData(pbmc))]

# make cell names unique
colnames(pbmc.covid19) <- make.unique(colnames(pbmc.covid19))

# normalize
normcounts(pbmc.covid19) <- normalize.raw.counts(counts(pbmc.covid19))
logcounts(pbmc.covid19) <- as(log2(normcounts(pbmc.covid19) + 1), "CsparseMatrix")

# batch-effects correction
# BPPARAM <- MulticoreParam(workers = 20)
# data("segList", package = "scMerge")
# scMerge_supervised <- scMerge(sce_combine = pbmc.covid19,
#                               ctl = segList$human$human_scSEG,
#                               assay_name = "scMerge_supervised",
#                               batch_name = 'sampleID',
#                               cell_type = colData(pbmc.covid19)[['majorType']],
#                               BPPARAM = BPPARAM)
# assay(pbmc.covid19, "scMerge_logcounts") <- assay(scMerge_supervised, "scMerge_supervised")


# save the dataset
cat('Dimension of the final dataset: ', dim(pbmc.covid19), '\n')
saveRDS(pbmc.covid19, file.path(save.path, "pbmc_covid.rds"))




