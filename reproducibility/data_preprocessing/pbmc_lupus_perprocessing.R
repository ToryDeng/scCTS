suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scMerge))
suppressPackageStartupMessages(library(tidyverse))
source("~/scCTS/R/utils.R")

load.path <- "/data/scRNA-seq/SLE_GSE96583/raw/GSE96583_RAW"
save.path <- "~/scCTS/reproducibility/data_preprocessing"

# load count matrices
sample.ids <- c('GSM2560245','GSM2560246', 'GSM2560247','GSM2560248','GSM2560249')
sample.list <- list()
for (sample.id in sample.ids){
  sample.data.path <- file.path(load.path, sample.id)
  sample.list[[sample.id]] <- Read10X(sample.data.path)
  cat('Accession number: ', sample.id,
      '       Dimension: ', dim(sample.list[[sample.id]]), '\n')
}


# genes in batch1 are the same, and genes in batch2 are the same,
# but genes in different batches are not the same
all(rownames(sample.list$GSM2560245) == rownames(sample.list$GSM2560246) &
      rownames(sample.list$GSM2560246) == rownames(sample.list$GSM2560247))
all(rownames(sample.list$GSM2560248) == rownames(sample.list$GSM2560249))


# concatenate expression matrices
batch1.expression <- cbind(sample.list$GSM2560245, sample.list$GSM2560246, sample.list$GSM2560247)
batch2.expression <- cbind(sample.list$GSM2560248, sample.list$GSM2560249)

# duplicated genes
cat('Number of duplicated genes in batch1: ', sum(duplicated(rownames(batch1.expression))), '\n')
cat('Number of duplicated genes in batch2: ', sum(duplicated(rownames(batch2.expression))), '\n')

# read cell-level information
batch1.cell.info <- read_tsv("~/cts_gene/data/GSE96583/GSE96583_batch1.total.tsne.df.tsv",
                             skip = 1, col_names=c('cell.id', 'tsne1', 'tsne2', 'batch', 'cluster', 'multiplets', 'cell.type', 'ind'),
                             col_types = cols(ind = col_character())) %>%
  column_to_rownames('cell.id')

batch2.cell.info <- read_tsv("~/cts_gene/data/GSE96583/GSE96583_batch2.total.tsne.df.tsv",
                             skip = 1, col_names=c('cell.id', 'tsne1',	'tsne2', 'ind', 'stim', 'cluster', 'cell.type', 'multiplets'),
                             col_types = cols(ind = col_character())) %>%
  column_to_rownames('cell.id')

# a suffix "1" in mismatched cells in `rownames(batch.cell.info)`
# batch1.is_mismatched_cells <- colnames(batch1.expression) != rownames(batch1.cell.info)
# colnames(batch1.expression)[batch1.is_mismatched_cells] <- str_c(colnames(batch1.expression)[batch1.is_mismatched_cells], 1)
# sum(colnames(batch1.expression) != rownames(batch1.cell.info))

# batch2.is_mismatched_cells <- colnames(batch2.expression) != rownames(batch2.cell.info)
# colnames(batch2.expression)[batch2.is_mismatched_cells] <- str_c(colnames(batch2.expression)[batch2.is_mismatched_cells], 1)
# sum(colnames(batch2.expression) != rownames(batch2.cell.info))



batch1.cell.info %<>%
  dplyr::filter(multiplets == 'singlet',
                cell.type != 'Megakaryocytes',
                rownames(.) == colnames(batch1.expression)) %>%
  tidyr::drop_na(c('cell.type', 'ind')) %>%
  dplyr::filter() %>%
  dplyr::rename(subject = ind, celltype=cell.type) %>%
  dplyr::select(c('subject', 'celltype'))

batch2.cell.info %<>%
  dplyr::filter(multiplets == 'singlet',
                cell.type != 'Megakaryocytes',
                rownames(.) == colnames(batch2.expression)) %>%
  tidyr::drop_na(c('cell.type', 'ind')) %>%
  dplyr::rename(celltype=cell.type) %>%
  tidyr::unite('subject', c('ind', 'stim')) %>%
  dplyr::select(c('subject', 'celltype'))


# create SingleCellExperiment objects
sce.batch1 <- SingleCellExperiment(
  assays = list(counts=batch1.expression[,rownames(batch1.cell.info)]),
  colData = batch1.cell.info
)

sce.batch2 <- SingleCellExperiment(
  assays = list(counts=batch2.expression[,rownames(batch2.cell.info)]),
  colData = batch2.cell.info
)

# create sce list for subjects
sce.subjects <- list()
for (sub in unique(colData(sce.batch1)[['subject']])){
  sce.subjects[[sub]] <- sce.batch1[,colData(sce.batch1)[['subject']] == sub]
}
for (sub in unique(colData(sce.batch2)[['subject']])){
  sce.subjects[[sub]] <- sce.batch2[,colData(sce.batch2)[['subject']] == sub]
}


# clean memory
# rm(list=setdiff(ls(), c('sce.subjects', 'normalize.raw.counts')))

# filter low expressed genes
for (sub in names(sce.subjects)){
  cat('Current subject: ', sub, '\n')
  sce <-sce.subjects[[sub]]
  sce.subjects[[sub]] <- sce[nexprs(sce, byrow=TRUE) > 2,]
}


# merge subjects
pbmc.lupus <- sce_cbind(sce.subjects, cut_off_batch=0, cut_off_overall=0,
                  exprs=c('counts'), colData_names=c('subject', 'celltype'))
colData(pbmc.lupus) <- colData(pbmc.lupus)[c('subject', 'celltype')]
# normalize
normcounts(pbmc.lupus) <- normalize.raw.counts(counts(pbmc.lupus))
logcounts(pbmc.lupus) <- as(log2(normcounts(pbmc.lupus) + 1), "CsparseMatrix")
# make cell names unique
colnames(pbmc.lupus) <- make.unique(colnames(pbmc.lupus))
# save the dataset
cat('Dimension of the final dataset: ', dim(pbmc.lupus), '\n')
saveRDS(pbmc.lupus, file.path(save.path, 'pbmc_lupus.rds'))



