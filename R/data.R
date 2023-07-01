#' A simulated scRNA-seq dataset
#'
#' The code for generating it:
#' \code{set.seed(123)}
#' \code{nGenes <- 2000}
#' \code{nCells <- 10000}
#' \code{rawCounts <- matrix(rnbinom(nGenes * nCells, mu=1, size=10),nGenes,nCells)}
#' \code{rownames(rawCounts) <- paste0("Gene", 1:nGenes)}
#' \code{colnames(rawCounts) <- paste0("Cell", 1:nCells)}
#' \code{celltypes <- paste0("celltype", sample.int(5, nCells, replace = T))}
#' \code{subjects <- paste0("subject", sample.int(10, nCells, replace = T))}
#' \code{cellAnnot <- data.frame(celltype=celltypes, subject=subjects)}
#' \code{sim.sce <- SingleCellExperiment(assays=list(counts=raw.counts), colData=cell.annot)}
#'
#'
#' @format ## `sim.sce`
#' This dataset contains 2k genes and 10k cells.
#'
#' @source
#' This dataset is generated in R.
#'
"sim.sce"
