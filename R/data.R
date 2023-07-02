#' A simulated scRNA-seq dataset
#'
#' The code for generating it:
#' ```
#' set.seed(123)
#' nGenes <- 2000
#' nCells <- 10000
#' ```
#'
#' \code{rawCounts <- matrix(rnbinom(nGenes * nCells, mu=1, size=10),nGenes,nCells)}
#' \code{rownames(rawCounts) <- paste0("Gene", 1:nGenes)}
#' \code{colnames(rawCounts) <- paste0("Cell", 1:nCells)}
#' \code{celltypes <- paste0("celltype", sample.int(5, nCells, replace = T))}
#' \code{subjects <- paste0("subject", sample.int(10, nCells, replace = T))}
#' \code{cellAnnot <- data.frame(celltype=celltypes, subject=subjects)}
#' \code{sim.sce <- SingleCellExperiment(assays=list(counts=raw.counts), colData=cell.annot)}
#'
#'
#' @format
#' This dataset contains 2k genes and 10k cells.
#'
#' @usage
#' library(popuDE)
#' data("sim.sce")
"sim.sce"
