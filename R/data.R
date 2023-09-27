#' A simulated scRNA-seq dataset
#'
#' The counts are sampled from negative binomial distribution.
#'
#'
#' @docType data
#'
#' @format
#' This dataset contains 2k genes and 10k cells in a \code{SingleCellExperiment} object.
#'
#' @usage data("sim.sce")
#'
#' @details
#'
#' The code for generating it:
#' ```{r}
#' set.seed(123)
#' nGenes <- 2000
#' nCells <- 10000
#' rawCounts <- matrix(rnbinom(nGenes * nCells, mu=1, size=10),nGenes,nCells)
#' rownames(rawCounts) <- paste0("Gene", 1:nGenes)
#' colnames(rawCounts) <- paste0("Cell", 1:nCells)
#' celltypes <- paste0("celltype", sample.int(5, nCells, replace = T))
#' subjects <- paste0("subject", sample.int(10, nCells, replace = T))
#' cellAnnot <- data.frame(celltype=celltypes, subject=subjects)
#' sim.sce <- SingleCellExperiment(assays=list(counts=rawCounts), colData=cellAnnot)
#' ```
#' @md
"sim.sce"
