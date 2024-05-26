

#' Normalize raw counts This function normalize (not log-normalize!) the raw
#' counts using a target sum of 1e4.
#'
#' @param X the raw count matrix (genes by cells).
#'
#' @return a matrix containing the normalized counts.
#' @importFrom Matrix t colSums
#' @importFrom cli cli_alert_info
#'
normalize.raw.counts <- function(X){
  cli_text("Normalizing the raw counts...")
  X <- t(t(X) / colSums(X)) * 1e4  # target sum is set to 1e4
  return(X)
}



create.null.matrix <- function(columns, rows){
  matrix.tmp <- matrix(NA, nrow=length(rows), ncol=length(columns))
  colnames(matrix.tmp) <- columns # column is unique cell types
  rownames(matrix.tmp) <- rows # row is gene names
  return(matrix.tmp)
}


create.res.list <- function(items, columns, rows){
  matrix.tmp <- create.null.matrix(columns, rows)
  subject.res <- list()
  for (item in items){
    subject.res[[item]] <- matrix.tmp
  }
  return(subject.res)
}


#' Get the expression matrix
#'
#' @param sce A SingleCellExperiment object
#' @param use.raw Whether use raw counts data. Default is \code{FALSE}. If
#'   \code{normalize=FALSE}, the function returns the raw counts directly. If
#'   \code{normalize=TRUE}, the function returns the normcounts.
#' @param use.norm.rep Which representation in \code{assayNames(sce)} will be
#'   used.
#' @param normalize Whether normalize the data. Default is \code{FALSE}. Only
#'   valid when \code{use.raw=TRUE}.
#'
#' @return A dense expression matrix (genes by cells).
#'
#' @import SummarizedExperiment
#' @importFrom methods is
get.expression.matrix <- function(
    sce,
    use.raw=FALSE,
    use.norm.rep=NULL,
    normalize=FALSE
    ){
  stopifnot(is(sce, 'SingleCellExperiment'))
  stopifnot(all(is.logical(c(use.raw, normalize))))
  stopifnot(is.null(use.norm.rep) | is.character(use.norm.rep))

  if (use.raw){  # retrieve the raw counts (and normalize them)
    stopifnot(is.null(use.norm.rep))
    Y <- counts(sce)
    cli_alert_info("Using the raw counts")
    if (normalize){
      Y <- normalize.raw.counts(Y)
    }
  }else{  # directly use the stored (log-)normalized counts
    if (is.null(use.norm.rep)){ # if not provide assay name
      stopifnot('normcounts' %in% assayNames(sce))
      Y <- normcounts(sce)
      cli_text("Using the assay 'normcounts'")
    }else{
      stopifnot(use.norm.rep %in% assayNames(sce))
      Y <- assay(sce, use.norm.rep)
      cli_text("Using the assay {.val {use.norm.rep}}")
    }
  }
  return(as.matrix(Y))
}


#' Inverse the log-transformation
#'
#' @param X a matrix containing the normalized counts (genes by cells).
#' @param log.base the base of log-transformation.
#'
#' @return the raw count matrix (genes by cells).
#'
#' @importFrom cli cli_text
#'
inverse.log <- function(X, log.base=c('natural', 10, 2)){
  stopifnot(log.base %in% c('natural', 10, 2))
  cli_text("Reverting the logcounts (base: {.val {log.base}})...")
  if (log.base == 'natural'){
    X <- expm1(X)
  }else{
    X <- log.base ^ X - 1
  }
  return(X)
}



#' Set parallel computation
#'
#' @param numCores the number of cores to use. Default is \code{NULL}, which
#'   means using all but one possible cores.
#'
#' @return the number of actually used cores.
#' @import doParallel
#' @importFrom parallel detectCores
#' @importFrom cli cli_alert_info

set.parallel.computation <- function(numCores = NULL){
  stopifnot(is.null(numCores) | is.numeric(numCores))

  if (is.null(numCores)){
    numCores <- detectCores() - 1
  }
  registerDoParallel(numCores) # registered cores number for parallel computation
  cli_alert_info("Using {.val {numCores}} workers for computation.")
  return(numCores)
}


store.test.res <- function(test.res, sub.res, names.map, current.ct){
  for (name in names(names.map)){
    name.in.sub.res = names.map[[name]]
    sub.res[[name.in.sub.res]][,current.ct] = test.res[,name]
  }
  return(sub.res)
}



check.raw.counts <- function(expr){
  return(ifelse(all(as.matrix(expr) %% 1 == 0), TRUE, FALSE))
}


#' Perform a stratified sampling on the cells from specified cell types
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param sub.rep The name of subject labels.
#' @param ct.rep The name of cell type labels.
#' @param fraction The sampling fraction between 0 and 1.
#' @param subset.cts Cell types that require sampling. Must be a subset of \code{unique(colData(sce)[[ct.rep]])}.
#'
#'
#' @return the sampling \code{SingleCellExperiment} object.
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom cli cli_text
#'
stratified.sampling <- function(sce, sub.rep, ct.rep, fraction, subset.cts){
  celltypes <- colData(sce)[[ct.rep]]
  stopifnot(subset.cts %in% unique(celltypes))
  stopifnot(between(fraction, 0, 1))
  cli_text("Dataset size before sampling: {.var {dim(sce)}}")
  cli_text("Cell types to sample cells: {.var {subset.cts}}")
  cli_text("Number of cells in specified types before sampling:
           {.var {table(celltypes)[subset.cts]}}")

  subset_cells <- colData(sce) %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    dplyr::filter(get(ct.rep) %in% subset.cts) %>%
    dplyr::group_by(get(sub.rep), get(ct.rep)) %>%
    dplyr::sample_frac(size=fraction) %>%
    dplyr::pull(get("rowname"))
  cells.to.select <- !(celltypes %in% subset.cts) | colnames(sce) %in% subset_cells
  sce <- sce[, cells.to.select]

  cli_text("Dataset size after sampling: {.var {dim(sce)}}")
  cli_text("Number of cells in specified types before sampling:
         {.var {table(colData(sce)[[ct.rep]])[subset.cts]}}")
  return(sce)
}


get.variable.name <- function(x){
  return(deparse1(substitute(x)))
}


