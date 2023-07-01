

#' Normalize raw counts
#' This function normalize (not log-normalize!) the raw counts using a target sum of 1e4.
#'
#' @param X the raw count matrix (genes by cells).
#'
#' @return a matrix containing the normalized counts.
#' @importFrom Matrix t colSums

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
#' @param use.raw Whether use raw counts data. Default is `FALSE`. If `normalize=FALSE`, the function
#' returns the raw counts directly. If `normalize=TRUE`, the function returns the
#' normcounts.
#' @param use.norm.rep Which representation in \code{assayNames(sce)} will be used.
#' @param normalize Whether normalize the data. Default is `FALSE`. Only valid when `use.raw=TRUE`.
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

  if (use.raw){
    Y <- counts(sce)
    cli_text("Using the raw counts")
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
#' @param numCores the number of cores to use. Default is \code{NULL}, which means
#' using all but one possible cores.
#'
#' @return the number of actually used cores.
#' @import doParallel
#' @importFrom parallel detectCores
#' @importFrom cli cli_text

set.parallel.computation <- function(numCores = NULL){
  stopifnot(is.null(numCores) | is.numeric(numCores))

  if (is.null(numCores)){
    numCores <- detectCores() - 1
  }
  registerDoParallel(numCores) # registered cores number for parallel computation
  cli_text("Using {.val {numCores}} workers for computation.")
  return(numCores)
}


store.test.res <- function(test.res, sub.res, names.map, current.ct){
  for (name in names(names.map)){
    name.in.sub.res = names.map[[name]]
    sub.res[[name.in.sub.res]][,current.ct] = test.res[,name]
  }
  return(sub.res)
}









