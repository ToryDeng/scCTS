
#' A wrapper function of the classic methods to find over-expressed DE genes
#'
#' @param sce A \code{SingleCellExperiment} object. Should contain normalized count matrix,
#' subject and cell type info.
#' @param use.raw Whether to use the raw count matrix, which is normalized internally.
#' @param use.norm.rep Which count matrix should be used. Default is the matrix
#' accessed by \code{normcounts()}.
#' @param normalize Whether to normalize the raw counts. Only relevant when \code{use.raw=TRUE}.
#' @param log.input Whether the input expression matrix is log-transformed.
#' @param log.base The base of log-transformation.
#' @param subject.rep Column name for the subject info.
#' @param celltype.rep Column name for the cell type info.
#' @param per.subject If \code{TRUE}, perform DE analysis for each sample.
#' Otherwise perform DE analysis for all cells.
#' @param method  A string indicating which DE method to use.
#' @param numCores Number of cores to use. Default is \code{NULL}.
#' @param python.path The path to the python. Need to install \pkg{nsforest} in advance.
#'
#' @return A list of 3-dim arrays. Each array is corresponding to a stat. The first dim
#' is for genes, the second dim is for cell types, and the last dim is for subjects.
#' @import SingleCellExperiment
#' @importFrom abind abind
#' @importFrom stringr str_glue
#' @export
#'
#' @examples
#' # load the simulated data
#' data(sim.sce)
#'
#' # run the Wilcoxon test per subject
#' res <- runBaselineMethod(sim.sce, use.raw = TRUE, normalize=TRUE,
#'                          subject.rep='subject', celltype.rep='celltype',
#'                          per.subject=TRUE, method='wilcox', numCores=2)
#'
runBaselineMethod <- function(
    sce,
    use.raw=FALSE,
    use.norm.rep=NULL,
    normalize=FALSE,
    log.input=FALSE,
    log.base=2,
    subject.rep='subject',
    celltype.rep='celltype',
    per.subject=TRUE,
    method=c('wilcox', 'twelch', 'DEseq2', 'NSforest'),
    numCores=NULL,
    python.path=NULL
){
  # check arguments
  stopifnot(all(c(subject.rep, celltype.rep) %in% names(colData(sce))),
            is(per.subject, "logical"))
  # set parallel computation
  numCores.used <- set.parallel.computation(numCores)
  # get normalized count matrix
  Y <- get.expression.matrix(sce, use.raw, use.norm.rep, normalize)
  # get (unique) cell types
  celltypes <- colData(sce)[[celltype.rep]] # convert to char list

  if (per.subject){
    cli_h1("Subject-level {.emph {method}} method")

    subjects <- colData(sce)[[subject.rep]] # convert to char list
    unique.subjects <- sort(unique(subjects))

    sub.res.list <- list()
    cli_progress_bar("Analyzing each subject", total = length(unique.subjects), type = "tasks")
    for (sub in unique.subjects){
      # cli_h2("Current subject: {sub}")
      sub.Y <- Y[,subjects == sub]
      sub.cts <- celltypes[subjects == sub]
      sub.result = switch(
        method,
        "wilcox" = BaselineMethod.wilcox(sub.Y, sub.cts, numCores.used),
        "twelch" = BaselineMethod.twelch(sub.Y, sub.cts, numCores.used, log.input, log.base),
        "DEseq2" = BaselineMethod.DEseq2(sub.Y, sub.cts, numCores.used),
        "NSforest" = BaselineMethod.NSforest(sub.Y, sub.cts, python.path),
        stop(str_glue("No method matched for {method}"))
      )
      # set the name of the last dim of each array as subject name, and then store
      sub.res.list[[sub]] <- lapply(sub.result, function(arr){dimnames(arr)[3] <- sub;arr})
      cli_progress_update()
    }
    # combine results for multiple subjects
    DE.res <- list()
    for (name in names(sub.res.list[[1]])){
      DE.res[[name]] <- abind(lapply(sub.res.list, function(lst) lst[[name]]))
    }
    return(DE.res)
  }else{
    cli_h1("Population-level {.emph {method}} method")
    cli_text("{.emph Note: this mode needs batch-effects correction in advance.}")
    subjects <- colData(sce)[[subject.rep]]
    all.result = switch(
      method,
      "wilcox" = BaselineMethod.wilcox(Y, celltypes, numCores.used),
      "twelch" = BaselineMethod.twelch(Y, celltypes, numCores.used, log.input, log.base),
      "DEseq2" = BaselineMethod.DEseq2(Y, celltypes, numCores.used),
      "NSforest" = BaselineMethod.NSforest(Y, celltypes, python.path),
      stop(str_glue("No method matched for {method}"))
    )
    all.result <- lapply(all.result, function(arr){dimnames(arr)[3] <- "all";arr})
    return(all.result)
  }
}



#' Naive Wilcoxon test
#'
#' @param expr A gene by cell matrix storing the expression values
#' @param celltypes A vector indicating cell types of each cell
#' @param nCores.used The number of cores actually used
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#' is genes, the second dim is cell types, and the last dim is a single subject.
#' @importFrom matrixTests row_wilcoxon_twosample
#' @importFrom foreach %dopar% foreach
#' @importFrom data.table rbindlist
#' @importFrom iterators iter
#' @importFrom parallel splitIndices
#' @importFrom stats p.adjust
#'
BaselineMethod.wilcox <- function(expr, celltypes, nCores.used){
  # both lognorm and norm data are acceptable for wilcoxon test
  unique.celltypes <- sort(unique(celltypes))

  # create 3-dim empty arrays to store results for a single subject (the last dim is 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(unique.celltypes), 1),
                     dimnames = list(rownames(expr), unique.celltypes))
  wilcox.res <- list('wilcox.stat_info' = array.tmp,
                     'wilcox.pval_info' = array.tmp,
                     'wilcox.fdr_info' = array.tmp)

  # compare each cell type with other cells
  for(ucelltype in unique.celltypes){
    # cli_text("Comparing {.val {ucelltype}} with other cells...")
    is.current.celltype <- (celltypes == ucelltype)
    not.current.celltype <- !is.current.celltype
    # split genes into chunks
    idx <- NULL  # to prevent "no visible binding for global variable" when checking
    chunks.res.list <- foreach(idx = iter(splitIndices(nrow(expr), nCores.used))) %dopar%{
      chunk.res <- row_wilcoxon_twosample(
        expr[idx,is.current.celltype], expr[idx,not.current.celltype], alternative = 'greater'
      )
      chunk.res[,c('statistic','pvalue')]
    }
    ct.res <- rbindlist(chunks.res.list)
    wilcox.res$wilcox.stat_info[,ucelltype,1] <- ct.res[['statistic']]
    wilcox.res$wilcox.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    wilcox.res$wilcox.fdr_info[,ucelltype,1] <- p.adjust(ct.res[['pvalue']], 'fdr')
  }
  return(wilcox.res)
}

#' t-test
#'
#' @param expr A gene by cell matrix storing the expression values
#' @param celltypes A vector indicating cell types of each cell
#' @param nCores.used The number of cores actually used
#' @param log.input Whether the input expression matrix is log-transformed.
#' @param log.base The base of log-transformation.
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#' is genes, the second dim is cell types, and the last dim is a single subject.
#' @importFrom matrixTests row_t_welch
#' @importFrom foreach %dopar% foreach
#' @importFrom data.table rbindlist
#' @importFrom iterators iter
#' @importFrom parallel splitIndices
#' @importFrom stats p.adjust
#'
#'
BaselineMethod.twelch <- function(expr, celltypes, nCores.used, log.input, log.base){
  if (log.input){expr <- inverse.log(expr, log.base)}
  unique.celltypes <- sort(unique(celltypes))
  # create 3-dim empty arrays to store results for a single subject (i.e, the last dim is 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(unique.celltypes), 1),
                     dimnames = list(rownames(expr), unique.celltypes))
  twelch.res <- list('twelch.stat_info' = array.tmp,
                     'twelch.pval_info' = array.tmp,
                     'twelch.fdr_info' = array.tmp)

  # compare each cell type with other cells
  for(ucelltype in unique.celltypes){
    # cli_text("Comparing {.val {ucelltype}} with other cells...")
    is.current.celltype <- (celltypes == ucelltype)
    not.current.celltype <- !is.current.celltype
    # split genes into chunks
    idx <- NULL  # to prevent "no visible binding for global variable" when checking
    chunks.res.list <- foreach(idx = iter(splitIndices(nrow(expr), nCores.used))) %dopar%{
      chunk.res <- row_t_welch(
        expr[idx,is.current.celltype], expr[idx,not.current.celltype], alternative = 'greater'
      )
      chunk.res[,c('statistic','pvalue')]
    }
    ct.res <- rbindlist(chunks.res.list)
    twelch.res$twelch.stat_info[,ucelltype,1] <- ct.res[['statistic']]
    twelch.res$twelch.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    twelch.res$twelch.fdr_info[,ucelltype,1] <- p.adjust(ct.res[['pvalue']], 'fdr')
  }
  return(twelch.res)
}



#' ZINB-WaVE + DESeq2
#'
#' @param expr A gene by cell matrix storing the expression values. Should be raw counts.
#' @param celltypes A vector indicating cell types of each cell
#' @param nCores.used The number of cores actually used
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#' is genes, the second dim is cell types, and the last dim is a single subject.
#' @importFrom DESeq2 DESeqDataSet DESeq results
#' @importFrom BiocParallel MulticoreParam
#' @importFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom zinbwave zinbwave
#'
BaselineMethod.DEseq2 <- function(expr, celltypes, nCores.used){
  stopifnot(all(expr %% 1 == 0))
  BPPARAM <- MulticoreParam(nCores.used)

  ucelltypes <- unique(celltypes)
  cleaned.celltypes <- make.names(celltypes)
  cleaned.ucelltypes <- make.names(ucelltypes)
  # create 3-dim empty arrays to store results for a single subject (i.e, the last dim is 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(ucelltypes), 1),
                     dimnames = list(rownames(expr), ucelltypes))
  DESeq2.res <- list('DESeq2.stat_info' = array.tmp,
                     'DESeq2.pval_info' = array.tmp,
                     'DESeq2.fdr_info' = array.tmp,
                     'DESeq2.log2FC_info' = array.tmp)

  # build a SummarizedExperiment
  core <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(expr)),
                               colData=data.frame(celltype = as.factor(cleaned.celltypes)))
  # ZINB-WaVE
  core_zinb <- zinbwave(core, K = 2, epsilon=1000, observationalWeights = TRUE, BPPARAM=BPPARAM)
  # DESeq2
  for (cleaned.uct in cleaned.ucelltypes){
    # build two groups
    is.current.celltype <- (cleaned.celltypes == cleaned.uct)
    colData(core_zinb)['group'] <- as.factor(ifelse(is.current.celltype, cleaned.uct, "others"))

    dds <- DESeqDataSet(core_zinb, design = ~ group)
    dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6, parallel=T, BPPARAM=BPPARAM)
    ct.res <- results(object = dds, contrast = c("group", cleaned.uct, "others"),
                              alpha = 0.05, pAdjustMethod ='fdr', altHypothesis="greater",
                              parallel=T, BPPARAM=BPPARAM)
    # store results
    ucelltype <- ucelltypes[cleaned.ucelltypes == cleaned.uct]
    DESeq2.res$DESeq2.stat_info[,ucelltype,1] <- ct.res[['stat']]
    DESeq2.res$DESeq2.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    DESeq2.res$DESeq2.fdr_info[,ucelltype,1] <- ct.res[['padj']]
    DESeq2.res$DESeq2.log2FC_info[,ucelltype,1] <- ct.res[['log2FoldChange']]
  }
  return(DESeq2.res)
}



#' NS-Forest
#'
#' @param expr A gene by cell matrix storing the expression values. Should be log-normalized counts.
#' @param celltypes A vector indicating cell types of each cell
#' @param python.path The path to the python. Need to install \pkg{nsforest} in advance.
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#' is genes, the second dim is cell types, and the last dim is a single subject.
#' \strong{Note}: the length of the first dim may not equal the total number of genes,
#' since \emph{negative markers} defined in the \href{https://doi.org/10.1101/gr.275569.121}{paper}
#' are filtered. Besides some values may be \code{NA}s due to the inconsistency between
#' \emph{positive markers} of each cell type.
#'
#' @importFrom reticulate use_python import r_to_py
#' @importFrom utils read.csv
#'
BaselineMethod.NSforest <- function(expr, celltypes, python.path){
  # import python packages
  use_python(python.path)
  nsforest <- import("nsforest")
  ad <- import("anndata")

  ucelltypes <- unique(celltypes)
  X <- r_to_py(t(expr))
  obs <- r_to_py(data.frame(celltype = celltypes))
  # create the anndata object
  adata <- ad$AnnData(X=X, obs=obs)
  adata$var_names <- r_to_py(rownames(expr))
  # run NSForest for all cell types
  # the detailed results are saved in ./NSForest_outputs/NSForest_supplementary.csv
  nsforest$NSForest(adata, cluster_header='celltype',
                    n_top_genes=adata$n_vars, n_binary_genes=adata$n_vars)
  supp <- read.csv("./NSForest_outputs/NSForest_supplementary.csv")
  # create 3-dim empty arrays to store results for a single subject (i.e, the last dim is 1)
  unique_binary_genes <- unique(supp$binary_genes)
  unique_celltypes <- unique(supp$clusterName)
  array.tmp <- array(NA, dim=c(length(unique_binary_genes), length(unique_celltypes), 1),
                     dimnames = list(unique_binary_genes, unique_celltypes))
  NSforest.res <- list('NSforest.rf_feature_importance' = array.tmp,
                       'NSforest.binary_score' = array.tmp,
                       'NSforest.cluster_median' = array.tmp)
  for(i in 1:nrow(supp)){
    g <- supp$binary_genes[i]
    ct <- supp$clusterName[i]
    NSforest.res$NSforest.rf_feature_importance[g, ct, 1] <- supp$rf_feature_importance[i]
    NSforest.res$NSforest.binary_score[g, ct, 1] <- supp$binary_score[i]
    NSforest.res$NSforest.cluster_median[g, ct, 1] <- supp$cluster_median[i]
  }
  return(NSforest.res)
}








