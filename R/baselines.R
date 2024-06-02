
#' A wrapper function of classic differential expression tests and
#' feature/marker selection methods.
#'
#' Currently, differential expression (DE) tests implemented are:
#' \enumerate{
#'   \item Wilcoxon test
#'   \item t-test
#'   \item ZINB-WaVE + DESeq2
#' }
#' These tests are aimed to find over-expressed DE genes.\cr Feature/marker
#' selection methods implemented are:
#' \enumerate{
#'   \item NS-Forest
#'   \item FEAST
#'   \item scGeneFit
#' }
#' For each cell type, these methods select a predefined number of markers.
#'
#' @details
#' For ZINB-WaVE + DESeq2, if \code{per.subject=TRUE}, the subject labels and cell
#' type labels will both be included in the formula. Otherwise only the subject
#' labels will be included in the formula.
#'
#' NS-Forest and scGeneFit are python packages. You need to install them through
#' \code{pip install nsforest scGeneFit} and provide the python path before you
#' call this function.
#'
#' Values in \code{celltype.ngenes} must be integers. You can specify it by
#' \code{celltype.ngenes = list(celltype1 = 3000L)}.
#'
#' @param sce A \code{SingleCellExperiment} object. Should contain normalized
#'   count matrix, subject and cell type info.
#' @param use.norm.rep Which count matrix should be used. Default is the matrix
#'   accessed by \code{normcounts()}.
#' @param subject.rep Column name for the subject info.
#' @param celltype.rep Column name for the cell type info.
#' @param per.subject If \code{TRUE}, calculate for each sample. Otherwise
#'   ignore the subject labels and calculate for all subjects.
#' @param method  A string indicating which DE method to use.
#' @param celltype.ngenes A named list. The names are cell types, and the values
#'   are number of genes to be selected for that cell type.
#' @param numCores Number of cores to use. Default is \code{NULL}, which means
#'   using all but one of the CPU cores.
#' @param python.path The path to the python.
#'
#' @return \itemize{
#'   \item For Wilcoxon and t tests, if \code{per.subject=TRUE}, a list of
#'   3-dim arrays. Each array is corresponding to a type of information (e.g.,
#'   P-values). The first dim is corresponding to genes, the second dim is
#'   corresponding to cell types, and the last dim is corresponding to subjects.
#'   If \code{per.subject=FALSE}, the last dim has a size of 1 and a name "all".
#'   \item For ZINB-WaVE + DESeq2, a list of 3-dim arrays. The first two dims are
#'   the same as Wilcoxon and t tests, but the last dim always has length 1.
#'   \item For feature/marker selection methods, if \code{per.subject=TRUE}, a
#'   list of sublists. Each sublist contains cell-type-specific selection
#'   results of each subject. If \code{per.subject=FALSE}, a list of
#'   cell-type-specific selection results regarding all subjects as a whole. }
#' @import SingleCellExperiment
#' @importFrom abind abind
#' @importFrom stringr str_glue
#' @importFrom S4Vectors metadata
#' @export
#'
#' @examples
#' # load the simulated data
#' data(sim.sce)
#'
#' # normalize
#' sim.sce <- scater::logNormCounts(sim.sce, log=FALSE)
#'
#' # run the Wilcoxon test per subject
#' res <- runBaselineMethod(sim.sce,
#'                          subject.rep='subject', celltype.rep='celltype',
#'                          per.subject=TRUE, method='wilcox', numCores=2)
#'
#' @references Risso, Davide, et al. "A general and flexible method for signal
#'   extraction from single-cell RNA-seq data." Nature communications 9.1
#'   (2018): 284.
#' @references Love, Michael, Simon Anders, and Wolfgang Huber. "Differential
#'   analysis of count data–the DESeq2 package." Genome Biol 15.550 (2014):
#'   10-1186.
#' @references Aevermann, Brian, et al. "A machine learning method for the
#'   discovery of minimum marker gene combinations for cell type identification
#'   from single-cell RNA sequencing." Genome research 31.10 (2021): 1767-1780.
#' @references Su, Kenong, Tianwei Yu, and Hao Wu. "Accurate feature selection
#'   improves single-cell RNA-seq cell clustering." Briefings in bioinformatics
#'   22.5 (2021): bbab034.
#' @references Dumitrascu, Bianca, et al. "Optimal marker gene selection for
#'   cell type discrimination in single cell analyses." Nature communications
#'   12.1 (2021): 1186.
#'
runBaselineMethod <- function(
    sce,
    use.norm.rep=NULL,
    subject.rep='subject',
    celltype.rep='celltype',
    per.subject=TRUE,
    method=c('wilcox', 'twelch', 'DEseq2', 'NSforest', 'FEAST', 'scGeneFit'),
    celltype.ngenes=NULL,
    numCores=NULL,
    python.path=NULL
){
  # check arguments
  stopifnot(all(c(subject.rep, celltype.rep) %in% names(colData(sce))),
            is(per.subject, "logical"))
  # set parallel computation
  numCores.used <- set.parallel.computation(numCores)

  # check python path
  if (method %in% baselines.python() & is.null(python.path)){
    cli_abort("{.var {method}} requires a python path.")
  }
  if (!(method %in% baselines.python()) & !is.null(python.path)){
    cli_alert_warning("Invalid argument python.path={.file {python.path}} for {.var {method}}.")
  }
  # get (normalized) count matrix
  if (method %in% baselines.counts()){
    Y <- get.expression.matrix(sce, use.raw=TRUE)
  }else{
    Y <- get.expression.matrix(sce, use.raw=FALSE, use.norm.rep=use.norm.rep)
  }
  # get cell types
  celltypes <- colData(sce)[[celltype.rep]] # convert to char list
  # check the list of numbers of selected markers for each cell type
  if (method %in% baselines.fixnumber()){
    # marker numbers should be given
    if (is.null(celltype.ngenes)){
      cli_abort("{.var {method}} requires predefined numbers of markers for each cell type.")
    }
    # cell type names should be equal
    if (!setequal(names(celltype.ngenes), celltypes)){
      cli_abort("Differences between: {.var {celltype.ngenes}} and {.var {celltypes}}.")
    }
    # marker numbers should be integers
    if (!is.integer(unlist(celltype.ngenes, use.names=FALSE))){
      cli_abort("Predefined marker numbers should be all integers, got {.var {celltype.ngenes}}.")
    }
  }
  if (!(method %in% baselines.fixnumber()) & !is.null(celltype.ngenes)){
    cli_alert_warning("Invalid argument celltype.ngenes={.var {celltype.ngenes}} for {.var {method}}.")
  }
  # prepare cache dir, currently only for ZINB-WaVe + DESeq2
  name.exists <- exists("dataset_name", where=metadata(sce))
  dataset.name <- ifelse(name.exists, metadata(sce)$dataset_name, deparse1(substitute(sce)))  # get variable's name
  cache.path <- file.path("cache", dataset.name)

  start.time <- Sys.time()
  if (per.subject){
    cli_h1("Subject-level {.emph {method}} method")
    cli_alert_info("Start at {.var {start.time}}")

    subjects <- colData(sce)[[subject.rep]] # convert to char list
    unique.subjects <- sort(unique(subjects))

    if (method != "DEseq2"){
      sub.res.list <- list()
      cli_progress_bar("Analyzing each subject", total = length(unique.subjects), type = "tasks")
      for (sub in unique.subjects){
        # cli_h2("Current subject: {sub}")
        sub.Y <- Y[,subjects == sub]
        sub.cts <- celltypes[subjects == sub]
        sub.result = switch(
          method,
          "wilcox" = BaselineMethod.wilcox(sub.Y, sub.cts, numCores.used),
          "twelch" = BaselineMethod.twelch(sub.Y, sub.cts, numCores.used),
          # "DEseq2" = BaselineMethod.DEseq2(sub.Y, sub.cts, numCores.used, cache.path),
          "NSforest" = BaselineMethod.NSforest(sub.Y, sub.cts, celltype.ngenes, python.path, numCores.used),
          "FEAST" = BaselineMethod.FEAST(sub.Y, sub.cts, celltype.ngenes, numCores.used),
          "scGeneFit" = BaselineMethod.scGeneFit(sub.Y, sub.cts, celltype.ngenes, python.path),
          stop(str_glue("No method matched for {method}"))
        )
        if (!(method %in% baselines.fixnumber())){
          # for DE tests, set the name of the last dim of each array as subject
          # name, and then store
          sub.res.list[[sub]] <- lapply(sub.result, function(arr){dimnames(arr)[3] <- sub;arr})
        }else{
          sub.res.list[[sub]] <- sub.result
        }
        cli_progress_update()
      }

      # post-analysis processing: combine results for multiple subjects
      if (!(method %in% baselines.fixnumber())){
        DE.res <- list()
        for (name in names(sub.res.list[[1]])){
          DE.res[[name]] <- abind(lapply(sub.res.list, function(lst) lst[[name]]))
        }
        final.res <- DE.res
    }else{
      final.res <- sub.res.list
    }
  }else{  # method == "DEseq2"
    final.res <- BaselineMethod.DEseq2(Y, celltypes, subjects, numCores.used, cache.path)
    final.res <- lapply(final.res, function(arr){dimnames(arr)[3] <- "per.subject";arr})
  }
  }else{
    cli_h1("Population-level {.emph {method}} method")
    cli_alert_warning("{.emph Note: this mode needs batch-effects correction in advance.}")
    cli_alert_info("Start at {.var {start.time}}")

    subjects <- colData(sce)[[subject.rep]]
    all.result = switch(
      method,
      "wilcox" = BaselineMethod.wilcox(Y, celltypes, numCores.used),
      "twelch" = BaselineMethod.twelch(Y, celltypes, numCores.used),
      "DEseq2" = BaselineMethod.DEseq2(Y, celltypes, NULL, numCores.used, cache.path),
      "NSforest" = BaselineMethod.NSforest(Y, celltypes, celltype.ngenes, python.path, numCores.used),
      "FEAST" = BaselineMethod.FEAST(Y, celltypes, celltype.ngenes, numCores.used),
      "scGeneFit" = BaselineMethod.scGeneFit(Y, celltypes, celltype.ngenes, python.path),
      stop(str_glue("No method matched for {method}"))
    )
    if (!(method %in% baselines.fixnumber())){
      all.result <- lapply(all.result, function(arr){dimnames(arr)[3] <- "all";arr})
    }
    final.res <- all.result
  }
  end.time <- Sys.time()
  diff.time <- as.numeric(end.time - start.time, units = "secs")
  cli_alert_info("Ends at {.var {end.time}}  Totoal: {.var {diff.time}} seconds ({prettyunits::pretty_sec(diff.time)})")

  return(final.res)
}


# returns methods implemented in python
baselines.python <- function(){
  return(c("NSforest", "scGeneFit"))
}

# returns methods that require a predefined number of markers
baselines.fixnumber <- function(){
  return(c("NSforest", "FEAST", "scGeneFit"))
}

# returns methods that require raw counts as inputs
baselines.counts <- function(){
  return(c("DEseq2"))
}



#' Naive Wilcoxon test
#'
#' @param expr A gene by cell matrix storing the expression values
#' @param celltypes A vector indicating cell types of each cell
#' @param nCores.used The number of cores actually used
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#'   is genes, the second dim is cell types, and the last dim is subjects.
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
#' @inheritParams BaselineMethod.wilcox
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#'   is genes, the second dim is cell types, and the last dim is subjects.
#' @importFrom matrixTests row_t_welch
#' @importFrom foreach %dopar% foreach
#' @importFrom data.table rbindlist
#' @importFrom iterators iter
#' @importFrom parallel splitIndices
#' @importFrom stats p.adjust
#'
#'
BaselineMethod.twelch <- function(expr, celltypes, nCores.used){

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
#' @inheritParams BaselineMethod.wilcox
#' @param subjects Subject-level labels
#' @param cache.path Path to the cached RDS file
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first
#'   dim is genes, the second dim is cell types, and the last dim is a single
#'   subject.
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stringr str_glue
#' @importFrom S4Vectors SimpleList
#' @import SingleCellExperiment cli
#'
BaselineMethod.DEseq2 <- function(expr, celltypes, subjects=NULL, nCores.used=NULL, cache.path=NULL){
  for (pkg in c("zinbwave", "DESeq2")){
    if(!requireNamespace(pkg)){
      cli_abort("This function requires the {.pkg {pkg}} package.")
    }
  }
  BPPARAM <- MulticoreParam(nCores.used)

  ucelltypes <- unique(celltypes)
  cleaned.celltypes <- make.names(celltypes)
  cleaned.ucelltypes <- make.names(ucelltypes)
  # create 3-dim empty arrays to store results for a single subject (i.e, the last dim has length 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(ucelltypes), 1),
                     dimnames = list(rownames(expr), ucelltypes))
  DESeq2.res <- list('DESeq2.stat_info' = array.tmp,
                     'DESeq2.pval_info' = array.tmp,
                     'DESeq2.fdr_info' = array.tmp,
                     'DESeq2.log2FC_info' = array.tmp)

  zinbwave.path <- file.path(cache.path, "zinbwave.rds")
  if (file.exists(zinbwave.path)){  # directly load the cache
    zinb <- readRDS(zinbwave.path)
    cli_alert_success("Cached ZINB-WaVe results loaded from {.file {zinbwave.path}}.")
  }else{  # run zinbwave and save cache
    # build a SingleCellExperiment object, according to the values of `subjects`
    if (is.null(subjects)){
      cli_alert_info("Subject-level labels not included.")
      colD <- data.frame(celltype = factor(cleaned.celltypes))
      core <- SingleCellExperiment(assays=list(counts=expr), colData=colD)
    }else{
      cli_alert_info("Subject-level labels included.")
      colD <- data.frame(celltype = factor(cleaned.celltypes), subject = factor(make.names(subjects)))
      core <- SingleCellExperiment(assays=list(counts=expr), colData=colD)
    }

    # ZINB-WaVE, specify `K = 0` to only compute observational weights
    zinb <- zinbwave::zinbwave(core, K=0, observationalWeights=TRUE, BPPARAM=BPPARAM, epsilon=1e12)
    mode(assay(zinb)) <- "integer"  # to prevent the message "converting counts to integer mode"
    if (!dir.exists(cache.path)){dir.create(cache.path, recursive=TRUE)}
    saveRDS(zinb, zinbwave.path)
    cli_alert_success("ZINB-WaVe results cached at {.file {zinbwave.path}}.")
  }

  # DESeq2
  cli_progress_bar("Running DESeq2", total = length(cleaned.ucelltypes))
  for (cleaned.uct in cleaned.ucelltypes){
    # prepare cache path
    ct.path <- file.path(cache.path, str_glue("{cleaned.uct}.rds"))
    # check whther the cell type cache exists
    if (file.exists(ct.path)){
      ct.res <- readRDS(ct.path)
    }else{
      # if the cache doesn't exist, build two cell type groups: current cell type vs all others
      is.current.celltype <- (cleaned.celltypes == cleaned.uct)
      colData(zinb)['group'] <- as.factor(ifelse(is.current.celltype, cleaned.uct, "others"))
      tryCatch({
        if (is.null(subjects)){
          dds <- DESeq2::DESeqDataSet(zinb, design = ~ group)
        }else{  # add subjects and their interactions with cell types (group)
          dds <- DESeq2::DESeqDataSet(zinb, design = ~ group + subject + (group * subject))
        }

        dds <- DESeq2::DESeq(
          dds, sfType="poscounts", useT=TRUE, minmu=1e-6, minRep=Inf, fitType='local',
          parallel=T, BPPARAM=BPPARAM, quiet=TRUE
        )
        ct.res <- DESeq2::results(object = dds, contrast = c("group", cleaned.uct, "others"),
                                  alpha = 0.05, pAdjustMethod ='fdr', altHypothesis="greater",
                                  parallel=T, BPPARAM=BPPARAM)
      }, error=function(e){print(e);cli_alert_warning("Error occured when processing {.var {cleaned.uct}}, continue...")})
      saveRDS(ct.res, ct.path)
    }

    # store results
    ucelltype <- ucelltypes[cleaned.ucelltypes == cleaned.uct]
    DESeq2.res$DESeq2.stat_info[,ucelltype,1] <- ct.res[['stat']]
    DESeq2.res$DESeq2.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    DESeq2.res$DESeq2.fdr_info[,ucelltype,1] <- ct.res[['padj']]
    DESeq2.res$DESeq2.log2FC_info[,ucelltype,1] <- ct.res[['log2FoldChange']]
    # update progress bar
    cli_progress_update(status = str_glue("last finished: {cleaned.uct}"))
  }
  cli_progress_done()
  return(DESeq2.res)
}





#' NS-Forest
#'
#' Running this function will create a folder ./NSForest_outputs/ in the current
#' working directory. For each cell type, the genes should be sorted first by
#' binary scores and then by feature importances from the random forest model.
#' The number of actually selected features may be less than the given number of
#' selected features, since \emph{negative markers} defined in the
#' \href{https://doi.org/10.1101/gr.275569.121}{paper} are filtered.
#'
#' @inheritParams BaselineMethod.wilcox
#' @param celltype.ngenes A named list. The names are cell types, and the values
#'   are number of features selected for that cell type.
#' @param python.path The path to the python.
#'
#' @return A named list. Names are unique cell types. Values are selected
#'   features for that cell type.
#'
#' @import dplyr
#' @importFrom readr read_csv
#' @importFrom purrr map2
#'
BaselineMethod.NSforest <- function(expr, celltypes, celltype.ngenes, python.path, nCores.used=NULL){
  # import python packages
  pkg <- "reticulate"
  if(!requireNamespace(pkg)){
    cli_abort("This function requires the {.pkg {pkg}} package.")
  }
  reticulate::use_python(python.path)
  nsforest <- reticulate::import("nsforest")
  ad <- reticulate::import("anndata")

  ucelltypes <- unique(celltypes)
  X <- reticulate::r_to_py(t(as.data.frame(expr)))
  obs <- reticulate::r_to_py(data.frame(celltype=celltypes, row.names=colnames(expr)))
  # create the anndata object
  adata <- ad$AnnData(X=X, obs=obs)
  adata$var_names <- reticulate::r_to_py(rownames(expr))
  # run NSForest for all cell types
  # the detailed results are saved in ./NSForest_outputs/NSForest_supplementary.csv
  reticulate::py_capture_output(nsforest$NSForest(
    adata, cluster_header='celltype', n_top_genes=adata$n_vars, n_binary_genes=adata$n_vars, n_jobs=as.integer(nCores.used)
  ))
  supp <- read_csv("./NSForest_outputs/NSForest_supplementary.csv", show_col_types=FALSE)

  process_cluster <- function(cluster_name, size, data) {
    filtered_data <- data %>%
      filter(get("clusterName") == cluster_name) %>%
      arrange(desc(get("binary_score")), desc(get("rf_feature_importance"))) %>%
      slice_head(n = size)  # use get() to prevent the no visible binding issue

    return(filtered_data$binary_genes)
  }

  NSforest.res <- map2(.x = names(celltype.ngenes),
                       .y = celltype.ngenes,
                       .f = function(name, size) process_cluster(name, size, supp))
  names(NSforest.res) <- names(celltype.ngenes)
  return(NSforest.res)
}


#' FEAST
#'
#' @inheritParams BaselineMethod.wilcox
#' @param celltype.ngenes A named list. Names are cell types, and the values
#'   are number of features selected for that cell type.
#'
#' @return A list of each cell type's highly variable genes. The number of HVGs
#'   are equal to the given number.
#' @importFrom purrr quietly
#'
BaselineMethod.FEAST <- function(expr, celltypes, celltype.ngenes, nCores.used){
  pkg <- "FEAST"
  if(!requireNamespace(pkg)){
    cli_abort("This function requires the {.pkg {pkg}} package.")
  }
  unique.celltypes <- sort(unique(celltypes))
  FEAST.res <- list()
  for(ucelltype in unique.celltypes){
    ngenes <- celltype.ngenes[[ucelltype]]
    tryCatch({# quietly returns a list of all messages and original return values
      idxs <- quietly(FEAST::FEAST_fast)(expr[,celltypes == ucelltype], nProc=nCores.used)$result
      FEAST.res[[ucelltype]] <- rownames(expr)[idxs[1:ngenes]]
    }, error=function(e){print(e);cli_alert_warning("Error occured when processing {.var {ucelltype}}, continue...")})
  }
  return(FEAST.res)
}


#' scGeneFit
#'
#' @inheritParams BaselineMethod.NSforest
#' @return A named list. Names are unique cell types. Values are selected
#'   features for that cell type.
#'
BaselineMethod.scGeneFit <- function(expr, celltypes, celltype.ngenes, python.path){
  # import python packages
  pkg <- "reticulate"
  if(!requireNamespace(pkg)){
    cli_abort("This function requires the {.pkg {pkg}} package.")
  }
  reticulate::use_python(python.path)
  scGeneFit <- reticulate::import("scGeneFit")

  unique.celltypes <- sort(unique(celltypes))
  scGeneFit.res <- list()
  for(ucelltype in unique.celltypes){
    ngenes <- celltype.ngenes[[ucelltype]]
    X <- reticulate::r_to_py(t(expr))
    cell.labels <- reticulate::r_to_py(ifelse(celltypes == ucelltype, ucelltype, "others"))
    # + 1 since in python indices start from 0
    reticulate::py_capture_output(
      idxs <- scGeneFit$functions$get_markers(X, cell.labels, num_markers=ngenes) + 1
    )
    scGeneFit.res[[ucelltype]] <- rownames(expr)[idxs]
  }
  return(scGeneFit.res)
}

