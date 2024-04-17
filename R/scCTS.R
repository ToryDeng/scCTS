#' The main function of scCTS. It performs DE analysis on normalized (not
#' log-normalized) count matrix.
#'
#' @param sce A \code{SingleCellExperiment} object. Should contain normalized
#'   count matrix, subject and cell type info.
#' @param use.raw Whether to use the raw counts. If \code{TRUE}, the raw counts
#'   are then normalized by \pkg{scCTS}.
#' @param use.norm.rep Which count matrix should be used. Default is the matrix
#'   accessed by \code{normcounts()}.
#' @param subject.rep The name of column that stores subject labels of cells in
#'   \code{colData} slot.
#' @param celltype.rep The name of column that stores cell type labels in
#'   \code{colData} slot.
#' @param log.input Whether the input expression matrix is log-transformed.
#' @param log.base The base of log-transformation.
#' @param effect_thres Threshold for filtering genes with negative mean
#'   (\eqn{m_{gk}} in Eq 3 in manuscript)
#' @param maxiter Maximum iteration number
#' @param tol EM stop control
#' @param numCores Number of cores for parallel computation.
#' @param min.cutoff Remove samples with extreme small log fold change for
#'   robust estimation. Default is quantile 0.05.
#' @param max.cutoff Remove samples with extreme large log fold change for
#'   robust estimation. Default is quantile 0.95.
#' @param numCores Number of cores to use. Default is the number of all possible
#'   cores minus 1.
#' @param verbose Whether to print details when the function is running.
#'
#' @return A list containing lists for every cell type. Each list contains
#'   posterior probability and estimates of parameters for a cell type.
#' @import SingleCellExperiment cli
#' @importFrom foreach %dopar% foreach
#' @export
#'
#' @examples
#' # load the simulated data
#' data(sim.sce)
#'
#' # run scCTS
#' res <- scCTS(sim.sce, use.raw = TRUE,
#'              subject.rep='subject', celltype.rep='celltype',
#'              numCores=2)
#'
scCTS <- function(
    sce,
    use.raw=FALSE,
    use.norm.rep=NULL,
    subject.rep='subject',
    celltype.rep='celltype',
    log.input=FALSE,
    log.base=2,
    effect_thres= 0.01,
    maxiter=1000,
    tol=1e-3,
    min.cutoff = 0.05,
    max.cutoff = 0.95,
    numCores=NULL,
    verbose=FALSE
){

  # check arguments
  stopifnot(all(c(subject.rep, celltype.rep) %in% names(colData(sce))))
  stopifnot(all(is.numeric(c(effect_thres, maxiter, tol, min.cutoff, max.cutoff))))

  start.time <- Sys.time()
  cli_h1("scCTS analysis")
  cli_alert_info("Start at {.var {start.time}}")

  # get normalized count matrix
  Ynorm <- get.expression.matrix(sce, use.raw, use.norm.rep, normalize=TRUE)
  if (log.input){
    Ynorm <- inverse.log(Ynorm, log.base)
  }
  # set parallel computation
  set.parallel.computation(numCores)
  # get subjects and cell types
  subjects <- colData(sce)[[subject.rep]]
  celltypes <- colData(sce)[[celltype.rep]]
  unique.subjects <- sort(unique(subjects))

  # prepare info for EM estimation
  subject.idx <- NULL # to prevent "no visible binding for global variable" when checking
  data.info <- foreach(subject.idx = 1:length(unique.subjects)) %dopar% {
  # data.info <- list()
  # for (subject.idx in 1:length(unique.subjects)){
    current.sub <- unique.subjects[subject.idx]
    cli_text("Collecting info on subject {.val {current.sub}}...")
    target_samples <- which(subjects == current.sub)
    Ynorm.tmp <- Ynorm[,target_samples]  # expressions of a single subject
    celltypes.tmp <- celltypes[target_samples]  # cell types of a single subject
    data.info.tmp <- data.info.collect(Ynorm.tmp, celltypes.tmp)
    cli_text("Finish info collectionon on subject {.val {current.sub}}.")
    data.info.tmp
    # data.info[[current.sub]] <- data.info.tmp
  }
  names(data.info) <- unique.subjects
  data.info <- data.info.reorganize(data.info)

  # perform EM estimation
  res <- EM_estimate(data.info, effect_thres=effect_thres, maxiter=maxiter,
                     tol=tol, min.cutoff=min.cutoff, max.cutoff=max.cutoff,
                     verbose=verbose)
  end.time <- Sys.time()
  diff.time <- as.numeric(end.time - start.time, units = "secs")
  cli_text("scCTS completes at {.var {end.time}}  Totoal: {.var {diff.time}} seconds ({prettyunits::pretty_sec(diff.time)})")
  return(res)
}



