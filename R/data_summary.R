
#' data.info.collect
#'
#' calculate effect size and corresponding variance.
#'
#' *Note!*
#' This function is not clean. Some of the calculated statistics are not used
#' in current manuscript.
#' Such statistics will be marked with `IGNORE!`
#'
#'
#' @param Y gene expression matrix (normalized to 10k)
#' @param celltype cell type info of cells
#'
#' @return All necessary sample level summary statistics
#' @import DelayedArray
#' @importFrom stats var
#' @importFrom matrixTests row_wilcoxon_twosample
#'

data.info.collect <- function(Y, celltype){
  gene.num <- nrow(Y)  # gene number
  unique.celltypes <- sort(unique(celltype)) # names of cell types
  celltype.num <- length(unique.celltypes) # number of cell types, NOT number of cells
  gene.names <- rownames(Y) # gene names

  res <- list() # store result
  matrix.tmp <-   matrix(NA,nrow = gene.num, ncol = celltype.num) # create a matrix
  colnames(matrix.tmp) <- unique.celltypes # column is cell type
  rownames(matrix.tmp) <- gene.names # row is genes

  # create null variables to store different statistics
  res[['expr']] = res[['expr_var']] = res[['expr_se']] =
    res[['log_expr']] = res[['log_expr_se']] =
    res[['expr_remain']] = res[['expr_remain_var']] = res[['expr_remain_se']] =
    res[['log_expr_remain']] = res[['log_expr_remain_se']] =
    res[['effect']] = res[['effect_se']] =
    res[['effect_2']] = res[['effect_2_se']] =
    res[['wilcox.stat']] = res[['wilcox.pval']] =
    res[['nonzero_rate']] = res[['nonzero_exprmean']] =
    res[['nonzero_log_expr']] <- matrix.tmp

  # do the loop for each cell type
  for(ct.ix in unique.celltypes){
    cell.ix <- (celltype == ct.ix) # index for target cell type
    expr.ct.tmp <- Y[,cell.ix]  # gene expression of cells from target cell type
    expr.remain.tmp <- Y[,!cell.ix]  # gene expression of cells from non-targer cell types

    res[['expr']][,ct.ix] <- DelayedArray::rowMeans(expr.ct.tmp)  # mean expression of target cell type
    res[['expr_var']][,ct.ix] <- apply(expr.ct.tmp,1,var) # variance of expression for target cell type
    res[['expr_se']][,ct.ix] <- sqrt( res[['expr_var']][,ct.ix] / sum(cell.ix)) # sample error for target cell type

    res[['log_expr']][,ct.ix] <- log2(res[['expr']][,ct.ix] + 1) # log2 transformed gene expression + psudo count 1
    res[['log_expr_se']][,ct.ix] <-
      sqrt(res[['expr_var']][,ct.ix]/(sum(cell.ix)*((res[['expr']][,ct.ix]+1)*log(2))^2))

    res[['expr_remain']][,ct.ix] <- DelayedArray::rowMeans(expr.remain.tmp) # IGNORE! mean expression of non-target cells
    res[['expr_remain_var']][,ct.ix] <- apply(expr.remain.tmp,1,var) # IGNORE! variance of expression for non-target cell types
    res[['expr_remain_se']][,ct.ix] <- sqrt( res[['expr_remain_var']][,ct.ix] / sum(!cell.ix)) # GNORE! sample error for non-target cell type

    res[['log_expr_remain']][,ct.ix] <- log2(res[['expr_remain']][,ct.ix] + 1) # IGNORE! log2 transformed expression + psduo count 1
    res[['log_expr_remain_se']][,ct.ix] <-
      sqrt(res[['expr_remain_var']][,ct.ix]/(sum(!cell.ix)*((res[['expr_remain']][,ct.ix]+1)*log(2))^2)) # GNORE!


    res[['effect_2']][,ct.ix] <-  res[['log_expr']][,ct.ix] - res[['log_expr_remain']][,ct.ix] # GNORE!
    res[['effect_2_se']][,ct.ix] <- sqrt( res[['log_expr_remain_se']][,ct.ix]^2 +  # IGNORE!
                                            res[['log_expr_se']][,ct.ix]^2 )


    # wilcoxon rank sum test one vs. others (a common used/traditional way)
    test.res <- row_wilcoxon_twosample(Y[,cell.ix], Y[,!cell.ix], alternative = 'greater')

    res[['wilcox.stat']][,ct.ix] = test.res[,'statistic']
    res[['wilcox.pval']][,ct.ix] = test.res[,'pvalue']


    nonzero.ix.tmp <- expr.ct.tmp != 0

    res[['nonzero_rate']][,ct.ix] <- apply(nonzero.ix.tmp,1,mean)  # IGNORE!drop out rate, not used here, but for previous data explore
    res[['nonzero_exprmean']][,ct.ix] <- res[['expr']][,ct.ix] / res[['nonzero_rate']][,ct.ix] # GNORE! mean of non-zero expression
    res[['nonzero_log_expr']][,ct.ix] <- log2(res[['nonzero_exprmean']][,ct.ix]+1) # IGNORE!

  }

  ### Loop again, for each cell type calculate effect size (a.k.a log2 fold change defined in Eq 2 in manuscript)
  ### and its corresponding variance defined in Eq 6 in Supplementary file
  ### (looks like Eq6 is not correct, it should be summation + not minus - between the two terms) ####
  for(ct.ix in unique.celltypes){
    res[['effect']][,ct.ix] <- res[['log_expr']][,ct.ix] -
      log2(rowMeans(res[['expr']][,! unique.celltypes %in% ct.ix],na.rm = T) +1)

    se.part1 <- res[['log_expr_se']][,ct.ix]^2
    se.part2 <- 1/(log(2)*(rowSums(res[['expr']][,! unique.celltypes %in% ct.ix],na.rm = T) +
                             (celltype.num-1)))^2 * rowSums(res[['expr_se']][,! unique.celltypes %in% ct.ix]^2)

    res[['effect_se']][,ct.ix] <- sqrt(se.part1 + se.part2) ### I took a square here, in following analysis we need variance
  }


  res[['nonzero_exprmean']][is.na(res[['nonzero_exprmean']])] <- 0 # just avoid na
  res[['nonzero_log_expr']][is.na(res[['nonzero_log_expr']])] <- 0 # just avoid na

  return(res)
}


#' Re-organize the data generated by \code{data.info.collect}
#' Convert the sublists containing necessary stats into arrays
#'
#' @param marker.info The data generated by \code{data.info.collect}.
#'
#' @return A list containing converted assays
#' @importFrom cli cli_text
#' @importFrom stats p.adjust
#'

data.info.reorganize <- function(marker.info){
  cli_text("Re-organizing data info...")
  nGenes <- nrow(marker.info[[1]]$expr)
  nCT <- ncol(marker.info[[1]]$expr)
  nInds <- length(marker.info)

  array.tmp <- array(NA, dim=c(nGenes, nCT, nInds))
  dimnames(array.tmp) <- list(rownames(marker.info[[1]]$expr),
                              colnames(marker.info[[1]]$expr),
                              names(marker.info))

  expr_info = expr_var_info =
    log_expr_info = log_expr_se_info =
    expr_remain_info = expr_remain_var_info =
    log_expr_remain_info = log_expr_remain_se_info =
    effect_info = effect_se_info =
    effect_2_info = effect_2_se_info =
    wilcox.stat_info = wilcox.pval_info = wilcox.fdr_info =
    nonzero_rate_info = nonzero_exprmean_info =
    nonzero_log_expr_info <- array.tmp

  subjects.tmp <- names(marker.info)
  for(sub.ix in subjects.tmp){

    expr_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$expr)

    expr_var_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$expr_var)

    log_expr_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$log_expr)

    log_expr_se_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$log_expr_se)

    expr_remain_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$expr_remain)

    expr_remain_var_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$expr_remain_var)

    log_expr_remain_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$log_expr_remain)

    log_expr_remain_se_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$log_expr_remain_se)


    effect_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$effect)

    effect_se_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$effect_se)

    effect_2_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$effect_2)

    effect_2_se_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$effect_2_se)

    wilcox.stat_info[,,sub.ix] <-
      as.matrix(marker.info[[sub.ix]]$wilcox.stat)

    wilcox.pval_info[,,sub.ix] <- as.matrix(marker.info[[sub.ix]]$wilcox.pval)

    nonzero_rate_info[,,sub.ix] <- as.matrix(marker.info[[sub.ix]]$nonzero_rate)
    nonzero_exprmean_info[,,sub.ix] <- as.matrix(marker.info[[sub.ix]]$nonzero_exprmean)
    nonzero_log_expr_info[,,sub.ix] <- as.matrix(marker.info[[sub.ix]]$nonzero_log_expr)

    for(cell.ix in 1:nCT){
      wilcox.fdr_info[, cell.ix, sub.ix] <-
        p.adjust(c(wilcox.pval_info[, cell.ix, sub.ix]), 'fdr')
    }
  }

  res <- list('expr_info'=expr_info,
              'expr_var_info' = expr_var_info,
              'log_expr_info' = log_expr_info,
              'log_expr_se_info' = log_expr_se_info,
              'expr_remain_info'=expr_remain_info,
              'expr_remain_var_info' = expr_remain_var_info,
              'log_expr_remain_info' = log_expr_remain_info,
              'log_expr_remain_se_info' = log_expr_remain_se_info,
              'effect_info' = effect_info,
              'effect_se_info' = effect_se_info,
              'effect_2_info' = effect_2_info,
              'effect_2_se_info' = effect_2_se_info,
              'wilcox.stat_info' = wilcox.stat_info,
              'wilcox.pval_info' = wilcox.pval_info,
              'wilcox.fdr_info' = wilcox.fdr_info,
              'nonzero_rate_info' = nonzero_rate_info,
              'nonzero_exprmean_info' = nonzero_exprmean_info,
              'nonzero_log_expr_info' = nonzero_log_expr_info)
  cli_text("Finish data info re-organization.")
  return(res)
}
