## ----knitr-options, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scCTS)

## -----------------------------------------------------------------------------
data("sim.sce")
sim.sce

## -----------------------------------------------------------------------------
library(scater)

sim.sce <- logNormCounts(sim.sce, log=FALSE)

## -----------------------------------------------------------------------------
scCTS.res <- scCTS(sim.sce, subject.rep='subject', celltype.rep='celltype', numCores=2)

## -----------------------------------------------------------------------------
head(scCTS.res$celltype1$pp.d1)

## -----------------------------------------------------------------------------
bs.res <- runBaselineMethod(sim.sce, per.subject=TRUE,
                            celltype.rep='celltype', subject.rep='subject', 
                            method='wilcox', numCores=2)

## -----------------------------------------------------------------------------
head(bs.res$wilcox.fdr_info[,,'subject1'])

## -----------------------------------------------------------------------------
sessionInfo()

