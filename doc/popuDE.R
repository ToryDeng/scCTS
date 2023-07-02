## ----knitr-options, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(popuDE)

## -----------------------------------------------------------------------------
data("sim.sce")
sim.sce

## -----------------------------------------------------------------------------
popuDE.res <- popuDE(sim.sce,use.raw = TRUE, subject.rep='subject', celltype.rep='celltype', numCores=2, effect_thres=0, tol=1e-5)

## -----------------------------------------------------------------------------
head(popuDE.res$celltype1$pp.d1)

## -----------------------------------------------------------------------------
bs.res <- runBaselineMethod(sim.sce,use.raw=TRUE, per.subject=TRUE,  
                            celltype.rep='celltype', subject.rep='subject', 
                            method='wilcox', numCores=2)

## -----------------------------------------------------------------------------
head(bs.res$wilcox.fdr_info[,,'subject1'])

## -----------------------------------------------------------------------------
sessionInfo()

