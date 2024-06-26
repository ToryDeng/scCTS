% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baselines.R
\name{BaselineMethod.DEseq2}
\alias{BaselineMethod.DEseq2}
\title{ZINB-WaVE + DESeq2}
\usage{
BaselineMethod.DEseq2(
  expr,
  celltypes,
  subjects = NULL,
  nCores.used = NULL,
  cache.path = NULL
)
}
\arguments{
\item{expr}{A gene by cell matrix storing the expression values}

\item{celltypes}{A vector indicating cell types of each cell}

\item{subjects}{Subject-level labels}

\item{nCores.used}{The number of cores actually used}

\item{cache.path}{Path to the cached RDS file}
}
\value{
A list of 3-dim arrays. Each array corresponds to a stat. The first
  dim is genes, the second dim is cell types, and the last dim is a single
  subject.
}
\description{
ZINB-WaVE + DESeq2
}
