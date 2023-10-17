# scCTS: identifying the cell type specific marker genes from population-level single-cell RNA-seq

`scCTS` is an R package for the statistical modeling of the 
gene differential expression (DE) in scRNA-seq data. It identifies cell-type specific genes (markers) that consistently appear in historical population-level scRNA-seq (scRNA-seq) data. `scCTS` is built on top of the R package [`SingleCellExperiment`](https://bioconductor.org/packages/devel/bioc/html/SingleCellExperiment.html) and supports parallel computation.

Except from our proposed method, `scCTS` also provides a common interface for classic DE methods such as the Wilcoxon rank-sum test, t-test and [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).


## Installation
You can install `scCTS` from [GitHub](https://github.com/luxiao10/scCTS) using the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package:

```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github('ToryDeng/scCTS', dependencies=T, build_vignettes = T)
library(scCTS)
```

## Getting started
Here we give a simple example to demonstrate how to run `scCTS`. Once the package is installed, you can load the simulated dataset included in the package:

```R
data(sim.sce)
```

`sim.sce` is a `SingleCellExperiment` object with 2,000 genes and 10,000 cells.
Next, you can run `scCTS` with a single line of code:

```R
res <- scCTS(sim.sce, use.raw = TRUE, subject.rep='subject', celltype.rep='celltype', numCores=2)
```

Some explanations about the parameters:
- **use.raw:** Specifies whether to use the raw counts. If set to `TRUE`, the raw counts will be normalized by `scCTS`.
- **subject.rep:** The name of the column that stores subject labels of cells in the `colData` slot.
- **celltype.rep:** The name of the column that stores cell type labels in the `colData` slot.
- **numCores:** Number of cores for parallel computation.


In the [tested environment](#tested-environment), the code finishes running within a minute. The return value `res` is a list containing lists for each cell type. Each list contains
posterior probabilities of genes and parameter estimations for a particular cell type. For example, you can extract the posterior probabilities of genes to show DE in `celltype1` using the following code:

```R
res$celltype1$pp.d1
```

For more details about how to run `scCTS` and classic DE methods, please refer to `vignette("scCTS")`.


## Reproducibility
If you want to reproduce results shown in the paper, please refer to the directory [reproducibility/](reproducibility/) in this repo.


## Tested environment
### Environment 1
- CPU: AMD Ryzen Threadripper 3990X 64-Core Processor
- Memory: 256GB
- System: Ubuntu 20.04.6 LTS
- R version: 4.3.0

### Environment 2
- CPU: Intel(R) Xeon(R) Gold 6240R CPU @ 2.40GHz
- Memory: 256GB
- System: Ubuntu 22.04.3 LTS
- R version: 4.3.1

## Citation
Coming soon.
