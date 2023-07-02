# popuDE: Population-level Differetial Expression Analysis for scRNA-seq data
`popuDE` is an R package for the statistical modeling of the 
gene differential expression (DE) in scRNA-seq data. It identifies cell-type specific genes (markers) that consistently appear in historical population-level scRNA-seq (scRNA-seq) data. `popuDE` is built on top of the R package [`SingleCellExperiment`](https://bioconductor.org/packages/devel/bioc/html/SingleCellExperiment.html) and supports parallel computation.

Except from our proposed method, `popuDE` also provides a common interface for classic DE methods such as the Wilcoxon test, t-test and [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).


## Installation
You can install `popuDE` from [GitHub](https://github.com/luxiao10/cts_gene) using the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package:
```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github('luxiao10/cts_gene', dependencies=T, build_vignettes = T)
library(popuDE)
```




## Getting started
Here we give a simple example to demonstrate how to run `popuDE`. Once the package is installed, you can load the simulated dataset included in the package:
```R
data(sim.sce)
```
`sim.sce` is a `SingleCellExperiment` object with 2,000 genes and 10,000 cells.
Next, you can run `popuDE` with a single line of code:
```R
res <- popuDE(sim.sce, use.raw = TRUE, subject.rep='subject', celltype.rep='celltype', numCores=2)
```
Some explanations about the parameters:
- **use.raw:** Specifies whether to use the raw counts. If set to `TRUE`, the raw counts will be normalized by `popuDE`.
- **subject.rep:** The name of the column that stores subject labels of cells in the `colData` slot.
- **celltype.rep:** The name of the column that stores cell type labels in the `colData` slot.
- **numCores:** Number of cores for parallel computation.

In the [tested environment](##Tested environment), the code finishes running within a minute. The return value `res` is a list containing lists for each cell type. Each list contains
posterior probabilities of genes and parameter estimations for a particular cell type. For example, you can extract the posterior probabilities of genes to show DE in `celltype1` using the following code:
```R
res$celltype1$pp.d1
```
For more details about how to run `popuDE` and classic DE methods, please refer to `vignette("popuDE")`.


## Reproducibility
If you want to reproduce results shown in the paper, please refer to the directory [reproducibility/](reproducibility/) in this repo.


## Tested environment
- CPU: AMD Ryzen Threadripper 3990X 64-Core Processor
- Memory: 256GB
- System: Ubuntu 20.04.6 LTS
- R version: 4.3.0

## Citation
Coming soon.
