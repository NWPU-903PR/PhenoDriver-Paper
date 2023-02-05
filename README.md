# PhenoDriver: Interpretable framework for studying personalized phenotype-associated driver genes in breast cancer

## Introduction
In this repository, we provide the R scripts to reproduce the results and figures and figures of the paper.

## Dependency
Install dependency packages:
```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("org.Hs.eg.db", "DESeq2", "AnnotationDbi", "clusterProfiler")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("NWPU-903PR/PhenoDriverR")
```

## Directory Tree
```
├─data                             [the code of biRFR framework]
| ├─BRCA_NT_Rawcounts.csv          [the gene expression matrix of normal samples for breast cancer]
| ├─BRCA_TP_Rawcounts.csv          [the gene expression matrix of Tumor patients for breast cancer]
| ├─MC3_BRCA.maf.gz                [Breast cancer MC3 mutation data after compress]
| ├─NCBI2Reactome.txt              [NCBI gene ID to Reactome pathway terms mapping (v79)]
| └─SignalTransductionNet.txt      [Signal Transduction Network (STN) edge list file]
├─main.R                           [R script for detecting personalized driver genes]
```

