# CLIPER: CLusterIng-enabled Peak-to-gEne Regression

CLIPER is a Bayesian framework for **peak-to-gene association** in single-cell multi-omic data.  
It enables clustering of correlated peaks, reduces dimensionality, and provides interpretable **cis-regulatory modules** with effect size estimates.

---

## Workflow

<p align="center">
  <img src="figures/workflow.png" alt="CLIPER workflow" width="700">
</p>

---

## Install required R packages:
```{r}
install.packages(c("data.table", "dplyr", "ggplot2", "GenomicRanges", "glmnet"))
```

## Prepare input data
- Peaks (X): Accessible regions in chr:start-end format
- Gene expression (y): Normalized counts (per metacell or cell)

## Output
- Peak–gene associations
- Cluster labels for peaks
- Estimated effect sizes
