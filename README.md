# SWORD

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/SWORD)](https://CRAN.R-project.org/package=SWORD)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**SWORD** (**S**upport vector machine **W**eighted **O**blique **R**andom
**D**ecision forests) fits random forests of oblique regression trees. At each
internal node a Weighted-SVM classifier with a linear kernel defines the
oblique splitting hyperplane as a sparse linear combination of predictors,
rather than a single axis-aligned threshold.

The base learner is **TORS** (**T**ree **O**blique for **R**egression with
weighted **S**VM), introduced in Carta & Frigau (2025).

## Overview

The main functions are:

- `TORS()` — fit a single TORS oblique regression tree.
- `SWORD()` — fit a SWORD forest of TORS trees with bootstrap sampling and
  out-of-bag (OOB) error estimation.
- `predict()` — S3 prediction method for both `tors_flat` and `sword_flat`
  objects; handles factor/character predictors automatically.
- `VI_SWORD()` — Oblique Impurity-Weighted Variable Importance (OIW-VI).
- `pdp_sword()` — partial dependence plots for a single predictor.
- `plot_tors_text()` — text representation of an oblique tree.
- `plot()` — base-R tree graph (`tors_flat`) or OOB convergence curve
  (`sword_flat`).
- `plot_oob_fit()` — OOB predicted vs. observed scatter plot.
- `plot_vi_sword()` — horizontal variable importance barplot.

## Installation

You can install the development version of SWORD from
[GitHub](https://github.com/acarta88/SWORD-rcpp) with:

``` r
# install.packages("remotes")
remotes::install_github("acarta88/SWORD-rcpp")
```

## Usage

We use the **Boston** housing dataset (`MASS` package): 506 observations,
13 numeric predictors, response `medv` (median house value).

``` r
library(SWORD)
data(Boston, package = "MASS")

set.seed(42)
n     <- nrow(Boston)
train <- sample(n, floor(0.8 * n))
test  <- setdiff(seq_len(n), train)

X_train <- Boston[train, setdiff(names(Boston), "medv")]
y_train <- Boston$medv[train]
X_test  <- Boston[test,  setdiff(names(Boston), "medv")]
y_test  <- Boston$medv[test]
```

### Fit a single TORS tree

``` r
set.seed(1)
tree <- TORS(X_train, y_train, nmin = 15, cp = 0.02)
tree
#> ---------- TORS (Tree Oblique for Regression with weighted SVM)
#> 
#>            N observations: 404
#>              N predictors: 13
#>                   N nodes: 29 (15 internal, 14 leaves)
#>                 Max depth: 6
#>          Avg obs in leaf : 28.9
#>          Min obs in leaf : 5
#>     Top-correlated feats: 2
#>              Correlation: Pearson
#>            Weight scheme: scale
#> 
#> -----------------------------------------
```

`plot_tors_text()` prints the tree structure with the oblique hyperplane
equation at each internal node:

``` r
plot_tors_text(tree, top_k = 2)
#> [1] root  n=404  mean=22.53
#>   [2] L  0.47*lstat + 0.88*rm < 0  n=221  mean=17.62
#>     [4] L  0.61*lstat + 0.79*dis < 0  n=116  mean=13.89
#>       ...
#>   [3] R  0.47*lstat + 0.88*rm >= 0  n=183  mean=28.44
#>       ...
```

### Fit a SWORD forest

``` r
set.seed(42)
forest <- SWORD(X_train, y_train, m = 50, nmin = 10, cp = 0.01,
                OOB = TRUE, verbose = FALSE)
forest
#> ---------- SWORD (Support vector machine Weighted Oblique Random Decision forest)
#> 
#>            N observations: 404
#>              N predictors: 13
#>                   N trees: 50
#>      Predictors per split: 13
#>      Top-correlated feats: random [2, 13]
#>       Avg leaves per tree: 18.420
#>          Min obs in leaf : 2
#>            OOB stat type: RMSE / R²
#>           OOB stat value: 3.2415 / 0.8724
#>            Weight scheme: scale
#> 
#> -----------------------------------------
```

### Predict on the test set

``` r
y_hat <- predict(forest, X_test)

cat(sprintf("Test RMSE : %.3f\n", sqrt(mean((y_hat - y_test)^2))))
cat(sprintf("Test MAE  : %.3f\n", mean(abs(y_hat - y_test))))
cat(sprintf("Test R²   : %.3f\n", cor(y_hat, y_test)^2))
#> Test RMSE : 3.247
#> Test MAE  : 2.198
#> Test R²   : 0.871
```

Both `TORS()` and `SWORD()` also accept an R formula:

``` r
df      <- cbind(X_train, medv = y_train)
forest2 <- SWORD(medv ~ ., data = df, m = 50, verbose = FALSE)
predict(forest2, Boston[test, ])   # response column is silently dropped
```

### Variable importance

`VI_SWORD()` computes the **OIW-VI** (Oblique Impurity-Weighted Variable
Importance): the deviance reduction at each node is distributed among
predictors proportionally to their normalised absolute Weighted-SVM
coefficients, then averaged across all trees.

``` r
vi <- VI_SWORD(forest)
print(round(vi, 3))
#>   lstat      rm ptratio     nox     dis    crim     age     tax   indus
#>   0.281   0.218   0.118   0.089   0.078   0.063   0.051   0.038   0.032
#>    chas      zn     rad       b
#>   0.017   0.009   0.004   0.002

plot_vi_sword(forest, top_n = 13)
```

### Partial dependence

`pdp_sword()` marginalises forest predictions over the full training set
at a grid of values for a single predictor.

``` r
par(mfrow = c(1, 2))
pdp_sword(forest, X_train, "lstat", main = "PDP: lstat")
pdp_sword(forest, X_train, "rm",    main = "PDP: rm")
```

### OOB diagnostics

``` r
plot(forest, which = "oob")   # OOB RMSE vs number of trees
plot_oob_fit(forest, y_train) # OOB predicted vs. observed
```

## Getting help

If you encounter a bug, please file an issue with a minimal reproducible
example on [GitHub](https://github.com/acarta88/SWORD-rcpp/issues).

## References

Carta, A. and Frigau, L. (2025). Tree oblique for regression with weighted
support vector machine. *Computational Statistics*, **40**, 5257–5291.
<https://doi.org/10.1007/s00180-025-01647-w>

Yang, X., Song, Q. and Cao, A. (2005). Weighted support vector machine for
data classification. In *Proceedings of the 2005 IEEE International Joint
Conference on Neural Networks (IJCNN)*.
<https://doi.org/10.1109/IJCNN.2005.1555965>

Xu, T. et al. (2024). *WeightSVM: Subject/Instance Weighted Support Vector
Machines*. R package version 1.7-16.
<https://doi.org/10.32614/CRAN.package.WeightSVM>
