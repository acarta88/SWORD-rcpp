# SWORD 0.1.0

* First CRAN release.
* Implements TORS (Tree Oblique for Regression with weighted SVM) and the
  SWORD random forest of TORS trees.
* C++ coordinate-descent solver (Rcpp) for fast weighted-SVM splits.
* OIW-VI variable importance, partial dependence plots, and interactive tree
  visualisation via visNetwork.
* Parallel tree fitting via the future/furrr framework (optional).
