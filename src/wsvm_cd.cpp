// =============================================================================
// wsvm_cd.cpp
// Weighted C-SVM with linear kernel — cyclic coordinate descent on the dual.
//
// Key idea: bias augmentation x̃_i = [x_i; 1] eliminates the equality
// constraint Σ α_i y_i = 0, giving an unconstrained box-only dual:
//
//   max_α  Σ α_i - (1/2) Σ_ij α_i α_j y_i y_j <x̃_i, x̃_j>
//   s.t.   0 ≤ α_i ≤ C · w_i   for all i
//
// Coordinate ascent update for α_i:
//   G_i    = 1 - y_i · (w^T x_i + b)          (dual gradient)
//   Q_ii   = ||x_i||² + 1                       (diagonal of Q matrix)
//   α_i   ← clamp(α_i + G_i / Q_ii, 0, C·w_i)
//
// Primal maintained explicitly: w += Δα_i · y_i · x_i,  b += Δα_i · y_i
// Cost: O(p) per variable update — avoids O(n·p) kernel recomputation.
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>   // std::max, std::min

using namespace Rcpp;

//' Weighted linear C-SVM via coordinate descent (bias augmentation)
//'
//' @param x        Numeric matrix (n x p) of scaled features.
//' @param y_bin    Integer vector of class labels in \{-1, +1\}.
//' @param w_obs    Numeric vector of per-observation weights (positive).
//' @param C        SVM cost parameter (default 1.0).
//' @param max_iter Maximum number of full cyclic passes (default 5000).
//' @param tol      Convergence threshold: stop when max|delta_alpha| < tol (default 1e-3).
//'
//' @return A named list with elements \code{w} (primal weights, length p) and
//'   \code{b} (scalar bias); decision = \code{w^T x_scaled + b}.
//' @references
//' Yang, X., Song, Q. and Cao, A. (2005). Weighted support vector machine for
//' data classification. In \emph{Proceedings of the 2005 IEEE International
//' Joint Conference on Neural Networks (IJCNN)}.
//' \doi{10.1109/IJCNN.2005.1555965}
//'
//' Xu, T. et al. (2024). \emph{WeightSVM: Subject/Instance Weighted Support
//' Vector Machines}. R package version 1.7-16.
//' \doi{10.32614/CRAN.package.WeightSVM}
//' @export
// [[Rcpp::export]]
List wsvm_cd_cpp(NumericMatrix x,
                 IntegerVector y_bin,
                 NumericVector w_obs,
                 double C        = 1.0,
                 int    max_iter = 5000,
                 double tol      = 1e-3) {

  const int n = x.nrow();
  const int p = x.ncol();

  std::vector<double> alpha(n, 0.0);
  std::vector<double> w(p,  0.0);
  double b = 0.0;

  // Precompute Q_ii = ||x_i||^2 + 1  (bias dimension contributes 1^2 = 1)
  std::vector<double> Qii(n);
  for (int i = 0; i < n; i++) {
    double s = 1.0;
    for (int j = 0; j < p; j++) s += x(i, j) * x(i, j);
    Qii[i] = s;
  }

  // Cyclic coordinate ascent
  for (int iter = 0; iter < max_iter; iter++) {
    double max_delta = 0.0;

    for (int i = 0; i < n; i++) {
      const int yi = y_bin[i];

      double f = b;
      for (int j = 0; j < p; j++) f += w[j] * x(i, j);

      const double G     = 1.0 - static_cast<double>(yi) * f;
      const double Ci    = C * w_obs[i];
      const double a_old = alpha[i];
      const double a_new = std::max(0.0, std::min(a_old + G / Qii[i], Ci));
      const double delta = a_new - a_old;

      if (std::abs(delta) < 1e-15) continue;

      alpha[i] = a_new;

      const double yi_d = static_cast<double>(yi) * delta;
      for (int j = 0; j < p; j++) w[j] += yi_d * x(i, j);
      b += yi_d;

      if (std::abs(delta) > max_delta) max_delta = std::abs(delta);
    }

    if (max_delta < tol) break;
  }

  return List::create(Named("w") = NumericVector(w.begin(), w.end()),
                      Named("b") = b);
}
