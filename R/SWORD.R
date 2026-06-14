# ==============================================================================
# SWORD.R
# SWORD: Support vector machine Weighted Oblique Random Decision forests
#
# Self-contained R implementation backed by a C++ coordinate-descent
# Weighted-SVM solver (wsvm_cd_cpp, compiled from src/wsvm_cd.cpp via R CMD INSTALL).
#
# Public API:
#   TORS()      \u2014 fit a TORS oblique regression tree  (returns "tors_flat")
#   SWORD()     \u2014 fit a SWORD forest of TORS trees    (returns "sword_flat")
#   predict()   \u2014 S3 dispatch for tors_flat / sword_flat
#   print()     \u2014 S3 dispatch
#   summary()   \u2014 S3 dispatch
#   VI_SWORD()  \u2014 Oblique Impurity-Weighted Variable Importance (OIW-VI)
#
# Authors: Carta A., Frigau L.
# ==============================================================================

#' @useDynLib SWORD, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor median mad sd var quantile model.frame model.response delete.response na.omit na.pass reformulate model.matrix terms predict setNames
#' @importFrom graphics par axis points legend abline barplot segments
#' @importFrom grDevices adjustcolor
NULL


# ==============================================================================
# INTERNAL HELPERS  (.sword_* prefix \u2014 not exported)
# ==============================================================================

.sword_validate <- function(Covariates, y) {
  if (!is.data.frame(Covariates) && !is.matrix(Covariates))
    stop("'Covariates' must be a data.frame or matrix.")
  if (!is.numeric(y))
    stop("'y' must be a numeric vector.")
  if (length(y) != nrow(Covariates))
    stop("length(y) (", length(y), ") != nrow(Covariates) (", nrow(Covariates), ").")
  if (length(y) < 2L)
    stop("At least 2 observations are required.")
  na_y  <- sum(!is.finite(y))
  # Count non-finite cells per column: numeric/logical columns use is.finite,
  # factor/character columns (one-hot encoded later) use is.na so that valid
  # categorical predictors are not flagged as non-finite.
  if (is.data.frame(Covariates)) {
    na_x <- sum(vapply(Covariates, function(col) {
      if (is.numeric(col)) sum(!is.finite(col))
      else                 sum(is.na(col))
    }, numeric(1L)))
  } else {  # matrix
    na_x <- if (is.numeric(Covariates)) sum(!is.finite(Covariates))
            else                        sum(is.na(Covariates))
  }
  if (na_y > 0L || na_x > 0L)
    stop("Data contains ", na_y, " non-finite value(s) in y and ",
         na_x, " in Covariates. Remove them with na.omit() or complete.cases() first.")
  if (stats::var(y) == 0)
    stop("'y' has zero variance; the response must not be constant.")
  invisible(NULL)
}

.sword_deviance <- function(y) {
  if (length(y) <= 1L) return(0)
  var(y) * (length(y) - 1L)
}

.sword_top_n <- function(scores, n_top) {
  order(scores, decreasing = TRUE)[seq_len(n_top)]
}

.sword_fix_sv_colnames <- function(svm_obj, mat) {
  if (!is.null(svm_obj$SV) && is.null(colnames(svm_obj$SV)))
    colnames(svm_obj$SV) <- colnames(mat)
  svm_obj
}

.sword_coef_unscaled <- function(tree_SVM, Covariates) {

  w         <- t(tree_SVM$coefs) %*% tree_SVM$SV
  intercept <- tree_SVM$rho

  m <- rep(0.0, ncol(Covariates))
  names(m) <- colnames(Covariates)

  if (dim(w)[2L] == 1L) {
    names(w) <- colnames(w)[1L]
    m[names(w)] <- w
    sd_cov       <- apply(Covariates, 2L, sd, na.rm = TRUE)
    m_new        <- as.matrix(m) / sd_cov
    int_new      <- (-intercept -
                       as.matrix(w) %*%
                       as.matrix(colMeans(Covariates[, names(w), drop = FALSE]) /
                                   sd_cov[names(w)]))
    w_unscaled   <- t(m_new)
    int_unscaled <- -as.numeric(int_new)

  } else if (all(!tree_SVM[["scaled"]])) {
    m[colnames(w)] <- w
    w_unscaled     <- t(as.matrix(m))
    int_unscaled   <- -as.numeric(intercept)

  } else {
    m[colnames(w)] <- w
    sd_cov         <- apply(Covariates, 2L, sd, na.rm = TRUE)
    m_new          <- as.matrix(m) / sd_cov
    int_new        <- (-intercept -
                         as.matrix(w) %*%
                         as.matrix(colMeans(Covariates[, colnames(w), drop = FALSE]) /
                                     sd_cov[colnames(w)]))
    w_unscaled     <- t(m_new)
    int_unscaled   <- -as.numeric(int_new)
  }

  coefficients <- c(w_unscaled, int_unscaled)
  # Replace all non-finite values (NA, NaN, +/-Inf): a degenerate node with a
  # constant predictor yields sd = 0 and hence Inf after back-transformation,
  # which is.na() would miss and would propagate as NA into the split test.
  coefficients[!is.finite(coefficients)] <- 0
  names(coefficients) <- c(colnames(w_unscaled), "Int")
  coefficients
}

.sword_coef_scaled <- function(tree_SVM, Covariates) {

  w         <- t(tree_SVM$coefs) %*% tree_SVM$SV
  intercept <- tree_SVM$rho

  m <- rep(0.0, ncol(Covariates))
  names(m) <- colnames(Covariates)

  if (dim(w)[2L] == 1L) {
    names(w) <- colnames(w)[1L]
    m[names(w)] <- w
  } else {
    m[colnames(w)] <- w
  }

  w_scaled   <- t(as.matrix(m))
  int_scaled <- -as.numeric(intercept)

  coefficients <- c(w_scaled, int_scaled)
  coefficients[!is.finite(coefficients)] <- 0
  names(coefficients) <- c(colnames(w_scaled), "Int")
  coefficients
}

.sword_pred_coef <- function(data, coeffs) {
  int_idx   <- which(names(coeffs) == "Int")
  intercept <- coeffs[[int_idx]]
  w         <- coeffs[-int_idx]
  drop(as.matrix(data) %*% w - intercept) < 0
}


# ==============================================================================
# CATEGORICAL ENCODING HELPERS
# ==============================================================================

# Converts factor/character/logical columns to dummy variables via model.matrix.
# Returns list(mat, enc_frm, contrasts, fac_lvls); enc_frm is NULL when no
# encoding is needed (all-numeric input \u2192 mat = as.matrix(X)).
.sword_encode_covariates <- function(X) {
  is_fac <- vapply(X, function(col) is.factor(col) || is.character(col), logical(1L))
  is_log <- vapply(X, is.logical, logical(1L))

  if (!any(is_fac) && !any(is_log))
    return(list(mat = as.matrix(X), enc_frm = NULL, contrasts = NULL,
                fac_lvls = NULL, enc_assign = NULL, enc_var_nms = NULL))

  for (nm in names(X)[is_log]) X[[nm]] <- as.integer(X[[nm]])
  for (nm in names(X)[is_fac]) X[[nm]] <- as.factor(X[[nm]])

  fac_lvls <- lapply(X[is_fac], levels)

  # Backtick-escape names so model.matrix treats them as column names, not calls
  safe_nms <- paste0("`", names(X), "`")
  enc_frm  <- stats::reformulate(safe_nms, intercept = FALSE)
  mm       <- stats::model.matrix(enc_frm, data = X)

  list(mat        = mm,
       enc_frm    = enc_frm,
       contrasts  = attr(mm, "contrasts"),
       fac_lvls   = fac_lvls,
       enc_assign  = attr(mm, "assign"),
       enc_var_nms = names(X))
}

# Applies the encoding stored at fit time to new data, ensuring identical
# factor levels and dummy columns.
.sword_apply_encoding <- function(X, enc) {
  for (nm in names(enc$fac_lvls)) {
    if (nm %in% names(X))
      X[[nm]] <- factor(as.character(X[[nm]]), levels = enc$fac_lvls[[nm]])
  }
  for (nm in names(X))
    if (is.logical(X[[nm]])) X[[nm]] <- as.integer(X[[nm]])

  # Build the design matrix without dropping rows: an unseen factor level
  # becomes NA after the factor() re-level above; with na.action = na.pass the
  # row is retained and its dummy columns are NA, which we then set to 0 so the
  # observation maps to the (all-zero) reference encoding instead of being lost.
  old_na <- getOption("na.action")
  options(na.action = "na.pass")
  on.exit(options(na.action = old_na), add = TRUE)

  mm <- stats::model.matrix(enc$enc_frm, data = X, contrasts.arg = enc$contrasts)
  mm[is.na(mm)] <- 0
  as.data.frame(mm)
}


# ==============================================================================
# C++ SOLVER HELPERS
# ==============================================================================

# Builds a WeightSVM-compatible pseudo-object from the C++ solver output.
# The rest of the pipeline (.sword_coef_unscaled, .sword_coef_scaled) expects:
#   w    <- t(coefs) %*% SV   \u2192  scaled weights
#   rho  <- tree_SVM$rho      \u2192  bias (sign convention: score = w^T x - rho)
#   scaled                    \u2192  triggers back-transformation
.make_pseudo_svm <- function(w_s, b_s, feat_names) {
  sv_mat <- matrix(w_s, nrow = 1L, ncol = length(w_s),
                   dimnames = list(NULL, feat_names))
  list(
    coefs  = matrix(1.0, 1L, 1L),
    SV     = sv_mat,
    rho    = -b_s,
    scaled = rep(TRUE, length(w_s))
  )
}

# Axis-aligned fallback when the C++ solver returns all-zero weights.
# Splits on the feature most correlated with y at its median.
.axis_aligned_fallback <- function(x_mat, y, mu, sigma, feat_names) {
  abs_cors <- abs(drop(cor(y, x_mat)))
  abs_cors[!is.finite(abs_cors)] <- 0
  best_j   <- which.max(abs_cors)
  med_j    <- median(x_mat[, best_j])
  w_s      <- numeric(length(feat_names))
  w_s[best_j] <- 1.0
  b_s      <- -(med_j - mu[best_j]) / sigma[best_j]
  .make_pseudo_svm(w_s, b_s, feat_names)
}

# Dispatch to C++ solver; for nu-classification applies the LIBSVM nu->C
# per-class transformation so the same coordinate-descent solver handles both.
.call_wsvm <- function(x_scaled, y_t, weights,
                       type_of_svm, cost_c, cost_nu, tolerance) {
  weights[!is.finite(weights) | weights <= 0] <- 1e-10
  if (type_of_svm == "nu-classification") {
    n_pos <- sum(y_t ==  1L)
    n_neg <- sum(y_t == -1L)
    C_pos <- if (n_pos > 0L) 1.0 / (cost_nu * n_pos) else 1.0
    C_neg <- if (n_neg > 0L) 1.0 / (cost_nu * n_neg) else 1.0
    eff_w <- ifelse(y_t == 1L, weights * C_pos, weights * C_neg)
    wsvm_cd_cpp(x_scaled, y_t, eff_w, C = 1.0, max_iter = 5000L, tol = tolerance)
  } else {
    wsvm_cd_cpp(x_scaled, y_t, weights, C = cost_c, max_iter = 5000L, tol = tolerance)
  }
}


# ==============================================================================
# .sword_best_split  \u2014 C++ coordinate-descent solver
# ==============================================================================
.sword_best_split <- function(x, y, n_perc, n_top_cor,
                               rf_var, rand_ntopcor, relation,
                               weight_scheme = c("scale", "robust"),
                               type_of_svm   = c("C-classification", "nu-classification"),
                               cost_c    = 1,
                               cost_nu   = 0.5,
                               tolerance = 0.001) {

  x <- as.matrix(x)

  if (!is.null(rf_var) && rf_var >= n_top_cor)
    x <- x[, sample(ncol(x), ceiling(rf_var)), drop = FALSE]

  non_const <- x[, apply(x, 2L, function(col) max(col) != min(col)), drop = FALSE]

  if (is.null(dim(non_const)) || ncol(non_const) == 0L) {

    selected <- x

  } else if (relation == "Pearson") {

    abs_cors <- abs(drop(cor(y, non_const)))
    names(abs_cors) <- colnames(non_const)
    abs_cors[!is.finite(abs_cors)] <- 0

    thresh_pass <- which(abs_cors >= 1)
    if (length(thresh_pass) >= 1L) {
      selected <- non_const[, thresh_pass, drop = FALSE]
    } else {
      n_k <- if (rand_ntopcor) sample(seq(2L, rf_var), 1L) else n_top_cor
      selected <- non_const[,
        .sword_top_n(abs_cors, min(n_k, ncol(non_const))),
        drop = FALSE]
    }

  } else {  # MI

    if (!requireNamespace("infotheo", quietly = TRUE))
      stop("Package 'infotheo' is required for relation = 'MI'. ",
           "Install with: install.packages('infotheo')")
    disc_x <- infotheo::discretize(non_const)
    disc_y <- infotheo::discretize(data.frame(y))
    mi_scores <- vapply(seq_len(ncol(disc_x)), function(i)
      infotheo::mutinformation(disc_x[, i, drop = FALSE], disc_y),
      numeric(1L))
    n_k <- if (rand_ntopcor) sample(seq(2L, rf_var), 1L) else n_top_cor
    selected <- non_const[,
      .sword_top_n(mi_scores, min(n_k, ncol(non_const))),
      drop = FALSE]
  }

  x_mat      <- selected
  feat_names <- colnames(x_mat)
  n_obs      <- nrow(x_mat)

  # Scale x (same convention as WeightSVM's internal scale = TRUE)
  mu    <- colMeans(x_mat)
  sigma <- apply(x_mat, 2L, sd)
  sigma[sigma == 0 | !is.finite(sigma)] <- 1.0
  x_scaled <- sweep(sweep(x_mat, 2L, mu, "-"), 2L, sigma, "/")

  # Observation weights
  if (weight_scheme == "scale") {
    weights <- as.numeric(abs(scale(y)))
    weights[weights == 0 | !is.finite(weights)] <- 1e-10
  } else {
    med_y <- median(y)
    s     <- mad(y)
    if (!is.finite(s) || s == 0) s <- sd(y)
    if (!is.finite(s) || s == 0) s <- 1
    weights <- abs(y - med_y) / s
    weights[weights == 0 | !is.finite(weights)] <- 1e-10
  }

  if (n_perc == 1L) {

    t   <- median(y, na.rm = TRUE)
    y_t <- as.integer(y <= t) * 2L - 1L
    if (length(unique(y_t)) == 1L)
      y_t <- as.integer(y < t) * 2L - 1L

    res <- .call_wsvm(x_scaled, y_t, weights,
                      type_of_svm, cost_c, cost_nu, tolerance)
    if (all(res$w == 0))
      return(.axis_aligned_fallback(x_mat, y, mu, sigma, feat_names))
    return(.make_pseudo_svm(res$w, res$b, feat_names))
  }

  # n_perc > 1: scan candidate thresholds, pick min deviance
  thr_all  <- unique(quantile(y,
    seq(0, 1, max(1 / (n_perc + 1), 1 / (n_obs + 1))),
    na.rm = TRUE))
  thr_cand <- thr_all[-length(thr_all)]

  deviances <- vapply(thr_cand, function(t) {
    y_t <- as.integer(y <= t) * 2L - 1L
    if (length(unique(y_t)) == 1L) return(.sword_deviance(y))
    res   <- .call_wsvm(x_scaled, y_t, weights,
                        type_of_svm, cost_c, cost_nu, tolerance)
    score <- drop(x_scaled %*% res$w) + res$b
    .sword_deviance(y[score < 0]) + .sword_deviance(y[score >= 0])
  }, numeric(1L))

  t_best <- thr_cand[which.min(deviances)]
  y_t    <- as.integer(y <= t_best) * 2L - 1L
  res    <- .call_wsvm(x_scaled, y_t, weights,
                       type_of_svm, cost_c, cost_nu, tolerance)
  if (all(res$w == 0))
    return(.axis_aligned_fallback(x_mat, y, mu, sigma, feat_names))
  .make_pseudo_svm(res$w, res$b, feat_names)
}


# ==============================================================================
# TREE INTERNALS  (unchanged from flat version)
# ==============================================================================

.sword_node_split <- function(x_sub, y_sub,
                              n_perc, n_top_cor,
                              rf_var, rand_ntopcor, relation,
                              weight_scheme, type_of_svm, cost_c, cost_nu,
                              tolerance = 0.001) {

  tree_SVM      <- .sword_best_split(x_sub, y_sub,
                                     n_perc = n_perc, n_top_cor = n_top_cor,
                                     rf_var = rf_var, rand_ntopcor = rand_ntopcor,
                                     relation = relation,
                                     weight_scheme = weight_scheme,
                                     type_of_svm = type_of_svm,
                                     cost_c = cost_c, cost_nu = cost_nu,
                                     tolerance = tolerance)
  coef_unscaled <- .sword_coef_unscaled(tree_SVM, x_sub)
  go_left       <- .sword_pred_coef(x_sub, coef_unscaled)

  list(
    go_left       = go_left,
    n_left        = sum( go_left),
    n_right       = sum(!go_left),
    tree_SVM      = tree_SVM,
    coef_unscaled = coef_unscaled
  )
}

.sword_vi_tree <- function(tree) {
  feat_names <- tree$feature_names
  int_col    <- which(colnames(tree$scaled_coeffs) == "Int")
  internal   <- which(!tree$is_leaf)

  if (length(internal) == 0L)
    return(setNames(rep(0.0, length(feat_names)), feat_names))

  coef_mat  <- abs(tree$scaled_coeffs[internal, -int_col, drop = FALSE])
  row_sums  <- rowSums(coef_mat)
  row_sums[row_sums == 0] <- 1
  norm_mat  <- coef_mat / row_sums

  delta_dev <- tree$dev_presplit[internal] - tree$dev[internal]
  delta_dev[is.na(delta_dev)] <- 0

  vi <- colSums(norm_mat * delta_dev, na.rm = TRUE)
  vi[feat_names]
}

.sword_grow_flat <- function(x_mat, y,
                              nmin, minleaf, cp, n_perc, n_top_cor,
                              original_deviance,
                              rf_var, rand_ntopcor, relation,
                              weight_scheme, type_of_svm, cost_c, cost_nu,
                              tolerance) {

  feature_names <- colnames(x_mat)
  all_names     <- c(feature_names, "Int")
  p_total       <- length(all_names)

  max_nodes <- 2L * ceiling(nrow(x_mat) / max(minleaf, 1L)) + 10L

  is_leaf_v   <- logical(max_nodes)
  left_ch     <- integer(max_nodes)
  right_ch    <- integer(max_nodes)
  leaf_mean_v <- rep(NA_real_, max_nodes)
  dev_pre_v   <- rep(NA_real_, max_nodes)
  dev_v       <- numeric(max_nodes)
  n_obs_v     <- integer(max_nodes)
  depth_v     <- integer(max_nodes)
  coeffs_m    <- matrix(0.0, nrow = max_nodes, ncol = p_total,
                         dimnames = list(NULL, all_names))
  sc_m        <- matrix(0.0, nrow = max_nodes, ncol = p_total,
                         dimnames = list(NULL, all_names))

  node_count <- 0L
  stack      <- vector("list", 64L)
  stack_top  <- 1L
  stack[[1L]] <- list(seq_len(nrow(x_mat)), 0L, 0L, TRUE)

  while (stack_top > 0L) {

    frame     <- stack[[stack_top]]
    stack_top <- stack_top - 1L

    row_idx   <- frame[[1L]]
    depth     <- frame[[2L]]
    parent_id <- frame[[3L]]
    is_left_c <- frame[[4L]]

    node_count <- node_count + 1L
    node_id    <- node_count

    if (parent_id > 0L) {
      if (is_left_c) left_ch[parent_id]  <- node_id
      else           right_ch[parent_id] <- node_id
    }

    y_sub  <- y[row_idx]
    n_node <- length(row_idx)
    n_obs_v[node_id] <- n_node
    depth_v[node_id] <- depth

    if (length(unique(y_sub)) == 1L || n_node < nmin) {
      is_leaf_v[node_id]   <- TRUE
      leaf_mean_v[node_id] <- mean(y_sub)
      dev_v[node_id]       <- .sword_deviance(y_sub)
      next
    }

    x_sub <- x_mat[row_idx, , drop = FALSE]
    nd    <- .sword_node_split(x_sub, y_sub,
                               n_perc = n_perc, n_top_cor = n_top_cor,
                               rf_var = rf_var, rand_ntopcor = rand_ntopcor,
                               relation = relation, weight_scheme = weight_scheme,
                               type_of_svm = type_of_svm, cost_c = cost_c,
                               cost_nu = cost_nu, tolerance = tolerance)

    y_left   <- y_sub[ nd$go_left]
    y_right  <- y_sub[!nd$go_left]
    dev_post <- .sword_deviance(y_left) + .sword_deviance(y_right)
    cp_node  <- dev_post / original_deviance

    # Defensive: a degenerate split (e.g. all-constant predictors at the node)
    # can yield NA membership; treat such a node as a terminal leaf rather than
    # letting NA reach the comparison below.
    no_valid_split <- anyNA(nd$go_left) || is.na(nd$n_left) ||
                      is.na(nd$n_right) || is.na(cp_node)

    if (no_valid_split ||
        nd$n_left < minleaf || nd$n_right < minleaf || cp_node < cp) {
      is_leaf_v[node_id]   <- TRUE
      leaf_mean_v[node_id] <- mean(y_sub)
      dev_v[node_id]       <- .sword_deviance(y_sub)
      next
    }

    co <- nd$coef_unscaled
    sc <- .sword_coef_scaled(nd$tree_SVM, x_sub)

    is_leaf_v[node_id]  <- FALSE
    dev_v[node_id]      <- dev_post
    dev_pre_v[node_id]  <- .sword_deviance(y_sub)
    coeffs_m[node_id, ] <- co
    sc_m[node_id, ]     <- sc

    left_idx  <- row_idx[ nd$go_left]
    right_idx <- row_idx[!nd$go_left]

    stack_top <- stack_top + 1L
    if (stack_top > length(stack)) length(stack) <- length(stack) + 32L
    stack[[stack_top]] <- list(right_idx, depth + 1L, node_id, FALSE)

    stack_top <- stack_top + 1L
    if (stack_top > length(stack)) length(stack) <- length(stack) + 32L
    stack[[stack_top]] <- list(left_idx,  depth + 1L, node_id, TRUE)
  }

  N <- node_count
  structure(
    list(
      n_nodes       = N,
      feature_names = feature_names,
      is_leaf       = is_leaf_v[seq_len(N)],
      left_child    = left_ch[seq_len(N)],
      right_child   = right_ch[seq_len(N)],
      leaf_mean     = leaf_mean_v[seq_len(N)],
      coeffs        = coeffs_m[seq_len(N), , drop = FALSE],
      scaled_coeffs = sc_m[seq_len(N), , drop = FALSE],
      dev_presplit  = dev_pre_v[seq_len(N)],
      dev           = dev_v[seq_len(N)],
      n_obs         = n_obs_v[seq_len(N)],
      depth         = depth_v[seq_len(N)]
    ),
    class = "tors_flat"
  )
}


# ==============================================================================
# PUBLIC API
# ==============================================================================

#' Fit a TORS oblique regression tree
#'
#' Grows one Tree Oblique for Regression with weighted SVM (TORS). At each
#' internal node a Weighted-SVM classifier with a linear kernel defines the
#' oblique splitting hyperplane as a sparse linear combination of predictors.
#' The tree is stored in a flat (parallel vector) representation for memory
#' efficiency and fast prediction.
#'
#' @param Covariates   data.frame of predictors (no response column), or a
#'   \code{formula} for the formula interface.
#' @param y            numeric response vector (length = nrow(Covariates)).
#' @param data         data.frame; required when \code{Covariates} is a formula.
#' @param nmin         minimum node size to attempt a split (default 5).
#' @param minleaf      minimum leaf size after a split (default round(nmin/3)).
#' @param cp           complexity parameter: minimum relative post-split deviance
#'                     (default 0.01; use 0 to grow the full tree).
#' @param n_perc       number of candidate split thresholds to evaluate per node
#'   (default 1). When \code{n_perc = 1} the response is dichotomised at its
#'   median (the recommended setting). When \code{n_perc > 1} the algorithm
#'   evaluates that many evenly-spaced quantile thresholds and picks the one
#'   that minimises node deviance, at the cost of higher computation.
#' @param n_top_cor     number of top-correlated features per split (default 2).
#' @param rf_var       features randomly sub-sampled per split (default ncol(Covariates)).
#' @param rand_ntopcor if TRUE, randomly draw the cardinality of top features (default FALSE).
#' @param relation     feature selection criterion: \code{"Pearson"} or \code{"MI"}.
#' @param weight_scheme observation weighting: \code{"scale"} or \code{"robust"}.
#' @param type_of_svm  SVM type: \code{"C-classification"} uses the cost parameter
#'   \code{cost_c}; \code{"nu-classification"} uses \code{cost_nu} (fraction of
#'   margin errors/support vectors, in (0, 1]).
#' @param cost_c       SVM cost parameter C; used when
#'   \code{type_of_svm = "C-classification"} (default 1).
#' @param cost_nu      nu parameter; used when
#'   \code{type_of_svm = "nu-classification"} (default 0.5).
#' @param tolerance    convergence tolerance for the C++ solver (default 0.001).
#'
#' @return An object of class \code{"tors_flat"}.
#' @references
#' Carta, A. and Frigau, L. (2025). Tree oblique for regression with weighted
#' support vector machine. \emph{Computational Statistics}, \bold{40}, 5257--5291.
#' \doi{10.1007/s00180-025-01647-w}
#'
#' Yang, X., Song, Q. and Cao, A. (2005). Weighted support vector machine for
#' data classification. In \emph{Proceedings of the 2005 IEEE International
#' Joint Conference on Neural Networks (IJCNN)}.
#' \doi{10.1109/IJCNN.2005.1555965}
#'
#' Xu, T. et al. (2024). \emph{WeightSVM: Subject/Instance Weighted Support
#' Vector Machines}. R package version 1.7-16.
#' \doi{10.32614/CRAN.package.WeightSVM}
#' @seealso \code{\link{SWORD}}, \code{\link{predict.tors_flat}}
#' @examples
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' tree <- TORS(X, y, nmin = 10, cp = 0.01)
#' print(tree)
#'
#' # formula interface
#' tree2 <- TORS(y ~ ., data = cbind(X, y = y), nmin = 10, cp = 0.01)
#' @export
TORS <- function(Covariates, y = NULL,
                 data          = NULL,
                 nmin          = 5,
                 minleaf       = round(nmin / 3),
                 cp            = 0.01,
                 n_perc        = 1,
                 n_top_cor      = 2,
                 rf_var        = NULL,
                 rand_ntopcor  = FALSE,
                 relation      = "Pearson",
                 weight_scheme = "scale",
                 type_of_svm   = "C-classification",
                 cost_c        = 1,
                 cost_nu       = 0.5,
                 tolerance     = 0.001) {

  # Formula interface: TORS(y ~ x1 + x2, data = df)
  terms_obj <- NULL
  if (inherits(Covariates, "formula")) {
    if (is.null(data))
      stop("'data' must be provided when using the formula interface.")
    mf         <- stats::model.frame(Covariates, data = data, na.action = stats::na.omit)
    terms_obj  <- attr(mf, "terms")
    y          <- stats::model.response(mf)
    Covariates <- mf[, -1L, drop = FALSE]
  }

  .sword_validate(Covariates, y)
  enc   <- .sword_encode_covariates(Covariates)
  x_mat <- enc$mat
  if (is.null(rf_var)) rf_var <- ncol(x_mat)

  tree <- .sword_grow_flat(
    x_mat             = x_mat,
    y                 = y,
    nmin              = nmin,
    minleaf           = minleaf,
    cp                = cp,
    n_perc            = n_perc,
    n_top_cor          = n_top_cor,
    original_deviance = .sword_deviance(y),
    rf_var            = rf_var,
    rand_ntopcor      = rand_ntopcor,
    relation          = relation,
    weight_scheme     = weight_scheme,
    type_of_svm       = type_of_svm,
    cost_c            = cost_c,
    cost_nu           = cost_nu,
    tolerance         = tolerance
  )
  tree$call_params <- list(
    n             = nrow(x_mat),
    p             = ncol(x_mat),
    nmin          = nmin,
    minleaf       = minleaf,
    n_top_cor      = n_top_cor,
    rf_var        = rf_var,
    rand_ntopcor  = rand_ntopcor,
    relation      = relation,
    weight_scheme = weight_scheme
  )
  tree$terms    <- terms_obj
  tree$encoding <- if (is.null(enc$enc_frm)) NULL else enc[c("enc_frm", "contrasts", "fac_lvls", "enc_assign", "enc_var_nms")]
  tree
}


#' Fit a SWORD forest
#'
#' Fits SWORD (Support vector machine Weighted Oblique Random Decision forests),
#' a random forest of TORS oblique regression trees. Each tree is trained on a
#' bootstrap sample; diversity is further increased by random predictor subspacing
#' and optional randomisation of the number of features per split. Out-of-bag
#' (OOB) predictions and metrics are computed when \code{oob = TRUE}. Parallel
#' fitting is supported via \code{furrr} when the caller sets a \code{future}
#' plan before invoking \code{SWORD}.
#'
#' @param Covariates   data.frame of predictors, or a \code{formula} for the
#'   formula interface.
#' @param y            numeric response vector.
#' @param data         data.frame; required when \code{Covariates} is a formula.
#' @param nmin         minimum node size (default 5).
#' @param minleaf      minimum leaf size (default 2).
#' @param cp           complexity parameter (default 0). Unlike \code{\link{TORS}},
#'   forests use \code{cp = 0} (full trees) by design: bootstrap sampling provides
#'   regularisation and unpruned trees maximise diversity.
#' @param n_top_cor     top-correlated features per split (default 2).
#' @param m            number of trees (default 100).
#' @param rf_var       features sub-sampled per split (default ncol(Covariates)).
#' @param rand_ntopcor randomly draw top-feature cardinality (default TRUE).
#' @param relation     \code{"Pearson"} or \code{"MI"}.
#' @param weight_scheme \code{"scale"} or \code{"robust"}.
#' @param type_of_svm  SVM type: \code{"C-classification"} uses \code{cost_c};
#'   \code{"nu-classification"} uses \code{cost_nu} (default \code{"C-classification"}).
#' @param cost_c       SVM cost parameter C; used when
#'   \code{type_of_svm = "C-classification"} (default 1).
#' @param cost_nu      nu parameter; used when
#'   \code{type_of_svm = "nu-classification"} (default 0.5).
#' @param tolerance    convergence tolerance for the C++ solver (default 0.001).
#' @param seed_bs      integer seed offset for bootstrap samples (default 25).
#' @param parallel     if TRUE, fit trees in parallel via \code{furrr::future_map}
#'                     (default FALSE). Set a \code{future} plan before calling.
#' @param chunk        if TRUE, process trees in chunks (default FALSE).
#' @param n_chunks     number of chunks; NULL = ceiling(m / n_workers) (default NULL).
#' @param n_workers    number of workers used only to derive the default
#'   \code{n_chunks} when \code{chunk = TRUE} and \code{n_chunks = NULL}
#'   (default 1). The actual parallelism is controlled by the \code{future}
#'   plan set by the caller.
#' @param oob          if TRUE, compute OOB predictions and metrics (default TRUE).
#' @param verbose      if TRUE, print tree index during sequential fitting (default FALSE).
#' @param timeout_tree max seconds per tree; NULL = no limit (default 300). Trees
#'   that exceed the limit are skipped and excluded from aggregation; the count
#'   of skipped trees is stored in \code{object$n_skipped} and printed by
#'   \code{\link{print.sword_flat}}.
#'
#' @details
#' \strong{Split threshold.} SWORD always uses \code{n_perc = 1} (median
#' dichotomisation) for all trees. This is a deliberate design choice: median
#' splits are faster and, combined with bootstrap sampling and random-\eqn{\gamma}
#' selection, provide sufficient diversity without the overhead of a quantile
#' search. Users who require multi-threshold splits should use \code{\link{TORS}}
#' directly with \code{n_perc > 1}.
#'
#' \strong{Reproducibility.} SWORD controls reproducibility through \code{seed_bs}
#' rather than the caller's ambient RNG. Internally, bootstrap sample \eqn{i} is
#' drawn with \code{set.seed(i + seed_bs)}, and the caller's \code{.Random.seed}
#' is restored on exit. This means two calls with the same \code{seed_bs} always
#' produce identical forests regardless of any outer \code{set.seed()}. To obtain
#' a different forest, change \code{seed_bs}.
#'
#' \strong{Missing values.} The formula interface applies \code{na.omit} to the
#' model frame. The data.frame interface expects complete cases; rows with
#' \code{NA} in \code{Covariates} or \code{y} trigger an error with a count of
#' the affected rows. Use \code{\link[stats]{na.omit}} or
#' \code{\link[stats]{complete.cases}} to pre-clean the data if needed.
#'
#' @return An object of class \code{"sword_flat"}.
#' @references
#' Carta, A. and Frigau, L. (2025). Tree oblique for regression with weighted
#' support vector machine. \emph{Computational Statistics}, \bold{40}, 5257--5291.
#' \doi{10.1007/s00180-025-01647-w}
#'
#' Yang, X., Song, Q. and Cao, A. (2005). Weighted support vector machine for
#' data classification. In \emph{Proceedings of the 2005 IEEE International
#' Joint Conference on Neural Networks (IJCNN)}.
#' \doi{10.1109/IJCNN.2005.1555965}
#'
#' Xu, T. et al. (2024). \emph{WeightSVM: Subject/Instance Weighted Support
#' Vector Machines}. R package version 1.7-16.
#' \doi{10.32614/CRAN.package.WeightSVM}
#' @seealso \code{\link{TORS}}, \code{\link{predict.sword_flat}}, \code{\link{VI_SWORD}}
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = TRUE, verbose = FALSE)
#' print(forest)
#'
#' # formula interface
#' forest2 <- SWORD(y ~ ., data = cbind(X, y = y), m = 10, verbose = FALSE)
#' }
#' @export
SWORD <- function(
    Covariates,
    y             = NULL,
    data          = NULL,
    nmin          = 5,
    minleaf       = 2,
    cp            = 0.0,
    n_top_cor      = 2,
    m             = 100,
    rf_var        = NULL,
    rand_ntopcor  = TRUE,
    relation      = "Pearson",
    weight_scheme = "scale",
    type_of_svm   = "C-classification",
    cost_c        = 1,
    cost_nu       = 0.5,
    tolerance     = 0.001,
    seed_bs       = 25,
    parallel      = FALSE,
    chunk         = FALSE,
    n_chunks      = NULL,
    n_workers     = 1,
    oob           = TRUE,
    verbose       = FALSE,
    timeout_tree  = 300
) {

  # Formula interface: SWORD(y ~ x1 + x2, data = df)
  terms_obj <- NULL
  if (inherits(Covariates, "formula")) {
    if (is.null(data))
      stop("'data' must be provided when using the formula interface.")
    mf         <- stats::model.frame(Covariates, data = data, na.action = stats::na.omit)
    terms_obj  <- attr(mf, "terms")
    y          <- stats::model.response(mf)
    Covariates <- mf[, -1L, drop = FALSE]
  }

  .sword_validate(Covariates, y)
  enc           <- .sword_encode_covariates(Covariates)
  Covariates    <- as.data.frame(enc$mat)
  feature_names <- colnames(Covariates)
  if (is.null(rf_var)) rf_var <- ncol(Covariates)
  n <- nrow(Covariates)

  call_params <- list(
    n             = n,
    p             = length(feature_names),
    m             = m,
    rf_var        = rf_var,
    n_top_cor      = n_top_cor,
    rand_ntopcor  = rand_ntopcor,
    minleaf       = minleaf,
    relation      = relation,
    weight_scheme = weight_scheme
  )

  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    .GlobalEnv$.Random.seed else NULL
  on.exit({
    if (is.null(rng_state)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", rng_state, envir = .GlobalEnv)
    }
  }, add = TRUE)

  index_list <- lapply(seq_len(m), function(i) {
    set.seed(i + seed_bs)
    boot_idx <- sample(n, n, replace = TRUE)
    list(boot = boot_idx,
         oob  = setdiff(seq_len(n), boot_idx),
         seed = i + seed_bs)
  })

  build_tree <- function(idx) {
    set.seed(idx$seed)

    baggingData   <- Covariates[idx$boot, , drop = FALSE]
    oobData       <- Covariates[idx$oob,  , drop = FALSE]
    y_bagging     <- y[idx$boot]

    msg <- NULL
    t0  <- proc.time()
    tree <- tryCatch({
      if (!is.null(timeout_tree)) setTimeLimit(elapsed = timeout_tree, transient = TRUE)
      res <- TORS(
        baggingData, y_bagging,
        nmin = nmin, minleaf = minleaf, cp = cp,
        n_perc = 1L, n_top_cor = n_top_cor,
        rf_var = rf_var, rand_ntopcor = rand_ntopcor,
        relation = relation, weight_scheme = weight_scheme,
        type_of_svm = type_of_svm, cost_c = cost_c, cost_nu = cost_nu,
        tolerance = tolerance
      )
      if (!is.null(timeout_tree)) setTimeLimit(elapsed = Inf, transient = FALSE)
      res
    }, error = function(e) {
      setTimeLimit(elapsed = Inf, transient = FALSE)
      msg <<- if (grepl("time|elapsed", conditionMessage(e), ignore.case = TRUE))
        sprintf("  [TIMEOUT] seed=%d exceeded %gs \u2014 skipped", idx$seed, timeout_tree)
      else
        sprintf("  [ERROR] seed=%d: %s", idx$seed, conditionMessage(e))
      NULL
    })

    t_tree <- as.numeric((proc.time() - t0)[["elapsed"]])

    if (is.null(tree))
      return(list(tree = NULL, oob_preds = rep(NA_real_, n),
                  time_tree = t_tree, time_oob = 0, skipped = TRUE, msg = msg))

    oob_preds <- rep(NA_real_, n)
    t_oob     <- 0
    if (oob && length(idx$oob) > 0L) {
      t1                 <- proc.time()
      oob_preds[idx$oob] <- predict(tree, oobData)
      t_oob              <- as.numeric((proc.time() - t1)[["elapsed"]])
    }

    list(tree = tree, oob_preds = oob_preds,
         time_tree = t_tree, time_oob = t_oob, skipped = FALSE, msg = NULL)
  }

  # ---------------------------------------------------------------------------
  # PARALLEL branch
  # ---------------------------------------------------------------------------
  if (parallel) {

    if (!requireNamespace("furrr", quietly = TRUE))
      stop("Package 'furrr' is required for parallel fitting. ",
           "Install with: install.packages('furrr')")
    if (!requireNamespace("progressr", quietly = TRUE))
      stop("Package 'progressr' is required for parallel fitting. ",
           "Install with: install.packages('progressr')")

    furrr_opts <- furrr::furrr_options(
      # seed = NULL: disables future's RNG-misuse warning without changing the
      # worker's RNG kind.  Reproducibility is managed explicitly by
      # set.seed(idx$seed) inside build_tree, which uses the worker's default
      # Mersenne-Twister — identical to the sequential branch.
      # (seed = TRUE would switch workers to L'Ecuyer-CMRG, breaking seq == par;
      #  seed = FALSE leaves the check active and emits warnings.)
      seed        = NULL,
      globals     = TRUE,
      scheduling  = Inf,
      chunk_size  = NULL
    )

    if (chunk) {
      if (is.null(n_chunks)) n_chunks <- ceiling(m / n_workers)
      chunks  <- split(index_list,
                       ceiling(seq_along(index_list) / ceiling(m / n_chunks)))
      results <- list()
      progressr::with_progress({
        for (ch in chunks) {
          results <- c(results, furrr::future_map(
            ch, build_tree,
            .options  = furrr_opts,
            .progress = TRUE
          ))
        }
      })
    } else {
      progressr::with_progress({
        results <- furrr::future_map(
          index_list, build_tree,
          .options  = furrr_opts,
          .progress = TRUE
        )
      })
    }

    skipped   <- vapply(results, function(r) isTRUE(r$skipped), logical(1L))
    n_skipped <- sum(skipped)

    # Print messages captured inside workers (message() is unreliable in workers)
    for (msg_txt in Filter(Negate(is.null), lapply(results, `[[`, "msg")))
      message(msg_txt)
    if (n_skipped > 0L) {
    message("  Skipped trees: ", n_skipped, " / ", m)
    if (n_skipped / m >= 0.1)
      warning(n_skipped, " / ", m, " trees were skipped (timeout or error). ",
              "Consider increasing 'timeout_tree' or reducing data complexity.",
              call. = FALSE)
  }

    treeList  <- lapply(results[!skipped], `[[`, "tree")
    time_tree <- vapply(results, `[[`, numeric(1L), "time_tree")
    time_oob  <- vapply(results, `[[`, numeric(1L), "time_oob")

    enc_out <- if (is.null(enc$enc_frm)) NULL else
      enc[c("enc_frm", "contrasts", "fac_lvls", "enc_assign", "enc_var_nms")]

    if (oob && any(!skipped)) {
      # Build n x m oob matrix (all trees, skipped trees are all-NA columns)
      oob_matrix_full <- do.call(cbind, lapply(results, `[[`, "oob_preds"))

      # Overall oob predictions: average over all trees where obs was oob
      oob_counts <- rowSums(!is.na(oob_matrix_full))
      oob_sums   <- rowSums(oob_matrix_full, na.rm = TRUE)
      final_oob  <- ifelse(oob_counts > 0L, oob_sums / oob_counts, NA_real_)

      # Convergence curve: use only valid (non-skipped) columns, in submission order
      valid_cols <- colSums(!is.na(oob_matrix_full)) > 0L
      oob_matrix <- oob_matrix_full[, valid_cols, drop = FALSE]
      m_valid    <- ncol(oob_matrix)
      oob_errors <- vapply(seq_len(m_valid), function(k) {
        partial  <- oob_matrix[, seq_len(k), drop = FALSE]
        row_mean <- rowSums(partial, na.rm = TRUE) / rowSums(!is.na(partial))
        mean((row_mean - y)^2, na.rm = TRUE)
      }, numeric(1L))

      mse_oob <- mean((final_oob - y)^2, na.rm = TRUE)
      return(structure(list(
        trees               = treeList,
        feature_names       = feature_names,
        call_params         = call_params,
        terms               = terms_obj,
        encoding            = enc_out,
        y_train             = y,
        oob_predictions     = final_oob,
        MSE                 = mse_oob,
        RMSE                = sqrt(mse_oob),
        MAE                 = mean(abs(final_oob - y), na.rm = TRUE),
        Rsquared            = cor(final_oob, y, use = "complete.obs")^2,
        oob_matrix          = oob_matrix,
        oob_errors_per_iter = oob_errors,
        time_tree           = time_tree,
        time_oob            = time_oob,
        n_skipped           = n_skipped
      ), class = "sword_flat"))
    }

    return(structure(list(
      trees         = treeList,
      feature_names = feature_names,
      call_params   = call_params,
      terms         = terms_obj,
      encoding      = enc_out,
      time_tree     = time_tree,
      time_oob      = time_oob,
      n_skipped     = n_skipped
    ), class = "sword_flat"))
  }

  # ---------------------------------------------------------------------------
  # SEQUENTIAL branch
  # ---------------------------------------------------------------------------
  treeList  <- vector("list", m)
  time_tree <- numeric(m)
  time_oob  <- numeric(m)
  n_skipped <- 0L

  oob_sum    <- rep(0.0, n)
  oob_counts <- rep(0L,  n)
  oob_matrix <- matrix(NA_real_, nrow = n, ncol = m)

  for (i in seq_len(m)) {
    if (verbose) message(i)
    idx      <- index_list[[i]]
    boot_idx <- idx$boot
    oob_idx  <- idx$oob

    baggingData    <- Covariates[boot_idx, , drop = FALSE]
    oobData        <- Covariates[oob_idx,  , drop = FALSE]
    y_bagging      <- y[boot_idx]

    set.seed(idx$seed)   # mirror build_tree: each tree starts from the same RNG state
    t0 <- proc.time()
    tree <- tryCatch({
      if (!is.null(timeout_tree)) setTimeLimit(elapsed = timeout_tree, transient = TRUE)
      res <- TORS(
        baggingData, y_bagging,
        nmin = nmin, minleaf = minleaf, cp = cp,
        n_perc = 1L, n_top_cor = n_top_cor,
        rf_var = rf_var, rand_ntopcor = rand_ntopcor,
        relation = relation, weight_scheme = weight_scheme,
        type_of_svm = type_of_svm, cost_c = cost_c, cost_nu = cost_nu,
        tolerance = tolerance
      )
      if (!is.null(timeout_tree)) setTimeLimit(elapsed = Inf, transient = FALSE)
      res
    }, error = function(e) {
      setTimeLimit(elapsed = Inf, transient = FALSE)
      if (grepl("time|elapsed", conditionMessage(e), ignore.case = TRUE))
        message("  [TIMEOUT] Tree ", i, " exceeded ", timeout_tree, "s \u2014 skipped")
      else
        message("  [ERROR] Tree ", i, ": ", conditionMessage(e))
      NULL
    })

    time_tree[i] <- as.numeric((proc.time() - t0)[["elapsed"]])

    if (is.null(tree)) { n_skipped <- n_skipped + 1L; next }

    treeList[[i]] <- tree

    if (oob && length(oob_idx) > 0L) {
      t1              <- proc.time()
      oob_preds       <- predict(tree, oobData)
      time_oob[i]     <- as.numeric((proc.time() - t1)[["elapsed"]])
      oob_sum[oob_idx]       <- oob_sum[oob_idx]    + oob_preds
      oob_counts[oob_idx]    <- oob_counts[oob_idx] + 1L
      oob_matrix[oob_idx, i] <- oob_preds
    } else {
      time_oob[i] <- 0
    }
  }

  if (n_skipped > 0L) {
    message("  Skipped trees: ", n_skipped, " / ", m)
    if (n_skipped / m >= 0.1)
      warning(n_skipped, " / ", m, " trees were skipped (timeout or error). ",
              "Consider increasing 'timeout_tree' or reducing data complexity.",
              call. = FALSE)
  }

  treeList   <- Filter(Negate(is.null), treeList)
  valid_cols <- colSums(!is.na(oob_matrix)) > 0L
  oob_matrix <- oob_matrix[, valid_cols, drop = FALSE]

  base_out <- list(
    trees         = treeList,
    feature_names = feature_names,
    call_params   = call_params,
    terms         = terms_obj,
    encoding      = if (is.null(enc$enc_frm)) NULL else enc[c("enc_frm", "contrasts", "fac_lvls", "enc_assign", "enc_var_nms")],
    time_tree     = time_tree,
    time_oob      = time_oob,
    n_skipped     = n_skipped
  )

  if (oob && ncol(oob_matrix) > 0L) {
    final_oob <- ifelse(oob_counts > 0L, oob_sum / oob_counts, NA_real_)

    m_valid    <- ncol(oob_matrix)
    oob_errors <- vapply(seq_len(m_valid), function(k) {
      partial  <- oob_matrix[, seq_len(k), drop = FALSE]
      row_mean <- rowSums(partial, na.rm = TRUE) / rowSums(!is.na(partial))
      mean((row_mean - y)^2, na.rm = TRUE)
    }, numeric(1L))

    mse_oob <- mean((final_oob - y)^2, na.rm = TRUE)

    return(structure(
      c(base_out, list(
        y_train             = y,
        oob_predictions     = final_oob,
        MSE                 = mse_oob,
        RMSE                = sqrt(mse_oob),
        MAE                 = mean(abs(final_oob - y), na.rm = TRUE),
        Rsquared            = cor(final_oob, y, use = "complete.obs")^2,
        oob_matrix          = oob_matrix,
        oob_errors_per_iter = oob_errors
      )),
      class = "sword_flat"
    ))
  }

  structure(base_out, class = "sword_flat")
}


# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Predict from a fitted TORS tree
#'
#' @param object \code{tors_flat} object from \code{\link{TORS}}.
#' @param newdata data.frame of predictors.
#' @param ... ignored.
#' @return Numeric vector of predictions.
#' @examples
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' tree  <- TORS(X, y, nmin = 10, cp = 0.01)
#' preds <- predict(tree, X[1:5, ])
#' @export
predict.tors_flat <- function(object, newdata, ...) {

  # When trained with formula, apply the same preprocessing (NA removal,
  # factor levels). Also strips the response column if present in newdata.
  if (!is.null(object$terms)) {
    newdata <- stats::model.frame(
      stats::delete.response(object$terms),
      data      = newdata,
      na.action = stats::na.pass
    )
  }

  if (!is.null(object$encoding))
    newdata <- .sword_apply_encoding(newdata, object$encoding)

  # Columns not in feature_names are silently dropped; missing ones error.
  missing_cols <- setdiff(object$feature_names, colnames(newdata))
  if (length(missing_cols) > 0L)
    stop("newdata missing columns: ", paste0(missing_cols, collapse = ", "))

  int_col  <- which(colnames(object$coeffs) == "Int")
  feat_col <- seq_len(ncol(object$coeffs))[-int_col]

  data_mat <- as.matrix(newdata[, object$feature_names, drop = FALSE])
  w_mat    <- object$coeffs[, feat_col, drop = FALSE]
  int_vec  <- object$coeffs[, int_col]

  node <- rep(1L, nrow(data_mat))

  repeat {
    at_internal <- !object$is_leaf[node]
    if (!any(at_internal)) break

    by_node <- split(which(at_internal), node[at_internal])

    for (nd_char in names(by_node)) {
      nd      <- as.integer(nd_char)
      obs     <- by_node[[nd_char]]
      scores  <- drop(data_mat[obs, , drop = FALSE] %*% w_mat[nd, ]) - int_vec[nd]
      go_left <- scores < 0
      node[obs[ go_left]] <- object$left_child[nd]
      node[obs[!go_left]] <- object$right_child[nd]
    }
  }

  object$leaf_mean[node]
}

#' @rdname TORS
#' @param x \code{tors_flat} object (for \code{print}).
#' @param ... ignored.
#' @export
print.tors_flat <- function(x, ...) {
  n_leaves   <- sum(x$is_leaf)
  n_internal <- x$n_nodes - n_leaves
  cp         <- x$call_params
  obs_leaf   <- x$n_obs[x$is_leaf]

  lbl <- 25
  cat("----------", "TORS (Tree Oblique for Regression with weighted SVM)\n\n")
  if (!is.null(cp)) {
    cat(sprintf("%*s: %d\n", lbl, "N observations",      cp$n))
    cat(sprintf("%*s: %d\n", lbl, "N predictors",        cp$p))
  }
  cat(sprintf("%*s: %d (%d internal, %d leaves)\n",
              lbl, "N nodes", x$n_nodes, n_internal, n_leaves))
  cat(sprintf("%*s: %d\n",   lbl, "Max depth",           max(x$depth, na.rm = TRUE)))
  if (length(obs_leaf) > 0L) {
    cat(sprintf("%*s: %.1f\n", lbl, "Avg obs in leaf",   mean(obs_leaf)))
    cat(sprintf("%*s: %d\n",   lbl, "Min obs in leaf",   min(obs_leaf)))
  }
  if (!is.null(cp)) {
    topcor_str <- if (isTRUE(cp$rand_ntopcor))
      sprintf("random [2, %d]", cp$rf_var) else as.character(cp$n_top_cor)
    cat(sprintf("%*s: %s\n", lbl, "Top-correlated feats", topcor_str))
    cat(sprintf("%*s: %s\n", lbl, "Correlation",          cp$relation))
    cat(sprintf("%*s: %s\n", lbl, "Weight scheme",        cp$weight_scheme))
  }
  cat("\n-----------------------------------------\n")
  invisible(x)
}

#' @rdname TORS
#' @param object \code{tors_flat} object (for \code{summary}).
#' @export
summary.tors_flat <- function(object, ...) {
  cp         <- object$call_params
  n_leaves   <- sum(object$is_leaf)
  n_internal <- object$n_nodes - n_leaves
  internal   <- which(!object$is_leaf)
  obs_leaf   <- object$n_obs[object$is_leaf]

  lbl <- 25
  cat("----------", "TORS (Tree Oblique for Regression with weighted SVM) \u2014 summary\n\n")

  if (!is.null(cp)) {
    cat(sprintf("%*s: %d\n", lbl, "N observations",      cp$n))
    cat(sprintf("%*s: %d\n", lbl, "N predictors",        cp$p))
    topcor_str <- if (isTRUE(cp$rand_ntopcor))
      sprintf("random [2, %d]", cp$rf_var) else as.character(cp$n_top_cor)
    cat(sprintf("%*s: %s\n", lbl, "Top-correlated feats", topcor_str))
    cat(sprintf("%*s: %s\n", lbl, "Correlation",          cp$relation))
    cat(sprintf("%*s: %s\n", lbl, "Weight scheme",        cp$weight_scheme))
  }

  cat("\n")
  cat(sprintf("%*s: %d (%d internal, %d leaves)\n",
              lbl, "N nodes", object$n_nodes, n_internal, n_leaves))
  cat(sprintf("%*s: %d\n",   lbl, "Max depth",           max(object$depth, na.rm = TRUE)))
  if (length(obs_leaf) > 0L)
    cat(sprintf("%*s: mean=%.1f  min=%d  max=%d\n",
                lbl, "Obs in leaf", mean(obs_leaf), min(obs_leaf), max(obs_leaf)))

  if (length(internal) > 0L) {
    cat("\n")
    delta <- object$dev_presplit[internal] - object$dev[internal]
    cat(sprintf("%*s: %.4f\n", lbl, "Total deviance gain", sum(delta, na.rm = TRUE)))
    cat(sprintf("%*s: %.4f\n", lbl, "Avg gain per split",  mean(delta, na.rm = TRUE)))
  }

  cat("\n-----------------------------------------\n")
  invisible(object)
}

#' Predict from a fitted SWORD forest
#'
#' @param object   \code{sword_flat} object from \code{\link{SWORD}}.
#' @param newdata  data.frame of predictors.
#' @param aggregate Logical; if \code{TRUE} (default) return the ensemble mean.
#'   If \code{FALSE} return an \eqn{n \times m} matrix of per-tree predictions.
#' @param interval Character; \code{"none"} (default) or \code{"prediction"}.
#'   When \code{"prediction"}, the function returns a three-column matrix with
#'   columns \code{fit}, \code{lwr}, and \code{upr} built from the empirical
#'   quantiles of the per-tree predictions.  Ignored when \code{aggregate = FALSE}.
#' @param level Confidence level for the prediction interval (default \code{0.95}).
#'   Ignored unless \code{interval = "prediction"}.
#' @param ... ignored.
#' @return A numeric vector (default), an \eqn{n \times m} matrix
#'   (\code{aggregate = FALSE}), or a three-column matrix
#'   (\code{interval = "prediction"}).
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = FALSE, verbose = FALSE)
#' preds  <- predict(forest, X[1:5, ])
#' # Per-tree predictions matrix
#' mat    <- predict(forest, X[1:5, ], aggregate = FALSE)
#' # 95% prediction intervals
#' pi95   <- predict(forest, X[1:5, ], interval = "prediction")
#' }
#' @export
predict.sword_flat <- function(object, newdata,
                               aggregate = TRUE,
                               interval  = c("none", "prediction"),
                               level     = 0.95,
                               ...) {
  interval <- match.arg(interval)

  # When trained with formula, strip response and align factor levels.
  if (!is.null(object$terms)) {
    newdata <- stats::model.frame(
      stats::delete.response(object$terms),
      data      = newdata,
      na.action = stats::na.pass
    )
  }

  if (!is.null(object$encoding))
    newdata <- .sword_apply_encoding(newdata, object$encoding)

  # Columns not in feature_names are silently dropped; missing ones error.
  missing_cols <- setdiff(object$feature_names, colnames(newdata))
  if (length(missing_cols) > 0L)
    stop("newdata missing columns: ", paste0(missing_cols, collapse = ", "))

  pred_mat <- vapply(object$trees,
                     function(tree) predict(tree, newdata),
                     numeric(nrow(newdata)))
  # Force to matrix (vapply drops dimension when nrow(newdata) == 1).
  pred_mat <- matrix(pred_mat, nrow = nrow(newdata))

  if (!aggregate) return(pred_mat)

  fit <- rowMeans(pred_mat)

  if (interval == "none") return(fit)

  alpha <- (1 - level) / 2
  lwr   <- apply(pred_mat, 1L, stats::quantile, probs = alpha,       names = FALSE)
  upr   <- apply(pred_mat, 1L, stats::quantile, probs = 1 - alpha,   names = FALSE)
  cbind(fit = fit, lwr = lwr, upr = upr)
}

#' @rdname SWORD
#' @param x \code{sword_flat} object (for \code{print}).
#' @param ... ignored.
#' @export
print.sword_flat <- function(x, ...) {
  cp      <- x$call_params
  n_trees <- length(x$trees)

  avg_leaves <- if (n_trees > 0L)
    mean(vapply(x$trees, function(t) sum(t$is_leaf), integer(1L)))
  else NA_real_

  lbl <- 25
  cat("----------", "SWORD (Support vector machine Weighted Oblique Random Decision forest)\n\n")
  if (!is.null(cp)) {
    cat(sprintf("%*s: %d\n",   lbl, "N observations",       cp$n))
    cat(sprintf("%*s: %d\n",   lbl, "N predictors",         cp$p))
  }
  cat(sprintf("%*s: %d\n",   lbl, "N trees",               n_trees))
  if (!is.null(cp)) {
    cat(sprintf("%*s: %d\n",   lbl, "Predictors per split", cp$rf_var))
    topcor_str <- if (isTRUE(cp$rand_ntopcor))
      sprintf("random [2, %d]", cp$rf_var) else as.character(cp$n_top_cor)
    cat(sprintf("%*s: %s\n",   lbl, "Top-correlated feats", topcor_str))
  }
  if (!is.na(avg_leaves))
    cat(sprintf("%*s: %.3f\n", lbl, "Avg leaves per tree",  avg_leaves))
  if (!is.null(cp))
    cat(sprintf("%*s: %d\n",   lbl, "Min obs in leaf",      cp$minleaf))
  if (!is.null(x$RMSE)) {
    cat(sprintf("%*s: %s\n",   lbl, "OOB stat type",        "RMSE / R\u00b2"))
    cat(sprintf("%*s: %.4f / %.4f\n", lbl, "OOB stat value", x$RMSE, x$Rsquared))
  }
  if (!is.null(cp)) {
    cat(sprintf("%*s: %s\n",   lbl, "Weight scheme",        cp$weight_scheme))
    cat(sprintf("%*s: %s\n",   lbl, "Correlation",          cp$relation))
  }
  if (!is.null(x$n_skipped) && x$n_skipped > 0L)
    cat(sprintf("%*s: %d\n",   lbl, "Skipped trees",        x$n_skipped))
  cat("\n-----------------------------------------\n")
  invisible(x)
}

#' @rdname SWORD
#' @param object \code{sword_flat} object (for \code{summary}).
#' @export
summary.sword_flat <- function(object, ...) {
  cp      <- object$call_params
  n_trees <- length(object$trees)

  n_nodes  <- vapply(object$trees, `[[`, integer(1L), "n_nodes")
  n_leaves <- vapply(object$trees, function(t) sum(t$is_leaf), integer(1L))
  depths   <- vapply(object$trees, function(t) max(t$depth, na.rm = TRUE), integer(1L))

  lbl <- 27
  cat("----------", "SWORD (Support vector machine Weighted Oblique Random Decision forest) \u2014 summary\n\n")

  if (!is.null(cp)) {
    cat(sprintf("%*s: %d\n", lbl, "N observations",       cp$n))
    cat(sprintf("%*s: %d\n", lbl, "N predictors",         cp$p))
  }
  cat(sprintf("%*s: %d\n",   lbl, "N trees",              n_trees))
  if (!is.null(cp)) {
    cat(sprintf("%*s: %d\n", lbl, "Predictors per split", cp$rf_var))
    topcor_str <- if (isTRUE(cp$rand_ntopcor))
      sprintf("random [2, %d]", cp$rf_var) else as.character(cp$n_top_cor)
    cat(sprintf("%*s: %s\n", lbl, "Top-correlated feats", topcor_str))
    cat(sprintf("%*s: %s\n", lbl, "Correlation",          cp$relation))
    cat(sprintf("%*s: %s\n", lbl, "Weight scheme",        cp$weight_scheme))
  }

  if (n_trees > 0L) {
    cat("\n")
    cat(sprintf("%*s: mean=%.1f  min=%d  max=%d\n",
                lbl, "Nodes per tree",  mean(n_nodes),  min(n_nodes),  max(n_nodes)))
    cat(sprintf("%*s: mean=%.1f  min=%d  max=%d\n",
                lbl, "Leaves per tree", mean(n_leaves), min(n_leaves), max(n_leaves)))
    cat(sprintf("%*s: mean=%.1f  min=%d  max=%d\n",
                lbl, "Max depth",       mean(depths),   min(depths),   max(depths)))
  }

  if (!is.null(object$RMSE)) {
    cat("\n")
    cat(sprintf("%*s: %.4f\n", lbl, "OOB RMSE",  object$RMSE))
    cat(sprintf("%*s: %.4f\n", lbl, "OOB MAE",   object$MAE))
    cat(sprintf("%*s: %.4f\n", lbl, "OOB R\u00b2",    object$Rsquared))
  }

  if (!is.null(object$n_skipped) && object$n_skipped > 0L) {
    cat("\n")
    cat(sprintf("%*s: %d\n",   lbl, "Skipped trees", object$n_skipped))
  }

  cat("\n-----------------------------------------\n")
  invisible(object)
}


# ==============================================================================
# nobs / fitted / residuals
# ==============================================================================

#' Number of observations used to fit a TORS tree
#'
#' @param object \code{tors_flat} object from \code{\link{TORS}}.
#' @param ... ignored.
#' @return Integer scalar.
#' @examples
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(100 * 3), 100, 3,
#'                 dimnames = list(NULL, paste0("x", 1:3))))
#' y <- X$x1 + rnorm(100)
#' tree <- TORS(X, y)
#' nobs(tree)
#' @importFrom stats nobs
#' @export
nobs.tors_flat <- function(object, ...) object$call_params$n

#' Number of observations used to fit a SWORD forest
#'
#' @param object \code{sword_flat} object from \code{\link{SWORD}}.
#' @param ... ignored.
#' @return Integer scalar.
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(100 * 3), 100, 3,
#'                 dimnames = list(NULL, paste0("x", 1:3))))
#' y <- X$x1 + rnorm(100)
#' forest <- SWORD(X, y, m = 10, verbose = FALSE)
#' nobs(forest)
#' }
#' @importFrom stats nobs
#' @export
nobs.sword_flat <- function(object, ...) object$call_params$n

#' Out-of-bag fitted values from a SWORD forest
#'
#' Returns the out-of-bag (OOB) predictions accumulated during forest fitting.
#' Only available when the forest was fitted with \code{oob = TRUE}.
#'
#' @param object \code{sword_flat} object from \code{\link{SWORD}}.
#' @param ... ignored.
#' @return Named numeric vector of OOB predictions (one per training observation),
#'   or \code{NULL} with a message if \code{oob = FALSE} was used.
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(100 * 3), 100, 3,
#'                 dimnames = list(NULL, paste0("x", 1:3))))
#' y <- X$x1 + rnorm(100)
#' forest <- SWORD(X, y, m = 10, oob = TRUE, verbose = FALSE)
#' head(fitted(forest))
#' }
#' @importFrom stats fitted
#' @export
fitted.sword_flat <- function(object, ...) {
  if (is.null(object$oob_predictions))
    message("OOB predictions not available: refit with oob = TRUE.")
  object$oob_predictions
}

#' Out-of-bag residuals from a SWORD forest
#'
#' Returns \code{y_train - oob_predictions} for each training observation.
#' Only available when the forest was fitted with \code{oob = TRUE}.
#'
#' @param object \code{sword_flat} object from \code{\link{SWORD}}.
#' @param ... ignored.
#' @return Numeric vector of OOB residuals (one per training observation),
#'   or \code{NULL} with a message if \code{oob = FALSE} was used.
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(100 * 3), 100, 3,
#'                 dimnames = list(NULL, paste0("x", 1:3))))
#' y <- X$x1 + rnorm(100)
#' forest <- SWORD(X, y, m = 10, oob = TRUE, verbose = FALSE)
#' head(residuals(forest))
#' }
#' @importFrom stats residuals
#' @export
residuals.sword_flat <- function(object, ...) {
  if (is.null(object$oob_predictions) || is.null(object$y_train)) {
    message("OOB residuals not available: refit with oob = TRUE.")
    return(NULL)
  }
  object$y_train - object$oob_predictions
}


# ==============================================================================
# VARIABLE IMPORTANCE
# ==============================================================================

#' Oblique Impurity-Weighted Variable Importance (OIW-VI)
#'
#' Computes the Oblique Impurity-Weighted Variable Importance (OIW-VI) for a
#' fitted SWORD forest. At each internal node the deviance reduction is
#' distributed among the predictors involved in the oblique split proportionally
#' to their normalised absolute coefficients in the standardised feature space.
#' These node-level scores are then averaged across all trees.
#'
#' When the forest was trained on data with factor or character columns,
#' dummy-column contributions are automatically summed back to the original
#' variable level (e.g. \code{cutFair + cutGood + ...} -> \code{cut}).
#' Set \code{raw = TRUE} to get the raw dummy-column importances instead.
#'
#' @param forest \code{sword_flat} object from \code{\link{SWORD}}.
#' @param raw    logical; if \code{TRUE} return importances at the dummy-column
#'   level instead of the original variable level (default \code{FALSE}).
#' @return Named numeric vector of importance scores, sorted descending.
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = FALSE, verbose = FALSE)
#' VI_SWORD(forest)
#' }
#' @export
VI_SWORD <- function(forest, raw = FALSE) {
  if (length(forest$trees) == 0L)
    stop("No valid trees in the forest (all trees were skipped).")
  vi_list <- lapply(forest$trees, .sword_vi_tree)
  vi_mat  <- do.call(rbind, vi_list)
  vi_raw  <- colMeans(vi_mat, na.rm = TRUE)

  enc <- forest$encoding
  if (raw || is.null(enc) || is.null(enc$enc_assign)) {
    return(sort(vi_raw, decreasing = TRUE))
  }

  # Sum dummy-column contributions back to the original variable
  orig_vars  <- enc$enc_var_nms[enc$enc_assign]
  vi_grouped <- tapply(vi_raw, orig_vars, sum)
  sort(vi_grouped, decreasing = TRUE)
}


# ==============================================================================
# VISUALIZATION  (internal helper + two exported functions)
# ==============================================================================

# Formats a named coefficient vector as a hyperplane equation string.
# e.g.  "0.31\u00b7lstat - 0.72\u00b7rm + 1.45 = 0"
# Zero coefficients are dropped; top_k keeps only the largest-magnitude terms.
.sword_eq_string <- function(betas, digits = 3L, top_k = NULL,
                              include_intercept = TRUE) {
  if (is.null(betas) || length(betas) == 0L)
    return("No split")

  b0 <- 0
  if ("Int" %in% names(betas)) {
    b0    <- betas[["Int"]]
    betas <- betas[names(betas) != "Int"]
  }

  betas <- betas[abs(betas) > 0]
  if (!length(betas)) {
    return(if (include_intercept && b0 != 0)
      paste0(format(round(b0, digits), trim = TRUE), " = 0")
      else "No split")
  }

  if (!is.null(top_k) && top_k < length(betas)) {
    ord   <- order(abs(betas), decreasing = TRUE)
    betas <- betas[ord][seq_len(top_k)]
  }

  terms <- paste0(
    ifelse(sign(betas) >= 0, "+ ", "- "),
    format(round(abs(betas), digits), trim = TRUE),
    "\u00b7", names(betas)
  )
  lhs <- sub("^\\+\\s", "", paste(terms, collapse = " "))

  if (include_intercept) {
    inter <- if (b0 == 0) "" else
      paste0(" ", ifelse(b0 >= 0, "+ ", "- "),
             format(round(abs(b0), digits), trim = TRUE))
    paste0(lhs, inter, " = 0")
  } else {
    paste0(lhs, " = c")
  }
}


#' Indented text representation of a TORS tree
#'
#' Prints the tree node by node with an indented layout. Each internal node
#' shows the hyperplane equation; each leaf shows the predicted mean and
#' sample size.
#'
#' @param tree      \code{tors_flat} object from \code{\link{TORS}}.
#' @param use_scaled logical; use scaled-space coefficients (default \code{FALSE}).
#' @param top_k    maximum number of terms per equation (default 3).
#' @param digits   decimal places for coefficients (default 3).
#' @return Invisibly returns \code{tree}.
#' @seealso \code{\link{TORS}}, \code{\link{plot_tors_visnet}}
#' @examples
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' tree <- TORS(X, y, nmin = 10, cp = 0.01)
#' plot_tors_text(tree, top_k = 2)
#' @export
plot_tors_text <- function(tree, use_scaled = FALSE, top_k = 3L, digits = 3L) {

  cat("TORS \u2014 text representation\n")
  cat(strrep("-", 60L), "\n")

  stack <- list(list(nd = 1L, indent = 0L, branch = ""))
  while (length(stack) > 0L) {
    frame  <- stack[[length(stack)]]
    stack  <- stack[-length(stack)]

    nd     <- frame$nd
    indent <- frame$indent
    branch <- frame$branch
    prefix <- if (indent == 0L) "" else
      paste0(strrep("  ", indent - 1L), branch)

    if (tree$is_leaf[nd]) {
      cat(sprintf("%sLeaf  n=%d  mean=%.4f\n",
                  prefix, tree$n_obs[nd], tree$leaf_mean[nd]))
    } else {
      coef_row <- if (use_scaled) tree$scaled_coeffs[nd, ] else tree$coeffs[nd, ]
      eq       <- .sword_eq_string(coef_row, digits = digits, top_k = top_k)
      cat(sprintf("%sNode %d  n=%d  [ %s ]\n", prefix, nd, tree$n_obs[nd], eq))
      # Push right before left so left is popped (and printed) first
      stack <- c(stack,
        list(list(nd = tree$right_child[nd], indent = indent + 1L, branch = "R-- ")),
        list(list(nd = tree$left_child[nd],  indent = indent + 1L, branch = "L-- "))
      )
    }
  }

  invisible(tree)
}


#' Interactive zoomable visualisation of a TORS tree
#'
#' Renders an interactive HTML widget of the tree. Supports zoom (mouse wheel),
#' pan (click-drag), and hover tooltips with the full hyperplane equation.
#' Requires the \pkg{visNetwork} package.
#'
#' @param tree      \code{tors_flat} object from \code{\link{TORS}}.
#' @param use_scaled logical; use scaled-space coefficients for labels
#'                  (default \code{TRUE}).
#' @param top_k    maximum number of equation terms shown on each node
#'                  (default 4). Full equation is always in the hover tooltip.
#' @param digits   decimal places for coefficients (default 3).
#' @return A \code{visNetwork} HTML widget (rendered as a side-effect in
#'   RStudio / a browser).
#' @seealso \code{\link{TORS}}, \code{\link{plot_tors_text}}
#' @examples
#' \donttest{
#' if (requireNamespace("visNetwork", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                   dimnames = list(NULL, paste0("x", 1:5))))
#'   y <- X$x1 * 2 - X$x2 + rnorm(150)
#'   tree <- TORS(X, y, nmin = 10, cp = 0.01)
#'   plot_tors_visnet(tree)
#' }
#' }
#' @export
plot_tors_visnet <- function(tree, use_scaled = TRUE, top_k = 4L, digits = 3L) {

  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("Package 'visNetwork' is required. Install with: install.packages('visNetwork')")

  N        <- tree$n_nodes
  labels   <- character(N)
  tooltips <- character(N)

  for (nd in seq_len(N)) {
    if (tree$is_leaf[nd]) {
      labels[nd]   <- sprintf("mean=%.3f\nn=%d", tree$leaf_mean[nd], tree$n_obs[nd])
      tooltips[nd] <- sprintf("Leaf | mean=%.4f | n=%d | dev=%.3f",
                               tree$leaf_mean[nd], tree$n_obs[nd], tree$dev[nd])
    } else {
      coef_row  <- if (use_scaled) tree$scaled_coeffs[nd, ] else tree$coeffs[nd, ]
      eq_short  <- .sword_eq_string(coef_row, digits = digits, top_k = top_k)
      eq_full   <- .sword_eq_string(coef_row, digits = digits, top_k = NULL)
      d_imp     <- if (!is.na(tree$dev_presplit[nd]))
        tree$dev_presplit[nd] - tree$dev[nd] else NA_real_
      labels[nd]   <- sprintf("n=%d\n%s", tree$n_obs[nd], eq_short)
      tooltips[nd] <- sprintf("Node %d | n=%d | \u0394imp=%.3f<br>%s",
                               nd, tree$n_obs[nd],
                               ifelse(is.na(d_imp), 0, d_imp), eq_full)
    }
  }

  nodes <- data.frame(
    id    = seq_len(N),
    label = labels,
    title = tooltips,
    color = ifelse(tree$is_leaf, "#a8d8a8", "#a8c8f8")
  )

  from_v <- integer(0L)
  to_v   <- integer(0L)
  for (nd in seq_len(N)) {
    if (!tree$is_leaf[nd]) {
      from_v <- c(from_v, nd, nd)
      to_v   <- c(to_v,   tree$left_child[nd], tree$right_child[nd])
    }
  }
  edges <- data.frame(from = from_v, to = to_v)

  visNetwork::visNetwork(nodes, edges, width = "100%", height = "700px") |>
    visNetwork::visHierarchicalLayout(direction = "UD", sortMethod = "directed") |>
    visNetwork::visNodes(shape = "box",
                         font  = list(size = 11L, face = "monospace")) |>
    visNetwork::visEdges(arrows = "to", color = list(color = "gray60")) |>
    visNetwork::visOptions(highlightNearest = TRUE) |>
    visNetwork::visInteraction(navigationButtons = TRUE,
                               zoomView          = TRUE,
                               dragView          = TRUE)
}


# ==============================================================================
# ADDITIONAL PLOT FUNCTIONS
# ==============================================================================

# Base-R tree structure plot (internal fallback when visNetwork is not available)
.plot_tors_base <- function(tree) {
  N     <- tree$n_nodes
  x_pos <- numeric(N)
  y_pos <- -tree$depth

  leaf_counter <- 0L
  assign_x <- function(nd) {
    if (tree$is_leaf[nd]) {
      leaf_counter <<- leaf_counter + 1L
      x_pos[nd]    <<- leaf_counter
    } else {
      assign_x(tree$left_child[nd])
      assign_x(tree$right_child[nd])
      x_pos[nd] <<- (x_pos[tree$left_child[nd]] + x_pos[tree$right_child[nd]]) / 2
    }
  }
  assign_x(1L)

  max_depth <- max(tree$depth, na.rm = TRUE)
  n_leaves  <- sum(tree$is_leaf)
  cex_vec   <- 0.6 + tree$n_obs / max(tree$n_obs) * 1.8

  old_mar <- par(mar = c(2, 3, 3, 1))
  on.exit(par(old_mar), add = TRUE)

  plot(x_pos, y_pos,
       type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Depth",
       main = sprintf("TORS  (%d nodes, %d leaves)", N, n_leaves),
       bty  = "n",
       xlim = c(0.5, leaf_counter + 0.5),
       ylim = c(-max_depth - 0.5, 0.5))

  axis(2, at = seq(0, -max_depth), labels = seq(0, max_depth), las = 1)

  for (nd in seq_len(N)) {
    if (!tree$is_leaf[nd]) {
      lc <- tree$left_child[nd]
      rc <- tree$right_child[nd]
      segments(x_pos[nd], y_pos[nd], x_pos[lc], y_pos[lc], col = "gray70", lwd = 0.8)
      segments(x_pos[nd], y_pos[nd], x_pos[rc], y_pos[rc], col = "gray70", lwd = 0.8)
    }
  }

  points(x_pos, y_pos,
         pch = ifelse(tree$is_leaf, 16L, 15L),
         col = ifelse(tree$is_leaf, "#2ecc71", "#3498db"),
         cex = cex_vec)

  legend("topright",
         legend = c("Internal node", "Leaf"),
         pch    = c(15L, 16L),
         col    = c("#3498db", "#2ecc71"),
         pt.cex = 1.2, bty = "n")

  invisible(tree)
}

# Internal: OOB MSE convergence line plot
.plot_oob_convergence <- function(forest) {
  errs <- forest$oob_errors_per_iter
  if (is.null(errs) || length(errs) == 0L)
    stop("No OOB errors found. Fit SWORD with oob = TRUE.")

  m      <- length(errs)
  i_best <- which.min(errs)

  plot(seq_len(m), errs,
       type = "l", lwd = 2L, col = "#3498db",
       xlab = "Number of trees", ylab = "OOB MSE",
       main = "SWORD \u2014 OOB convergence",
       bty  = "l")
  abline(h = min(errs), lty = 2L, col = "gray60")
  points(i_best, errs[i_best], pch = 19L, col = "#e74c3c", cex = 1.4)
  legend("topright",
         legend = sprintf("min MSE = %.4f  (tree %d)", errs[i_best], i_best),
         bty = "n", text.col = "#e74c3c")
  invisible(forest)
}


#' Plot a fitted TORS tree
#'
#' Dispatches to \code{\link{plot_tors_visnet}} when \pkg{visNetwork} is
#' installed; otherwise draws a base-R structure plot (node size proportional
#' to observation count).
#'
#' @param x   \code{tors_flat} object from \code{\link{TORS}}.
#' @param use_scaled logical; use scaled-space coefficients (default \code{TRUE}).
#' @param top_k    max equation terms shown per node (default 4).
#' @param digits   decimal places for coefficients (default 3).
#' @param ...  ignored.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{plot_tors_visnet}}, \code{\link{plot_tors_text}}
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' tree <- TORS(X, y, nmin = 10, cp = 0.01)
#' plot(tree)
#' }
#' @export
plot.tors_flat <- function(x, use_scaled = TRUE, top_k = 4L, digits = 3L, ...) {
  if (requireNamespace("visNetwork", quietly = TRUE)) {
    plot_tors_visnet(x, use_scaled = use_scaled, top_k = top_k, digits = digits)
  } else {
    .plot_tors_base(x)
  }
}


#' Plot method for a SWORD forest
#'
#' @param x     \code{sword_flat} object from \code{\link{SWORD}}.
#' @param which \code{"oob"} for OOB MSE convergence curve (requires
#'   \code{oob = TRUE} at fit time); \code{"vi"} for variable importance
#'   barplot.
#' @param top_n maximum variables shown when \code{which = "vi"} (default 20).
#' @param ...   further arguments passed to \code{\link{plot_vi_sword}} when
#'   \code{which = "vi"}.
#' @return Invisibly returns \code{x}.
#' @seealso \code{\link{plot_vi_sword}}, \code{\link{plot_oob_fit}},
#'   \code{\link{pdp_sword}}
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = TRUE, verbose = FALSE)
#' plot(forest, which = "oob")
#' plot(forest, which = "vi")
#' }
#' @export
plot.sword_flat <- function(x, which = c("oob", "vi"), top_n = 20L, ...) {
  which <- match.arg(which)
  if (which == "oob") .plot_oob_convergence(x)
  else                plot_vi_sword(x, top_n = top_n, ...)
  invisible(x)
}


#' Variable importance barplot for a SWORD forest
#'
#' Calls \code{\link{VI_SWORD}} and draws a horizontal barplot of
#' deviance-reduction importance scores, sorted descending.
#'
#' @param forest \code{sword_flat} object from \code{\link{SWORD}}.
#' @param top_n  maximum number of variables to display (default 20).
#' @param col    bar fill colour (default \code{"#3498db"}).
#' @param raw    logical; if \code{TRUE} show importances at the dummy-column
#'   level instead of the original variable level (default \code{FALSE}).
#' @param main   plot title; \code{NULL} generates a default title.
#' @param ...    further arguments passed to \code{barplot}.
#' @return Invisibly returns the named importance vector (full length, before
#'   truncation to \code{top_n}).
#' @seealso \code{\link{VI_SWORD}}, \code{\link{plot.sword_flat}}
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = FALSE, verbose = FALSE)
#' plot_vi_sword(forest)
#' }
#' @export
plot_vi_sword <- function(forest, top_n = 20L, col = "#3498db",
                          raw = FALSE, main = NULL, ...) {
  vi <- VI_SWORD(forest, raw = raw)
  vi <- vi[vi > 0]
  if (length(vi) == 0L) {
    message("All variable importances are zero.")
    return(invisible(vi))
  }
  vi_show <- if (length(vi) > top_n) vi[seq_len(top_n)] else vi
  vi_plot  <- rev(vi_show)

  if (is.null(main))
    main <- sprintf("SWORD \u2014 Variable Importance  (top %d)", length(vi_plot))

  max_name <- max(nchar(names(vi_plot)))
  old_mar  <- par(mar = c(4, max_name * 0.55 + 1, 3, 1))
  on.exit(par(old_mar), add = TRUE)

  barplot(vi_plot,
          horiz  = TRUE,
          las    = 1L,
          col    = col,
          border = NA,
          xlab   = "Mean deviance reduction",
          main   = main,
          ...)
  invisible(vi)
}


#' OOB predicted vs observed scatter plot
#'
#' Plots out-of-bag predictions against the true response values with an
#' identity line and RMSE / R\eqn{^2} annotation. Requires that the forest was
#' fit with \code{oob = TRUE}.
#'
#' @param forest \code{sword_flat} object fit with \code{oob = TRUE}.
#' @param y      numeric vector of true response values in the same order as
#'               the training data passed to \code{\link{SWORD}}.
#' @param main   plot title (default \code{"SWORD - OOB fit"}).
#' @param ...    further graphical arguments passed to \code{plot}.
#' @return Invisibly returns a list with elements \code{rmse}, \code{mae}, and
#'   \code{r2}.
#' @seealso \code{\link{SWORD}}, \code{\link{plot.sword_flat}}
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = TRUE, verbose = FALSE)
#' plot_oob_fit(forest, y)
#' }
#' @export
plot_oob_fit <- function(forest, y, main = "SWORD \u2014 OOB fit", ...) {
  yhat <- forest$oob_predictions
  if (is.null(yhat))
    stop("No OOB predictions found. Fit SWORD with oob = TRUE.")

  keep <- !is.na(yhat)
  y_k  <- y[keep]
  yh_k <- yhat[keep]
  rmse <- sqrt(mean((yh_k - y_k)^2))
  mae  <- mean(abs(yh_k - y_k))
  r2   <- cor(yh_k, y_k)^2

  lims <- range(c(y_k, yh_k))
  plot(y_k, yh_k,
       xlab = "Observed", ylab = "OOB predicted",
       main = main,
       pch  = 16L, col = adjustcolor("#3498db", 0.5),
       xlim = lims, ylim = lims,
       bty  = "l", ...)
  abline(0, 1, col = "gray50", lty = 2L)
  legend("topleft",
         legend = sprintf("RMSE=%.3f   MAE=%.3f   R\u00b2=%.3f", rmse, mae, r2),
         bty = "n")

  invisible(list(rmse = rmse, mae = mae, r2 = r2))
}


#' Partial dependence plot for a SWORD forest
#'
#' Computes the marginal effect of a single predictor on forest predictions by
#' fixing the variable to a grid of values and averaging predictions over the
#' full data matrix (Friedman 2001). A rug of observed values is added by
#' default.
#'
#' @param forest \code{sword_flat} object from \code{\link{SWORD}}.
#' @param X      data.frame of predictor values (same columns used to fit the
#'               forest; typically the training set).
#' @param var    character string; name of the variable to profile.
#' @param n_grid number of equally-spaced grid points over the variable's
#'               observed range (default 50).
#' @param rug    logical; add a rug of observed values (default \code{TRUE}).
#' @param main   plot title; \code{NULL} generates a default title.
#' @param xlab   x-axis label (default: the variable name).
#' @param ylab   y-axis label (default \code{"Partial dependence"}).
#' @param ...    further graphical arguments passed to \code{plot}.
#' @return Invisibly returns a data.frame with columns \code{x} (grid values,
#'   numeric for continuous variables or character for discrete/factor variables)
#'   and \code{yhat} (partial dependence values).
#' @seealso \code{\link{SWORD}}, \code{\link{plot.sword_flat}}
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(150 * 5), 150, 5,
#'                 dimnames = list(NULL, paste0("x", 1:5))))
#' y <- X$x1 * 2 - X$x2 + rnorm(150)
#' forest <- SWORD(X, y, m = 10, oob = FALSE, verbose = FALSE)
#' pdp_sword(forest, X, "x1")
#' }
#' @export
pdp_sword <- function(forest, X, var, n_grid = 50L, rug = TRUE,
                      main = NULL, xlab = var,
                      ylab = "Partial dependence", ...) {
  if (!var %in% names(X))
    stop("Variable '", var, "' not found in X.")

  col_vals    <- X[[var]]
  uniq_vals   <- if (is.factor(col_vals)) levels(col_vals)
                 else sort(unique(na.omit(col_vals)))
  is_discrete <- is.factor(col_vals) || is.character(col_vals) ||
                 length(uniq_vals) <= 10L

  grid_vals <- if (is_discrete) uniq_vals
               else seq(min(col_vals, na.rm = TRUE),
                        max(col_vals, na.rm = TRUE),
                        length.out = n_grid)

  if (is.null(main))
    main <- sprintf("SWORD \u2014 Partial dependence: %s", var)

  pd <- vapply(grid_vals, function(v) {
    X_tmp <- X
    if (is.factor(col_vals)) {
      X_tmp[[var]] <- factor(rep(as.character(v), nrow(X)),
                             levels = levels(col_vals))
    } else {
      X_tmp[[var]] <- v
    }
    mean(predict(forest, X_tmp))
  }, numeric(1L))

  if (is_discrete) {
    old_mar <- par(mar = c(5, 4, 3, 1))
    on.exit(par(old_mar), add = TRUE)
    barplot(pd,
            names.arg = as.character(grid_vals),
            col    = "#3498db",
            border = NA,
            xlab   = xlab,
            ylab   = ylab,
            main   = main,
            ...)
  } else {
    plot(grid_vals, pd,
         type = "l", lwd = 2L, col = "#3498db",
         xlab = xlab, ylab = ylab, main = main,
         bty  = "l", ...)
    if (rug) rug(col_vals, col = adjustcolor("gray40", 0.4))
  }

  invisible(data.frame(x = grid_vals, yhat = pd))
}
