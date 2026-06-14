## Edge-case and correctness tests for SWORD / TORS

set.seed(42)
X_base <- data.frame(matrix(rnorm(120 * 5), 120, 5,
                             dimnames = list(NULL, paste0("x", 1:5))))
y_base <- X_base$x1 * 2 - X_base$x2 + rnorm(120)

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("TORS errors on non-numeric y", {
  expect_error(TORS(X_base, as.character(y_base)), "must be a numeric")
})

test_that("TORS errors on length mismatch", {
  expect_error(TORS(X_base, y_base[-1L]), "length\\(y\\)")
})

test_that("TORS errors on constant y", {
  expect_error(TORS(X_base, rep(1, nrow(X_base))), "zero variance")
})

test_that("TORS errors on NA in y", {
  y_na <- y_base; y_na[1L] <- NA
  expect_error(TORS(X_base, y_na), "non-finite")
})

test_that("TORS errors on NA in X", {
  X_na <- X_base; X_na[1L, 1L] <- NA
  expect_error(TORS(X_na, y_base), "non-finite")
})

test_that("SWORD errors on mismatched dimensions", {
  expect_error(SWORD(X_base, y_base[-1L], m = 5), "length\\(y\\)")
})

# ---------------------------------------------------------------------------
# Single predictor (p = 1)
# ---------------------------------------------------------------------------

test_that("TORS works with p = 1", {
  X1 <- data.frame(x1 = X_base$x1)
  y1 <- X_base$x1 * 3 + rnorm(120)
  tree <- TORS(X1, y1, nmin = 10, cp = 0.05)
  expect_s3_class(tree, "tors_flat")
  preds <- predict(tree, X1)
  expect_length(preds, 120L)
})

test_that("SWORD works with p = 1", {
  X1 <- data.frame(x1 = X_base$x1)
  y1 <- X1$x1 * 3 + rnorm(120)
  forest <- SWORD(X1, y1, m = 5, oob = FALSE, verbose = FALSE)
  expect_s3_class(forest, "sword_flat")
  expect_length(predict(forest, X1), 120L)
})

# ---------------------------------------------------------------------------
# Unseen factor levels at predict() time
# ---------------------------------------------------------------------------

test_that("predict handles unseen factor levels gracefully", {
  X_fac <- data.frame(
    num = rnorm(100),
    cat = factor(sample(c("A", "B"), 100, replace = TRUE))
  )
  y_fac <- X_fac$num + as.integer(X_fac$cat == "A") + rnorm(100)
  forest <- SWORD(X_fac, y_fac, m = 5, oob = FALSE, verbose = FALSE)

  # New data with an unseen level "C" — should not error, maps to zero dummy
  X_new <- data.frame(num = rnorm(10), cat = factor(c("A", "B", "C", "C",
                                                       "A", "B", "C", "A",
                                                       "B", "C")))
  expect_no_error(preds <- predict(forest, X_new))
  expect_length(preds, 10L)
  expect_true(all(is.finite(preds)))
})

# ---------------------------------------------------------------------------
# Extra columns in newdata are silently dropped
# ---------------------------------------------------------------------------

test_that("predict.sword_flat drops extra columns silently", {
  forest <- SWORD(X_base, y_base, m = 5, oob = FALSE, verbose = FALSE)
  X_extra <- cbind(X_base, response = y_base, junk = 99)
  expect_no_error(preds <- predict(forest, X_extra))
  expect_length(preds, nrow(X_base))
})

# ---------------------------------------------------------------------------
# nobs / fitted / residuals
# ---------------------------------------------------------------------------

test_that("nobs returns training sample size", {
  tree   <- TORS(X_base, y_base, nmin = 10, cp = 0.05)
  forest <- SWORD(X_base, y_base, m = 5, oob = TRUE, verbose = FALSE)
  expect_equal(nobs(tree),   nrow(X_base))
  expect_equal(nobs(forest), nrow(X_base))
})

test_that("fitted returns oob predictions when oob = TRUE", {
  # m large enough that every observation is out-of-bag for at least one tree
  # (otherwise never-OOB observations get a legitimate NA, as in other RF pkgs).
  forest <- SWORD(X_base, y_base, m = 60, oob = TRUE, verbose = FALSE)
  f <- fitted(forest)
  expect_length(f, nrow(X_base))
  expect_true(all(is.finite(f)))
})

test_that("fitted returns NULL with message when oob = FALSE", {
  forest <- SWORD(X_base, y_base, m = 5, oob = FALSE, verbose = FALSE)
  expect_message(f <- fitted(forest), "oob")
  expect_null(f)
})

test_that("residuals returns y_train - oob_predictions", {
  forest <- SWORD(X_base, y_base, m = 10, oob = TRUE, verbose = FALSE)
  r <- residuals(forest)
  expect_length(r, nrow(X_base))
  expect_equal(r, forest$y_train - forest$oob_predictions)
})

# ---------------------------------------------------------------------------
# VI correctness: signal variable ranks above noise
# ---------------------------------------------------------------------------

test_that("VI_SWORD ranks signal variable above noise", {
  set.seed(7)
  n <- 200
  signal <- rnorm(n)
  noise  <- matrix(rnorm(n * 9), n, 9)
  X_vi   <- as.data.frame(cbind(signal = signal, noise))
  y_vi   <- signal * 3 + rnorm(n, sd = 0.5)

  forest_vi <- SWORD(X_vi, y_vi, m = 50, oob = FALSE, verbose = FALSE)
  vi        <- VI_SWORD(forest_vi)
  expect_equal(names(vi)[1L], "signal",
               info = "Signal variable should rank first in VI_SWORD()")
})

# ---------------------------------------------------------------------------
# m = 1 edge case
# ---------------------------------------------------------------------------

test_that("SWORD works with m = 1", {
  forest <- SWORD(X_base, y_base, m = 1, oob = TRUE, verbose = FALSE)
  expect_s3_class(forest, "sword_flat")
  expect_length(forest$trees, 1L)
})

# ---------------------------------------------------------------------------
# H3: axis-aligned fallback (via internal helper)
# ---------------------------------------------------------------------------

test_that(".axis_aligned_fallback produces a valid pseudo-SVM object", {
  set.seed(1)
  n <- 20; p <- 3
  x_mat    <- matrix(rnorm(n * p), n, p,
                     dimnames = list(NULL, paste0("f", 1:p)))
  y        <- x_mat[, 1] * 2 + rnorm(n)
  mu       <- colMeans(x_mat)
  sigma    <- apply(x_mat, 2, sd)
  sigma[sigma == 0] <- 1

  fb <- SWORD:::.axis_aligned_fallback(x_mat, y, mu, sigma, colnames(x_mat))

  # Must be a list with the expected pseudo-SVM fields
  expect_type(fb, "list")
  expect_true(all(c("coefs", "SV", "rho", "scaled") %in% names(fb)))

  # SV should be a 1-row matrix with one non-zero entry
  expect_equal(nrow(fb$SV), 1L)
  expect_equal(ncol(fb$SV), p)
  expect_equal(sum(fb$SV != 0), 1L)

  # Most correlated feature should be f1; split at its median => left iff x < median
  med_f1   <- median(x_mat[, "f1"])
  coef_u   <- SWORD:::.sword_coef_unscaled(fb, as.data.frame(x_mat))
  go_left  <- SWORD:::.sword_pred_coef(as.data.frame(x_mat), coef_u)
  expected <- x_mat[, "f1"] < med_f1
  expect_equal(go_left, expected)
})
