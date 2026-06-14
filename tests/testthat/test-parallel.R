skip_if_not_installed("future")
skip_if_not_installed("furrr")
skip_if_not_installed("progressr")

# Small dataset — keep tests fast
set.seed(1)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)

# Helper: run expr under a 2-worker multisession plan, restore sequential after
with_par <- function(expr) {
  future::plan("multisession", workers = 2L)
  on.exit(future::plan("sequential"), add = TRUE)
  force(expr)
}

# Fit reference forests once — reused across multiple tests
forest_seq <- SWORD(X, y, m = 5L, oob = TRUE,  verbose = FALSE, parallel = FALSE)
forest_par <- with_par(SWORD(X, y, m = 5L, oob = TRUE,  verbose = FALSE, parallel = TRUE))


# ==============================================================================
# A  Sequential == Parallel  (same seeds → identical trees → identical results)
# ==============================================================================

test_that("A1: parallel predict matches sequential predict", {
  expect_equal(predict(forest_par, X), predict(forest_seq, X))
})

test_that("A2: parallel oob predictions match sequential", {
  expect_equal(forest_par$oob_predictions, forest_seq$oob_predictions)
})

test_that("A3: parallel oob metrics match sequential", {
  expect_equal(forest_par$RMSE,     forest_seq$RMSE)
  expect_equal(forest_par$MSE,      forest_seq$MSE)
  expect_equal(forest_par$MAE,      forest_seq$MAE)
  expect_equal(forest_par$Rsquared, forest_seq$Rsquared)
})

test_that("A4: parallel oob_errors_per_iter matches sequential", {
  expect_equal(forest_par$oob_errors_per_iter, forest_seq$oob_errors_per_iter)
})

test_that("A5: parallel is reproducible — two identical runs give identical results", {
  forest_par2 <- with_par(
    SWORD(X, y, m = 5L, oob = TRUE, verbose = FALSE, parallel = TRUE)
  )
  expect_equal(predict(forest_par, X), predict(forest_par2, X))
  expect_equal(forest_par$oob_predictions, forest_par2$oob_predictions)
})


# ==============================================================================
# B  Regression tests — one test per fixed bug
# ==============================================================================

# B1 — terms bug: formula-fitted forest fitted in parallel must be predictable
test_that("B1: predict works after parallel formula fit (regression: missing terms field)", {
  df         <- cbind(X, y = y)
  forest_frm <- with_par(
    SWORD(y ~ ., data = df, m = 4L, oob = FALSE, verbose = FALSE, parallel = TRUE)
  )
  expect_s3_class(forest_frm, "sword_flat")
  preds <- predict(forest_frm, df)
  expect_length(preds, nrow(df))
  expect_true(all(is.finite(preds)))
})

# B2 — oob matrix dimensions
test_that("B2: oob_matrix has n rows and one column per valid tree in parallel mode", {
  expect_equal(nrow(forest_par$oob_matrix), nrow(X))
  # No trees were skipped → ncol must equal the number of trees requested
  expect_equal(ncol(forest_par$oob_matrix), length(forest_par$trees))
})

# B3 — oob predictions coverage
# With m bootstrap samples, P(obs never oob) ≈ (1-1/e)^m ≈ 0.63^m.
# For m = 5 that is ~10 %, so some NAs are expected and correct.
# We only require that the majority of observations are covered.
test_that("B3: oob_predictions length is correct and most observations are covered (parallel)", {
  expect_length(forest_par$oob_predictions, nrow(X))
  n_covered <- sum(!is.na(forest_par$oob_predictions))
  expect_gt(n_covered, 0L)
  expect_gt(n_covered / nrow(X), 0.7)   # ≥70 % covered with m=5 is expected
  # Covered predictions must be finite
  expect_true(all(is.finite(
    forest_par$oob_predictions[!is.na(forest_par$oob_predictions)]
  )))
})

# B4 — convergence curve length and values
test_that("B4: oob_errors_per_iter length == ncol(oob_matrix) in parallel mode", {
  expect_length(
    forest_par$oob_errors_per_iter,
    ncol(forest_par$oob_matrix)
  )
})

test_that("B4b: oob_errors_per_iter values are finite and non-negative in parallel mode", {
  errs <- forest_par$oob_errors_per_iter
  expect_true(all(is.finite(errs)))
  expect_true(all(errs >= 0))
})

test_that("B4c: oob_errors_per_iter[m_valid] equals MSE from final_oob (parallel)", {
  last_err <- tail(forest_par$oob_errors_per_iter, 1L)
  expect_equal(last_err, forest_par$MSE, tolerance = 1e-10)
})

# B5 — chunk path: same predictions and oob as non-chunked parallel
test_that("B5: chunk = TRUE parallel gives identical results to non-chunked parallel", {
  forest_chunk <- with_par(
    SWORD(X, y, m = 5L, oob = TRUE, verbose = FALSE,
          parallel = TRUE, chunk = TRUE, n_chunks = 2L)
  )
  expect_equal(predict(forest_chunk, X), predict(forest_par, X))
  expect_equal(forest_chunk$oob_predictions, forest_par$oob_predictions)
  expect_equal(forest_chunk$oob_errors_per_iter, forest_par$oob_errors_per_iter)
})

# B6 — oob = FALSE in parallel (no oob fields, no crash)
test_that("B6: SWORD parallel with oob = FALSE returns valid forest without oob fields", {
  forest_no_oob <- with_par(
    SWORD(X, y, m = 4L, oob = FALSE, verbose = FALSE, parallel = TRUE)
  )
  expect_s3_class(forest_no_oob, "sword_flat")
  expect_length(forest_no_oob$trees, 4L)
  expect_null(forest_no_oob$oob_predictions)
  expect_null(forest_no_oob$RMSE)
  preds <- predict(forest_no_oob, X)
  expect_length(preds, nrow(X))
  expect_true(all(is.finite(preds)))
})
