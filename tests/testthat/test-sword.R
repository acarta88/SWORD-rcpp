set.seed(1)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)

test_that("SWORD fits without error", {
  forest <- SWORD(X, y, m = 5, oob = FALSE, verbose = FALSE)
  expect_s3_class(forest, "sword_flat")
})

test_that("SWORD returns the requested number of trees", {
  forest <- SWORD(X, y, m = 4, oob = FALSE, verbose = FALSE)
  expect_length(forest$trees, 4L)
})

test_that("predict.sword_flat returns finite numeric vector", {
  forest <- SWORD(X, y, m = 5, oob = FALSE, verbose = FALSE)
  preds  <- predict(forest, X)
  expect_type(preds, "double")
  expect_length(preds, nrow(X))
  expect_true(all(is.finite(preds)))
})

test_that("SWORD oob = TRUE computes predictions and RMSE", {
  forest <- SWORD(X, y, m = 5, oob = TRUE, verbose = FALSE)
  expect_false(is.null(forest$oob_predictions))
  expect_length(forest$oob_predictions, nrow(X))
  expect_true(is.finite(forest$RMSE))
  expect_true(is.finite(forest$Rsquared))
})

test_that("SWORD formula interface works", {
  df     <- cbind(X, y = y)
  forest <- SWORD(y ~ ., data = df, m = 4, oob = FALSE, verbose = FALSE)
  expect_s3_class(forest, "sword_flat")
  preds <- predict(forest, df)
  expect_length(preds, nrow(df))
})

test_that("predict.sword_flat works with a single-row newdata (regression: rowMeans on vector)", {
  forest <- SWORD(X, y, m = 5, oob = FALSE, verbose = FALSE)
  p1 <- predict(forest, X[1L, , drop = FALSE])
  pn <- predict(forest, X)
  expect_length(p1, 1L)
  expect_true(is.finite(p1))
  expect_equal(p1, pn[[1L]])   # consistent with the full predict
})

test_that("SWORD with nu-classification runs without error", {
  forest <- SWORD(X, y, m = 3, oob = FALSE, verbose = FALSE,
                  type_of_svm = "nu-classification", cost_nu = 0.5)
  expect_s3_class(forest, "sword_flat")
  expect_true(all(is.finite(predict(forest, X))))
})

# ---------------------------------------------------------------------------
# H2: predict.sword_flat — aggregate = FALSE / interval = "prediction"
# ---------------------------------------------------------------------------

test_that("predict with aggregate = FALSE returns n x m matrix", {
  forest <- SWORD(X, y, m = 5, oob = FALSE, verbose = FALSE)
  mat    <- predict(forest, X, aggregate = FALSE)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(nrow(X), 5L))
  expect_true(all(is.finite(mat)))
})

test_that("predict with interval = 'prediction' returns 3-column matrix", {
  forest <- SWORD(X, y, m = 5, oob = FALSE, verbose = FALSE)
  pi95   <- predict(forest, X, interval = "prediction")
  expect_true(is.matrix(pi95))
  expect_equal(dim(pi95), c(nrow(X), 3L))
  expect_equal(colnames(pi95), c("fit", "lwr", "upr"))
  expect_true(all(is.finite(pi95)))
  expect_true(all(pi95[, "lwr"] <= pi95[, "fit"]))
  expect_true(all(pi95[, "fit"] <= pi95[, "upr"]))
})

test_that("predict interval level changes interval width", {
  forest  <- SWORD(X, y, m = 10, oob = FALSE, verbose = FALSE)
  pi50    <- predict(forest, X, interval = "prediction", level = 0.50)
  pi95    <- predict(forest, X, interval = "prediction", level = 0.95)
  # 95% intervals should be at least as wide as 50% intervals
  width50 <- pi50[, "upr"] - pi50[, "lwr"]
  width95 <- pi95[, "upr"] - pi95[, "lwr"]
  expect_true(all(width95 >= width50 - 1e-10))
})

test_that("predict aggregate = FALSE single row returns 1 x m matrix", {
  forest <- SWORD(X, y, m = 5, oob = FALSE, verbose = FALSE)
  mat1   <- predict(forest, X[1L, , drop = FALSE], aggregate = FALSE)
  expect_equal(dim(mat1), c(1L, 5L))
})
