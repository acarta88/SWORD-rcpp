set.seed(1)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)

test_that("SWORD fits without error", {
  forest <- SWORD(X, y, m = 5, OOB = FALSE, verbose = FALSE)
  expect_s3_class(forest, "sword_flat")
})

test_that("SWORD returns the requested number of trees", {
  forest <- SWORD(X, y, m = 4, OOB = FALSE, verbose = FALSE)
  expect_length(forest$trees, 4L)
})

test_that("predict.sword_flat returns finite numeric vector", {
  forest <- SWORD(X, y, m = 5, OOB = FALSE, verbose = FALSE)
  preds  <- predict(forest, X)
  expect_type(preds, "double")
  expect_length(preds, nrow(X))
  expect_true(all(is.finite(preds)))
})

test_that("SWORD OOB = TRUE computes predictions and RMSE", {
  forest <- SWORD(X, y, m = 5, OOB = TRUE, verbose = FALSE)
  expect_false(is.null(forest$oob_predictions))
  expect_length(forest$oob_predictions, nrow(X))
  expect_true(is.finite(forest$RMSE))
  expect_true(is.finite(forest$Rsquared))
})

test_that("SWORD formula interface works", {
  df     <- cbind(X, y = y)
  forest <- SWORD(y ~ ., data = df, m = 4, OOB = FALSE, verbose = FALSE)
  expect_s3_class(forest, "sword_flat")
  preds <- predict(forest, df)
  expect_length(preds, nrow(df))
})

test_that("SWORD with nu-classification runs without error", {
  forest <- SWORD(X, y, m = 3, OOB = FALSE, verbose = FALSE,
                  type_of_svm = "nu-classification", cost_nu = 0.5)
  expect_s3_class(forest, "sword_flat")
  expect_true(all(is.finite(predict(forest, X))))
})
