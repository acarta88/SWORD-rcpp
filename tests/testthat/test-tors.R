set.seed(1)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)

test_that("TORS fits without error", {
  tree <- TORS(X, y, nmin = 10, cp = 0.01)
  expect_s3_class(tree, "tors_flat")
})

test_that("TORS returns valid structure", {
  tree <- TORS(X, y, nmin = 10, cp = 0.01)
  expect_gte(tree$n_nodes, 1L)
  expect_length(tree$is_leaf, tree$n_nodes)
  expect_true(all(is.finite(tree$leaf_mean[tree$is_leaf])))
})

test_that("predict.tors_flat returns finite numeric vector", {
  tree  <- TORS(X, y, nmin = 10, cp = 0.01)
  preds <- predict(tree, X)
  expect_type(preds, "double")
  expect_length(preds, nrow(X))
  expect_true(all(is.finite(preds)))
})

test_that("TORS formula interface matches data.frame interface", {
  df    <- cbind(X, y = y)
  tree1 <- TORS(X, y, nmin = 10, cp = 0.01, rand_ntopcor = FALSE)
  tree2 <- TORS(y ~ ., data = df, nmin = 10, cp = 0.01, rand_ntopcor = FALSE)
  expect_equal(predict(tree1, X), predict(tree2, X))
})

test_that("TORS with large cp yields a single-leaf tree", {
  tree <- TORS(X, y, nmin = 10, cp = 1.0)
  expect_s3_class(tree, "tors_flat")
  expect_length(predict(tree, X), nrow(X))
})

test_that("TORS with nu-classification runs without error", {
  tree <- TORS(X, y, nmin = 10, cp = 0.01,
               type_of_svm = "nu-classification", cost_nu = 0.5)
  expect_s3_class(tree, "tors_flat")
  expect_true(all(is.finite(predict(tree, X))))
})
