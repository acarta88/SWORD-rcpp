set.seed(1)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)
forest <- SWORD(X, y, m = 5, OOB = FALSE, verbose = FALSE)

test_that("VI_SWORD returns a named numeric vector", {
  vi <- VI_SWORD(forest)
  expect_type(vi, "double")
  expect_named(vi)
  expect_length(vi, ncol(X))
})

test_that("VI_SWORD scores are non-negative", {
  vi <- VI_SWORD(forest)
  expect_true(all(vi >= 0))
})

test_that("VI_SWORD scores sum to a positive value", {
  vi <- VI_SWORD(forest)
  expect_gt(sum(vi), 0)
})

test_that("VI_SWORD is sorted descending", {
  vi <- VI_SWORD(forest)
  expect_equal(vi, sort(vi, decreasing = TRUE))
})
