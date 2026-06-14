set.seed(42)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)

tree   <- TORS(X, y, nmin = 10, cp = 0.01)
forest <- SWORD(X, y, m = 5, oob = TRUE, verbose = FALSE)

# --- print.tors_flat ----------------------------------------------------------

test_that("print.tors_flat runs without error", {
  expect_no_error(print(tree))
})

test_that("print.tors_flat returns x invisibly", {
  expect_identical(print(tree), tree)
})

# --- summary.tors_flat --------------------------------------------------------

test_that("summary.tors_flat runs without error", {
  expect_no_error(summary(tree))
})

# --- print.sword_flat ---------------------------------------------------------

test_that("print.sword_flat runs without error", {
  expect_no_error(print(forest))
})

test_that("print.sword_flat returns x invisibly", {
  expect_identical(print(forest), forest)
})

# --- summary.sword_flat -------------------------------------------------------

test_that("summary.sword_flat runs without error", {
  expect_no_error(summary(forest))
})

# --- print/summary on stump (single leaf) -------------------------------------

test_that("print.tors_flat works on stump (cp = 1)", {
  stump <- TORS(X, y, nmin = 10, cp = 1.0)
  expect_no_error(print(stump))
  expect_no_error(summary(stump))
})

# --- predict error paths ------------------------------------------------------

test_that("predict.tors_flat errors on missing columns", {
  expect_error(predict(tree, X[, 1:2, drop = FALSE]), "missing columns")
})

test_that("predict.sword_flat errors on missing columns", {
  expect_error(predict(forest, X[, 1:2, drop = FALSE]), "missing columns")
})

# --- VI_SWORD on forest with factor encoding (grouped VI) ---------------------

test_that("VI_SWORD returns grouped VI when forest has factor predictors", {
  set.seed(7)
  Xf <- data.frame(
    x1  = rnorm(80),
    x2  = rnorm(80),
    grp = factor(sample(c("a", "b", "c"), 80, replace = TRUE))
  )
  yf <- Xf$x1 - Xf$x2 + rnorm(80)
  ff <- SWORD(Xf, yf, m = 5, oob = FALSE, verbose = FALSE)
  vi <- VI_SWORD(ff)
  expect_true("grp" %in% names(vi))
})

test_that("VI_SWORD raw = TRUE skips grouping", {
  set.seed(7)
  Xf <- data.frame(
    x1  = rnorm(80),
    grp = factor(sample(c("a", "b"), 80, replace = TRUE))
  )
  yf <- Xf$x1 + rnorm(80)
  ff <- SWORD(Xf, yf, m = 5, oob = FALSE, verbose = FALSE)
  vi_raw <- VI_SWORD(ff, raw = TRUE)
  expect_false("grp" %in% names(vi_raw))
})
