set.seed(42)
X <- data.frame(matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4))))
y <- X$x1 * 2 - X$x2 + rnorm(80)

tree   <- TORS(X, y, nmin = 10, cp = 0.01)
forest <- SWORD(X, y, m = 5, OOB = TRUE, verbose = FALSE)

# --- plot_tors_text -----------------------------------------------------------

test_that("plot_tors_text runs without error", {
  expect_no_error(plot_tors_text(tree))
})

test_that("plot_tors_text with use_scaled = TRUE runs without error", {
  expect_no_error(plot_tors_text(tree, use_scaled = TRUE))
})

test_that("plot_tors_text on stump runs without error", {
  stump <- TORS(X, y, nmin = 10, cp = 1.0)
  expect_no_error(plot_tors_text(stump))
})

# --- plot.tors_flat (base fallback) -------------------------------------------

test_that("plot.tors_flat runs without error", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot(tree))
})

# --- plot.sword_flat ----------------------------------------------------------

test_that("plot.sword_flat which = 'oob' runs without error", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot(forest, which = "oob"))
})

test_that("plot.sword_flat which = 'vi' runs without error", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot(forest, which = "vi"))
})

# --- plot_vi_sword ------------------------------------------------------------

test_that("plot_vi_sword runs without error", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot_vi_sword(forest))
})

test_that("plot_vi_sword top_n truncates correctly", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot_vi_sword(forest, top_n = 2L))
})

# --- plot_oob_fit -------------------------------------------------------------

test_that("plot_oob_fit runs without error", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(plot_oob_fit(forest, y))
})

test_that("plot_oob_fit errors when OOB = FALSE", {
  f_no_oob <- SWORD(X, y, m = 3, OOB = FALSE, verbose = FALSE)
  expect_error(plot_oob_fit(f_no_oob, y), "OOB = TRUE")
})

# --- pdp_sword ----------------------------------------------------------------

test_that("pdp_sword runs without error on continuous variable", {
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(pdp_sword(forest, X, "x1"))
})

test_that("pdp_sword errors on unknown variable", {
  expect_error(pdp_sword(forest, X, "zzz"), "not found")
})

test_that("pdp_sword runs on discrete variable (few unique values)", {
  set.seed(3)
  Xd <- data.frame(
    x1  = rnorm(80),
    grp = sample(1:4, 80, replace = TRUE)
  )
  yd <- Xd$x1 + rnorm(80)
  fd <- SWORD(Xd, yd, m = 4, OOB = FALSE, verbose = FALSE)
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  on.exit({ dev.off(); unlink(tf) }, add = TRUE)
  expect_no_error(pdp_sword(fd, Xd, "grp"))
})
