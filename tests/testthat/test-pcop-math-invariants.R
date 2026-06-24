test_that("backend is invariant to coordinate translation", {
  set.seed(14)
  n <- 120L
  t <- seq(-2, 2, length.out = n)
  x <- cbind(t, t^2 + rnorm(n, sd = 0.03))
  offset <- c(0, 5000)
  shifted <- sweep(x, 2, offset, "+")

  base <- Rpcop:::pcop_backend(x, 0.3, 1.5)
  shifted_fit <- Rpcop:::pcop_backend(shifted, 0.3, 1.5)

  expect_identical(dim(shifted_fit), dim(base))
  shifted_fit[, 6:7] <- sweep(shifted_fit[, 6:7, drop = FALSE], 2, offset, "-")
  expect_equal(shifted_fit, base, tolerance = 1e-5)
})

test_that("large finite coordinates do not crash the public wrapper", {
  set.seed(12)
  n <- 120L
  t <- seq(-2, 2, length.out = n)
  x <- cbind(t, t^2 + rnorm(n, sd = 0.03))
  shifted <- sweep(x, 2, c(20000, -15000), "+")

  fit <- pcop(shifted, plot.true = FALSE)
  expect_s3_class(fit, "pcop")
  expect_true(all(is.finite(as.matrix(fit$pcop.f1))))
  expect_false(fit$parameters$scaled_backend)
})

test_that("finalized densities integrate to one on the returned arc-length grid", {
  set.seed(15)
  n <- 120L
  t <- seq(-2, 2, length.out = n)
  x <- cbind(t, t^2 + rnorm(n, sd = 0.03))

  out <- Rpcop:::pcop_backend(x, 0.3, 1.5)
  param <- out[, 2]
  dens <- out[, 3]
  expect_true(length(param) > 1L)
  expect_true(all(diff(param) > 0))

  width <- numeric(length(param))
  width[1] <- param[2] - param[1]
  width[length(param)] <- param[length(param)] - param[length(param) - 1L]
  if (length(param) > 2L) {
    width[2:(length(param) - 1L)] <-
      (param[3:length(param)] - param[1:(length(param) - 2L)]) / 2
  }
  expect_equal(sum(dens * width), 1, tolerance = 1e-5)
})

test_that("one-dimensional input returns a valid trivial curve", {
  x <- matrix(seq_len(60), ncol = 1)

  fit <- pcop(x, plot.true = FALSE)
  expect_s3_class(fit, "pcop")
  expect_identical(dim(fit$pcop.f1), c(1L, 6L))
  expect_equal(fit$pcop.f1$pop1, mean(x))
  expect_true(is.matrix(fit$pcop.f2$s))
  expect_identical(ncol(fit$pcop.f2$s), 1L)
})

test_that("degenerate multidimensional inputs do not expose INF variance sentinel", {
  constant <- matrix(1, nrow = 100L, ncol = 2L)
  constant_fit <- pcop(constant, plot.true = FALSE)

  expect_s3_class(constant_fit, "pcop")
  expect_identical(nrow(constant_fit$pcop.f1), 1L)
  expect_true(is.finite(constant_fit$pcop.f1$orth.var))
  expect_equal(constant_fit$pcop.f1$orth.var, 0)

  collinear <- cbind(seq_len(100L), seq_len(100L))
  collinear_fit <- pcop(collinear, plot.true = FALSE)

  expect_s3_class(collinear_fit, "pcop")
  expect_identical(nrow(collinear_fit$pcop.f1), 1L)
  expect_true(is.finite(collinear_fit$pcop.f1$orth.var))
  expect_false(isTRUE(all.equal(collinear_fit$pcop.f1$orth.var, 9999)))
})
