#' srr_stats_pcop_stress_tests
#'
#' @srrstats {G5.3} Stress tests assert finite backend and wrapper outputs for
#' representative valid point clouds.
#' @srrstats {G5.6} Synthetic data include Gaussian, low-rank, clustered,
#' outlier, and anisotropic point clouds with fixed seeds.
#' @srrstats {G5.9} Fixed-seed stress fixtures document behavior under noisy and
#' anisotropic inputs.
#' @srrstats {UL7.1} Stress fixtures include anisotropic and low-rank point
#' clouds to cover scale-sensitive and near-collinear inputs.
#' @srrstats {UL7.3} Wrapper tests assert that submitted dimnames are retained in
#' `input_names`.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {UL7.2} `pcop()` does not produce group-size labels.
#' @srrstatsNA {UL7.4} `pcop()` does not implement prediction for new data.
#' @srrstatsNA {UL7.5} `pcop()` does not implement batch processing.
#' @srrstatsNA {UL7.5a} See UL7.5.
#' @noRd
NULL

test_that("hard backend stress cases remain finite and non-empty", {
  cases <- stress_backend_cases()
  expect_identical(length(cases), if (stress_run_slow_cases()) 5L else 4L)

  for (case in cases) {
    msg <- sprintf("backend case: %s", case$id)
    out <- tryCatch(
      Rpcop:::pcop_backend(case$x, case$Cd, case$Ch),
      error = identity
    )
    expect_false(inherits(out, "error"), info = msg)
    if (inherits(out, "error")) {
      next
    }
    expect_true(is.matrix(out), info = msg)
    expect_true(nrow(out) > 0L, info = msg)
    expect_true(ncol(out) > 0L, info = msg)
    expect_true(all(is.finite(out)), info = msg)
  }
})

test_that("hard wrapper stress cases return well-formed outputs", {
  cases <- stress_wrapper_cases()
  expect_identical(length(cases), 4L)

  for (case in cases) {
    msg <- sprintf("wrapper case: %s", case$id)
    fit <- tryCatch(
      pcop(case$x, Ch = case$Ch, Cd = case$Cd, plot.true = FALSE),
      error = identity
    )
    expect_false(inherits(fit, "error"), info = msg)
    if (inherits(fit, "error")) {
      next
    }

    expect_true(is.list(fit), info = msg)
    expect_true(all(c("pcop.f1", "pcop.f2") %in% names(fit)), info = msg)
    expect_s3_class(fit, "pcop")
    expect_equal(fit$parameters$Ch, case$Ch, info = msg)
    expect_equal(fit$parameters$Cd, case$Cd, info = msg)
    expect_identical(fit$parameters$dimension, case$p, info = msg)

    expect_true(is.data.frame(fit$pcop.f1), info = msg)
    expect_true(nrow(fit$pcop.f1) > 0L, info = msg)
    expect_identical(colnames(fit$pcop.f1), stress_expected_cols(case$p), info = msg)
    expect_true(all(is.finite(as.matrix(fit$pcop.f1))), info = msg)

    expect_true(is.list(fit$pcop.f2), info = msg)
    expect_true(is.matrix(fit$pcop.f2$s), info = msg)
    expect_true(nrow(fit$pcop.f2$s) > 0L, info = msg)
    expect_true(all(is.finite(fit$pcop.f2$s)), info = msg)
  }
})

test_that("pcop methods expose compact metadata and plotting", {
  case <- stress_wrapper_cases()[[1L]]
  rownames(case$x) <- sprintf("row%03d", seq_len(nrow(case$x)))
  colnames(case$x) <- c("x", "y")

  fit <- pcop(case$x, Ch = case$Ch, Cd = case$Cd, plot.true = FALSE)
  expect_s3_class(fit, "pcop")
  expect_identical(fit$input_names, dimnames(case$x))

  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary.pcop")
  expect_identical(summary_fit$dimension, case$p)
  expect_identical(summary_fit$parameters$dimension, case$p)

  expect_output(print(fit), "Principal curve of oriented points")
  expect_output(print(summary_fit), "Principal curve of oriented points summary")

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit))
})
