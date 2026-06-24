#' srr_stats_pcop_validation_tests
#'
#' @srrstats {G5.2} Validation tests cover representative user-facing errors.
#' @srrstats {G5.2b} Each validation branch added in `pcop()` has a direct
#' expectation against its error message.
#' @srrstats {G5.8} Edge-condition tests cover invalid types, non-finite values,
#' empty input, tuning-parameter ranges, and plotting preconditions.
#' @srrstats {UL7.0} Inappropriate input data types are rejected with expected
#' error messages.
#' @noRd
NULL

test_that("pcop validates x before calling native code", {
  expect_error(
    pcop(data.frame(a = c("a", "b"), b = c("c", "d"))),
    "`x` must be a numeric matrix or data frame.",
    fixed = TRUE
  )
  expect_error(
    pcop(matrix(c(1, NA, 2, 3), ncol = 2)),
    "`x` must contain only finite values.",
    fixed = TRUE
  )
  expect_error(
    pcop(matrix(numeric(), nrow = 0, ncol = 2)),
    "`x` must contain at least two rows.",
    fixed = TRUE
  )
})

test_that("pcop validates tuning parameters and plotting option", {
  x <- matrix(seq_len(200), ncol = 2)

  expect_error(
    pcop(x, Ch = NA_real_),
    "`Ch` must be a finite numeric scalar.",
    fixed = TRUE
  )
  expect_error(
    pcop(x, Ch = 0.4),
    "`Ch` must be between 0.5 and 1.5.",
    fixed = TRUE
  )
  expect_error(
    pcop(x, Cd = c(0.3, 0.4)),
    "`Cd` must be a finite numeric scalar.",
    fixed = TRUE
  )
  expect_error(
    pcop(x, Cd = 0.6),
    "`Cd` must be between 0.25 and 0.5.",
    fixed = TRUE
  )
  expect_error(
    pcop(x, plot.true = NA),
    "`plot.true` must be TRUE or FALSE.",
    fixed = TRUE
  )
  expect_error(
    pcop(matrix(seq_len(100), ncol = 1), plot.true = TRUE),
    "`plot.true = TRUE` requires at least two columns in `x`.",
    fixed = TRUE
  )
})
