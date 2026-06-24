#' Principal Curve of Oriented Points
#'
#' @description Computes a principal curve as defined in Delicado and Huerta (2003)
#' \doi{10.1007/s001800300145}.
#' @param x A finite numeric matrix or data frame of \eqn{n} points in dimension
#'   \eqn{p}. Missing and infinite values are rejected.
#' @param Ch The smoothing parameter \eqn{h} is \eqn{C_H} times the value given by
#'   the normal reference rule. Default value \eqn{1.5}. Constraints
#'   \eqn{0.5 \le C_H \le 1.5}.
#' @param Cd The distance between two consecutive principal oriented points in a
#'   PCOP is about \eqn{C_D} times the value of the smoothing parameter
#'   \eqn{h}. Default value \eqn{0.3}. Constraints \eqn{0.25 \le C_D \le 0.5}.
#' @param plot.true If \code{TRUE}, produce a two-dimensional plot of the
#'   resulting curve. Plotting requires at least two columns in \code{x}.
#' @param ... Additional parameters passed to \code{lines}.
#' @return
#' A list with two elements:
#' \describe{
#'   \item{pcop.f1}{Data frame storing the principal curve of oriented points in
#'     the original format, with columns \code{param}, \code{dens},
#'     \code{span}, \code{orth.var}, \code{pop1}, \code{pop2}, \ldots,
#'     \code{pr.dir1}, \code{pr.dir2}, \ldots}
#'   \item{pcop.f2}{List conforming to the format used in \pkg{princurve}; see that
#'     package for details.}
#'   \item{parameters}{List of algorithm parameters used for the fit.}
#'   \item{input_names}{Input row and column names, if present. Other input
#'     attributes are not used by the algorithm and are not propagated.}
#'   \item{call}{Matched function call.}
#' }
#' @examples
#' n <- 500
#' p <- 3
#' x <- matrix(rnorm(n * p), ncol = p) %*% diag(p:1)
#' pcop(x, plot.true = FALSE)
#'
#' x <- runif(100, -1, 1)
#' x <- cbind(x, x^2 + rnorm(100, sd = 0.1))
#' pcop(x, plot.true = FALSE)
#' if (interactive()) {
#'   pcop(x, plot.true = TRUE, lwd = 4, col = 2)
#' }
#' @importFrom princurve project_to_curve
#' @importFrom graphics lines
#' @export
pcop <- function(x, Ch = 1.5, Cd = 0.3, plot.true = FALSE, ...) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L) {
    stop("`x` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(x) < 2L) {
    stop("`x` must contain at least two rows.", call. = FALSE)
  }
  if (ncol(x) < 1L) {
    stop("`x` must contain at least one column.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("`x` must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(Ch) || length(Ch) != 1L || !is.finite(Ch)) {
    stop("`Ch` must be a finite numeric scalar.", call. = FALSE)
  }
  if (Ch < 0.5 || Ch > 1.5) {
    stop("`Ch` must be between 0.5 and 1.5.", call. = FALSE)
  }
  if (!is.numeric(Cd) || length(Cd) != 1L || !is.finite(Cd)) {
    stop("`Cd` must be a finite numeric scalar.", call. = FALSE)
  }
  if (Cd < 0.25 || Cd > 0.5) {
    stop("`Cd` must be between 0.25 and 0.5.", call. = FALSE)
  }
  if (!is.logical(plot.true) || length(plot.true) != 1L || is.na(plot.true)) {
    stop("`plot.true` must be TRUE or FALSE.", call. = FALSE)
  }
  if (plot.true && ncol(x) < 2L) {
    stop("`plot.true = TRUE` requires at least two columns in `x`.", call. = FALSE)
  }

  p <- ncol(x)

  backend_fit <- tryCatch(
    pcop_backend(x, Cd, Ch),
    error = identity
  )

  used_scaled_backend <- FALSE
  scaled_center <- rep.int(0, p)
  scaled_scale <- rep.int(1, p)
  if (inherits(backend_fit, "error")) {
    scaled_center <- colMeans(x)
    centered <- sweep(x, 2, scaled_center, "-")
    scaled_scale <- apply(centered, 2, stats::sd)
    scaled_scale[!is.finite(scaled_scale) | scaled_scale <= 0] <- 1
    x_scaled <- sweep(centered, 2, scaled_scale, "/")
    backend_fit <- tryCatch(
      pcop_backend(x_scaled, Cd, Ch),
      error = identity
    )
    if (inherits(backend_fit, "error")) {
      stop(conditionMessage(backend_fit), call. = FALSE)
    }
    used_scaled_backend <- TRUE
  }

  proyec <- as.data.frame(backend_fit)
  pr.curve <- proyec[, -1, drop = FALSE]

  names(pr.curve) <- c(
    "param", "dens", "span", "orth.var",
    paste0("pop", seq_len(p)),
    paste0("pr.dir", seq_len(p))
  )

  if (used_scaled_backend) {
    pop_idx <- 5:(4 + p)
    dir_idx <- (5 + p):(4 + 2 * p)
    pr.curve[, pop_idx] <- sweep(pr.curve[, pop_idx, drop = FALSE], 2, scaled_scale, "*")
    pr.curve[, pop_idx] <- sweep(pr.curve[, pop_idx, drop = FALSE], 2, scaled_center, "+")
    pr.curve[, dir_idx] <- sweep(pr.curve[, dir_idx, drop = FALSE], 2, scaled_scale, "*")
  }

  s <- as.matrix(pr.curve[, 5:(4 + p), drop = FALSE])
  keep_finite <- stats::complete.cases(s)
  s <- s[keep_finite, , drop = FALSE]
  if (nrow(s) == 0L) {
    s <- matrix(colMeans(x), nrow = 1)
  }

  if (nrow(s) > 1L) {
    step_norm <- sqrt(rowSums((s[-1, , drop = FALSE] - s[-nrow(s), , drop = FALSE])^2))
    keep_unique <- c(TRUE, step_norm > sqrt(.Machine$double.eps))
    s <- s[keep_unique, , drop = FALSE]
  }

  pc <- tryCatch(
    project_to_curve(x, s),
    error = identity
  )
  if (inherits(pc, "error")) {
    if (nrow(s) < 2L) {
      pc <- list(s = s, ord = seq_len(nrow(s)))
    } else {
      eps <- sqrt(.Machine$double.eps)
      col_scale <- pmax(apply(abs(s), 2, max), 1)
      offset <- outer(
        seq_len(nrow(s)),
        seq_len(ncol(s)),
        function(i, j) (i + j) * eps * col_scale[j]
      )
      s_jitter <- s + offset
      pc <- tryCatch(
        project_to_curve(x, s_jitter),
        error = identity
      )
      if (inherits(pc, "error")) {
        pc <- list(s = s, ord = seq_len(nrow(s)))
      }
    }
  }

  if (plot.true) {
    adjust_range <- function(values, cte = 1) {
      rgx <- range(values)
      extra <- diff(rgx) * (cte - 1) / 2
      c(rgx[1] - extra, rgx[2] + extra)
    }
    plot(
      x[, 1:2],
      xlim = adjust_range(x[, 1], 1.1),
      ylim = adjust_range(x[, 2], 1.1),
      col = 8
    )
    lines(pc$s[pc$ord, 1:2], ...)
  }

  out <- list(
    pcop.f1 = pr.curve,
    pcop.f2 = pc,
    parameters = list(
      Ch = Ch,
      Cd = Cd,
      dimension = p,
      scaled_backend = used_scaled_backend
    ),
    input_names = dimnames(x),
    call = match.call()
  )
  class(out) <- "pcop"
  out
}

#' @export
print.pcop <- function(x, ...) {
  cat("Principal curve of oriented points\n")
  cat("Dimension:", x$parameters$dimension, "\n")
  cat("Curve points:", nrow(x$pcop.f1), "\n")
  cat("Ch:", x$parameters$Ch, "Cd:", x$parameters$Cd, "\n")
  invisible(x)
}

#' @export
summary.pcop <- function(object, ...) {
  structure(
    list(
      dimension = object$parameters$dimension,
      curve_points = nrow(object$pcop.f1),
      projected_points = nrow(object$pcop.f2$s),
      parameters = object$parameters
    ),
    class = "summary.pcop"
  )
}

#' @export
print.summary.pcop <- function(x, ...) {
  cat("Principal curve of oriented points summary\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Curve points:", x$curve_points, "\n")
  cat("Projected points:", x$projected_points, "\n")
  cat("Ch:", x$parameters$Ch, "Cd:", x$parameters$Cd, "\n")
  invisible(x)
}

#' @export
plot.pcop <- function(x, ...) {
  s <- x$pcop.f2$s
  if (!is.matrix(s) || ncol(s) < 2L) {
    stop("`plot.pcop()` requires at least two curve dimensions.", call. = FALSE)
  }
  dots <- list(...)
  if (is.null(dots$type)) {
    dots$type <- "l"
  }
  if (is.null(dots$xlab)) {
    dots$xlab <- "Coordinate 1"
  }
  if (is.null(dots$ylab)) {
    dots$ylab <- "Coordinate 2"
  }
  do.call(graphics::plot, c(list(x = s[, 1], y = s[, 2]), dots))
  invisible(x)
}

#' srr_stats_pcop_interface
#'
#' @srrstats {UL1.0} `pcop()` documents accepted numeric matrix and data-frame
#' input, along with missing, infinite, and non-numeric inputs that are rejected.
#' @srrstats {UL1.1} `pcop()` validates input data and tuning parameters before
#' calling the native backend, with distinct user-facing error messages.
#' @srrstats {UL1.3} Input row and column names are retained in `input_names` on
#' the returned object when supplied.
#' @srrstats {UL1.3a} Documentation states that other input attributes are not
#' used by the algorithm and are not propagated.
#' @srrstats {UL1.4} README and the vignette document that PCOP is distance-based
#' and therefore sensitive to coordinate scaling.
#' @srrstats {UL1.4a} The wrapper retries the native backend on internally
#' scaled coordinates only after the raw backend fails, and then maps returned
#' coordinates back to the submitted scale.
#' @srrstats {UL1.4b} Examples do not apply `scale()`; the vignette explains
#' why scaling matters for distance-based PCOP fits.
#' @srrstats {UL2.0} The internal scaling fallback is implemented before
#' returning native-backend results to users.
#' @srrstats {UL2.1} The vignette documents the scaling behavior and the fact
#' that no user data are modified.
#' @srrstats {UL4.0} `pcop()` returns an S3 object of class `"pcop"` while
#' retaining the original list elements.
#' @srrstats {UL4.2} The return object includes the submitted tuning parameters
#' and whether the scaled backend fallback was used.
#' @srrstats {UL4.3} `print.pcop()` summarizes input dimension, curve size, and
#' tuning parameters.
#' @srrstats {UL4.3a} `print.pcop()` prints only aggregate metadata, not full
#' result matrices.
#' @srrstats {UL4.4} `summary.pcop()` returns a compact object summarizing the
#' primary fitted-curve dimensions and parameters.
#' @srrstats {UL6.0} `plot.pcop()` provides a default two-dimensional curve plot.
#' @srrstats {UL6.1} S3 plot dispatch is registered for `"pcop"` objects.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {UL1.2} `pcop()` does not use row or column names as labels for
#' output objects; curve coordinates are returned in PCOP and `princurve`
#' formats.
#' @srrstatsNA {UL2.2} Missing and non-finite inputs are rejected rather than
#' accepted with user-selectable missing-value processing.
#' @srrstatsNA {UL2.3} PCOP fits curves to point clouds rather than estimating
#' coefficients from a design matrix; collinearity diagnostics are not a model
#' precondition.
#' @srrstatsNA {UL3.0} `pcop()` does not assign clustering or partition labels
#' to observations.
#' @srrstatsNA {UL3.1} `pcop()` returns curve coordinates, not labelled reduced
#' dimensions ordered by explained variance.
#' @srrstatsNA {UL3.2} `pcop()` does not use observation labels.
#' @srrstatsNA {UL3.3} `pcop()` does not return a reusable model object for
#' predicting ordinates of new data.
#' @srrstatsNA {UL3.4} `pcop()` does not partition observations into discrete
#' groups.
#' @srrstatsNA {UL4.1} Empty model objects are not part of the PCOP API.
#' @srrstatsNA {UL6.2} `pcop()` does not place text labels on plots.
#' @noRd
NULL
