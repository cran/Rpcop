# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Computes a principal curve as defined in Delicado (2001). DO NOT use this function unless you know what you are doing. Use `pcop()` instead.
#'
#' @param x    See `pcop()`
#' @param c_h  See `pcop()`
#' @param c_d  See `pcop()`
#' @return A numeric matrix to be parsed by `pcop()`.
pcop_backend <- function(x, c_d, c_h) {
    .Call('_Rpcop_pcop_backend', PACKAGE = 'Rpcop', x, c_d, c_h)
}

