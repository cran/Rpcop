stress_n_for_dim <- function(p) {
  50L * as.integer(p)
}

stress_run_slow_cases <- function() {
  tolower(Sys.getenv("RPCOP_RUN_SLOW_STRESS", unset = "false")) %in%
    c("1", "true", "yes")
}

stress_gen_gaussian <- function(n, p, seed) {
  set.seed(seed)
  matrix(rnorm(n * p), ncol = p)
}

stress_gen_anisotropic <- function(n, p, seed) {
  z <- stress_gen_gaussian(n, p, seed)
  scales <- exp(seq(0, log(1e3), length.out = p))
  sweep(z, 2, scales, "*")
}

stress_gen_rank2 <- function(n, p, seed) {
  set.seed(seed)
  u <- matrix(rnorm(n * 2), ncol = 2)
  loadings <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
  signal <- u %*% t(loadings)
  noise <- matrix(rnorm(n * p, sd = 1e-4), ncol = p)
  signal + noise
}

stress_gen_clusters_outliers <- function(n, p, seed) {
  set.seed(seed)
  k <- 4
  centers <- matrix(rnorm(k * p, sd = 3), nrow = k, ncol = p)
  cluster <- sample.int(k, n, replace = TRUE, prob = c(0.40, 0.30, 0.20, 0.10))
  x <- centers[cluster, , drop = FALSE] + matrix(rnorm(n * p, sd = 0.30), ncol = p)
  out_n <- max(2L, floor(0.02 * n))
  out_idx <- sample.int(n, out_n)
  x[out_idx, ] <- x[out_idx, ] + matrix(rnorm(out_n * p, sd = 15), ncol = p)
  x
}

stress_gen_heavy_tail <- function(n, p, seed) {
  set.seed(seed)
  matrix(rt(n * p, df = 3), ncol = p)
}

stress_gen_duplicates <- function(n, p, seed) {
  set.seed(seed)
  pool_n <- max(20L, ceiling(0.12 * n))
  pool <- matrix(rnorm(pool_n * p), ncol = p)
  idx <- sample.int(pool_n, n, replace = TRUE)
  x <- pool[idx, , drop = FALSE]
  jit_n <- max(1L, floor(0.10 * n))
  jit_idx <- sample.int(n, jit_n)
  x[jit_idx, ] <- x[jit_idx, ] + matrix(rnorm(jit_n * p, sd = 1e-5), ncol = p)
  x
}

stress_gen_scaled <- function(n, p, seed, scale_factor) {
  scale_factor * stress_gen_gaussian(n, p, seed)
}

stress_expected_cols <- function(p) {
  c(
    "param",
    "dens",
    "span",
    "orth.var",
    paste0("pop", seq_len(p)),
    paste0("pr.dir", seq_len(p))
  )
}

stress_backend_cases <- function(include_slow = stress_run_slow_cases()) {
  cases <- list(
    list(
      id = "gaussian_p2",
      p = 2L,
      x = stress_gen_gaussian(stress_n_for_dim(2L), 2L, seed = 1002L),
      Cd = 0.30,
      Ch = 1.00
    ),
    list(
      id = "anisotropic_p3",
      p = 3L,
      x = stress_gen_anisotropic(stress_n_for_dim(3L), 3L, seed = 2003L),
      Cd = 0.25,
      Ch = 1.50
    ),
    list(
      id = "rank2_p5",
      p = 5L,
      x = stress_gen_rank2(stress_n_for_dim(5L), 5L, seed = 3005L),
      Cd = 0.50,
      Ch = 0.50
    ),
    list(
      id = "clusters_outliers_p5",
      p = 5L,
      x = stress_gen_clusters_outliers(stress_n_for_dim(5L), 5L, seed = 4010L),
      Cd = 0.30,
      Ch = 1.50
    )
  )
  if (include_slow) {
    cases[[length(cases) + 1L]] <- list(
      id = "clusters_outliers_p8",
      p = 8L,
      x = stress_gen_clusters_outliers(stress_n_for_dim(8L), 8L, seed = 4010L),
      Cd = 0.30,
      Ch = 1.50
    )
  }
  cases
}

stress_wrapper_cases <- function() {
  cases <- list(
    list(
      id = "wrapper_parabola_p2",
      p = 2L,
      x = local({
        set.seed(5001L)
        n <- 450L
        t <- runif(n, -2, 2)
        cbind(t, t^2 + rnorm(n, sd = 0.07))
      }),
      Cd = 0.30,
      Ch = 1.50
    ),
    list(
      id = "wrapper_scurve_p3",
      p = 3L,
      x = local({
        set.seed(5002L)
        n <- 420L
        t <- runif(n, -pi, pi)
        cbind(sin(t), t / pi, sign(t) * (1 - cos(t)) + rnorm(n, sd = 0.05))
      }),
      Cd = 0.25,
      Ch = 1.00
    ),
    list(
      id = "wrapper_rank2_p5",
      p = 5L,
      x = stress_gen_rank2(stress_n_for_dim(5L), 5L, seed = 5003L),
      Cd = 0.50,
      Ch = 0.50
    )
  )
  cases[[length(cases) + 1L]] <- list(
    id = "wrapper_anisotropic_p8",
    p = 8L,
    x = stress_gen_anisotropic(stress_n_for_dim(8L), 8L, seed = 5004L),
    Cd = 0.30,
    Ch = 1.00
  )
  cases
}
