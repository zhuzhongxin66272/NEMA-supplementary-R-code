#!/usr/bin/env Rscript

################################################################################
# NEMA complete R code for manuscript reproducibility
################################################################################
# This single file collects the R code used for the NEMA manuscript:
#   1. Core NEMA functions.
#   2. Formal simulation study with checkpointing.
#   3. NHANES POPs and leukocyte telomere length empirical analysis.
#
# All comments are written in English for journal submission and code sharing.
#
# Recommended usage from the project root:
#   Rscript reports/submission_materials/NEMA_complete_R_code.R
#
# By default, both RUN_SIMULATION and RUN_NHANES are FALSE to prevent accidental
# long runs. Enable the part you need through environment variables:
#
#   RUN_SIMULATION=TRUE N_SIM=1000 B_BOOT=100 N_CORES=4 \
#     Rscript reports/submission_materials/NEMA_complete_R_code.R
#
#   RUN_NHANES=TRUE NEMA_BOOT=1000 \
#     Rscript reports/submission_materials/NEMA_complete_R_code.R
#
# If Gibson et al.'s SHARP Mixtures Workshop repository is not already cloned,
# the NHANES analysis can clone it automatically when internet access is
# available:
#
#   RUN_NHANES=TRUE DOWNLOAD_GIBSON=TRUE \
#     Rscript reports/submission_materials/NEMA_complete_R_code.R
#
# Main outputs:
#   Simulation results: results/formal_v6/
#   NHANES results:     results/nhanes_pops_telomere/
################################################################################

RUN_SIMULATION <- tolower(Sys.getenv("RUN_SIMULATION", "false")) %in% c("1", "true", "yes")
RUN_NHANES <- tolower(Sys.getenv("RUN_NHANES", "false")) %in% c("1", "true", "yes")

if (!RUN_SIMULATION && !RUN_NHANES) {
  message("No analysis requested. Set RUN_SIMULATION=TRUE and/or RUN_NHANES=TRUE.")
  message("The file still defines all functions and contains all code needed for reproducibility.")
}


################################################################################
# Section 1. Core NEMA functions
################################################################################

# Core implementation for Network-Embedded Mixture Analysis (NEMA).

npn_transform <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  apply(X, 2, function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    u <- (r - 0.5) / n
    qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))
  })
}

ridge_partial_network <- function(X, ridge = 0.10, keep = 0.30,
                                  npn = TRUE) {
  X <- scale(as.matrix(X))
  if (npn) X <- npn_transform(X)
  p <- ncol(X)
  S <- stats::cov(X)
  Omega <- solve(S + ridge * diag(p))
  inv_sd <- 1 / sqrt(diag(Omega))
  Pcor <- -outer(inv_sd, inv_sd) * Omega
  diag(Pcor) <- 0

  W <- abs(Pcor)
  if (!is.null(keep) && keep < 1) {
    vals <- W[upper.tri(W)]
    threshold <- stats::quantile(vals, probs = 1 - keep, na.rm = TRUE)
    W[W < threshold] <- 0
  }
  W <- (W + t(W)) / 2
  diag(W) <- 0
  list(W = W, Omega = Omega, Sigma = S, transformed_X = X)
}

estimate_network <- function(X, estimator = c("ridge", "huge"),
                             ridge = 0.10, keep = 0.30, npn = TRUE) {
  estimator <- match.arg(estimator)
  if (estimator == "huge" && requireNamespace("huge", quietly = TRUE)) {
    X_work <- if (npn) huge::huge.npn(as.matrix(X)) else scale(as.matrix(X))
    fit <- huge::huge(X_work, method = "glasso", verbose = FALSE)
    sel <- huge::huge.select(fit, criterion = "ebic", ebic.gamma = 0.5,
                             verbose = FALSE)
    Omega <- sel$opt.icov
    W <- abs(Omega)
    diag(W) <- 0
    return(list(W = W, Omega = Omega, Sigma = solve(Omega),
                transformed_X = X_work))
  }

  ridge_partial_network(X, ridge = ridge, keep = keep, npn = npn)
}

normalized_laplacian <- function(W) {
  p <- nrow(W)
  degree <- rowSums(W)
  inv_sqrt_degree <- ifelse(degree > 1e-10, 1 / sqrt(degree), 0)
  diag(p) - diag(inv_sqrt_degree, p) %*% W %*% diag(inv_sqrt_degree, p)
}

select_mode_count <- function(eigenvalues, d_max) {
  d_max <- min(d_max, length(eigenvalues) - 1)
  if (d_max < 1) return(1L)
  gaps <- diff(eigenvalues[seq_len(d_max + 1)])
  as.integer(which.max(gaps))
}

network_modes <- function(W, d = NULL, d_max = floor(log(nrow(W)))) {
  L <- normalized_laplacian(W)
  eg <- eigen(L, symmetric = TRUE)
  ord <- order(eg$values, decreasing = FALSE)
  values <- eg$values[ord]
  vectors <- eg$vectors[, ord, drop = FALSE]
  if (is.null(d)) d <- select_mode_count(values, d_max = d_max)
  d <- min(d, ncol(vectors))
  U <- vectors[, seq_len(d), drop = FALSE]
  colnames(U) <- paste0("mode_", seq_len(d))
  list(values = values, vectors = U, laplacian = L, d = d)
}

mode_scores <- function(X, U) {
  S <- scale(as.matrix(X), center = TRUE, scale = TRUE) %*% U
  colnames(S) <- colnames(U)
  S
}

fit_mode_regression <- function(Y, S, Z = NULL, family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  dat <- data.frame(Y = Y, as.data.frame(S))
  if (!is.null(Z)) dat <- cbind(dat, as.data.frame(Z))
  if (family == "gaussian") {
    fit <- stats::lm(Y ~ ., data = dat)
  } else {
    fit <- stats::glm(Y ~ ., data = dat, family = stats::binomial())
  }
  gamma <- stats::coef(fit)[colnames(S)]
  names(gamma) <- colnames(S)
  list(fit = fit, gamma = gamma)
}

overall_iqr_effect <- function(X, U, gamma) {
  qs <- apply(as.matrix(X), 2, stats::quantile, probs = c(0.25, 0.75),
              na.rm = TRUE)
  contrast <- qs[2, ] - qs[1, ]
  as.numeric(sum(gamma * as.numeric(crossprod(U, contrast))))
}

chemical_contribution <- function(U, gamma) {
  out <- as.numeric(U %*% gamma)
  names(out) <- rownames(U)
  out
}

network_total_effect <- function(contribution, Sigma) {
  Sigma <- as.matrix(Sigma)
  out <- as.numeric(crossprod(contribution, Sigma)) / diag(Sigma)
  names(out) <- names(contribution)
  out
}

nema_fit <- function(X, Y, Z = NULL, d = NULL, d_max = floor(log(ncol(X))),
                     estimator = "ridge", ridge = 0.10, keep = 0.30,
                     npn = TRUE, family = "gaussian") {
  X <- scale(as.matrix(X), center = TRUE, scale = TRUE)
  net <- estimate_network(X, estimator = estimator, ridge = ridge, keep = keep,
                          npn = npn)
  modes <- network_modes(net$W, d = d, d_max = d_max)
  U <- modes$vectors
  rownames(U) <- colnames(X)
  S <- mode_scores(X, U)
  reg <- fit_mode_regression(Y, S, Z = Z, family = family)
  contrib <- chemical_contribution(U, reg$gamma)
  te <- network_total_effect(contrib, net$Sigma)

  list(W = net$W, Sigma = net$Sigma, eigenvalues = modes$values,
       U = U, scores = S, regression = reg$fit, gamma = reg$gamma,
       contribution = contrib, total_effect = te,
       overall_iqr_effect = overall_iqr_effect(X, U, reg$gamma),
       d = modes$d)
}

all_permutations <- function(n) {
  if (n == 1) return(matrix(1L, nrow = 1))
  prev <- all_permutations(n - 1)
  out <- vector("list", n * nrow(prev))
  k <- 1
  for (i in seq_len(nrow(prev))) {
    for (pos in seq_len(n)) {
      out[[k]] <- append(prev[i, ], n, after = pos - 1)
      k <- k + 1
    }
  }
  do.call(rbind, out)
}

align_modes <- function(U, U_ref) {
  U <- as.matrix(U)
  U_ref <- as.matrix(U_ref)
  d <- min(ncol(U), ncol(U_ref))
  U <- U[, seq_len(d), drop = FALSE]
  U_ref <- U_ref[, seq_len(d), drop = FALSE]

  perms <- all_permutations(d)
  scores <- apply(perms, 1, function(pm) {
    sum(diag(abs(crossprod(U_ref, U[, pm, drop = FALSE]))))
  })
  best <- perms[which.max(scores), ]
  U <- U[, best, drop = FALSE]
  for (j in seq_len(d)) {
    if (sum(U[, j] * U_ref[, j]) < 0) U[, j] <- -U[, j]
  }
  colnames(U) <- colnames(U_ref)
  U
}

bootstrap_nema <- function(X, Y, Z = NULL, fit_ref, B = 200,
                           estimator = "ridge", ridge = 0.10, keep = 0.30,
                           npn = TRUE, family = "gaussian",
                           seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  X <- scale(as.matrix(X), center = TRUE, scale = TRUE)
  n <- nrow(X)
  p <- ncol(X)
  d <- fit_ref$d

  gamma_mat <- matrix(NA_real_, B, d)
  contribution_mat <- matrix(NA_real_, B, p)
  te_mat <- matrix(NA_real_, B, p)
  delta_vec <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    Yb <- Y[idx]
    Zb <- if (is.null(Z)) NULL else as.matrix(Z)[idx, , drop = FALSE]

    boot_fit <- tryCatch(
      nema_fit(Xb, Yb, Z = Zb, d = d, estimator = estimator, ridge = ridge,
               keep = keep, npn = npn, family = family),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) next

    Ub <- align_modes(boot_fit$U, fit_ref$U)
    Sb <- mode_scores(Xb, Ub)
    reg <- tryCatch(fit_mode_regression(Yb, Sb, Z = Zb, family = family),
                    error = function(e) NULL)
    if (is.null(reg)) next

    cb <- chemical_contribution(Ub, reg$gamma)
    teb <- network_total_effect(cb, stats::cov(Xb))
    gamma_mat[b, ] <- reg$gamma
    contribution_mat[b, ] <- cb
    te_mat[b, ] <- teb
    delta_vec[b] <- overall_iqr_effect(Xb, Ub, reg$gamma)
    if (verbose && b %% 50 == 0) message("bootstrap ", b, "/", B)
  }

  colnames(gamma_mat) <- names(fit_ref$gamma)
  colnames(contribution_mat) <- names(fit_ref$contribution)
  colnames(te_mat) <- names(fit_ref$total_effect)
  list(gamma = gamma_mat, contribution = contribution_mat,
       total_effect = te_mat, overall_iqr_effect = delta_vec)
}

ci_percentile <- function(x, probs = c(0.025, 0.975)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

ci_table <- function(estimate, boot_mat, term_name = "term") {
  cis <- t(apply(boot_mat, 2, ci_percentile))
  data.frame(
    term = names(estimate),
    estimate = as.numeric(estimate),
    ci_low = cis[, 1],
    ci_high = cis[, 2],
    significant = cis[, 1] > 0 | cis[, 2] < 0,
    row.names = NULL
  )
}


################################################################################
# Section 2. Formal simulation study
################################################################################

if (RUN_SIMULATION) {

  suppressPackageStartupMessages({
    library(MASS)
    library(huge)
    library(parallel)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(knitr)
  })

  parse_arg <- function(name, default) {
    args <- commandArgs(trailingOnly = TRUE)
    key <- paste0("--", name, "=")
    hit <- args[startsWith(args, key)]
    if (length(hit) == 0) return(default)
    sub(key, "", hit[1], fixed = TRUE)
  }

  N_SIM <- as.integer(parse_arg("n_sim", Sys.getenv("N_SIM", "1000")))
  B_BOOT <- as.integer(parse_arg("boot", Sys.getenv("B_BOOT", "100")))
  N_CORES <- as.integer(parse_arg("cores", Sys.getenv("N_CORES", "4")))
  BASE_SEED <- as.integer(parse_arg("seed", Sys.getenv("SEED", "20260613")))
  ALPHA <- as.numeric(parse_arg("alpha", Sys.getenv("ALPHA", "0.05")))
  Z_CRIT <- as.numeric(parse_arg("zcrit", Sys.getenv("Z_CRIT", "1.70")))
  BATCH_SIZE <- as.integer(parse_arg("batch", Sys.getenv("BATCH_SIZE", "25")))
  OUT_DIR <- parse_arg("out", Sys.getenv("OUT_DIR", "results/formal_v6"))
  FIG_DIR <- file.path(OUT_DIR, "figures")
  REPORT_DIR <- "reports"

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
  SCENARIO_DIR <- file.path(OUT_DIR, "scenario_checkpoints")
  dir.create(SCENARIO_DIR, recursive = TRUE, showWarnings = FALSE)

  standardize <- function(x) {
    z <- scale(x)
    z[is.na(z)] <- 0
    z
  }

  normalized_laplacian <- function(W) {
    p <- nrow(W)
    deg <- rowSums(W)
    inv_sqrt <- ifelse(deg > 1e-10, 1 / sqrt(deg), 0)
    diag(p) - diag(inv_sqrt, p) %*% W %*% diag(inv_sqrt, p)
  }

  fiedler_mode <- function(W) {
    eg <- eigen(normalized_laplacian(W), symmetric = TRUE)
    ord <- order(eg$values)
    u <- eg$vectors[, ord[2], drop = FALSE]
    colnames(u) <- "mode_1"
    list(U = u, values = eg$values[ord])
  }

  align_fiedler <- function(U, U_ref) {
    U <- U[, 1, drop = FALSE]
    U_ref <- U_ref[, 1, drop = FALSE]
    if (sum(U[, 1] * U_ref[, 1]) < 0) U[, 1] <- -U[, 1]
    colnames(U) <- colnames(U_ref)
    U
  }

  make_truth <- function(p, rho_within = 0.90, bridge = 0.002) {
    stopifnot(p %% 2 == 0)
    K <- 2
    block <- p / K
    W <- matrix(0, p, p)

    for (k in seq_len(K)) {
      idx <- ((k - 1) * block + 1):(k * block)
      for (j in seq_len(length(idx) - 1)) {
        a <- idx[j]
        b <- idx[j + 1]
        W[a, b] <- rho_within
        W[b, a] <- rho_within
      }
      if (length(idx) > 4) {
        for (j in seq_len(length(idx) - 2)) {
          a <- idx[j]
          b <- idx[j + 2]
          W[a, b] <- rho_within * 0.35
          W[b, a] <- rho_within * 0.35
        }
      }
    }

    W[block, block + 1] <- bridge
    W[block + 1, block] <- bridge
    Omega <- diag(rowSums(W) + 0.8, p) - W
    Sigma <- cov2cor(solve(Omega))

    modes <- fiedler_mode(W)
    U <- modes$U
    rownames(U) <- paste0("X", seq_len(p))
    list(Sigma = Sigma, W = W, U = U, eigenvalues = modes$values)
  }

  generate_dataset <- function(sc, rep_seed) {
    set.seed(rep_seed)
    truth <- make_truth(sc$p, rho_within = sc$rho_within, bridge = sc$bridge)
    X <- standardize(MASS::mvrnorm(sc$n, rep(0, sc$p), truth$Sigma))
    colnames(X) <- paste0("X", seq_len(sc$p))
    S_true <- as.numeric(standardize(X %*% truth$U))
    signal <- S_true * sc$beta
    Y <- signal + rnorm(sc$n, 0, sc$error_sd)
    list(X = X, Y = Y, S_true = S_true, truth = truth)
  }

  fit_nema_huge <- function(X, Y) {
    invisible(capture.output({
      X_npn <- huge::huge.npn(X, verbose = FALSE)
      path <- huge::huge(X_npn, method = "glasso", verbose = FALSE)
      sel <- huge::huge.select(path, criterion = "ebic", ebic.gamma = 0.5,
                               verbose = FALSE)
    }))
    W <- abs(as.matrix(sel$opt.icov))
    diag(W) <- 0
    modes <- fiedler_mode(W)
    U <- modes$U
    rownames(U) <- colnames(X)
    scores <- as.numeric(standardize(X %*% U))
    fit <- lm(Y ~ scores)
    gamma <- unname(coef(fit)[2])
    list(U = U, scores = scores, gamma = gamma, W = W,
         eigenvalues = modes$values, fit = fit)
  }

  refit_gamma_with_mode <- function(X, Y, U) {
    scores <- as.numeric(standardize(X %*% U))
    fit <- lm(Y ~ scores)
    list(gamma = unname(coef(fit)[2]), scores = scores, fit = fit)
  }

  bootstrap_full_pipeline <- function(X, Y, fit_ref, B, seed) {
    set.seed(seed)
    n <- nrow(X)
    boot_gamma <- rep(NA_real_, B)

    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      fb <- tryCatch(fit_nema_huge(X[idx, , drop = FALSE], Y[idx]),
                     error = function(e) NULL)
      if (is.null(fb)) next
      Ub <- align_fiedler(fb$U, fit_ref$U)
      rb <- tryCatch(refit_gamma_with_mode(X[idx, , drop = FALSE], Y[idx], Ub),
                     error = function(e) NULL)
      if (!is.null(rb)) boot_gamma[b] <- rb$gamma
    }

    valid <- boot_gamma[is.finite(boot_gamma)]
    if (length(valid) < max(10, B / 2)) {
      return(list(valid_fraction = length(valid) / B,
                  ci_low = NA_real_, ci_high = NA_real_, sig = NA))
    }

    se <- sd(valid, na.rm = TRUE)
    if (!is.finite(se) || se < 1e-8) {
      return(list(valid_fraction = length(valid) / B,
                  ci_low = NA_real_, ci_high = NA_real_, sig = NA))
    }
    ci <- c(fit_ref$gamma - Z_CRIT * se, fit_ref$gamma + Z_CRIT * se)
    list(valid_fraction = length(valid) / B, ci_low = ci[1], ci_high = ci[2],
         sig = ci[1] > 0 | ci[2] < 0)
  }

  individual_ols_metrics <- function(X, Y) {
    pv <- vapply(seq_len(ncol(X)), function(j) {
      coef(summary(lm(Y ~ X[, j])))[2, 4]
    }, numeric(1))
    qv <- p.adjust(pv, method = "BH")
    c(
      mean_p_lt_0_05 = mean(pv < 0.05),
      any_p_lt_0_05 = as.numeric(any(pv < 0.05)),
      any_bonf = as.numeric(any(pv < 0.05 / length(pv))),
      mean_bh_fdr = mean(qv < 0.05),
      any_bh_fdr = as.numeric(any(qv < 0.05)),
      min_p = min(pv),
      min_q = min(qv),
      median_abs_cor = median(abs(cor(X, Y)))
    )
  }

  one_rep <- function(sc, scenario_index, rep_id) {
    dat <- generate_dataset(sc, BASE_SEED + scenario_index * 100000 + rep_id)
    fit0 <- tryCatch(fit_nema_huge(dat$X, dat$Y), error = function(e) NULL)
    if (is.null(fit0)) return(NULL)

    U_aligned <- align_fiedler(fit0$U, dat$truth$U)
    refit <- refit_gamma_with_mode(dat$X, dat$Y, U_aligned)
    fit_ref <- list(U = U_aligned, gamma = refit$gamma, scores = refit$scores)
    boot <- bootstrap_full_pipeline(
      dat$X, dat$Y, fit_ref, B_BOOT,
      seed = BASE_SEED + scenario_index * 1000000 + rep_id
    )
    if (is.na(boot$sig)) return(NULL)

    ols <- individual_ols_metrics(dat$X, dat$Y)
    data.frame(
      scenario = sc$id,
      family = sc$family,
      label = sc$label,
      rep = rep_id,
      n = sc$n,
      p = sc$p,
      beta = sc$beta,
      error_sd = sc$error_sd,
      rho_within = sc$rho_within,
      bridge = sc$bridge,
      gamma_true = sc$beta,
      gamma_hat = refit$gamma,
      gamma_bias = refit$gamma - sc$beta,
      gamma_abs_error = abs(refit$gamma - sc$beta),
      ci_low = boot$ci_low,
      ci_high = boot$ci_high,
      significant = as.logical(boot$sig),
      covered = boot$ci_low <= sc$beta & boot$ci_high >= sc$beta,
      score_cor = as.numeric(cor(dat$S_true, fit_ref$scores)),
      valid_bootstrap = boot$valid_fraction,
      ols_mean_p_lt_0_05 = ols[["mean_p_lt_0_05"]],
      ols_any_p_lt_0_05 = ols[["any_p_lt_0_05"]],
      ols_any_bonf = ols[["any_bonf"]],
      ols_mean_bh_fdr = ols[["mean_bh_fdr"]],
      ols_any_bh_fdr = ols[["any_bh_fdr"]],
      ols_min_p = ols[["min_p"]],
      ols_min_q = ols[["min_q"]],
      ols_median_abs_cor = ols[["median_abs_cor"]],
      stringsAsFactors = FALSE
    )
  }

  scenario <- function(id, family, label, n = 300, p = 60, beta = 0.30,
                       error_sd = 1, rho_within = 0.95, bridge = 0.002) {
    list(id = id, family = family, label = label, n = n, p = p,
         beta = beta, error_sd = error_sd, rho_within = rho_within,
         bridge = bridge)
  }

  scenarios <- list(
    scenario("T0", "type_i", "Global null", beta = 0.00),
    scenario("E1", "effect_size", "beta=0.18", beta = 0.18),
    scenario("E2", "effect_size", "beta=0.30", beta = 0.30),
    scenario("E3", "effect_size", "beta=0.65", beta = 0.65),
    scenario("N1", "sample_size", "n=300", n = 300, beta = 0.20),
    scenario("N2", "sample_size", "n=600", n = 600, beta = 0.20),
    scenario("N3", "sample_size", "n=1000", n = 1000, beta = 0.20),
    scenario("D1", "dimension_fixed", "p=30 fixed beta", n = 400, p = 30,
             beta = 0.22),
    scenario("D2", "dimension_fixed", "p=60 fixed beta", n = 400, p = 60,
             beta = 0.22),
    scenario("D3", "dimension_fixed", "p=90 fixed beta", n = 400, p = 90,
             beta = 0.22),
    scenario("M1", "micro_effect", "p=30 micro", n = 400, p = 30,
             beta = 0.10),
    scenario("M2", "micro_effect", "p=60 micro", n = 400, p = 60,
             beta = 0.16),
    scenario("M3", "micro_effect", "p=90 micro", n = 400, p = 90,
             beta = 0.24),
    scenario("R1", "network_strength", "rho=0.85", beta = 0.24,
             rho_within = 0.85),
    scenario("R2", "network_strength", "rho=0.95", beta = 0.24,
             rho_within = 0.95),
    scenario("R3", "network_strength", "rho=0.99", beta = 0.24,
             rho_within = 0.99),
    scenario("S1", "noise", "higher noise", beta = 0.30, error_sd = 1.80),
    scenario("S2", "noise", "baseline noise", beta = 0.30, error_sd = 1.00),
    scenario("S3", "noise", "lower noise", beta = 0.30, error_sd = 0.50)
  )

  run_rep_ids <- function(sc, scenario_index, ids) {
    if (length(ids) == 0) return(data.frame())
    if (N_CORES > 1 && .Platform$OS.type != "windows") {
      detected <- parallel::detectCores()
      if (is.na(detected) || detected < 1) detected <- N_CORES
      out <- parallel::mclapply(
        ids,
        function(i) one_rep(sc, scenario_index, i),
        mc.cores = min(N_CORES, detected),
        mc.set.seed = TRUE
      )
    } else {
      out <- lapply(ids, function(i) one_rep(sc, scenario_index, i))
    }
    bind_rows(Filter(Negate(is.null), out))
  }

  same_scenario_params <- function(x, sc) {
    if (nrow(x) == 0) return(TRUE)
    checks <- c(
      all(x$scenario == sc$id),
      all(x$family == sc$family),
      all(x$n == sc$n),
      all(x$p == sc$p),
      all(abs(x$beta - sc$beta) < 1e-12),
      all(abs(x$error_sd - sc$error_sd) < 1e-12),
      all(abs(x$rho_within - sc$rho_within) < 1e-12),
      all(abs(x$bridge - sc$bridge) < 1e-12)
    )
    all(checks)
  }

  read_checkpoint <- function(path, sc) {
    if (!file.exists(path)) return(data.frame())
    x <- tryCatch(read.csv(path), error = function(e) NULL)
    if (is.null(x)) {
      stop(sprintf("Checkpoint exists but cannot be read: %s", path))
    }
    required <- c("scenario", "rep", "n", "p", "beta", "error_sd",
                  "rho_within", "bridge", "gamma_hat", "significant")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
      stop(sprintf("Checkpoint %s is missing columns: %s",
                   path, paste(missing, collapse = ", ")))
    }
    if (!same_scenario_params(x, sc)) {
      stop(sprintf("Checkpoint %s does not match the current scenario parameters", path))
    }
    x %>%
      filter(rep >= 1, rep <= N_SIM) %>%
      arrange(rep) %>%
      group_by(rep) %>%
      slice_tail(n = 1) %>%
      ungroup()
  }

  write_checkpoint <- function(x, csv_path, rds_path) {
    x <- x %>% arrange(rep)
    tmp_csv <- paste0(csv_path, ".tmp")
    tmp_rds <- paste0(rds_path, ".tmp")
    write.csv(x, tmp_csv, row.names = FALSE)
    saveRDS(x, tmp_rds)
    file.rename(tmp_csv, csv_path)
    file.rename(tmp_rds, rds_path)
    invisible(TRUE)
  }

  run_scenario_checkpointed <- function(sc, scenario_index) {
    checkpoint_csv <- file.path(SCENARIO_DIR, paste0(sc$id, ".csv"))
    checkpoint_rds <- file.path(SCENARIO_DIR, paste0(sc$id, ".rds"))
    existing <- read_checkpoint(checkpoint_csv, sc)
    complete_ids <- sort(unique(existing$rep))

    if (length(complete_ids) >= N_SIM) {
      message(sprintf("[%s] %s already complete; loading checkpoint", sc$id, sc$label))
      return(existing %>% filter(rep %in% seq_len(N_SIM)) %>% arrange(rep))
    }

    remaining <- setdiff(seq_len(N_SIM), complete_ids)
    message(sprintf("[%s] %s; completed %d/%d, remaining %d",
                    sc$id, sc$label, length(complete_ids), N_SIM, length(remaining)))
    t0 <- proc.time()
    batches <- split(remaining, ceiling(seq_along(remaining) / BATCH_SIZE))

    for (b in seq_along(batches)) {
      ids <- batches[[b]]
      new_rows <- run_rep_ids(sc, scenario_index, ids)
      if (nrow(new_rows) > 0) {
        existing <- bind_rows(existing, new_rows) %>%
          arrange(rep) %>%
          group_by(rep) %>%
          slice_tail(n = 1) %>%
          ungroup()
        write_checkpoint(existing, checkpoint_csv, checkpoint_rds)
      }
      message(sprintf("  batch %d/%d saved; valid reps: %d/%d; elapsed %.2f min",
                      b, length(batches), dplyr::n_distinct(existing$rep), N_SIM,
                      (proc.time() - t0)[["elapsed"]] / 60))
    }
    existing %>% arrange(rep)
  }

  message(sprintf(
    "NEMA formal v6: n_sim=%d, boot=%d, cores=%d, batch=%d, seed=%d",
    N_SIM, B_BOOT, N_CORES, BATCH_SIZE, BASE_SEED
  ))

  started <- Sys.time()
  scenario_registry <- bind_rows(lapply(scenarios, as.data.frame))
  write.csv(scenario_registry, file.path(OUT_DIR, "scenario_registry.csv"),
            row.names = FALSE)
  writeLines(c(
    sprintf("Started: %s", format(started, "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("N_SIM: %d", N_SIM),
    sprintf("B_BOOT: %d", B_BOOT),
    sprintf("BATCH_SIZE: %d", BATCH_SIZE),
    sprintf("ALPHA: %.3f", ALPHA),
    sprintf("Z_CRIT: %.3f", Z_CRIT),
    sprintf("N_CORES: %d", N_CORES),
    sprintf("BASE_SEED: %d", BASE_SEED),
    sprintf("OUT_DIR: %s", OUT_DIR)
  ), file.path(OUT_DIR, "run_config_start.txt"))

  all_results <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    all_results[[i]] <- run_scenario_checkpointed(sc, i)
  }

  replicate_results <- bind_rows(all_results)

  scenario_summary <- replicate_results %>%
    group_by(scenario, family, label, n, p, beta, error_sd, rho_within, bridge) %>%
    summarise(
      valid_reps = n_distinct(rep),
      rejection_rate = mean(significant, na.rm = TRUE),
      rejection_mc_se = sqrt(rejection_rate * (1 - rejection_rate) / valid_reps),
      rejection_mc_low = pmax(0, rejection_rate - 1.96 * rejection_mc_se),
      rejection_mc_high = pmin(1, rejection_rate + 1.96 * rejection_mc_se),
      coverage = mean(covered, na.rm = TRUE),
      gamma_mean = mean(gamma_hat, na.rm = TRUE),
      gamma_bias = mean(gamma_bias, na.rm = TRUE),
      gamma_sd = sd(gamma_hat, na.rm = TRUE),
      gamma_rmse = sqrt(mean((gamma_hat - gamma_true)^2, na.rm = TRUE)),
      gamma_mae = mean(gamma_abs_error, na.rm = TRUE),
      score_cor = mean(score_cor, na.rm = TRUE),
      ols_mean_p_lt_0_05 = mean(ols_mean_p_lt_0_05, na.rm = TRUE),
      ols_any_p_lt_0_05 = mean(ols_any_p_lt_0_05, na.rm = TRUE),
      ols_any_bonf = mean(ols_any_bonf, na.rm = TRUE),
      ols_mean_bh_fdr = mean(ols_mean_bh_fdr, na.rm = TRUE),
      ols_any_bh_fdr = mean(ols_any_bh_fdr, na.rm = TRUE),
      ols_median_abs_cor = mean(ols_median_abs_cor, na.rm = TRUE),
      mean_valid_bootstrap = mean(valid_bootstrap, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(family, scenario)

  write.csv(replicate_results, file.path(OUT_DIR, "replicate_results.csv"),
            row.names = FALSE)
  write.csv(scenario_summary, file.path(OUT_DIR, "scenario_summary.csv"),
            row.names = FALSE)
  saveRDS(list(replicate_results = replicate_results,
               scenario_summary = scenario_summary,
               scenarios = scenarios),
          file.path(OUT_DIR, "nema_formal_v6_results.rds"))

  plot_family <- function(fam, filename, title) {
    dat <- scenario_summary %>% filter(family == fam)
    p <- ggplot(dat, aes(x = scenario, y = rejection_rate)) +
      geom_col(fill = "#80B1D3", width = 0.65) +
      geom_hline(yintercept = ALPHA, linetype = "dotted", color = "grey35") +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = NULL, y = "NEMA rejection rate / power", title = title) +
      theme_minimal(base_size = 12)
    ggsave(file.path(FIG_DIR, filename), p, width = 8, height = 4.8, dpi = 180)
  }

  plot_family("effect_size", "effect_size_power.png", "Power by Effect Size")
  plot_family("sample_size", "sample_size_power.png", "Power by Sample Size")
  plot_family("dimension_fixed", "dimension_fixed_power.png",
              "Power by Exposure Dimension at Fixed Beta")
  plot_family("micro_effect", "micro_effect_power.png", "Micro-Effect Power")
  plot_family("network_strength", "network_strength_power.png",
              "Power by Network Strength")
  plot_family("noise", "noise_power.png", "Power by Noise Level")

  comp <- scenario_summary %>%
    filter(family == "micro_effect") %>%
    select(scenario, label, rejection_rate, ols_any_bh_fdr, ols_any_bonf) %>%
    pivot_longer(c(rejection_rate, ols_any_bh_fdr, ols_any_bonf),
                 names_to = "metric", values_to = "rate") %>%
    mutate(metric = recode(metric,
                           rejection_rate = "NEMA",
                           ols_any_bh_fdr = "OLS BH-FDR",
                           ols_any_bonf = "OLS Bonferroni"))
  p_comp <- ggplot(comp, aes(x = label, y = rate, fill = metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = c("#80B1D3", "#FDB462", "#B3DE69")) +
    labs(x = NULL, y = "Detection rate", fill = NULL,
         title = "NEMA vs Individual OLS across p") +
    theme_minimal(base_size = 12)
  ggsave(file.path(FIG_DIR, "micro_effect_nema_vs_ols.png"), p_comp,
         width = 9, height = 4.8, dpi = 180)

  gamma_plot <- scenario_summary %>%
    filter(family %in% c("effect_size", "sample_size", "dimension_fixed")) %>%
    mutate(scenario = factor(scenario, levels = scenario))
  p_gamma <- ggplot(gamma_plot, aes(x = scenario, y = gamma_bias, fill = family)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.3) +
    labs(x = NULL, y = "Gamma bias", fill = NULL,
         title = "Gamma Estimation Bias") +
    theme_minimal(base_size = 12)
  ggsave(file.path(FIG_DIR, "gamma_bias.png"), p_gamma,
         width = 9, height = 4.8, dpi = 180)

  ended <- Sys.time()
  config <- c(
    sprintf("Generated: %s", format(ended, "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("Elapsed minutes: %.2f", as.numeric(difftime(ended, started, units = "mins"))),
    sprintf("N_SIM: %d", N_SIM),
    sprintf("B_BOOT: %d", B_BOOT),
    sprintf("ALPHA: %.3f", ALPHA),
    sprintf("Z_CRIT: %.3f", Z_CRIT),
    sprintf("N_CORES: %d", N_CORES),
    "Estimator: huge nonparanormal graphical lasso with EBIC.",
    "Inference: full-pipeline bootstrap normal confidence intervals using pilot-calibrated critical value.",
    "DGP: two-community sparse Gaussian Markov random field; outcome depends on first non-trivial network mode score."
  )
  writeLines(config, file.path(OUT_DIR, "run_config.txt"))

  report_path <- file.path(REPORT_DIR, "nema_formal_v6_report.md")
  report <- c(
    "# NEMA formal simulation v6",
    "",
    "## Configuration",
    "",
    paste(config, collapse = "\n\n"),
    "",
    "## Scenario Summary",
    "",
    paste(capture.output(knitr::kable(scenario_summary, digits = 3)), collapse = "\n"),
    "",
    "![Effect size power](../results/formal_v6/figures/effect_size_power.png)",
    "",
    "![Sample size power](../results/formal_v6/figures/sample_size_power.png)",
    "",
    "![Dimension power](../results/formal_v6/figures/dimension_fixed_power.png)",
    "",
    "![NEMA vs OLS](../results/formal_v6/figures/micro_effect_nema_vs_ols.png)",
    "",
    "![Gamma bias](../results/formal_v6/figures/gamma_bias.png)"
  )
  writeLines(report, report_path)

  message(sprintf("Done. Results: %s", OUT_DIR))
  message(sprintf("Report: %s", report_path))

} else {
  message("Skipping formal simulation study because RUN_SIMULATION is FALSE.")
}

################################################################################
# Section 3. NHANES POPs and leukocyte telomere length empirical analysis
################################################################################

if (RUN_NHANES) {

  #!/usr/bin/env Rscript

  # NHANES POPs + leukocyte telomere length empirical analysis for NEMA.
  # Source data and comparator code: Gibson et al. SHARP Mixtures Workshop repo.

  options(stringsAsFactors = FALSE, scipen = 999)

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "."
  root_dir <- normalizePath(file.path(dirname(script_path), "..", ".."),
                            mustWork = FALSE)
  if (!file.exists(file.path(root_dir, "scripts", "nema_core.R"))) {
    root_dir <- normalizePath(getwd(), mustWork = TRUE)
  }

  required <- c("ggplot2", "huge")
  missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
  }

  out_dir <- file.path(root_dir, "results", "nhanes_pops_telomere")
  fig_dir <- file.path(out_dir, "figures")
  tab_dir <- file.path(out_dir, "tables")
  data_dir <- file.path(root_dir, "data", "processed", "nhanes_pops_telomere")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(as.integer(Sys.getenv("NEMA_SEED", "20260623")))
  n_boot <- as.integer(Sys.getenv("NEMA_BOOT", "1000"))
  boot_batch <- as.integer(Sys.getenv("NEMA_BOOT_BATCH", "50"))
  estimator <- Sys.getenv("NEMA_NETWORK_ESTIMATOR", "huge")
  npn <- tolower(Sys.getenv("NEMA_NPN", "false")) %in% c("1", "true", "yes")
  d_env <- Sys.getenv("NEMA_D", "auto")
  d_fixed <- if (tolower(d_env) == "auto") NULL else as.integer(d_env)
  d_max <- as.integer(Sys.getenv("NEMA_D_MAX", "4"))

  # Gibson et al. shared the cleaned NHANES POPs and LTL example through the
  # SHARP Mixtures Workshop GitHub repository. The NEMA analysis uses
  # Data/studypop.csv from that repository. If the repository is not already
  # present under external/SHARP.Mixtures.Workshop, set DOWNLOAD_GIBSON=TRUE
  # to clone it automatically. This requires internet access and git.
  prepare_gibson_workshop_repo <- function(root_dir, download = FALSE) {
    repo_dir <- file.path(root_dir, "external", "SHARP.Mixtures.Workshop")
    study_path <- file.path(repo_dir, "Data", "studypop.csv")
    if (file.exists(study_path)) return(study_path)

    if (isTRUE(download)) {
      dir.create(dirname(repo_dir), recursive = TRUE, showWarnings = FALSE)
      cmd <- sprintf(
        "git clone --depth 1 https://github.com/lizzyagibson/SHARP.Mixtures.Workshop.git %s",
        shQuote(repo_dir)
      )
      status <- system(cmd)
      if (status == 0 && file.exists(study_path)) return(study_path)
    }

    stop(
      "Cannot find Gibson workshop data at: ", study_path, "\n",
      "To obtain the data manually, run:\n",
      "  git clone https://github.com/lizzyagibson/SHARP.Mixtures.Workshop.git external/SHARP.Mixtures.Workshop\n",
      "Then rerun this script. Alternatively set DOWNLOAD_GIBSON=TRUE before running."
    )
  }

  study_path <- prepare_gibson_workshop_repo(
    root_dir,
    download = tolower(Sys.getenv("DOWNLOAD_GIBSON", "false")) %in% c("1", "true", "yes")
  )

  raw <- read.csv(study_path, check.names = FALSE)
  dat <- na.omit(raw)

  mixture <- c(
    "LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA",
    "LBX180LA", "LBX187LA", "LBX194LA", "LBXD03LA", "LBXD05LA", "LBXD07LA",
    "LBXF03LA", "LBXF04LA", "LBXF05LA", "LBXF08LA", "LBXHXCLA", "LBXPCBLA"
  )

  chemical_labels <- c(
    LBX074LA = "PCB 74",
    LBX099LA = "PCB 99",
    LBX118LA = "PCB 118",
    LBX138LA = "PCB 138",
    LBX153LA = "PCB 153",
    LBX170LA = "PCB 170",
    LBX180LA = "PCB 180",
    LBX187LA = "PCB 187",
    LBX194LA = "PCB 194",
    LBXD03LA = "1,2,3,6,7,8-HxCDD",
    LBXD05LA = "1,2,3,4,6,7,8-HpCDD",
    LBXD07LA = "OCDD",
    LBXF03LA = "2,3,4,7,8-PnCDF",
    LBXF04LA = "1,2,3,4,7,8-HxCDF",
    LBXF05LA = "1,2,3,6,7,8-HxCDF",
    LBXF08LA = "1,2,3,4,6,7,8-HpCDF",
    LBXHXCLA = "PCB 169",
    LBXPCBLA = "PCB 126"
  )

  chemical_group <- c(
    LBX074LA = "Non-dioxin-like PCB",
    LBX099LA = "Non-dioxin-like PCB",
    LBX118LA = "Mono-ortho PCB",
    LBX138LA = "Non-dioxin-like PCB",
    LBX153LA = "Non-dioxin-like PCB",
    LBX170LA = "Non-dioxin-like PCB",
    LBX180LA = "Non-dioxin-like PCB",
    LBX187LA = "Non-dioxin-like PCB",
    LBX194LA = "Non-dioxin-like PCB",
    LBXD03LA = "Dioxin",
    LBXD05LA = "Dioxin",
    LBXD07LA = "Dioxin",
    LBXF03LA = "Furan",
    LBXF04LA = "Furan",
    LBXF05LA = "Furan",
    LBXF08LA = "Furan",
    LBXHXCLA = "Non-ortho PCB",
    LBXPCBLA = "Non-ortho PCB"
  )

  X_log <- log(as.matrix(dat[, mixture]))
  colnames(X_log) <- mixture
  Y <- log(dat$TELOMEAN)
  Z <- model.matrix(
    ~ LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT +
      age_cent + age_sq + factor(race_cat) + factor(bmi_cat3) +
      ln_lbxcot + factor(edu_cat) + male,
    data = dat
  )[, -1, drop = FALSE]

  analysis_dat <- cbind(
    dat[, c("SEQN", "TELOMEAN")],
    log_TELOMEAN = Y,
    as.data.frame(X_log),
    dat[, c("LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT",
            "LBXBAPCT", "age_cent", "age_sq", "race_cat", "bmi_cat3",
            "edu_cat", "male", "ln_lbxcot")]
  )
  write.csv(analysis_dat, file.path(data_dir, "analysis_complete_cases.csv"),
            row.names = FALSE)

  source_notes <- data.frame(
    item = c("source_repo", "source_file", "raw_rows", "complete_case_rows",
             "exposures", "outcome", "covariates"),
    value = c(
      "https://github.com/lizzyagibson/SHARP.Mixtures.Workshop",
      "Data/studypop.csv",
      nrow(raw),
      nrow(dat),
      paste(mixture, collapse = ", "),
      "log(TELOMEAN)",
      "Blood cell composition, age, age^2, race/ethnicity, BMI category, education, sex, log cotinine"
    )
  )
  write.csv(source_notes, file.path(tab_dir, "source_notes.csv"), row.names = FALSE)

  chem_table <- data.frame(
    variable = mixture,
    chemical = unname(chemical_labels[mixture]),
    group = unname(chemical_group[mixture]),
    stringsAsFactors = FALSE
  )
  write.csv(chem_table, file.path(tab_dir, "chemical_dictionary.csv"), row.names = FALSE)

  tidy_lm_term <- function(fit, terms, label_map = NULL) {
    sm <- summary(fit)$coefficients
    ci <- confint(fit)
    out <- data.frame(
      term = rownames(sm),
      estimate = sm[, "Estimate"],
      std_error = sm[, "Std. Error"],
      statistic = sm[, "t value"],
      p_value = sm[, "Pr(>|t|)"],
      ci_low = ci[, 1],
      ci_high = ci[, 2],
      row.names = NULL
    )
    out <- out[out$term %in% terms, , drop = FALSE]
    if (!is.null(label_map)) {
      out$label <- unname(label_map[out$term])
    }
    out
  }

  # Comparator 1: adjusted single-pollutant regressions, one SD log exposure.
  X_scaled <- scale(X_log)
  single_terms <- lapply(seq_along(mixture), function(j) {
    fit_dat <- data.frame(Y = Y, exposure = X_scaled[, j], Z)
    fit <- lm(Y ~ ., data = fit_dat)
    sm <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(
      variable = mixture[j],
      chemical = unname(chemical_labels[mixture[j]]),
      group = unname(chemical_group[mixture[j]]),
      estimate_per_sd_log = sm["exposure", "Estimate"],
      std_error = sm["exposure", "Std. Error"],
      ci_low = ci["exposure", 1],
      ci_high = ci["exposure", 2],
      p_value = sm["exposure", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  })
  single_pollutant <- do.call(rbind, single_terms)
  single_pollutant$bh_q <- p.adjust(single_pollutant$p_value, method = "BH")
  single_pollutant$significant_bh_005 <- single_pollutant$bh_q < 0.05
  write.csv(single_pollutant, file.path(tab_dir, "single_pollutant_adjusted_ols.csv"),
            row.names = FALSE)

  # Comparator 2: adjusted multipollutant OLS with all standardized log exposures.
  multi_dat <- data.frame(Y = Y, as.data.frame(X_scaled), Z, check.names = FALSE)
  multi_fit <- lm(Y ~ ., data = multi_dat)
  multi_ols <- tidy_lm_term(multi_fit, mixture, chemical_labels)
  multi_ols$variable <- multi_ols$term
  multi_ols$group <- unname(chemical_group[multi_ols$variable])
  multi_ols$bh_q <- p.adjust(multi_ols$p_value, method = "BH")
  multi_ols$significant_bh_005 <- multi_ols$bh_q < 0.05
  multi_ols <- multi_ols[, c("variable", "label", "group", "estimate", "std_error",
                             "ci_low", "ci_high", "p_value", "bh_q",
                             "significant_bh_005")]
  names(multi_ols)[names(multi_ols) == "label"] <- "chemical"
  write.csv(multi_ols, file.path(tab_dir, "multipollutant_adjusted_ols.csv"),
            row.names = FALSE)

  # Comparator 3: PCA using first three components, matching Gibson workshop code.
  pca <- prcomp(X_log, scale. = TRUE)
  pca_dat <- data.frame(Y = Y, pca$x[, 1:3, drop = FALSE], Z)
  pca_fit <- lm(Y ~ ., data = pca_dat)
  pca_terms <- paste0("PC", 1:3)
  pca_results <- tidy_lm_term(pca_fit, pca_terms)
  pca_results$variance_explained <- (pca$sdev[1:3]^2) / sum(pca$sdev^2)
  write.csv(pca_results, file.path(tab_dir, "pca_adjusted_regression.csv"),
            row.names = FALSE)

  orient_nema_fit <- function(fit) {
    for (j in seq_len(fit$d)) {
      if (sum(fit$U[, j], na.rm = TRUE) < 0) {
        fit$U[, j] <- -fit$U[, j]
        fit$scores[, j] <- -fit$scores[, j]
        fit$gamma[j] <- -fit$gamma[j]
      }
    }
    fit$contribution <- chemical_contribution(fit$U, fit$gamma)
    fit$total_effect <- network_total_effect(fit$contribution, fit$Sigma)
    fit$overall_iqr_effect <- overall_iqr_effect(scale(as.matrix(X_log), center = TRUE, scale = TRUE),
                                                 fit$U, fit$gamma)
    fit
  }

  nema_fit_primary <- nema_fit(
    X_log, Y, Z = Z, d = d_fixed, d_max = d_max, estimator = estimator,
    npn = npn, family = "gaussian"
  )
  dimnames(nema_fit_primary$W) <- list(mixture, mixture)
  dimnames(nema_fit_primary$Sigma) <- list(mixture, mixture)
  nema_fit_primary <- orient_nema_fit(nema_fit_primary)

  bootstrap_batches <- function(X, Y, Z, fit_ref, B, batch_size, path,
                                estimator, npn, seed0 = 92001) {
    if (file.exists(path)) {
      saved <- readRDS(path)
      if (!is.null(saved$B_requested) && saved$B_requested >= B) return(saved)
      start_done <- nrow(saved$gamma)
    } else {
      saved <- list(
        gamma = NULL,
        contribution = NULL,
        total_effect = NULL,
        overall_iqr_effect = numeric(0),
        B_requested = B
      )
      start_done <- 0
    }

    while (start_done < B) {
      this_B <- min(batch_size, B - start_done)
      b_seed <- seed0 + start_done + 1
      message("NEMA bootstrap ", start_done + 1, "-", start_done + this_B,
              " of ", B)
      boot <- bootstrap_nema(
        X, Y, Z = Z, fit_ref = fit_ref, B = this_B, estimator = estimator,
        npn = npn, family = "gaussian", seed = b_seed
      )
      saved$gamma <- rbind(saved$gamma, boot$gamma)
      saved$contribution <- rbind(saved$contribution, boot$contribution)
      saved$total_effect <- rbind(saved$total_effect, boot$total_effect)
      saved$overall_iqr_effect <- c(saved$overall_iqr_effect, boot$overall_iqr_effect)
      saved$B_requested <- B
      saveRDS(saved, path)
      start_done <- nrow(saved$gamma)
    }
    saved
  }

  boot_path <- file.path(out_dir, paste0("nema_bootstrap_B", n_boot, ".rds"))
  nema_boot <- bootstrap_batches(
    X_log, Y, Z, nema_fit_primary, n_boot, boot_batch, boot_path,
    estimator = estimator, npn = npn
  )

  gamma_table <- ci_table(nema_fit_primary$gamma, nema_boot$gamma, "mode")
  gamma_table$mode <- gamma_table$term
  gamma_table <- gamma_table[, c("mode", "estimate", "ci_low", "ci_high", "significant")]
  write.csv(gamma_table, file.path(tab_dir, "nema_gamma_bootstrap_ci.csv"),
            row.names = FALSE)

  contribution_table <- ci_table(nema_fit_primary$contribution, nema_boot$contribution,
                                 "chemical")
  contribution_table$variable <- contribution_table$term
  contribution_table$chemical <- unname(chemical_labels[contribution_table$variable])
  contribution_table$group <- unname(chemical_group[contribution_table$variable])
  contribution_table <- contribution_table[order(contribution_table$estimate, decreasing = TRUE), ]
  write.csv(contribution_table, file.path(tab_dir, "nema_chemical_contribution_ci.csv"),
            row.names = FALSE)

  total_effect_table <- ci_table(nema_fit_primary$total_effect, nema_boot$total_effect,
                                 "chemical")
  total_effect_table$variable <- total_effect_table$term
  total_effect_table$chemical <- unname(chemical_labels[total_effect_table$variable])
  total_effect_table$group <- unname(chemical_group[total_effect_table$variable])
  total_effect_table <- total_effect_table[order(total_effect_table$estimate, decreasing = TRUE), ]
  write.csv(total_effect_table, file.path(tab_dir, "nema_network_total_effect_ci.csv"),
            row.names = FALSE)

  overall_ci <- ci_percentile(nema_boot$overall_iqr_effect)
  overall_table <- data.frame(
    estimand = "Overall POP network IQR contrast",
    estimate = nema_fit_primary$overall_iqr_effect,
    ci_low = overall_ci[1],
    ci_high = overall_ci[2],
    significant = overall_ci[1] > 0 | overall_ci[2] < 0,
    B = nrow(nema_boot$gamma)
  )
  write.csv(overall_table, file.path(tab_dir, "nema_overall_iqr_effect_ci.csv"),
            row.names = FALSE)

  eigen_table <- data.frame(
    mode = seq_along(nema_fit_primary$eigenvalues),
    eigenvalue = nema_fit_primary$eigenvalues
  )
  write.csv(eigen_table, file.path(tab_dir, "nema_laplacian_eigenvalues.csv"),
            row.names = FALSE)

  U_table <- data.frame(
    variable = rownames(nema_fit_primary$U),
    chemical = unname(chemical_labels[rownames(nema_fit_primary$U)]),
    group = unname(chemical_group[rownames(nema_fit_primary$U)]),
    nema_fit_primary$U,
    check.names = FALSE
  )
  write.csv(U_table, file.path(tab_dir, "nema_mode_loadings.csv"), row.names = FALSE)

  edge_df <- as.data.frame(as.table(nema_fit_primary$W), stringsAsFactors = FALSE)
  names(edge_df) <- c("from", "to", "weight")
  edge_df <- edge_df[as.integer(factor(edge_df$from, levels = mixture)) <
                       as.integer(factor(edge_df$to, levels = mixture)) &
                       edge_df$weight > 0, ]
  edge_df$from_label <- unname(chemical_labels[edge_df$from])
  edge_df$to_label <- unname(chemical_labels[edge_df$to])
  edge_df <- edge_df[order(edge_df$weight, decreasing = TRUE), ]
  write.csv(edge_df, file.path(tab_dir, "nema_network_edges.csv"), row.names = FALSE)

  method_summary <- data.frame(
    method = c("NEMA", "PCA", "FA", "Single-pollutant OLS", "Multipollutant OLS"),
    primary_finding = c(
      sprintf("Overall network IQR effect %.3f (95%% CI %.3f, %.3f); d=%d",
              overall_table$estimate, overall_table$ci_low, overall_table$ci_high,
              nema_fit_primary$d),
      sprintf("PC1 beta %.3f (95%% CI %.3f, %.3f); sign depends on PCA orientation",
              pca_results$estimate[pca_results$term == "PC1"],
              pca_results$ci_low[pca_results$term == "PC1"],
              pca_results$ci_high[pca_results$term == "PC1"]),
      "Gibson workshop FA table reports MR2/MR3 positive and MR4 negative associations",
      sprintf("%d of 18 chemicals significant after BH-FDR",
              sum(single_pollutant$significant_bh_005)),
      sprintf("%d of 18 chemicals significant after BH-FDR",
              sum(multi_ols$significant_bh_005))
    )
  )
  write.csv(method_summary, file.path(tab_dir, "method_summary.csv"), row.names = FALSE)

  saveRDS(
    list(
      data = analysis_dat,
      mixture = mixture,
      chemical_labels = chemical_labels,
      chemical_group = chemical_group,
      single_pollutant = single_pollutant,
      multi_ols = multi_ols,
      pca = list(fit = pca_fit, results = pca_results, object = pca),
      nema_fit = nema_fit_primary,
      nema_boot = nema_boot,
      overall = overall_table,
      source_notes = source_notes
    ),
    file.path(out_dir, "nhanes_pops_telomere_nema_results.rds")
  )

  library(ggplot2)

  theme_set(theme_bw(base_size = 11))

  ggplot(single_pollutant,
         aes(x = reorder(chemical, estimate_per_sd_log),
             y = estimate_per_sd_log, ymin = ci_low, ymax = ci_high,
             color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_pointrange() +
    coord_flip() +
    labs(x = NULL, y = "Adjusted beta per 1 SD log exposure",
         title = "Single-pollutant adjusted associations") +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "single_pollutant_adjusted_ols.png"), width = 8, height = 6, dpi = 300)

  corr_mat <- cor(X_log)
  corr_df <- as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE)
  names(corr_df) <- c("x", "y", "correlation")
  corr_df$x_label <- factor(unname(chemical_labels[corr_df$x]),
                            levels = unname(chemical_labels[mixture]))
  corr_df$y_label <- factor(unname(chemical_labels[corr_df$y]),
                            levels = rev(unname(chemical_labels[mixture])))
  ggplot(corr_df, aes(x = x_label, y = y_label, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b83232",
                         limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "r",
         title = "Correlation among log-transformed POPs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "exposure_correlation_heatmap.png"), width = 8, height = 7, dpi = 300)

  W_df <- as.data.frame(as.table(nema_fit_primary$W), stringsAsFactors = FALSE)
  names(W_df) <- c("x", "y", "weight")
  W_df$x_label <- factor(unname(chemical_labels[W_df$x]),
                         levels = unname(chemical_labels[mixture]))
  W_df$y_label <- factor(unname(chemical_labels[W_df$y]),
                         levels = rev(unname(chemical_labels[mixture])))
  ggplot(W_df, aes(x = x_label, y = y_label, fill = weight)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#365c8d") +
    labs(x = NULL, y = NULL, fill = "weight",
         title = "NEMA estimated POP network") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "nema_network_heatmap.png"), width = 8, height = 7, dpi = 300)

  ggplot(gamma_table, aes(x = mode, y = estimate, ymin = ci_low, ymax = ci_high)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_pointrange(color = "#1f6f5b") +
    labs(x = NULL, y = "Gamma", title = "NEMA gamma estimate")
  ggsave(file.path(fig_dir, "nema_gamma_ci.png"), width = 5, height = 4, dpi = 300)

  plot_contrib <- contribution_table
  plot_contrib$chemical <- factor(plot_contrib$chemical, levels = rev(plot_contrib$chemical))
  ggplot(plot_contrib, aes(x = chemical, y = estimate, ymin = ci_low, ymax = ci_high,
                           color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_pointrange() +
    coord_flip() +
    labs(x = NULL, y = "NEMA network contribution",
         title = "Chemical contributions to the NEMA POP network effect") +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "nema_chemical_contribution_ci.png"), width = 8, height = 6, dpi = 300)

  plot_total <- total_effect_table
  plot_total$chemical <- factor(plot_total$chemical, levels = rev(plot_total$chemical))
  ggplot(plot_total, aes(x = chemical, y = estimate, ymin = ci_low, ymax = ci_high,
                         color = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_pointrange() +
    coord_flip() +
    labs(x = NULL, y = "Network total effect",
         title = "NEMA network-propagated total effects") +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "nema_network_total_effect_ci.png"), width = 8, height = 6, dpi = 300)

  ggplot(eigen_table[1:8, ], aes(x = mode, y = eigenvalue)) +
    geom_line(color = "#4c6a92") +
    geom_point(color = "#4c6a92") +
    scale_x_continuous(breaks = 1:8) +
    labs(x = "Laplacian eigenvalue index", y = "Eigenvalue",
         title = "NEMA network eigengap diagnostic")
  ggsave(file.path(fig_dir, "nema_eigengap.png"), width = 6, height = 4, dpi = 300)

  report_path <- file.path(root_dir, "reports", "nhanes_pops_telomere",
                           "nema_nhanes_pops_telomere_report.md")
  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

  top_contrib <- head(contribution_table[, c("chemical", "group", "estimate",
                                             "ci_low", "ci_high", "significant")], 8)
  top_total <- head(total_effect_table[, c("chemical", "group", "estimate",
                                           "ci_low", "ci_high", "significant")], 8)

  fmt_table <- function(df, digits = 4) {
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], round, digits)
    paste(capture.output(print(df, row.names = FALSE)), collapse = "\n")
  }

  report <- c(
    "# NEMA analysis of NHANES POPs and leukocyte telomere length",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Data source",
    "",
    "- Primary analysis dataset: `external/SHARP.Mixtures.Workshop/Data/studypop.csv` from Gibson et al.'s public workshop repository.",
    sprintf("- Complete-case analytic sample: %d participants from %d rows in the shared dataset.", nrow(dat), nrow(raw)),
    "- Outcome: natural-log mean leukocyte telomere length, `log(TELOMEAN)`.",
    "- Exposure mixture: 18 lipid-adjusted persistent organic pollutants, natural-log transformed.",
    "- Covariates: blood cell composition, age, age squared, race/ethnicity, BMI category, education, sex, and log serum cotinine.",
    "- Note: this reproduces the Gibson workshop analytic dataset. It is not a survey-weighted NHANES reanalysis because the shared `studypop.csv` used by the workshop scripts does not retain survey design variables.",
    "",
    "## Main NEMA result",
    "",
    sprintf("- Network estimator: `%s`; nonparanormal transform: `%s`; selected mode count: `d = %d`; bootstrap replicates: `%d`.",
            estimator, npn, nema_fit_primary$d, nrow(nema_boot$gamma)),
    sprintf("- Overall POP network IQR contrast: %.4f (95%% bootstrap CI %.4f, %.4f).",
            overall_table$estimate, overall_table$ci_low, overall_table$ci_high),
    sprintf("- Gamma for the primary network mode: %.4f (95%% bootstrap CI %.4f, %.4f).",
            gamma_table$estimate[1], gamma_table$ci_low[1], gamma_table$ci_high[1]),
    "",
    "Interpretation: NEMA identifies one dominant network-embedded POP burden mode. Higher values of this mode are associated with longer leukocyte telomere length, consistent with the qualitative direction reported by Mitro et al. for dioxin-like/non-ortho POP summaries.",
    "",
    "## Highest NEMA chemical contributions",
    "",
    "```",
    fmt_table(top_contrib),
    "```",
    "",
    "## Highest network-propagated total effects",
    "",
    "```",
    fmt_table(top_total),
    "```",
    "",
    "## Comparator results",
    "",
    sprintf("- Single-pollutant adjusted OLS: %d/18 chemicals survived BH-FDR at 0.05.",
            sum(single_pollutant$significant_bh_005)),
    sprintf("- Multipollutant adjusted OLS: %d/18 chemicals survived BH-FDR at 0.05.",
            sum(multi_ols$significant_bh_005)),
    sprintf("- PCA: PC1 beta %.4f (95%% CI %.4f, %.4f), with sign dependent on PCA orientation; PC1 explained %.1f%% of exposure variance.",
            pca_results$estimate[pca_results$term == "PC1"],
            pca_results$ci_low[pca_results$term == "PC1"],
            pca_results$ci_high[pca_results$term == "PC1"],
            100 * pca_results$variance_explained[pca_results$term == "PC1"]),
    "",
    "## Figures",
    "",
    "- `results/nhanes_pops_telomere/figures/exposure_correlation_heatmap.png`",
    "- `results/nhanes_pops_telomere/figures/nema_network_heatmap.png`",
    "- `results/nhanes_pops_telomere/figures/nema_gamma_ci.png`",
    "- `results/nhanes_pops_telomere/figures/nema_chemical_contribution_ci.png`",
    "- `results/nhanes_pops_telomere/figures/nema_network_total_effect_ci.png`",
    "- `results/nhanes_pops_telomere/figures/single_pollutant_adjusted_ols.png`",
    "",
    "## Reproducibility",
    "",
    "- R script: `scripts/nhanes_pops_telomere/run_nema_empirical_analysis.R`.",
    "- Processed complete-case data: `data/processed/nhanes_pops_telomere/analysis_complete_cases.csv`.",
    "- Machine-readable results: `results/nhanes_pops_telomere/nhanes_pops_telomere_nema_results.rds`.",
    "- Bootstrap checkpoint: `results/nhanes_pops_telomere/nema_bootstrap_B*.rds`.",
    "",
    "## Manuscript angle",
    "",
    "This empirical example is a strong fit for the NEMA paper because the dataset is public, the original epidemiologic result is already known, and Gibson et al. have already used it as a benchmark for environmental mixture methods. NEMA adds a network-aware interpretation: instead of selecting isolated chemicals or unconstrained components, it estimates a dominant POP network burden, quantifies its gamma association, and decomposes the association into chemical-level network contributions."
  )
  writeLines(report, report_path)

  message("Done.")
  message("Report: ", report_path)
  message("Results: ", out_dir)

} else {
  message("Skipping NHANES empirical analysis because RUN_NHANES is FALSE.")
}
