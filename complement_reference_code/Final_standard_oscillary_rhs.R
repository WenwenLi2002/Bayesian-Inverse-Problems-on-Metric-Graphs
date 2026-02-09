suppressPackageStartupMessages({
  library(rSPDE)
  library(MetricGraph)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------
# A) Output helpers
# ------------------------------------------------
make_outdir <- function(prefix = "inverse_mcmc_metricgraph") {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(getwd(), paste0(prefix, "_", stamp))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir
}

safe_save_rds <- function(obj, path) {
  ok <- tryCatch({
    saveRDS(obj, path)
    TRUE
  }, error = function(e) {
    message("Failed saveRDS: ", path, " | ", e$message)
    FALSE
  })
  invisible(ok)
}

safe_write_csv <- function(df, path) {
  ok <- tryCatch({
    write.csv(df, path, row.names = FALSE)
    TRUE
  }, error = function(e) {
    message("Failed write.csv: ", path, " | ", e$message)
    FALSE
  })
  invisible(ok)
}

# ------------------------------------------------
# Prior/SPDE hyperparameters
# ------------------------------------------------
range_prior <- 3
sigma_prior <- 0.05

alpha_base <- 1
alpha_spde <- 2 * alpha_base

# ------------------------------------------------
# Truth filtering: regenerate u if p-range too large
# ------------------------------------------------
P_RANGE_CAP <- 500           
MAX_TRUTH_TRIES <- 50L       

# ------------------------------------------------
# 1) Build metric graph + mesh
# ------------------------------------------------
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance = list(edge_edge = 1e-3, vertex_vertex = 1e-3, edge_vertex = 1e-3)
)
graph$build_mesh(h = 0.05)

meshV <- as.matrix(graph$mesh$V)
if (ncol(meshV) < 2) stop("Mesh node matrix has fewer than 2 columns!")
colnames(meshV) <- c("x", "y")
graph$mesh$V <- meshV

x <- meshV[, "x"]
y <- meshV[, "y"]

# ------------------------------------------------
# 2) FEM matrices and RHS: trig forcing on the embedded coordinates
#     fbar(z) = sin(2*pi*k*x01) * cos(2*pi*k*y01), where x01,y01 in [0,1]
#     then remove C-mean and L2(C)-normalize
# ------------------------------------------------
graph$compute_fem()
G <- graph$mesh$G
C <- graph$mesh$C

# ---- rescale embedding coordinates to [0,1] to control effective frequency
x01 <- (x - min(x)) / max(diff(range(x)), 1e-12)
y01 <- (y - min(y)) / max(diff(range(y)), 1e-12)

# ---- frequency parameter 
k_freq <- 2L

# ---- trig RHS (smooth, global, sign-changing)
fbar <- 10* sin(2*pi*k_freq*x01) * cos(2*pi*k_freq*y01)

# ---- remove C-mean (mass-weighted) to suppress near-constant mode
c_mean <- as.numeric(sum(C %*% fbar) / sum(C))
fbar <- fbar - c_mean

# ---- L2(C) normalize so forcing energy is comparable across RHS choices
norm_fbar <- sqrt(as.numeric(crossprod(fbar, C %*% fbar)))
if (is.finite(norm_fbar) && norm_fbar > 0) fbar <- fbar / norm_fbar

# ---- FEM load vector
f <- as.vector(C %*% fbar)

cat(sprintf(
  "Trig RHS: fbar=sin(2*pi*k*x01)*cos(2*pi*k*y01), k=%d | C-mean removed + L2(C) normalized\n",
  k_freq
))

# ------------------------------------------------
# 4) PDE solve 
# ------------------------------------------------
kappa <- 1

solve_pde <- function(u, jitter0 = 1e-10, max_jitter_tries = 6L) {
  n_mesh <- nrow(meshV)
  u <- as.numeric(u)
  if (length(u) != n_mesh || any(!is.finite(u))) return(rep(NA_real_, n_mesh))
  
  Du <- Matrix::Diagonal(n_mesh, exp(u))
  
  # robust symmetric assembly
  GD <- G %*% Du
  Gu <- (GD + Matrix::t(GD)) / 2
  
  L <- kappa^2 * C + Gu
  L <- Matrix::forceSymmetric(L)
  
  # scale-aware jitter
  diagL <- Matrix::diag(L)
  diag_scale <- mean(abs(diagL))
  if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
  
  for (k in 0:(max_jitter_tries - 1L)) {
    jit <- jitter0 * (10^k) * diag_scale
    Lk <- L + Matrix::Diagonal(n_mesh, jit)
    
    sol <- tryCatch({
      chol <- Matrix::Cholesky(Lk, LDL = FALSE, perm = TRUE)
      Matrix::solve(chol, f)
    }, error = function(e) NULL)
    
    if (!is.null(sol)) {
      sol <- as.numeric(sol)
      if (all(is.finite(sol))) return(sol)
    }
  }
  
  # last resort (may fail)
  sol <- tryCatch(as.numeric(Matrix::solve(L, f)), error = function(e) rep(NA_real_, n_mesh))
  sol
}

# ------------------------------------------------
# Helpers: pointwise mode, errors
# ------------------------------------------------
get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

rmse <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  sqrt(mean((est[ok] - truth[ok])^2))
}

# relative L2 error (Euclidean): sqrt(sum e^2 / sum truth^2)
rel_l2 <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  num <- sum((est[ok] - truth[ok])^2)
  den <- sum((truth[ok])^2)
  sqrt(num / max(den, 1e-12))
}

# range-normalized RMSE
rel_range <- function(est, truth, eps = 1e-12) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  rng <- diff(range(truth[ok], na.rm = TRUE))
  rmse(est, truth) / max(rng, eps)
}

# relative L2(Gamma) error using mass matrix C: ||v||^2 ≈ v^T C v
rel_l2_gamma <- function(est, truth, C, eps = 1e-12) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  e <- est[ok] - truth[ok]
  Ct <- C[ok, ok, drop = FALSE]
  num <- as.numeric(crossprod(e, Ct %*% e))
  den <- as.numeric(crossprod(truth[ok], Ct %*% truth[ok]))
  sqrt(num / max(den, eps))
}

# normalized RMSE for p
nrmse_range <- function(est, truth, eps = 1e-12) {
  rng <- diff(range(truth, na.rm = TRUE))
  sqrt(mean((est - truth)^2, na.rm = TRUE)) / max(rng, eps)
}
nrmse_sd <- function(est, truth, eps = 1e-12) {
  s <- sd(truth, na.rm = TRUE)
  sqrt(mean((est - truth)^2, na.rm = TRUE)) / max(s, eps)
}

# ------------------------------------------------
# 5) Fix observation locations ONCE (same obs_indices for all trials)
# ------------------------------------------------
set.seed(123)
obs_indices <- sample(1:nrow(meshV), 1998)

alpha_noise <- 0.02
beta_noise  <- 0.1

# ------------------------------------------------
# 6) pCN MCMC (ROBUST beta + stable adaptation)
# ------------------------------------------------
pCN_MCMC <- function(n_iter, beta0, u_init, neg_log_likelihood_fn,
                     T0 = 1, cooling_factor = 0.95, adapt_interval = 500,
                     show_iter_pb = FALSE, verbose = FALSE) {
  u_current <- as.numeric(u_init)
  n_dim <- length(u_current)
  samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
  acceptance <- numeric(n_iter)
  
  beta_min <- 0.01
  beta_max <- 0.99
  beta <- min(max(beta0, beta_min), beta_max)
  
  # stable adaptation hyperparameters
  target_acc <- 0.30
  kappa_adapt <- 0.10
  
  if (show_iter_pb) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
    pb_interval <- max(1, floor(n_iter / 1000))
  }
  
  # cache current potential
  phi_current <- neg_log_likelihood_fn(u_current, T = max(1, T0))
  if (!is.finite(phi_current)) phi_current <- Inf
  
  for (i in 1:n_iter) {
    T_current <- max(1, T0 * cooling_factor^(i / adapt_interval))
    
    # draw prior sample
    xi <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_spde,
                                 graph = graph, type = "mesh"))
    if (length(xi) != n_dim) xi <- xi[1:n_dim]
    if (any(!is.finite(xi))) xi[!is.finite(xi)] <- 0
    
    # clamp beta and propose
    beta <- min(max(beta, beta_min), beta_max)
    u_proposed <- sqrt(pmax(1 - beta^2, 0)) * u_current + beta * xi
    
    phi_proposed <- neg_log_likelihood_fn(u_proposed, T = T_current)
    if (!is.finite(phi_proposed)) phi_proposed <- Inf
    
    log_alpha <- -(phi_proposed - phi_current)
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      u_current <- u_proposed
      phi_current <- phi_proposed
      acceptance[i] <- 1
    } else {
      acceptance[i] <- 0
    }
    
    samples[i, ] <- u_current
    
    # adapt beta 
    if (i %% adapt_interval == 0) {
      idx <- (i - adapt_interval + 1):i
      current_accept_rate <- mean(acceptance[idx])
      beta <- beta * exp(kappa_adapt * (current_accept_rate - target_acc))
      beta <- min(max(beta, beta_min), beta_max)
      
      if (verbose) {
        cat("Iter:", i,
            "| acc:", sprintf("%.3f", current_accept_rate),
            "| beta:", sprintf("%.3f", beta),
            "| T:", sprintf("%.2f", T_current), "\n")
      }
    }
    
    if (show_iter_pb && (i %% pb_interval == 0 || i == n_iter)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  if (show_iter_pb) close(pb)
  
  list(samples = samples, acceptance = acceptance, beta_final = beta)
}

# ------------------------------------------------
# 7) SAVE bundles
# ------------------------------------------------
OUTDIR <- make_outdir(prefix = "inverse_mcmc_metricgraph")
cat("Saving everything under OUTDIR = ", OUTDIR, "\n")

save_global_bundle <- function(OUTDIR, n_iter, burn_in, beta0, gamma_smooth) {
  global_dir <- file.path(OUTDIR, "global")
  dir.create(global_dir, showWarnings = FALSE)
  
  global_bundle <- list(
    timestamp = Sys.time(),
    meshV = meshV,
    E = graph$mesh$E,
    C = C,
    G = G,
    kappa = kappa,
    src_sigma = src_sigma,
    source_center_idx = idx0,
    source_center_xy = c(x0 = x0, y0 = y0),
    # RHS (trig)
    fbar = fbar,
    f = f,
    rhs_type = "trig_sin_cos",
    k_freq = as.integer(k_freq),
    range_prior = range_prior,
    sigma_prior = sigma_prior,
    alpha_base = alpha_base,
    alpha_spde = alpha_spde,
    P_RANGE_CAP = P_RANGE_CAP,
    MAX_TRUTH_TRIES = MAX_TRUTH_TRIES,
    obs_indices = as.integer(obs_indices),
    alpha_noise = alpha_noise,
    beta_noise = beta_noise,
    n_iter = n_iter,
    burn_in = burn_in,
    beta0 = beta0,
    gamma_smooth = gamma_smooth,
    sessionInfo = sessionInfo()
  )
  
  safe_save_rds(global_bundle, file.path(global_dir, "global_bundle.rds"))
  capture.output(str(global_bundle), file = file.path(global_dir, "global_bundle_str.txt"))
}

save_one_trial <- function(OUTDIR, trial_id,
                           u_true, p_true,
                           obs_indices, observations,
                           u_post_mean, u_map_raw,
                           p_post_mean, p_map_raw,
                           a_true = NULL, a_post_mean = NULL, a_map_raw = NULL,
                           acc_rate_all, beta_final,
                           acceptance_vec = NULL,
                           sd_u_mean = NA_real_, sd_u_max = NA_real_,
                           posterior_samples_post_burnin = NULL,
                           save_thinned_samples = FALSE,
                           thin_every = 50L,
                           trial_settings = NULL) {
  
  trial_dir <- file.path(OUTDIR, sprintf("trial_%03d", trial_id))
  dir.create(trial_dir, showWarnings = FALSE)
  
  trial_bundle <- list(
    trial = trial_id,
    timestamp = Sys.time(),
    trial_settings = trial_settings,
    u_true = as.numeric(u_true),
    p_true = as.numeric(p_true),
    a_true = if (!is.null(a_true)) as.numeric(a_true) else NULL,
    obs_indices = as.integer(obs_indices),
    observations = as.numeric(observations),
    u_post_mean = as.numeric(u_post_mean),
    u_map_raw   = as.numeric(u_map_raw),
    p_post_mean = as.numeric(p_post_mean),
    p_map_raw   = as.numeric(p_map_raw),
    a_post_mean = if (!is.null(a_post_mean)) as.numeric(a_post_mean) else NULL,
    a_map_raw   = if (!is.null(a_map_raw)) as.numeric(a_map_raw) else NULL,
    acc_rate_all = as.numeric(acc_rate_all),
    beta_final   = as.numeric(beta_final),
    sd_u_mean = as.numeric(sd_u_mean),
    sd_u_max  = as.numeric(sd_u_max)
  )
  
  if (!is.null(acceptance_vec)) {
    trial_bundle$acceptance_vec <- as.numeric(acceptance_vec)
  }
  
  if (isTRUE(save_thinned_samples) && !is.null(posterior_samples_post_burnin)) {
    idx <- seq(1, nrow(posterior_samples_post_burnin), by = as.integer(thin_every))
    trial_bundle$posterior_samples_thin <- posterior_samples_post_burnin[idx, , drop = FALSE]
    trial_bundle$thin_every <- as.integer(thin_every)
  }
  
  safe_save_rds(trial_bundle, file.path(trial_dir, "trial_bundle.rds"))
  
  scale_stats <- data.frame(
    trial = trial_id,
    u_range = diff(range(u_true, na.rm = TRUE)),
    u_sd    = sd(u_true, na.rm = TRUE),
    p_range = diff(range(p_true, na.rm = TRUE)),
    p_sd    = sd(p_true, na.rm = TRUE)
  )
  safe_write_csv(scale_stats, file.path(trial_dir, "scale_stats.csv"))
}

# ------------------------------------------------
# 8) One-trial pipeline
# ------------------------------------------------
run_one_trial <- function(trial_id, OUTDIR,
                          n_iter = 60000,
                          burn_in = 30000,
                          beta0 = 0.3,
                          gamma_smooth = 0.1,
                          seed_truth = 1000,
                          seed_noise = 2000,
                          seed_init  = 3000,
                          seed_mcmc  = 4000,
                          inner_mcmc_pb = FALSE,
                          verbose_mcmc = FALSE,
                          save_thinned_samples = FALSE,
                          thin_every = 50L) {
  
  # (i) truth with filtering on range(p_true)
  truth_attempt <- 0L
  repeat {
    seed_truth_used <- seed_truth + trial_id + 100000L * truth_attempt
    set.seed(seed_truth_used)
    
    u_true <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_spde,
                                     graph = graph, type = "mesh"))
    p_true <- solve_pde(u_true)
    
    range_p <- diff(range(p_true, na.rm = TRUE))
    if (is.finite(range_p) && range_p <= P_RANGE_CAP) break
    
    truth_attempt <- truth_attempt + 1L
    if (truth_attempt >= MAX_TRUTH_TRIES) {
      warning(sprintf("trial %d: reached MAX_TRUTH_TRIES=%d but range_p still %.2f (>%.2f). Keeping last truth.",
                      trial_id, MAX_TRUTH_TRIES, range_p, P_RANGE_CAP))
      break
    }
  }
  
  a_true <- exp(u_true)
  
  # (ii) observations
  set.seed(seed_noise + trial_id)
  eps <- rnorm(length(obs_indices), mean = 0, sd = 1)
  observations <- p_true[obs_indices] + (alpha_noise * abs(p_true[obs_indices]) + beta_noise) * eps
  
  # (iii) neg log likelihood (return Inf if solve fails)
  neg_log_likelihood <- function(u, T = 1) {
    p_est <- solve_pde(u)
    if (any(!is.finite(p_est))) return(Inf)
    
    p_obs <- p_est[obs_indices]
    sd_vec <- alpha_noise * abs(p_obs) + beta_noise
    if (any(!is.finite(sd_vec)) || any(sd_vec <= 0)) return(Inf)
    
    like <- sum(log(sd_vec)) + 0.5 * sum(((p_obs - observations) / sd_vec)^2)
    like / max(T, 1e-12)
  }
  
  # (iv) init + MCMC
  set.seed(seed_init + trial_id)
  u_init <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_spde,
                                   graph = graph, type = "mesh"))
  if (length(u_init) != length(u_true)) u_init <- u_init[1:length(u_true)]
  
  set.seed(seed_mcmc + trial_id)
  mcmc_result <- pCN_MCMC(
    n_iter = n_iter,
    beta0 = beta0,
    u_init = u_init,
    neg_log_likelihood_fn = neg_log_likelihood,
    T0 = 1, cooling_factor = 0.95, adapt_interval = 500,
    show_iter_pb = inner_mcmc_pb,
    verbose = verbose_mcmc
  )
  
  posterior_samples <- mcmc_result$samples
  acceptance <- mcmc_result$acceptance
  beta_final <- mcmc_result$beta_final
  acc_rate_all <- mean(acceptance)
  
  posterior_samples_post_burnin <- posterior_samples[(burn_in + 1):n_iter, , drop = FALSE]
  
  # (v) estimators for u
  u_post_mean <- apply(posterior_samples_post_burnin, 2, mean)
  u_map_raw   <- apply(posterior_samples_post_burnin, 2, get_mode)
  
  # (v-a) a = exp(u) estimators
  a_post_mean <- exp(u_post_mean)
  a_map_raw   <- exp(u_map_raw)
  
  # (vi) plug-in p fields
  p_post_mean <- solve_pde(u_post_mean)
  p_map_raw   <- solve_pde(u_map_raw)
  
  # (vii) posterior sd summary for u
  sd_u <- apply(posterior_samples_post_burnin, 2, sd)
  sd_u_mean <- mean(sd_u, na.rm = TRUE)
  sd_u_max  <- max(sd_u, na.rm = TRUE)
  
  # (viii) metrics
  out <- list(
    trial = trial_id,
    acc_rate_all = acc_rate_all,
    beta_final = beta_final,
    
    # RMSE(u)
    rmse_u_mean = rmse(u_post_mean, u_true),
    rmse_u_map_raw = rmse(u_map_raw, u_true),
    
    # relative L2(u) (Euclidean)
    rermse_u_mean = rel_l2(u_post_mean, u_true),
    rermse_u_map_raw = rel_l2(u_map_raw, u_true),
    
    # range-normalized RMSE(u)
    relative_error_u_mean = rel_range(u_post_mean, u_true),
    relative_error_u_map_raw = rel_range(u_map_raw, u_true),
    
    # RMSE(a=exp(u))
    rmse_a_mean = rmse(a_post_mean, a_true),
    rmse_a_map_raw = rmse(a_map_raw, a_true),
    
    # relative L2(Gamma) error for a using mass matrix
    relL2G_a_mean = rel_l2_gamma(a_post_mean, a_true, C),
    relL2G_a_map_raw = rel_l2_gamma(a_map_raw, a_true, C),
    
    # RMSE(p)
    rmse_p_mean_all = rmse(p_post_mean, p_true),
    rmse_p_map_raw_all = rmse(p_map_raw, p_true),
    
    rmse_p_mean_obs = rmse(p_post_mean[obs_indices], p_true[obs_indices]),
    rmse_p_map_raw_obs = rmse(p_map_raw[obs_indices], p_true[obs_indices]),
    
    # relative L2(p) all nodes (Euclidean)
    rermse_p_mean_all = rel_l2(p_post_mean, p_true),
    rermse_p_map_raw_all = rel_l2(p_map_raw, p_true),
    
    # optional normalized RMSE (p)
    nrmse_range_p_mean_all = nrmse_range(p_post_mean, p_true),
    nrmse_range_p_map_raw_all = nrmse_range(p_map_raw, p_true),
    nrmse_sd_p_mean_all = nrmse_sd(p_post_mean, p_true),
    nrmse_sd_p_map_raw_all = nrmse_sd(p_map_raw, p_true),
    
    # posterior spread summary
    sd_u_mean = sd_u_mean,
    sd_u_max  = sd_u_max,
    
    # truth scale
    p_range_true = diff(range(p_true, na.rm = TRUE)),
    truth_attempt = truth_attempt
  )
  
  # (ix) SAVE per-trial bundle
  trial_settings <- list(
    n_iter = n_iter, burn_in = burn_in, beta0 = beta0, gamma_smooth = gamma_smooth,
    seed_truth_used = seed_truth_used,
    truth_attempt = truth_attempt,
    seed_noise = seed_noise + trial_id,
    seed_init  = seed_init  + trial_id,
    seed_mcmc  = seed_mcmc  + trial_id,
    P_RANGE_CAP = P_RANGE_CAP,
    MAX_TRUTH_TRIES = MAX_TRUTH_TRIES
  )
  
  save_one_trial(
    OUTDIR = OUTDIR, trial_id = trial_id,
    u_true = u_true, p_true = p_true,
    obs_indices = obs_indices, observations = observations,
    u_post_mean = u_post_mean, u_map_raw = u_map_raw,
    p_post_mean = p_post_mean, p_map_raw = p_map_raw,
    a_true = a_true, a_post_mean = a_post_mean, a_map_raw = a_map_raw,
    acc_rate_all = acc_rate_all, beta_final = beta_final,
    acceptance_vec = acceptance,
    sd_u_mean = sd_u_mean, sd_u_max = sd_u_max,
    posterior_samples_post_burnin = posterior_samples_post_burnin,
    save_thinned_samples = save_thinned_samples,
    thin_every = thin_every,
    trial_settings = trial_settings
  )
  
  as.data.frame(out)
}

# ------------------------------------------------
# 9) Outer loop settings
# ------------------------------------------------
n_trials <- 30
n_iter  <- 60000
burn_in <- 30000
beta0   <- 0.3
gamma_smooth <- 0.1

save_global_bundle(OUTDIR, n_iter = n_iter, burn_in = burn_in, beta0 = beta0, gamma_smooth = gamma_smooth)

metrics_list <- vector("list", n_trials)
pb_outer <- utils::txtProgressBar(min = 0, max = n_trials, style = 3)

for (t in 1:n_trials) {
  metrics_list[[t]] <- run_one_trial(
    trial_id = t, OUTDIR = OUTDIR,
    n_iter = n_iter, burn_in = burn_in, beta0 = beta0, gamma_smooth = gamma_smooth,
    inner_mcmc_pb = TRUE,
    verbose_mcmc = TRUE,
    save_thinned_samples = FALSE
  )
  utils::setTxtProgressBar(pb_outer, t)
}
close(pb_outer)

metrics_df <- do.call(rbind, metrics_list)

# ------------------------------------------------
# 10) Report stage: tables + single boxplots
# ------------------------------------------------
REPORT_DIR <- file.path(OUTDIR, "report_stage_full")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

safe_write_csv(metrics_df, file.path(REPORT_DIR, sprintf("metrics_trials_%d.csv", n_trials)))

num_cols <- setdiff(names(metrics_df), "trial")
summary_df <- data.frame(
  metric = num_cols,
  mean   = sapply(metrics_df[, num_cols, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
  sd     = sapply(metrics_df[, num_cols, drop = FALSE], function(z) sd(z, na.rm = TRUE))
)
safe_write_csv(summary_df, file.path(REPORT_DIR, sprintf("metrics_summary_%d.csv", n_trials)))

cat("\n================= Per-trial metrics (head) =================\n")
print(head(metrics_df))
cat("\n================= Summary mean ± sd =================\n")
print(summary_df)

# ----  core table
mean_sd <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
fmt_mean_sd <- function(m, s, digits = 4) sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), m, s)
get_ms <- function(df, col) {
  if (!col %in% names(df)) return(c(mean = NA_real_, sd = NA_real_))
  mean_sd(df[[col]])
}

ms <- list(
  rmse_u_mean       = get_ms(metrics_df, "rmse_u_mean"),
  rmse_u_map_raw    = get_ms(metrics_df, "rmse_u_map_raw"),
  rermse_u_mean     = get_ms(metrics_df, "rermse_u_mean"),
  rermse_u_map_raw  = get_ms(metrics_df, "rermse_u_map_raw"),
  relrange_u_mean   = get_ms(metrics_df, "relative_error_u_mean"),
  relrange_u_map_raw= get_ms(metrics_df, "relative_error_u_map_raw"),
  
  rmse_a_mean       = get_ms(metrics_df, "rmse_a_mean"),
  rmse_a_map_raw    = get_ms(metrics_df, "rmse_a_map_raw"),
  relL2G_a_mean     = get_ms(metrics_df, "relL2G_a_mean"),
  relL2G_a_map_raw  = get_ms(metrics_df, "relL2G_a_map_raw"),
  
  rmse_p_mean_all   = get_ms(metrics_df, "rmse_p_mean_all"),
  rmse_p_map_raw_all= get_ms(metrics_df, "rmse_p_map_raw_all"),
  rermse_p_mean_all = get_ms(metrics_df, "rermse_p_mean_all"),
  rermse_p_map_raw_all = get_ms(metrics_df, "rermse_p_map_raw_all")
)

paper_table <- data.frame(
  Metric = c("RMSE(u) over all nodes",
             "relative_error(u) = RMSE(u)/range(u_true)",
             "RERMSE(u) = relative L2(u) (Euclidean)",
             "RMSE(a=exp(u)) over all nodes",
             "RelErr_{L2(Gamma)}(a=exp(u)) (mass-matrix)",
             "RMSE(p) over all nodes",
             "RERMSE(p) = relative L2(p) (Euclidean)"),
  `Posterior mean` = c(
    fmt_mean_sd(ms$rmse_u_mean["mean"], ms$rmse_u_mean["sd"], digits = 4),
    fmt_mean_sd(ms$relrange_u_mean["mean"], ms$relrange_u_mean["sd"], digits = 4),
    fmt_mean_sd(ms$rermse_u_mean["mean"], ms$rermse_u_mean["sd"], digits = 4),
    fmt_mean_sd(ms$rmse_a_mean["mean"], ms$rmse_a_mean["sd"], digits = 4),
    fmt_mean_sd(ms$relL2G_a_mean["mean"], ms$relL2G_a_mean["sd"], digits = 4),
    fmt_mean_sd(ms$rmse_p_mean_all["mean"], ms$rmse_p_mean_all["sd"], digits = 4),
    fmt_mean_sd(ms$rermse_p_mean_all["mean"], ms$rermse_p_mean_all["sd"], digits = 4)
  ),
  `Raw MAP` = c(
    fmt_mean_sd(ms$rmse_u_map_raw["mean"], ms$rmse_u_map_raw["sd"], digits = 4),
    fmt_mean_sd(ms$relrange_u_map_raw["mean"], ms$relrange_u_map_raw["sd"], digits = 4),
    fmt_mean_sd(ms$rermse_u_map_raw["mean"], ms$rermse_u_map_raw["sd"], digits = 4),
    fmt_mean_sd(ms$rmse_a_map_raw["mean"], ms$rmse_a_map_raw["sd"], digits = 4),
    fmt_mean_sd(ms$relL2G_a_map_raw["mean"], ms$relL2G_a_map_raw["sd"], digits = 4),
    fmt_mean_sd(ms$rmse_p_map_raw_all["mean"], ms$rmse_p_map_raw_all["sd"], digits = 4),
    fmt_mean_sd(ms$rermse_p_map_raw_all["mean"], ms$rermse_p_map_raw_all["sd"], digits = 4)
  ),
  check.names = FALSE
)
safe_write_csv(paper_table, file.path(REPORT_DIR, "paper_summary_table.csv"))
cat("\n========== Paper core table (mean ± sd) ==========\n")
print(paper_table)

# ---- LaTeX table 
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}
fmt_mean_pm_latex <- function(m, s, digits = 4) {
  if (!is.finite(m) || !is.finite(s)) return("")
  sprintf("$%.*f \\\\pm %.*f$", digits, m, digits, s)
}

paper_table_tex <- paper_table
paper_table_tex$`Posterior mean` <- c(
  fmt_mean_pm_latex(ms$rmse_u_mean["mean"], ms$rmse_u_mean["sd"], 4),
  fmt_mean_pm_latex(ms$relrange_u_mean["mean"], ms$relrange_u_mean["sd"], 4),
  fmt_mean_pm_latex(ms$rermse_u_mean["mean"], ms$rermse_u_mean["sd"], 4),
  fmt_mean_pm_latex(ms$rmse_a_mean["mean"], ms$rmse_a_mean["sd"], 4),
  fmt_mean_pm_latex(ms$relL2G_a_mean["mean"], ms$relL2G_a_mean["sd"], 4),
  fmt_mean_pm_latex(ms$rmse_p_mean_all["mean"], ms$rmse_p_mean_all["sd"], 4),
  fmt_mean_pm_latex(ms$rermse_p_mean_all["mean"], ms$rermse_p_mean_all["sd"], 4)
)
paper_table_tex$`Raw MAP` <- c(
  fmt_mean_pm_latex(ms$rmse_u_map_raw["mean"], ms$rmse_u_map_raw["sd"], 4),
  fmt_mean_pm_latex(ms$relrange_u_map_raw["mean"], ms$relrange_u_map_raw["sd"], 4),
  fmt_mean_pm_latex(ms$rermse_u_map_raw["mean"], ms$rermse_u_map_raw["sd"], 4),
  fmt_mean_pm_latex(ms$rmse_a_map_raw["mean"], ms$rmse_a_map_raw["sd"], 4),
  fmt_mean_pm_latex(ms$relL2G_a_map_raw["mean"], ms$relL2G_a_map_raw["sd"], 4),
  fmt_mean_pm_latex(ms$rmse_p_map_raw_all["mean"], ms$rmse_p_map_raw_all["sd"], 4),
  fmt_mean_pm_latex(ms$rermse_p_map_raw_all["mean"], ms$rermse_p_map_raw_all["sd"], 4)
)

to_latex_booktabs <- function(df, caption, label) {
  coln <- latex_escape(names(df))
  body <- apply(df, 1, function(r) paste(r, collapse = " & "))
  paste0(
    "\\begin{table}[t]\n\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{lcc}\n\\toprule\n",
    paste(coln, collapse = " & "), " \\\\\n\\midrule\n",
    paste0(body, " \\\\ \n", collapse = ""),
    "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
  )
}

latex_tab <- to_latex_booktabs(
  data.frame(
    Metric = latex_escape(paper_table_tex$Metric),
    `Posterior mean` = latex_escape(paper_table_tex$`Posterior mean`),
    `Raw MAP` = latex_escape(paper_table_tex$`Raw MAP`),
    check.names = FALSE
  ),
  caption = "Error metrics aggregated over multiple independent trials (mean $\\pm$ standard deviation).",
  label = "tab:metrics_summary"
)
writeLines(latex_tab, con = file.path(REPORT_DIR, "paper_metrics_table.tex"))

# ------------------------------------------------
# 11) Single boxplots
# ------------------------------------------------
save_boxplot_2 <- function(df_long, outbase, title, ylab) {
  p <- ggplot(df_long, aes(x = estimator, y = value)) +
    geom_boxplot(outlier.size = 0.8) +
    theme_bw() +
    labs(title = title, x = "", y = ylab)
  print(p)
  ggsave(file.path(REPORT_DIR, paste0(outbase, ".png")), p, width = 6.5, height = 4.2, dpi = 220)
  ggsave(file.path(REPORT_DIR, paste0(outbase, ".pdf")), p, width = 6.5, height = 4.2)
  p
}

plots <- list()

plots$rmse_u <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "Posterior mean", value = metrics_df$rmse_u_mean),
    data.frame(trial = metrics_df$trial, estimator = "Raw MAP",        value = metrics_df$rmse_u_map_raw)
  ),
  "fig_rmse_u_boxplot", "RMSE(u) across trials", "RMSE(u)"
)

plots$relrange_u <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "Posterior mean", value = metrics_df$relative_error_u_mean),
    data.frame(trial = metrics_df$trial, estimator = "Raw MAP",        value = metrics_df$relative_error_u_map_raw)
  ),
  "fig_relative_error_u_range_boxplot", "relative_error(u)=RMSE/range across trials", "relative_error(u)"
)

plots$rermse_u <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "Posterior mean", value = metrics_df$rermse_u_mean),
    data.frame(trial = metrics_df$trial, estimator = "Raw MAP",        value = metrics_df$rermse_u_map_raw)
  ),
  "fig_rermse_u_boxplot", "RERMSE(u) (relative L2, Euclidean) across trials", "RERMSE(u)"
)

# NEW: a metrics
plots$rmse_a <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "Posterior mean", value = metrics_df$rmse_a_mean),
    data.frame(trial = metrics_df$trial, estimator = "Raw MAP",        value = metrics_df$rmse_a_map_raw)
  ),
  "fig_rmse_a_boxplot", "RMSE(a=exp(u)) across trials", "RMSE(a)"
)

plots$relL2G_a <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "Posterior mean", value = metrics_df$relL2G_a_mean),
    data.frame(trial = metrics_df$trial, estimator = "Raw MAP",        value = metrics_df$relL2G_a_map_raw)
  ),
  "fig_relL2G_a_boxplot", "RelErr_{L2(Gamma)}(a=exp(u)) across trials", "RelErr L2(Gamma)"
)

plots$rmse_p_all <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "p(mean u)",  value = metrics_df$rmse_p_mean_all),
    data.frame(trial = metrics_df$trial, estimator = "p(raw MAP)", value = metrics_df$rmse_p_map_raw_all)
  ),
  "fig_rmse_p_all_boxplot", "RMSE(p) over all nodes across trials", "RMSE(p)"
)

plots$rermse_p_all <- save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "p(mean u)",  value = metrics_df$rermse_p_mean_all),
    data.frame(trial = metrics_df$trial, estimator = "p(raw MAP)", value = metrics_df$rermse_p_map_raw_all)
  ),
  "fig_rermse_p_all_boxplot", "RERMSE(p) over all nodes across trials", "RERMSE(p)"
)

# Diagnostics: acceptance & beta
p_acc <- ggplot(metrics_df, aes(x = trial, y = acc_rate_all)) +
  geom_point() + geom_line() + theme_bw() +
  labs(title = "Overall acceptance rate per trial", x = "trial", y = "acceptance rate")
print(p_acc)
ggsave(file.path(REPORT_DIR, "fig_accept_rate.png"), p_acc, width = 6.5, height = 4.2, dpi = 220)
ggsave(file.path(REPORT_DIR, "fig_accept_rate.pdf"), p_acc, width = 6.5, height = 4.2)

p_beta <- ggplot(metrics_df, aes(x = trial, y = beta_final)) +
  geom_point() + geom_line() + theme_bw() +
  labs(title = "Final beta per trial", x = "trial", y = "beta_final")
print(p_beta)
ggsave(file.path(REPORT_DIR, "fig_beta_final.png"), p_beta, width = 6.5, height = 4.2, dpi = 220)
ggsave(file.path(REPORT_DIR, "fig_beta_final.pdf"), p_beta, width = 6.5, height = 4.2)

report_bundle <- list(
  timestamp = Sys.time(),
  outdir = OUTDIR,
  report_dir = REPORT_DIR,
  metrics_df = metrics_df,
  summary_df = summary_df,
  paper_table = paper_table,
  latex_table = latex_tab,
  plots = plots,
  diag_plots = list(p_acc = p_acc, p_beta = p_beta)
)
safe_save_rds(report_bundle, file.path(REPORT_DIR, "report_bundle.rds"))

cat("\nDONE.\n",
    "OUTDIR:\n  ", OUTDIR, "\n",
    "Key outputs:\n",
    "  - global/global_bundle.rds\n",
    "  - trial_###/trial_bundle.rds (per trial)\n",
    "  - report_stage_full/metrics_trials_*.csv\n",
    "  - report_stage_full/metrics_summary_*.csv\n",
    "  - report_stage_full/paper_summary_table.csv\n",
    "  - report_stage_full/paper_metrics_table.tex\n",
    "  - report_stage_full/fig_*.png/.pdf\n",
    "  - report_stage_full/report_bundle.rds\n", sep = "")
