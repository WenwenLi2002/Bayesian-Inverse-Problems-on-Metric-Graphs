suppressPackageStartupMessages({
  library(rSPDE)
  library(MetricGraph)
  library(Matrix)
  library(ggplot2)
  library(viridis)
  library(scales)
})

# ------------------------------------------------------------
# 0) Global parameters
# ------------------------------------------------------------
set.seed(123)

# Mesh / PDE
h_mesh  <- 0.2
kappa_L <- 1

# Prior/SPDE hyperparams 
alpha_prior <- 2
range_prior <- 3
sigma_prior <- 0.05

# Forward (fractional)
use_fractional <- TRUE
frac_order <- 1.5

# Observation model (heteroscedastic)
alpha_noise <- 0.03
beta_noise  <- 0.2
n_obs_max   <- 507
seed_obsidx <- 123   

# MCMC
n_iter        <- 8e4
burn_in_user  <- 3e4
beta_pcn_init <- 0.3
T0            <- 20
cooling       <- 0.95
adapt_every   <- 500

# Posterior var(p) (streaming + thinning)
compute_var_p <- TRUE
thin_p        <- 5

# Outer loop
n_trials <- 10

# Progress bars
inner_mcmc_pb <- TRUE   

# Saving (ON)
SAVE_ALL <- TRUE

# Base seeds (per-trial)
seed_truth_base <- 1000
seed_noise_base <- 2000
seed_init_base  <- 3000
seed_mcmc_base  <- 4000

# ------------------------------------------------------------
# Helper: safe exp to prevent overflow/underflow
# ------------------------------------------------------------
safe_exp <- function(u, clip = 20) {
  exp(pmax(pmin(as.numeric(u), clip), -clip))
}

# ------------------------------------------------------------
# 1) Build graph + mesh + FEM + RHS f 
# ------------------------------------------------------------
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance = list(edge_edge = 1e-3,
                   vertex_vertex = 1e-3,
                   edge_vertex   = 1e-3)
)
graph$build_mesh(h = h_mesh)

meshV <- as.matrix(graph$mesh$V)
if (nrow(meshV) < 2) stop("Mesh has too few vertices. Did you define a graph with edges?")
if (ncol(meshV) < 2) stop("Mesh vertex matrix has fewer than 2 columns!")
colnames(meshV) <- c("x", "y")
graph$mesh$V <- meshV
x <- meshV[, "x"]
y <- meshV[, "y"]

graph$compute_fem()
G <- graph$mesh$G
C <- graph$mesh$C
Cmat <- as(C, "dgCMatrix")

fbar <- x^2 - y^2
f <- as.numeric(Cmat %*% fbar)

n_nodes <- nrow(meshV)

# ------------------------------------------------------------
# 2) Assemble L(u) and solvers
# ------------------------------------------------------------
assemble_L <- function(u) {
  a  <- safe_exp(u, clip = 20)
  Du <- Diagonal(x = a)
  Gu <- (G %*% Du + Du %*% G) * 0.5
  L  <- kappa_L^2 * Cmat + Gu
  (L + t(L)) * 0.5
}

solve_pde_beta1 <- function(u) {
  L <- assemble_L(u)
  out <- tryCatch({
    chol <- Cholesky(L, LDL = FALSE, Imult = 0)
    as.numeric(solve(chol, f))
  }, error = function(e) rep(NA_real_, length(f)))
  out
}

solve_pde_fractional <- function(u, frac_order) {
  L <- assemble_L(u)
  out <- tryCatch({
    Ld  <- as.matrix(L)
    eig <- eigen(Ld, symmetric = TRUE)
    
    lam_floor <- 1e-12
    lam <- pmax(eig$values, lam_floor)
    
    vals <- lam^(-frac_order)
    V <- eig$vectors
    as.numeric(V %*% (vals * crossprod(V, f)))
  }, error = function(e) rep(NA_real_, length(f)))
  out
}

solve_pde <- function(u) {
  if (use_fractional) solve_pde_fractional(u, frac_order) else solve_pde_beta1(u)
}

# ------------------------------------------------------------
# 3) Fix observation locations 
# ------------------------------------------------------------
set.seed(seed_obsidx)
n_obs <- min(n_obs_max, n_nodes)
obs_indices <- sample.int(n_nodes, n_obs)

# ------------------------------------------------------------
# 4) Robust negative log-likelihood factory 
# ------------------------------------------------------------
make_neg_log_likelihood <- function(observations) {
  function(u, T = 1) {
    p_est <- solve_pde(u)
    if (any(!is.finite(p_est))) return(Inf)
    
    p_obs <- p_est[obs_indices]
    if (any(!is.finite(p_obs))) return(Inf)
    
    sd_v <- alpha_noise * abs(p_obs) + beta_noise
    sd_v <- pmax(sd_v, 1e-10)
    if (any(!is.finite(sd_v))) return(Inf)
    
    ll <- sum(log(sd_v)) + 0.5 * sum(((p_obs - observations) / sd_v)^2)
    if (!is.finite(ll)) return(Inf)
    
    ll / T
  }
}

# ------------------------------------------------------------
# 5) Robust pCN MCMC (supports per-trial likelihood)
# ------------------------------------------------------------
pCN_MCMC <- function(n_iter, beta_init, u_init, neg_log_likelihood_fn,
                     T0 = 20, cooling = 0.95, adapt_every = 500,
                     show_progress = TRUE) {
  
  u_curr <- as.numeric(u_init)
  n_dim  <- length(u_curr)
  
  samples <- matrix(NA_real_, nrow = n_iter, ncol = n_dim)
  acc     <- logical(n_iter)
  
  i_star <- ceiling(adapt_every * log(1 / T0) / log(cooling))
  beta   <- beta_init
  
  pb <- NULL
  pb_interval <- max(1L, floor(n_iter / 1000))
  if (show_progress) pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  
  for (i in seq_len(n_iter)) {
    
    T_curr <- max(1, T0 * cooling^(i / adapt_every))
    
    phi_c <- neg_log_likelihood_fn(u_curr, T = T_curr)
    
    xi <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_prior,
                                 graph = graph, type = "mesh"))
    if (length(xi) != n_dim) xi <- xi[seq_len(n_dim)]
    
    u_prop <- sqrt(1 - beta^2) * u_curr + beta * xi
    phi_p  <- neg_log_likelihood_fn(u_prop, T = T_curr)
    
    log_alpha <- phi_c - phi_p
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      u_curr <- u_prop
      acc[i] <- TRUE
    }
    samples[i, ] <- u_curr
    
    if (i %% adapt_every == 0 && i >= i_star) {
      rate <- mean(acc[(i - adapt_every + 1):i])
      beta <- beta * ifelse(rate < 0.36, 0.9, ifelse(rate > 0.44, 1.2, 1))
      beta <- min(max(beta, 1e-3), 0.999)
      cat(sprintf("\nIter %d | acc %.3f | beta %.3f | T %.2f\n", i, rate, beta, T_curr))
    }
    
    if (show_progress && (i %% pb_interval == 0L || i == n_iter)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  if (show_progress) close(pb)
  
  list(samples = samples, acceptance = acc, i_star = i_star, beta_final = beta)
}

# ------------------------------------------------------------
# 6) Estimators + errors
# ------------------------------------------------------------
get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

rmse <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  sqrt(mean((est[ok] - truth[ok])^2))
}

nrmse_range <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  rng <- diff(range(truth[ok], na.rm = TRUE))
  rmse(est[ok], truth[ok]) / rng
}

nrmse_sd <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  s <- sd(truth[ok], na.rm = TRUE)
  rmse(est[ok], truth[ok]) / s
}

# ------------------------------------------------------------
# 7) SAVE-ALL helpers
# ------------------------------------------------------------
make_outdir <- function(prefix = "inverse_mcmc_fractional") {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(getwd(), paste0(prefix, "_", stamp))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir
}

safe_save_rds <- function(obj, path) {
  tryCatch({ saveRDS(obj, path); TRUE },
           error = function(e) { message("Failed saveRDS: ", path, " | ", e$message); FALSE })
}

safe_write_csv <- function(df, path) {
  tryCatch({ write.csv(df, path, row.names = FALSE); TRUE },
           error = function(e) { message("Failed write.csv: ", path, " | ", e$message); FALSE })
}

save_global_bundle <- function(OUTDIR) {
  global_dir <- file.path(OUTDIR, "global")
  dir.create(global_dir, showWarnings = FALSE)
  
  global_bundle <- list(
    timestamp = Sys.time(),
    meshV = meshV,
    E     = graph$mesh$E,
    C = C, G = G,
    fbar = fbar, f = f,
    kappa_L = kappa_L,
    # prior
    alpha_prior = alpha_prior, range_prior = range_prior, sigma_prior = sigma_prior,
    # forward
    use_fractional = use_fractional, frac_order = frac_order,
    # obs
    obs_indices = obs_indices, n_obs = length(obs_indices),
    alpha_noise = alpha_noise, beta_noise = beta_noise,
    # mcmc
    n_iter = n_iter, burn_in_user = burn_in_user, beta_pcn_init = beta_pcn_init,
    T0 = T0, cooling = cooling, adapt_every = adapt_every,
    # var(p)
    compute_var_p = compute_var_p, thin_p = thin_p,
    sessionInfo = sessionInfo()
  )
  
  safe_save_rds(global_bundle, file.path(global_dir, "global_bundle.rds"))
  capture.output(str(global_bundle), file = file.path(global_dir, "global_bundle_str.txt"))
}

save_one_trial <- function(OUTDIR, trial_id,
                           u_true, p_true,
                           observations,
                           posterior_mean, u_map_raw,
                           p_mean, p_map_raw,
                           acc_rate_all, beta_final, i_star, burn_used,
                           var_u = NULL, sd_u = NULL,
                           var_p = NULL, sd_p = NULL,
                           acceptance_vec = NULL) {
  
  trial_dir <- file.path(OUTDIR, sprintf("trial_%03d", trial_id))
  dir.create(trial_dir, showWarnings = FALSE)
  
  trial_bundle <- list(
    trial = trial_id,
    timestamp = Sys.time(),
    # truths (you explicitly requested u_true per trial)
    u_true = as.numeric(u_true),
    p_true = as.numeric(p_true),
    # data
    obs_indices = as.integer(obs_indices),
    observations = as.numeric(observations),
    # estimators (u)
    u_post_mean = as.numeric(posterior_mean),
    u_map_raw   = as.numeric(u_map_raw),
    # plug-in p fields
    p_post_mean = as.numeric(p_mean),
    p_map_raw   = as.numeric(p_map_raw),
    # mcmc diag
    acc_rate_all = as.numeric(acc_rate_all),
    beta_final   = as.numeric(beta_final),
    i_star       = as.integer(i_star),
    burn_used    = as.integer(burn_used)
  )
  
  if (!is.null(acceptance_vec)) trial_bundle$acceptance_vec <- as.numeric(acceptance_vec)
  if (!is.null(var_u)) trial_bundle$var_u <- as.numeric(var_u)
  if (!is.null(sd_u))  trial_bundle$sd_u  <- as.numeric(sd_u)
  if (!is.null(var_p)) trial_bundle$var_p <- as.numeric(var_p)
  if (!is.null(sd_p))  trial_bundle$sd_p  <- as.numeric(sd_p)
  
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

save_report_bundle <- function(OUTDIR, metrics_df, summary_df, p_rmse, p_nrmse, p_acc, p_beta) {
  report_dir <- file.path(OUTDIR, "report_stage")
  dir.create(report_dir, showWarnings = FALSE)
  
  safe_write_csv(metrics_df, file.path(report_dir, "metrics_df.csv"))
  safe_write_csv(summary_df, file.path(report_dir, "summary_df.csv"))
  
  report_bundle <- list(
    timestamp = Sys.time(),
    metrics_df = metrics_df,
    summary_df = summary_df,
    plots = list(p_rmse = p_rmse, p_nrmse = p_nrmse, p_acc = p_acc, p_beta = p_beta)
  )
  safe_save_rds(report_bundle, file.path(report_dir, "report_bundle.rds"))
  
  ggsave(file.path(report_dir, "rmse_boxplot.png"),  p_rmse,  width = 10, height = 6, dpi = 200)
  ggsave(file.path(report_dir, "nrmse_boxplot.png"), p_nrmse, width = 10, height = 5, dpi = 200)
  ggsave(file.path(report_dir, "accept_rate.png"),   p_acc,   width = 9,  height = 4, dpi = 200)
  ggsave(file.path(report_dir, "beta_final.png"),    p_beta,  width = 9,  height = 4, dpi = 200)
}

# Create OUTDIR + save global once
if (SAVE_ALL) {
  OUTDIR <- make_outdir(prefix = "inverse_mcmc_fractional")
  message("Saving everything under: ", OUTDIR)
  save_global_bundle(OUTDIR)
}

# ------------------------------------------------------------
# 8) One-trial pipeline (fractional): truth -> data -> MCMC -> metrics 
# ------------------------------------------------------------
run_one_trial_fractional <- function(trial_id) {
  
  # (i) truth
  set.seed(seed_truth_base + trial_id)
  u_true <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_prior,
                                   graph = graph, type = "mesh"))
  if (length(u_true) != n_nodes) u_true <- u_true[seq_len(n_nodes)]
  p_true <- solve_pde(u_true)
  
  # (ii) data
  set.seed(seed_noise_base + trial_id)
  eps <- rnorm(length(obs_indices))
  sd_true <- alpha_noise * abs(p_true[obs_indices]) + beta_noise
  sd_true <- pmax(sd_true, 1e-10)
  observations <- p_true[obs_indices] + sd_true * eps
  
  # (iii) likelihood closure
  nll <- make_neg_log_likelihood(observations)
  
  # (iv) init + mcmc
  set.seed(seed_init_base + trial_id)
  u_init <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_prior,
                                   graph = graph, type = "mesh"))
  if (length(u_init) != n_nodes) u_init <- u_init[seq_len(n_nodes)]
  
  set.seed(seed_mcmc_base + trial_id)
  mcmc_res <- pCN_MCMC(
    n_iter = n_iter,
    beta_init = beta_pcn_init,
    u_init = u_init,
    neg_log_likelihood_fn = nll,
    T0 = T0, cooling = cooling, adapt_every = adapt_every,
    show_progress = inner_mcmc_pb
  )
  
  post_samps <- mcmc_res$samples
  acc_vec    <- mcmc_res$acceptance
  acc_rate   <- mean(acc_vec)
  beta_final <- mcmc_res$beta_final
  i_star     <- mcmc_res$i_star
  
  burn_used <- max(burn_in_user, i_star + 1000)
  if (burn_used >= n_iter - 10) stop("Burn-in too large relative to n_iter. Reduce burn_in_user or increase n_iter.")
  post <- post_samps[(burn_used + 1):n_iter, , drop = FALSE]
  
  # (v) estimators u
  posterior_mean <- colMeans(post)
  u_map_raw <- apply(post, 2, get_mode)
  
  # (vi) plug-in p fields
  p_mean    <- solve_pde(posterior_mean)
  p_map_raw <- solve_pde(u_map_raw)
  
  # (vii) var(u) 
  var_u <- apply(post, 2, var)
  sd_u  <- sqrt(var_u)
  
  # (viii) var(p) (streaming + thinning)
  var_p <- NULL; sd_p <- NULL
  if (isTRUE(compute_var_p)) {
    idx <- seq(1, nrow(post), by = thin_p)
    n_use <- length(idx)
    
    m <- rep(0, n_nodes)
    s <- rep(0, n_nodes)
    k <- 0L
    
    pb_p <- utils::txtProgressBar(min = 0, max = n_use, style = 3)
    for (jj in seq_along(idx)) {
      ii <- idx[jj]
      pj <- solve_pde(post[ii, ])
      if (any(!is.finite(pj))) { utils::setTxtProgressBar(pb_p, jj); next }
      k <- k + 1L
      if (k == 1L) {
        m <- pj
      } else {
        delta <- pj - m
        m <- m + delta / k
        s <- s + delta * (pj - m)
      }
      utils::setTxtProgressBar(pb_p, jj)
    }
    close(pb_p)
    
    if (k >= 2) {
      var_p <- s / (k - 1L)
      sd_p  <- sqrt(var_p)
    } else {
      warning(sprintf("Trial %d: too few finite p draws for var(p).", trial_id))
    }
  }
  
  # (ix) metrics row
  out <- data.frame(
    trial = trial_id,
    acc_rate_all = acc_rate,
    beta_final = beta_final,
    
    rmse_u_mean = rmse(posterior_mean, u_true),
    rmse_u_map  = rmse(u_map_raw, u_true),
    
    rmse_p_mean_all = rmse(p_mean, p_true),
    rmse_p_map_all  = rmse(p_map_raw, p_true),
    
    rmse_p_mean_obs = rmse(p_mean[obs_indices], p_true[obs_indices]),
    rmse_p_map_obs  = rmse(p_map_raw[obs_indices], p_true[obs_indices]),
    
    nrmse_range_p_mean_all = nrmse_range(p_mean, p_true),
    nrmse_range_p_map_all  = nrmse_range(p_map_raw, p_true),
    
    nrmse_sd_p_mean_all = nrmse_sd(p_mean, p_true),
    nrmse_sd_p_map_all  = nrmse_sd(p_map_raw, p_true),
    
    sd_u_mean = mean(sd_u, na.rm = TRUE),
    sd_u_max  = max(sd_u, na.rm = TRUE),
    
    burn_used = burn_used,
    i_star = i_star
  )
  
  # (x) save per trial (includes u_true as requested)
  if (SAVE_ALL) {
    save_one_trial(
      OUTDIR = OUTDIR, trial_id = trial_id,
      u_true = u_true, p_true = p_true,
      observations = observations,
      posterior_mean = posterior_mean, u_map_raw = u_map_raw,
      p_mean = p_mean, p_map_raw = p_map_raw,
      acc_rate_all = acc_rate, beta_final = beta_final,
      i_star = i_star, burn_used = burn_used,
      var_u = var_u, sd_u = sd_u,
      var_p = var_p, sd_p = sd_p,
      acceptance_vec = acc_vec
    )
  }
  
  out
}

# ------------------------------------------------------------
# 9) Outer loop: 10 trials with progress bar
# ------------------------------------------------------------
metrics_list <- vector("list", n_trials)
pb_outer <- utils::txtProgressBar(min = 0, max = n_trials, style = 3)

for (t in seq_len(n_trials)) {
  metrics_list[[t]] <- run_one_trial_fractional(trial_id = t)
  utils::setTxtProgressBar(pb_outer, t)
}
close(pb_outer)

metrics_df <- do.call(rbind, metrics_list)

# ------------------------------------------------------------
# 10) Summary table: mean ± sd over 10 trials
# ------------------------------------------------------------
num_cols <- setdiff(names(metrics_df), "trial")
summary_df <- data.frame(
  metric = num_cols,
  mean   = sapply(metrics_df[, num_cols, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
  sd     = sapply(metrics_df[, num_cols, drop = FALSE], function(z) sd(z, na.rm = TRUE))
)

cat("\n================= Per-trial metrics (first 6 rows) =================\n")
print(head(metrics_df))

cat("\n================= Summary over trials: mean ± sd =================\n")
print(summary_df)

# Save tables (also under OUTDIR/report_stage at the end)
if (SAVE_ALL) {
  safe_write_csv(metrics_df, file.path(OUTDIR, "metrics_trials.csv"))
  safe_write_csv(summary_df, file.path(OUTDIR, "metrics_summary.csv"))
}

# ------------------------------------------------------------
# 11) Plots: RMSE distributions (boxplots)
# ------------------------------------------------------------
rmse_long <- rbind(
  data.frame(trial = metrics_df$trial, quantity = "u", estimator = "Posterior mean", value = metrics_df$rmse_u_mean),
  data.frame(trial = metrics_df$trial, quantity = "u", estimator = "Raw MAP",        value = metrics_df$rmse_u_map),
  
  data.frame(trial = metrics_df$trial, quantity = "p (all nodes)", estimator = "p(mean u)",  value = metrics_df$rmse_p_mean_all),
  data.frame(trial = metrics_df$trial, quantity = "p (all nodes)", estimator = "p(raw MAP)", value = metrics_df$rmse_p_map_all),
  
  data.frame(trial = metrics_df$trial, quantity = "p (obs nodes)", estimator = "p(mean u)",  value = metrics_df$rmse_p_mean_obs),
  data.frame(trial = metrics_df$trial, quantity = "p (obs nodes)", estimator = "p(raw MAP)", value = metrics_df$rmse_p_map_obs)
)

p_rmse <- ggplot(rmse_long, aes(x = estimator, y = value)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ quantity, scales = "free_y") +
  theme_bw() +
  labs(title = sprintf("RMSE over %d fractional trials", n_trials), x = "", y = "RMSE")
print(p_rmse)

# ------------------------------------------------------------
# 12) Plots: Normalized RMSE 
# ------------------------------------------------------------
nrmse_long <- rbind(
  data.frame(trial = metrics_df$trial, norm = "NRMSE_range", estimator = "p(mean u)",  value = metrics_df$nrmse_range_p_mean_all),
  data.frame(trial = metrics_df$trial, norm = "NRMSE_range", estimator = "p(raw MAP)", value = metrics_df$nrmse_range_p_map_all),
  
  data.frame(trial = metrics_df$trial, norm = "NRMSE_sd", estimator = "p(mean u)",  value = metrics_df$nrmse_sd_p_mean_all),
  data.frame(trial = metrics_df$trial, norm = "NRMSE_sd", estimator = "p(raw MAP)", value = metrics_df$nrmse_sd_p_map_all)
)

p_nrmse <- ggplot(nrmse_long, aes(x = estimator, y = value)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ norm, scales = "free_y") +
  theme_bw() +
  labs(title = sprintf("Normalized RMSE for p over %d fractional trials", n_trials), x = "", y = "NRMSE")
print(p_nrmse)

# ------------------------------------------------------------
# 13) Plots: acceptance rate / final beta over trials
# ------------------------------------------------------------
p_acc <- ggplot(metrics_df, aes(x = trial, y = acc_rate_all)) +
  geom_point() + geom_line() +
  theme_bw() +
  labs(title = "Acceptance rate (overall) per trial", x = "trial", y = "acceptance rate")
print(p_acc)

p_beta <- ggplot(metrics_df, aes(x = trial, y = beta_final)) +
  geom_point() + geom_line() +
  theme_bw() +
  labs(title = "Final beta per trial", x = "trial", y = "beta_final")
print(p_beta)

# ------------------------------------------------------------
# 14) Save report bundle (tables + plots) under OUTDIR/report_stage
# ------------------------------------------------------------
if (SAVE_ALL) {
  save_report_bundle(OUTDIR, metrics_df, summary_df, p_rmse, p_nrmse, p_acc, p_beta)
  message("Report bundle saved under: ", file.path(OUTDIR, "report_stage"))
}

# Also save local copies (optional)
ggsave("rmse_boxplot_fractional.png",  p_rmse,  width = 10, height = 6, dpi = 200)
ggsave("nrmse_boxplot_fractional.png", p_nrmse, width = 10, height = 5, dpi = 200)
ggsave("accept_rate_fractional.png",   p_acc,   width = 9,  height = 4, dpi = 200)
ggsave("beta_final_fractional.png",    p_beta,  width = 9,  height = 4, dpi = 200)



