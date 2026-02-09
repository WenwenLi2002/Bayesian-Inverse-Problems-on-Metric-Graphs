suppressPackageStartupMessages({
  library(rSPDE)
  library(MetricGraph)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ----------------------------
# 0) USER SETTINGS
# ----------------------------
OUTDIR <- "/Users/wenwenli/inverse_mcmc_metricgraph_20260121_144015"
stopifnot(dir.exists(OUTDIR))
cat("Using OUTDIR =", OUTDIR, "\n")

DELTAS <- c(0, 0.1, 0.2, 0.5)

# Inference settings (consistent with the baseline)
n_iter  <- 60000
burn_in <- 30000
beta0   <- 0.3

# Progress 
inner_mcmc_pb <- TRUE
verbose_mcmc  <- TRUE

# save thinned samples 
save_thinned_samples <- FALSE
thin_every <- 50L

# Seeds:
# - Keep init & MCMC seeds the SAME across deltas for better comparability.
# - Use a delta-dependent seed ONLY for the perturbation noise added to y.
SEED_INIT_BASE  <- 3000L
SEED_MCMC_BASE  <- 4000L
SEED_PERT_BASE  <- 9000L

# ----------------------------
# 1) Load global bundle
# ----------------------------
GLOBAL_PATH <- file.path(OUTDIR, "global", "global_bundle.rds")
stopifnot(file.exists(GLOBAL_PATH))
global_bundle <- readRDS(GLOBAL_PATH)

# Pull essentials
meshV <- as.matrix(global_bundle$meshV)
C <- as(global_bundle$C, "dgCMatrix")
G <- as(global_bundle$G, "dgCMatrix")
f <- as.numeric(global_bundle$f)
kappa <- as.numeric(global_bundle$kappa)

range_prior <- as.numeric(global_bundle$range_prior)
sigma_prior <- as.numeric(global_bundle$sigma_prior)
alpha_spde  <- as.numeric(global_bundle$alpha_spde)

alpha_noise <- as.numeric(global_bundle$alpha_noise)
beta_noise  <- as.numeric(global_bundle$beta_noise)

cat("Loaded global bundle.\n",
    "  n_mesh =", nrow(meshV), "\n",
    "  n_f =", length(f), "\n",
    "  alpha_noise/beta_noise =", alpha_noise, "/", beta_noise, "\n",
    "  prior(range,sigma,alpha) =", range_prior, sigma_prior, alpha_spde, "\n")

# ----------------------------
# 2) Locate baseline trial bundles (delta=0)
# ----------------------------
trial_dirs_all <- list.dirs(OUTDIR, full.names = TRUE, recursive = FALSE)
trial_dirs0 <- trial_dirs_all[grepl("^trial_[0-9]{3}$", basename(trial_dirs_all))]
trial_paths0 <- file.path(trial_dirs0, "trial_bundle.rds")
trial_paths0 <- trial_paths0[file.exists(trial_paths0)]
stopifnot(length(trial_paths0) > 0)

trial_id_from_path <- function(path) {
  d <- basename(dirname(path))  # "trial_007"
  as.integer(sub("^trial_([0-9]{3})$", "\\1", d))
}
trial_ids0 <- sapply(trial_paths0, trial_id_from_path)
ord <- order(trial_ids0)
trial_paths0 <- trial_paths0[ord]
trial_ids0   <- trial_ids0[ord]
cat("Found baseline trials:", length(trial_paths0), "\n")

# ----------------------------
# 3) Build graph object for sample_spde 
#    Rebuild the same graph+mesh.
#    If mismatch occurs, we STOP to avoid wrong obs_indices mapping.
# ----------------------------
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance = list(edge_edge = 1e-3, vertex_vertex = 1e-3, edge_vertex = 1e-3)
)
graph$build_mesh(h = 0.05)
graph$compute_fem()

meshV_re <- as.matrix(graph$mesh$V)
if (nrow(meshV_re) != nrow(meshV)) {
  stop("Rebuilt mesh size != saved mesh size. Cannot safely rerun.\n",
       "Saved n_mesh=", nrow(meshV), ", rebuilt n_mesh=", nrow(meshV_re), "\n",
       "Suggestion: re-run with the exact same graph construction or save graph object next time.")
}

# check coordinates match roughly
coord_diff <- max(abs(meshV_re[,1:2] - meshV[,1:2]))
if (!is.finite(coord_diff) || coord_diff > 1e-6) {
  warning("Rebuilt mesh coordinates differ from saved mesh (max abs diff = ", coord_diff, ").\n",
          "If you suspect ordering differs, STOP and instead save/load the original graph object.")
} else {
  cat("Mesh coordinate check passed (max abs diff ~", coord_diff, ").\n")
}

# ----------------------------
# 4) PDE solver 
# ----------------------------
solve_pde <- function(u) {
  u <- as.numeric(u)
  n_mesh <- nrow(meshV)
  if (length(u) != n_mesh) u <- u[1:n_mesh]
  Du <- Matrix::Diagonal(n_mesh, exp(u))
  Gu <- (G %*% Du + Du %*% G) / 2
  L <- kappa^2 * C + Gu
  # Robust solve: if fail, return NA vector so likelihood becomes large
  out <- tryCatch({
    as.numeric(solve(L, f))
  }, error = function(e) {
    rep(NA_real_, n_mesh)
  })
  out
}

# ----------------------------
# 5) Helpers:  metrics, L2(Gamma) norms
# ----------------------------
get_mode <- function(x, kde_n = 256) {
  d <- density(x, n = kde_n)
  d$x[which.max(d$y)]
}

rmse <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  sqrt(mean((est[ok] - truth[ok])^2))
}

range_norm_rmse <- function(est, truth, eps = 1e-12) {
  truth <- as.numeric(truth)
  ok <- is.finite(truth)
  rng <- diff(range(truth[ok], na.rm = TRUE))
  rmse(est, truth) / max(rng, eps)
}

l2_gamma_norm <- function(v, Cmat) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  v2 <- v; v2[!ok] <- 0
  sqrt(as.numeric(Matrix::crossprod(v2, Cmat %*% v2)))
}

rel_l2_gamma <- function(est, truth, Cmat, eps = 1e-12) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  e <- est; t <- truth
  e[!ok] <- 0; t[!ok] <- 0
  num <- l2_gamma_norm(e - t, Cmat)
  den <- max(l2_gamma_norm(t, Cmat), eps)
  num / den
}

exp_clip <- function(u, clip = 20) {
  u <- as.numeric(u)
  exp(pmin(pmax(u, -clip), clip))
}

# ----------------------------
# 6) pCN MCMC 
# ----------------------------
pCN_MCMC <- function(n_iter, beta0, u_init, neg_log_likelihood_fn,
                     T0 = 20, cooling_factor = 0.95, adapt_interval = 500,
                     show_iter_pb = FALSE, verbose = FALSE) {
  u_current <- as.numeric(u_init)
  n_dim <- length(u_current)
  samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
  acceptance <- numeric(n_iter)
  
  beta <- beta0
  
  if (show_iter_pb) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
    pb_interval <- max(1, floor(n_iter / 1000))
  }
  
  # precompute once (dimension)
  for (i in 1:n_iter) {
    T_current <- max(1, T0 * cooling_factor^(i / adapt_interval))
    
    phi_current <- neg_log_likelihood_fn(u_current, T = T_current)
    
    xi <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_spde,
                                 graph = graph, type = "mesh"))
    if (length(xi) != n_dim) xi <- xi[1:n_dim]
    
    u_proposed <- sqrt(1 - beta^2) * u_current + beta * xi
    phi_proposed <- neg_log_likelihood_fn(u_proposed, T = T_current)
    
    log_alpha <- -(phi_proposed - phi_current)
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      u_current <- u_proposed
      acceptance[i] <- 1
    } else {
      acceptance[i] <- 0
    }
    
    samples[i, ] <- u_current
    
    if (i %% adapt_interval == 0) {
      current_accept_rate <- mean(acceptance[(i - adapt_interval + 1):i])
      if (current_accept_rate < 0.4 * 0.9) {
        beta <- max(beta * 0.9, 0.01)
      } else if (current_accept_rate > 0.4 * 1.1) {
        beta <- beta * 1.2
      }
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

# ----------------------------
# 7) Saving helpers
# ----------------------------
safe_save_rds <- function(obj, path) {
  tryCatch({ saveRDS(obj, path); TRUE }, error = function(e) {
    message("Failed saveRDS: ", path, " | ", e$message); FALSE
  })
}

safe_write_csv <- function(df, path) {
  tryCatch({ write.csv(df, path, row.names = FALSE); TRUE }, error = function(e) {
    message("Failed write.csv: ", path, " | ", e$message); FALSE
  })
}

delta_tag <- function(delta) {
  # 0.1 -> "delta_0p1"
  paste0("delta_", gsub("\\.", "p", sprintf("%.1f", delta)))
}

# ----------------------------
# 8) Run one trial for one perturbed delta
# ----------------------------
run_one_trial_with_observations <- function(trial_id, u_true, p_true, obs_indices,
                                            observations, n_iter, burn_in, beta0,
                                            seed_init, seed_mcmc,
                                            inner_mcmc_pb = TRUE, verbose_mcmc = TRUE) {
  obs_idx <- as.integer(obs_indices)
  obs_idx <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= length(p_true)]
  
  # neg log likelihood 
  neg_log_likelihood <- function(u, T = 1) {
    p_est <- solve_pde(u)
    if (any(!is.finite(p_est))) return(1e30)  
    p_obs <- p_est[obs_idx]
    sd_vec <- alpha_noise * abs(p_obs) + beta_noise
    if (any(!is.finite(sd_vec)) || any(sd_vec <= 0)) return(1e30)
    like <- sum(log(sd_vec)) + 0.5 * sum(((p_obs - observations) / sd_vec)^2)
    like / T
  }
  
  # init
  set.seed(seed_init)
  u_init <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_spde,
                                   graph = graph, type = "mesh"))
  if (length(u_init) != length(u_true)) u_init <- u_init[1:length(u_true)]
  
  # MCMC
  set.seed(seed_mcmc)
  mcmc_result <- pCN_MCMC(
    n_iter = n_iter, beta0 = beta0,
    u_init = u_init,
    neg_log_likelihood_fn = neg_log_likelihood,
    T0 = 20, cooling_factor = 0.95, adapt_interval = 500,
    show_iter_pb = inner_mcmc_pb,
    verbose = verbose_mcmc
  )
  
  posterior_samples <- mcmc_result$samples
  acceptance <- mcmc_result$acceptance
  beta_final <- mcmc_result$beta_final
  acc_rate_all <- mean(acceptance)
  
  posterior_samples_post <- posterior_samples[(burn_in + 1):n_iter, , drop = FALSE]
  
  # estimators for u
  u_post_mean <- colMeans(posterior_samples_post)
  u_map_raw   <- apply(posterior_samples_post, 2, get_mode)
  
  # plug-in p fields
  p_post_mean <- solve_pde(u_post_mean)
  p_map_raw   <- solve_pde(u_map_raw)
  
  # sd summaries
  sd_u <- apply(posterior_samples_post, 2, sd)
  sd_u_mean <- mean(sd_u, na.rm = TRUE)
  sd_u_max  <- max(sd_u, na.rm = TRUE)
  
  list(
    u_post_mean = u_post_mean,
    u_map_raw = u_map_raw,
    p_post_mean = p_post_mean,
    p_map_raw = p_map_raw,
    acceptance = acceptance,
    acc_rate_all = acc_rate_all,
    beta_final = beta_final,
    sd_u_mean = sd_u_mean,
    sd_u_max = sd_u_max,
    posterior_samples_post = posterior_samples_post
  )
}

# ----------------------------
# 9) Perturbation study runner (delta>0 only)
# ----------------------------
PERT_BASEDIR <- file.path(OUTDIR, "perturbation_study")
dir.create(PERT_BASEDIR, showWarnings = FALSE, recursive = TRUE)

# Save a small "study config"
study_config <- list(
  timestamp = Sys.time(),
  outdir = OUTDIR,
  deltas = DELTAS,
  n_iter = n_iter,
  burn_in = burn_in,
  beta0 = beta0,
  seed_init_base = SEED_INIT_BASE,
  seed_mcmc_base = SEED_MCMC_BASE,
  seed_pert_base = SEED_PERT_BASE,
  note = "Perturbation: y_delta = y0 + delta * sigma_true * N(0,1), sigma_true = alpha_noise*|p_true|+beta_noise"
)
safe_save_rds(study_config, file.path(PERT_BASEDIR, "perturbation_study_config.rds"))

# Run inference for each delta>0
deltas_run <- DELTAS[DELTAS > 0]
if (length(deltas_run) > 0) {
  for (delta in deltas_run) {
    dtag <- delta_tag(delta)
    DELTA_DIR <- file.path(PERT_BASEDIR, dtag)
    dir.create(DELTA_DIR, showWarnings = FALSE, recursive = TRUE)
    
    cat("\n=============================\n")
    cat("Running delta =", delta, " -> ", DELTA_DIR, "\n")
    cat("=============================\n")
    
    pb <- utils::txtProgressBar(min = 0, max = length(trial_paths0), style = 3)
    
    for (k in seq_along(trial_paths0)) {
      path0 <- trial_paths0[k]
      b0 <- readRDS(path0)
      tid <- if ("trial" %in% names(b0)) as.integer(b0$trial) else trial_ids0[k]
      
      # Required fields
      u_true <- as.numeric(b0$u_true)
      p_true <- as.numeric(b0$p_true)
      obs_idx <- as.integer(b0$obs_indices)
      y0 <- as.numeric(b0$observations)
      
      # Construct sigma_true at obs using p_true (data-generation scale)
      obs_idx2 <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= length(p_true)]
      sigma_true <- alpha_noise * abs(p_true[obs_idx2]) + beta_noise
      
      # Perturb y0 -> y_delta
      set.seed(SEED_PERT_BASE + 10000L * tid + as.integer(round(delta * 1000)))
      z <- rnorm(length(obs_idx2))
      y_delta <- y0
      # y0 length should match obs_idx;
      if (length(y_delta) != length(obs_idx2)) y_delta <- p_true[obs_idx2] + 0 * z  # fallback
      y_delta <- y_delta + delta * sigma_true * z
      
      # Run inference
      res <- run_one_trial_with_observations(
        trial_id = tid,
        u_true = u_true,
        p_true = p_true,
        obs_indices = obs_idx2,
        observations = y_delta,
        n_iter = n_iter,
        burn_in = burn_in,
        beta0 = beta0,
        seed_init = SEED_INIT_BASE + tid,
        seed_mcmc = SEED_MCMC_BASE + tid,
        inner_mcmc_pb = inner_mcmc_pb,
        verbose_mcmc = verbose_mcmc
      )
      
      # Save trial bundle (delta-specific)
      trial_dir <- file.path(DELTA_DIR, sprintf("trial_%03d", tid))
      dir.create(trial_dir, showWarnings = FALSE)
      
      trial_bundle_delta <- list(
        trial = tid,
        delta = delta,
        timestamp = Sys.time(),
        # Keep truth identical to baseline
        u_true = u_true,
        p_true = p_true,
        obs_indices = obs_idx2,
        # Save BOTH base & perturbed data
        observations_base = y0,
        observations = y_delta,
        perturb_seed = SEED_PERT_BASE + 10000L * tid + as.integer(round(delta * 1000)),
        perturb_z = z,
        sigma_true_obs = sigma_true,
        
        # Estimators
        u_post_mean = as.numeric(res$u_post_mean),
        u_map_raw   = as.numeric(res$u_map_raw),
        p_post_mean = as.numeric(res$p_post_mean),
        p_map_raw   = as.numeric(res$p_map_raw),
        
        # Diagnostics
        acc_rate_all = as.numeric(res$acc_rate_all),
        beta_final   = as.numeric(res$beta_final),
        acceptance_vec = as.numeric(res$acceptance),
        
        # Posterior spread summaries
        sd_u_mean = as.numeric(res$sd_u_mean),
        sd_u_max  = as.numeric(res$sd_u_max),
        
        # Truth scales
        u_range_true = diff(range(u_true, na.rm = TRUE)),
        p_range_true = diff(range(p_true, na.rm = TRUE))
      )
      
      if (isTRUE(save_thinned_samples)) {
        post <- res$posterior_samples_post
        idx <- seq(1, nrow(post), by = as.integer(thin_every))
        trial_bundle_delta$posterior_samples_thin <- post[idx, , drop = FALSE]
        trial_bundle_delta$thin_every <- as.integer(thin_every)
      }
      
      safe_save_rds(trial_bundle_delta, file.path(trial_dir, "trial_bundle.rds"))
      
      utils::setTxtProgressBar(pb, k)
    } # end trials
    close(pb)
  } # end deltas
}

# ============================================================
# 10) REPORT STAGE: load bundles for all deltas + plots
# ============================================================
REPORT_DIR <- file.path(OUTDIR, paste0("report_stage_perturbation_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("\nSaving combined report to:\n  ", REPORT_DIR, "\n")

# ---- helper: list trial_bundle.rds in a directory that has trial_###/
list_trial_bundles <- function(base_dir) {
  dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("^trial_[0-9]{3}$", basename(dirs))]
  paths <- file.path(dirs, "trial_bundle.rds")
  paths[file.exists(paths)]
}

# ---- read metrics from one trial bundle path
read_one_trial_metrics <- function(path, delta_value) {
  b <- readRDS(path)
  
  req <- c("u_true","p_true","u_post_mean","u_map_raw","p_post_mean","p_map_raw","obs_indices")
  miss <- setdiff(req, names(b))
  if (length(miss) > 0) return(NULL)
  
  tid <- if ("trial" %in% names(b)) as.integer(b$trial) else trial_id_from_path(path)
  
  u_true <- as.numeric(b$u_true)
  p_true <- as.numeric(b$p_true)
  u_mean <- as.numeric(b$u_post_mean)
  u_map  <- as.numeric(b$u_map_raw)
  p_mean <- as.numeric(b$p_post_mean)
  p_map  <- as.numeric(b$p_map_raw)
  
  obs_idx <- as.integer(b$obs_indices)
  obs_idx <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= length(p_true)]
  
  a_true <- exp_clip(u_true)
  a_mean <- exp_clip(u_mean)
  a_map  <- exp_clip(u_map)
  
  data.frame(
    delta = delta_value,
    trial = tid,
    trial_path = path,
    
    acc_rate_all = if ("acc_rate_all" %in% names(b)) as.numeric(b$acc_rate_all) else NA_real_,
    beta_final   = if ("beta_final"   %in% names(b)) as.numeric(b$beta_final)   else NA_real_,
    sd_u_mean    = if ("sd_u_mean"    %in% names(b)) as.numeric(b$sd_u_mean)    else NA_real_,
    sd_u_max     = if ("sd_u_max"     %in% names(b)) as.numeric(b$sd_u_max)     else NA_real_,
    
    rmse_u_mean    = rmse(u_mean, u_true),
    rmse_u_map_raw = rmse(u_map,  u_true),
    
    relative_error_u_mean    = range_norm_rmse(u_mean, u_true),
    relative_error_u_map_raw = range_norm_rmse(u_map,  u_true),
    
    rmse_p_mean_all    = rmse(p_mean, p_true),
    rmse_p_map_raw_all = rmse(p_map,  p_true),
    rmse_p_mean_obs    = if (length(obs_idx) > 0) rmse(p_mean[obs_idx], p_true[obs_idx]) else NA_real_,
    rmse_p_map_raw_obs = if (length(obs_idx) > 0) rmse(p_map[obs_idx],  p_true[obs_idx]) else NA_real_,
    
    relL2G_a_mean    = rel_l2_gamma(a_mean, a_true, C),
    relL2G_a_map_raw = rel_l2_gamma(a_map,  a_true, C),
    
    relL2G_p_mean_all    = rel_l2_gamma(p_mean, p_true, C),
    relL2G_p_map_raw_all = rel_l2_gamma(p_map,  p_true, C),
    
    u_range_true = diff(range(u_true, na.rm = TRUE)),
    p_range_true = diff(range(p_true, na.rm = TRUE))
  )
}

# ---- collect metrics for delta=0 from baseline OUTDIR/trial_###
paths0 <- list_trial_bundles(OUTDIR)
metrics0 <- lapply(paths0, read_one_trial_metrics, delta_value = 0)
metrics0 <- Filter(Negate(is.null), metrics0)
metrics0 <- do.call(rbind, metrics0)

# ---- collect metrics for each delta>0 from OUTDIR/perturbation_study/delta_*/
metrics_list <- list(metrics0)

for (delta in DELTAS[DELTAS > 0]) {
  dtag <- delta_tag(delta)
  DELTA_DIR <- file.path(PERT_BASEDIR, dtag)
  if (!dir.exists(DELTA_DIR)) next
  paths <- list_trial_bundles(DELTA_DIR)
  mm <- lapply(paths, read_one_trial_metrics, delta_value = delta)
  mm <- Filter(Negate(is.null), mm)
  if (length(mm) > 0) metrics_list[[dtag]] <- do.call(rbind, mm)
}

metrics_all <- do.call(rbind, metrics_list)
metrics_all <- metrics_all[order(metrics_all$delta, metrics_all$trial), ]

safe_write_csv(metrics_all, file.path(REPORT_DIR, "metrics_all_deltas.csv"))

# ---- summary by delta (mean ± sd across trials)
num_cols <- setdiff(names(metrics_all), c("delta","trial","trial_path"))
summary_by_delta <- do.call(rbind, lapply(split(metrics_all, metrics_all$delta), function(df) {
  data.frame(
    delta = unique(df$delta),
    metric = num_cols,
    mean = sapply(df[, num_cols, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
    sd   = sapply(df[, num_cols, drop = FALSE], function(z) sd(z, na.rm = TRUE))
  )
}))
safe_write_csv(summary_by_delta, file.path(REPORT_DIR, "metrics_summary_by_delta.csv"))

# ============================================================
# 11) WELL-POSEDNESS / STABILITY PLOTS (delta vs estimator changes)
# ============================================================

# Load baseline estimators for each trial (delta=0)
read_estimators <- function(path) {
  b <- readRDS(path)
  tid <- if ("trial" %in% names(b)) as.integer(b$trial) else trial_id_from_path(path)
  list(
    trial = tid,
    obs_indices = as.integer(b$obs_indices),
    y = if ("observations" %in% names(b)) as.numeric(b$observations) else if ("observations_base" %in% names(b)) as.numeric(b$observations_base) else NULL,
    u_mean = as.numeric(b$u_post_mean),
    u_map  = as.numeric(b$u_map_raw),
    p_mean = as.numeric(b$p_post_mean),
    p_map  = as.numeric(b$p_map_raw)
  )
}

base_paths <- list_trial_bundles(OUTDIR)
base_est <- lapply(base_paths, read_estimators)
base_est <- base_est[!sapply(base_est, is.null)]
base_by_trial <- setNames(base_est, sapply(base_est, function(x) as.character(x$trial)))

# For delta>0, compute change w.r.t baseline
change_rows <- list()

for (delta in DELTAS[DELTAS > 0]) {
  dtag <- delta_tag(delta)
  DELTA_DIR <- file.path(PERT_BASEDIR, dtag)
  if (!dir.exists(DELTA_DIR)) next
  
  paths_d <- list_trial_bundles(DELTA_DIR)
  for (p in paths_d) {
    bd <- read_estimators(p)
    tid <- bd$trial
    if (!as.character(tid) %in% names(base_by_trial)) next
    b0 <- base_by_trial[[as.character(tid)]]
    
    # data perturb size
    # use y stored in delta-bundle: observations (perturbed) and observations_base
    bb <- readRDS(p)
    y0 <- if (!is.null(bb$observations_base)) as.numeric(bb$observations_base) else b0$y
    yD <- if (!is.null(bb$observations)) as.numeric(bb$observations) else NULL
    data_rel <- NA_real_
    if (!is.null(y0) && !is.null(yD) && length(y0) == length(yD)) {
      data_rel <- sqrt(sum((yD - y0)^2)) / max(sqrt(sum(y0^2)), 1e-12)
    }
    
    # estimator changes in L2(Gamma) using mass matrix
    umean_rel <- rel_l2_gamma(bd$u_mean, b0$u_mean, C)
    umap_rel  <- rel_l2_gamma(bd$u_map,  b0$u_map,  C)
    pmean_rel <- rel_l2_gamma(bd$p_mean, b0$p_mean, C)
    pmap_rel  <- rel_l2_gamma(bd$p_map,  b0$p_map,  C)
    
    change_rows[[length(change_rows)+1]] <- data.frame(
      delta = delta,
      trial = tid,
      data_rel_perturb = data_rel,
      rel_change_u_mean_L2G = umean_rel,
      rel_change_u_map_L2G  = umap_rel,
      rel_change_p_mean_L2G = pmean_rel,
      rel_change_p_map_L2G  = pmap_rel
    )
  }
}

change_df <- if (length(change_rows) > 0) do.call(rbind, change_rows) else data.frame()
safe_write_csv(change_df, file.path(REPORT_DIR, "wellposedness_change_metrics.csv"))

# ---- Plot 1: boxplot of estimator changes vs delta
plot_box_delta <- function(df, ycol, title, ylab, outbase) {
  if (!nrow(df)) return(NULL)
  p <- ggplot(df, aes(x = factor(delta), y = .data[[ycol]])) +
    geom_boxplot(outlier.size = 0.8) +
    theme_bw() +
    labs(title = title, x = expression(delta), y = ylab)
  print(p)
  ggsave(file.path(REPORT_DIR, paste0(outbase, ".png")), p, width = 6.5, height = 4.2, dpi = 220)
  ggsave(file.path(REPORT_DIR, paste0(outbase, ".pdf")), p, width = 6.5, height = 4.2)
  p
}

plot_box_delta(change_df, "rel_change_u_mean_L2G",
               "Well-posedness: ||u_mean(delta) - u_mean(0)|| / ||u_mean(0)|| in L2(Γ)",
               "relative change (L2(Γ))",
               "fig_wellposed_u_mean_change_boxplot")

plot_box_delta(change_df, "rel_change_u_map_L2G",
               "Well-posedness: ||u_MAP(delta) - u_MAP(0)|| / ||u_MAP(0)|| in L2(Γ)",
               "relative change (L2(Γ))",
               "fig_wellposed_u_map_change_boxplot")

plot_box_delta(change_df, "rel_change_p_mean_L2G",
               "Well-posedness: ||p_mean(delta) - p_mean(0)|| / ||p_mean(0)|| in L2(Γ)",
               "relative change (L2(Γ))",
               "fig_wellposed_p_mean_change_boxplot")

plot_box_delta(change_df, "rel_change_p_map_L2G",
               "Well-posedness: ||p_MAP(delta) - p_MAP(0)|| / ||p_MAP(0)|| in L2(Γ)",
               "relative change (L2(Γ))",
               "fig_wellposed_p_map_change_boxplot")

# ---- Plot 2: scatter: estimator change vs data perturb size
if (nrow(change_df) > 0 && any(is.finite(change_df$data_rel_perturb))) {
  p_scatter <- ggplot(change_df, aes(x = data_rel_perturb, y = rel_change_u_mean_L2G, color = factor(delta))) +
    geom_point() +
    theme_bw() +
    labs(title = "Stability: estimator change vs data perturb size",
         x = "||y_delta - y_0|| / ||y_0|| (Euclidean)",
         y = "||u_mean(delta) - u_mean(0)|| / ||u_mean(0)|| (L2(Γ))",
         color = expression(delta))
  print(p_scatter)
  ggsave(file.path(REPORT_DIR, "fig_wellposed_scatter_u_mean.png"), p_scatter, width = 6.5, height = 4.2, dpi = 220)
  ggsave(file.path(REPORT_DIR, "fig_wellposed_scatter_u_mean.pdf"), p_scatter, width = 6.5, height = 4.2)
}

# ============================================================
# 12) Metric boxplots across deltas (the existing metrics)
# ============================================================
save_boxplot_metric_by_delta <- function(df, col_mean, col_map, outbase, title, ylab) {
  if (!(col_mean %in% names(df)) || !(col_map %in% names(df))) return(NULL)
  dd <- rbind(
    data.frame(delta = df$delta, estimator = "Posterior mean", value = df[[col_mean]]),
    data.frame(delta = df$delta, estimator = "Raw MAP",        value = df[[col_map]])
  )
  p <- ggplot(dd, aes(x = factor(delta), y = value)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~ estimator, nrow = 1) +
    theme_bw() +
    labs(title = title, x = expression(delta), y = ylab)
  print(p)
  ggsave(file.path(REPORT_DIR, paste0(outbase, ".png")), p, width = 8.5, height = 4.2, dpi = 220)
  ggsave(file.path(REPORT_DIR, paste0(outbase, ".pdf")), p, width = 8.5, height = 4.2)
  p
}

save_boxplot_metric_by_delta(metrics_all, "rmse_u_mean", "rmse_u_map_raw",
                             "fig_rmse_u_by_delta",
                             "RMSE(u) across deltas", "RMSE(u)")

save_boxplot_metric_by_delta(metrics_all, "relative_error_u_mean", "relative_error_u_map_raw",
                             "fig_relative_error_u_by_delta",
                             "relative_error_u = RMSE(u)/range(u_true) across deltas",
                             "RMSE(u)/range(u_true)")

save_boxplot_metric_by_delta(metrics_all, "rmse_p_mean_all", "rmse_p_map_raw_all",
                             "fig_rmse_p_all_by_delta",
                             "RMSE(p) (all nodes) across deltas", "RMSE(p)")

save_boxplot_metric_by_delta(metrics_all, "relL2G_a_mean", "relL2G_a_map_raw",
                             "fig_relL2G_a_by_delta",
                             "Relative L2(Γ) error of a=exp(u) across deltas",
                             "||a_est-a_true||/||a_true|| (L2(Γ))")

save_boxplot_metric_by_delta(metrics_all, "relL2G_p_mean_all", "relL2G_p_map_raw_all",
                             "fig_relL2G_p_by_delta",
                             "Relative L2(Γ) error of p across deltas",
                             "||p_est-p_true||/||p_true|| (L2(Γ))")

# ============================================================
# 13) Save report bundle
# ============================================================
report_bundle <- list(
  timestamp = Sys.time(),
  outdir = OUTDIR,
  report_dir = REPORT_DIR,
  deltas = DELTAS,
  metrics_all = metrics_all,
  summary_by_delta = summary_by_delta,
  change_df = change_df,
  config = study_config
)
safe_save_rds(report_bundle, file.path(REPORT_DIR, "report_bundle_perturbation.rds"))

cat("\nDONE.\n",
    "Report saved under:\n  ", REPORT_DIR, "\n",
    "Key files:\n",
    "  - metrics_all_deltas.csv\n",
    "  - metrics_summary_by_delta.csv\n",
    "  - wellposedness_change_metrics.csv\n",
    "  - fig_wellposed_*.(png/pdf)\n",
    "  - fig_*_by_delta.(png/pdf)\n",
    "  - report_bundle_perturbation.rds\n", sep = "")

