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
BASE_OUTDIR <- "/Users/wenwenli/inverse_mcmc_fractional_20260118_105241"
stopifnot(dir.exists(BASE_OUTDIR))
cat("Using BASE_OUTDIR =", BASE_OUTDIR, "\n")

# Use these baseline trials as the shared truths across all deltas
BASE_TRIAL_IDS <- 4:6  # trial_004, trial_005, trial_006

# New deltas to run (baseline delta=0 is read only)
DELTAS_NEW <- c(1.1, 1.3, 1.5)

# suffix to avoid overwriting old delta_* folders
RUN_SUFFIX <- "updated"


# MCMC settings for the perturb runs ( match the baseline )
n_iter  <- 80000
burn_in <- 30000
beta0   <- 0.3

# Seeds (keep init & MCMC seeds the SAME across deltas; perturb seed depends on delta)
SEED_INIT_BASE <- 3000L
SEED_MCMC_BASE <- 4000L
SEED_PERT_BASE <- 9000L

# Progress
inner_mcmc_pb <- TRUE
verbose_mcmc  <- TRUE

# Resume policy for delta>0 runs
RESUME <- TRUE           
OVERWRITE_EXISTING <- FALSE

# Mesh rebuild (match baseline)
H_MESH <- 0.2

# Fractional order fallback (read from global bundle if present)
FRAC_ORDER_FALLBACK <- 1.5

# write new delta folders directly under BASE_OUTDIR to match.
PERT_BASEDIR <- BASE_OUTDIR

# ----------------------------
# 1) Load global bundle (baseline)
# ----------------------------
GLOBAL_PATH <- file.path(BASE_OUTDIR, "global", "global_bundle.rds")
stopifnot(file.exists(GLOBAL_PATH))
global_bundle <- readRDS(GLOBAL_PATH)

meshV_saved <- as.matrix(global_bundle$meshV)
C <- as(global_bundle$C, "dgCMatrix")
G <- as(global_bundle$G, "dgCMatrix")
f <- as.numeric(global_bundle$f)

kappa_L <- if (!is.null(global_bundle$kappa_L)) as.numeric(global_bundle$kappa_L) else 1

# Prior params for sample_spde
alpha_prior <- if (!is.null(global_bundle$alpha_prior)) as.numeric(global_bundle$alpha_prior) else 2
range_prior <- if (!is.null(global_bundle$range_prior)) as.numeric(global_bundle$range_prior) else 3
sigma_prior <- if (!is.null(global_bundle$sigma_prior)) as.numeric(global_bundle$sigma_prior) else 0.05

# Noise model params (NO delta here)
alpha_noise <- if (!is.null(global_bundle$alpha_noise)) as.numeric(global_bundle$alpha_noise) else 0.03
beta_noise  <- if (!is.null(global_bundle$beta_noise))  as.numeric(global_bundle$beta_noise)  else 0.2

# Fractional settings
use_fractional <- if (!is.null(global_bundle$use_fractional)) isTRUE(global_bundle$use_fractional) else TRUE
frac_order <- if (!is.null(global_bundle$frac_order)) as.numeric(global_bundle$frac_order) else FRAC_ORDER_FALLBACK

cat("Loaded global bundle.\n",
    "  n_mesh =", nrow(meshV_saved), "\n",
    "  prior(alpha,range,sigma) =", alpha_prior, range_prior, sigma_prior, "\n",
    "  noise(alpha,beta) =", alpha_noise, beta_noise, "\n",
    "  use_fractional =", use_fractional, " frac_order =", frac_order, "\n")

# ----------------------------
# 2) Locate baseline trial bundles (delta=0)
# ----------------------------
trial_dir_name <- function(tid) sprintf("trial_%03d", tid)
trial_bundle_path0 <- function(tid) file.path(BASE_OUTDIR, trial_dir_name(tid), "trial_bundle.rds")

stopifnot(all(file.exists(sapply(BASE_TRIAL_IDS, trial_bundle_path0))))
cat("Found baseline trial bundles for:", paste(BASE_TRIAL_IDS, collapse = ","), "\n")

# ----------------------------
# 3) Rebuild graph for sample_spde ( match saved mesh ordering)
# ----------------------------
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance = list(edge_edge = 1e-3, vertex_vertex = 1e-3, edge_vertex = 1e-3)
)
graph$build_mesh(h = H_MESH)
graph$compute_fem()

meshV_re <- as.matrix(graph$mesh$V)
if (nrow(meshV_re) != nrow(meshV_saved)) {
  stop("Rebuilt mesh size != saved mesh size.\n",
       "Saved n_mesh=", nrow(meshV_saved), ", rebuilt n_mesh=", nrow(meshV_re), "\n",
       "Fix: set H_MESH to the baseline mesh size you used, or save/load the graph object.")
}

coord_diff <- max(abs(meshV_re[,1:2] - meshV_saved[,1:2]))
if (!is.finite(coord_diff) || coord_diff > 1e-6) {
  warning("Rebuilt mesh coordinates differ from saved mesh (max abs diff = ", coord_diff, ").\n",
          "If you suspect ordering differs, STOP and instead save/load the original graph object.")
} else {
  cat("Mesh coordinate check passed (max abs diff ~", coord_diff, ").\n")
}

n_nodes <- nrow(meshV_saved)

# ----------------------------
# 4) Forward solve: fractional PDE
# ----------------------------
safe_exp <- function(u, clip = 20) exp(pmax(pmin(as.numeric(u), clip), -clip))

assemble_L <- function(u) {
  u <- as.numeric(u)
  if (length(u) != n_nodes) u <- u[seq_len(n_nodes)]
  a  <- safe_exp(u, clip = 20)
  Du <- Matrix::Diagonal(n_nodes, a)
  Gu <- (G %*% Du + Du %*% G) * 0.5
  L  <- kappa_L^2 * C + Gu
  (L + t(L)) * 0.5
}

solve_pde_fractional <- function(u, frac_order) {
  L <- assemble_L(u)
  tryCatch({
    Ld  <- as.matrix(L)  # dense eig (baseline choice)
    eig <- eigen(Ld, symmetric = TRUE)
    
    lam_raw <- eig$values
    lam_floor <- max(1e-12, 1e-12 * max(abs(lam_raw)))
    lam <- pmax(lam_raw, lam_floor)
    
    vals <- lam^(-frac_order)
    V <- eig$vectors
    as.numeric(V %*% (vals * crossprod(V, f)))
  }, error = function(e) rep(NA_real_, n_nodes))
}

solve_pde <- function(u) {
  if (isTRUE(use_fractional)) solve_pde_fractional(u, frac_order) else {
    L <- assemble_L(u)
    tryCatch({
      chol <- Cholesky(L, LDL = FALSE, Imult = 0)
      as.numeric(solve(chol, f))
    }, error = function(e) rep(NA_real_, n_nodes))
  }
}

# ----------------------------
# 5) Helpers: metrics / norms
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

l2_gamma_norm <- function(v, Cmat) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  v2 <- v; v2[!ok] <- 0
  sqrt(as.numeric(Matrix::crossprod(v2, Cmat %*% v2)))
}

# embed an observation vector (length m) to full mesh vector (length n_nodes)
embed_obs_to_mesh <- function(obs_vec, obs_idx, n_nodes) {
  out <- rep(0, n_nodes)
  obs_idx <- as.integer(obs_idx)
  ok <- is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= n_nodes
  obs_idx <- obs_idx[ok]
  if (length(obs_idx) == 0) return(out)
  if (length(obs_vec) != length(obs_idx)) {
    stop("embed_obs_to_mesh: length mismatch: obs_vec vs obs_idx")
  }
  out[obs_idx] <- as.numeric(obs_vec)
  out
}

# ----------------------------
# 6) pCN MCMC (likelihood WITHOUT delta)
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
  
  for (i in 1:n_iter) {
    T_current <- max(1, T0 * cooling_factor^(i / adapt_interval))
    
    phi_current <- neg_log_likelihood_fn(u_current, T = T_current)
    
    xi <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_prior,
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
        beta <- min(beta * 1.2, 0.999)
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

run_one_trial_with_observations <- function(u_true, p_true, obs_indices, observations,
                                            n_iter, burn_in, beta0,
                                            seed_init, seed_mcmc,
                                            inner_mcmc_pb = TRUE, verbose_mcmc = TRUE) {
  
  obs_idx <- as.integer(obs_indices)
  obs_idx <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= length(p_true)]
  stopifnot(length(obs_idx) > 0)
  
  # Likelihood: sd(p_obs) = alpha_noise*|p_obs| + beta_noise  (NO delta)
  neg_log_likelihood <- function(u, T = 1) {
    p_est <- solve_pde(u)
    if (any(!is.finite(p_est))) return(1e30)
    p_obs <- p_est[obs_idx]
    sd_vec <- alpha_noise * abs(p_obs) + beta_noise
    if (any(!is.finite(sd_vec)) || any(sd_vec <= 0)) return(1e30)
    like <- sum(log(sd_vec)) + 0.5 * sum(((p_obs - observations) / sd_vec)^2)
    like / T
  }
  
  set.seed(seed_init)
  u_init <- as.numeric(sample_spde(range = range_prior, sigma = sigma_prior, alpha = alpha_prior,
                                   graph = graph, type = "mesh"))
  if (length(u_init) != length(u_true)) u_init <- u_init[1:length(u_true)]
  
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
  
  # estimators u
  u_post_mean <- colMeans(posterior_samples_post)
  u_map_raw   <- apply(posterior_samples_post, 2, get_mode)
  
  # plug-in p fields
  p_post_mean <- solve_pde(u_post_mean)
  p_map_raw   <- solve_pde(u_map_raw)
  
  sd_u <- apply(posterior_samples_post, 2, sd)
  sd_u_mean <- mean(sd_u, na.rm = TRUE)
  sd_u_max  <- max(sd_u, na.rm = TRUE)
  
  list(
    u_post_mean = u_post_mean,
    u_map_raw = u_map_raw,
    p_post_mean = p_post_mean,
    p_map_raw = p_map_raw,
    acc_rate_all = acc_rate_all,
    beta_final = beta_final,
    sd_u_mean = sd_u_mean,
    sd_u_max  = sd_u_max
  )
}

# ----------------------------
# 7) Folder helpers
# ----------------------------
delta_tag <- function(delta, suffix = NULL) {
  base <- paste0("delta_", gsub("\\.", "p", sprintf("%.1f", delta)))
  if (!is.null(suffix) && nzchar(suffix)) paste0(base, "_", suffix) else base
}

trial_dir_delta <- function(delta_dir, tid) file.path(delta_dir, sprintf("trial_%03d", tid))
trial_bundle_path_delta <- function(delta_dir, tid) file.path(trial_dir_delta(delta_dir, tid), "trial_bundle.rds")

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

# ----------------------------
# 8) Run delta > 0 perturbations (READ baseline truth/data)
# ----------------------------
cat("\n=== Running perturbation deltas (fractional) ===\n")
for (delta in DELTAS_NEW) {
  dtag <- delta_tag(delta, RUN_SUFFIX)
  DELTA_DIR <- file.path(PERT_BASEDIR, dtag)
  dir.create(DELTA_DIR, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n-----------------------------\n")
  cat("delta =", delta, " -> ", DELTA_DIR, "\n")
  cat("-----------------------------\n")
  
  pb <- utils::txtProgressBar(min = 0, max = length(BASE_TRIAL_IDS), style = 3)
  
  for (k in seq_along(BASE_TRIAL_IDS)) {
    tid <- BASE_TRIAL_IDS[k]
    
    out_path <- trial_bundle_path_delta(DELTA_DIR, tid)
    if (RESUME && file.exists(out_path) && !OVERWRITE_EXISTING) {
      message("delta=", delta, " trial=", tid, " already exists in ", dtag, "; skipping.")
      utils::setTxtProgressBar(pb, k)
      next
    }
    
    b0 <- readRDS(trial_bundle_path0(tid))
    
    u_true <- as.numeric(b0$u_true)
    p_true <- as.numeric(b0$p_true)
    obs_idx <- as.integer(b0$obs_indices)
    
    # baseline data y0
    y0 <- if (!is.null(b0$observations)) as.numeric(b0$observations) else stop("Baseline bundle missing observations.")
    
    obs_idx2 <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= length(p_true)]
    if (length(obs_idx2) == 0) stop("No valid obs indices in baseline trial ", tid)
    
    if (length(y0) != length(obs_idx2)) {
      stop("Baseline y0 length mismatch with obs_idx2 in trial ", tid,
           " | length(y0)=", length(y0), " length(obs_idx2)=", length(obs_idx2))
    }
    
    # sigma_true based on p_true at obs nodes (NO delta in likelihood)
    sigma_true <- alpha_noise * abs(p_true[obs_idx2]) + beta_noise
    
    # perturbation: y_delta = y0 + delta * sigma_true * z
    set.seed(SEED_PERT_BASE + 10000L * tid + as.integer(round(delta * 1000)))
    z <- rnorm(length(obs_idx2))
    y_delta <- y0 + delta * sigma_true * z
    
    # Run inference (likelihood unchanged)
    res <- run_one_trial_with_observations(
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
    
    # Save
    td <- trial_dir_delta(DELTA_DIR, tid)
    dir.create(td, showWarnings = FALSE, recursive = TRUE)
    
    tb <- list(
      trial = tid,
      delta = delta,
      run_suffix = RUN_SUFFIX,
      timestamp = Sys.time(),
      
      # shared truth/data from baseline
      u_true = u_true,
      p_true = p_true,
      obs_indices = obs_idx2,
      observations_base = y0,
      observations = y_delta,
      
      # perturbation bookkeeping
      perturb_seed = SEED_PERT_BASE + 10000L * tid + as.integer(round(delta * 1000)),
      perturb_z = z,
      sigma_true_obs = sigma_true,
      
      # estimators
      u_post_mean = as.numeric(res$u_post_mean),
      u_map_raw   = as.numeric(res$u_map_raw),
      p_post_mean = as.numeric(res$p_post_mean),
      p_map_raw   = as.numeric(res$p_map_raw),
      
      # diagnostics
      acc_rate_all = as.numeric(res$acc_rate_all),
      beta_final   = as.numeric(res$beta_final),
      sd_u_mean    = as.numeric(res$sd_u_mean),
      sd_u_max     = as.numeric(res$sd_u_max)
    )
    
    safe_save_rds(tb, out_path)
    
    utils::setTxtProgressBar(pb, k)
  }
  close(pb)
}

# ============================================================
# 9) REPORT STAGE: baseline (delta=0) + new deltas (from *_updated folders)
# ============================================================
REPORT_DIR <- file.path(
  BASE_OUTDIR,
  paste0("report_stage_fractional_perturb_", RUN_SUFFIX, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("\nSaving report to:\n  ", REPORT_DIR, "\n")

# ---- Read baseline estimators (delta=0) for the chosen trials
read_baseline_trial <- function(tid) {
  b <- readRDS(trial_bundle_path0(tid))
  list(
    trial = tid,
    u_true = as.numeric(b$u_true),
    p_true = as.numeric(b$p_true),
    obs_indices = as.integer(b$obs_indices),
    y0 = as.numeric(b$observations),
    u_mean0 = as.numeric(b$u_post_mean),
    u_map0  = as.numeric(b$u_map_raw),
    p_mean0 = as.numeric(b$p_post_mean),
    p_map0  = as.numeric(b$p_map_raw)
  )
}
base_list <- lapply(BASE_TRIAL_IDS, read_baseline_trial)
base_by_trial <- setNames(base_list, sapply(base_list, function(x) as.character(x$trial)))

# ---- Read delta trial bundle (delta>0) FROM *_updated directory
read_delta_trial <- function(delta, tid) {
  dtag <- delta_tag(delta, RUN_SUFFIX)
  path <- trial_bundle_path_delta(file.path(PERT_BASEDIR, dtag), tid)
  if (!file.exists(path)) return(NULL)
  b <- readRDS(path)
  list(
    trial = tid,
    delta = delta,
    obs_indices = as.integer(b$obs_indices),
    y0 = as.numeric(b$observations_base),
    yD = as.numeric(b$observations),
    sigma_true_obs = as.numeric(b$sigma_true_obs),
    
    u_mean = as.numeric(b$u_post_mean),
    u_map  = as.numeric(b$u_map_raw),
    p_mean = as.numeric(b$p_post_mean),
    p_map  = as.numeric(b$p_map_raw)
  )
}

# ---- Compute well-posedness change metrics
one_change_row <- function(delta, tid) {
  b0 <- base_by_trial[[as.character(tid)]]
  bd <- read_delta_trial(delta, tid)
  if (is.null(bd)) return(NULL)
  
  obs_idx <- b0$obs_indices
  obs_idx <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= n_nodes]
  if (length(obs_idx) == 0) return(NULL)
  
  # data perturb
  y0 <- bd$y0
  yD <- bd$yD
  if (length(y0) != length(yD)) return(NULL)
  
  dy <- yD - y0
  data_abs_eucl <- sqrt(sum(dy^2))
  data_rel_eucl <- data_abs_eucl / max(sqrt(sum(y0^2)), 1e-12)
  
  dy_full <- embed_obs_to_mesh(dy, obs_idx, n_nodes)
  y0_full <- embed_obs_to_mesh(y0, obs_idx, n_nodes)
  
  data_abs_L2G <- l2_gamma_norm(dy_full, C)
  data_rel_L2G <- data_abs_L2G / max(l2_gamma_norm(y0_full, C), 1e-12)
  
  # sigma_true L2G (for A_* metrics)
  sig <- bd$sigma_true_obs
  sig_full <- embed_obs_to_mesh(sig, obs_idx, n_nodes)
  sig_L2G <- l2_gamma_norm(sig_full, C)
  denom_A <- max(delta * sig_L2G, 1e-12)
  
  # a-fields
  a_mean0 <- safe_exp(b0$u_mean0)
  a_map0  <- safe_exp(b0$u_map0)
  a_mean  <- safe_exp(bd$u_mean)
  a_map   <- safe_exp(bd$u_map)
  
  # abs changes in L2(Gamma)
  abs_change_u_mean_L2G <- l2_gamma_norm(bd$u_mean - b0$u_mean0, C)
  abs_change_u_map_L2G  <- l2_gamma_norm(bd$u_map  - b0$u_map0,  C)
  
  abs_change_p_mean_L2G <- l2_gamma_norm(bd$p_mean - b0$p_mean0, C)
  abs_change_p_map_L2G  <- l2_gamma_norm(bd$p_map  - b0$p_map0,  C)
  
  abs_change_a_mean_L2G <- l2_gamma_norm(a_mean - a_mean0, C)
  abs_change_a_map_L2G  <- l2_gamma_norm(a_map  - a_map0,  C)
  
  # relative changes in L2(Gamma)
  rel_change_u_mean_L2G <- abs_change_u_mean_L2G / max(l2_gamma_norm(b0$u_mean0, C), 1e-12)
  rel_change_u_map_L2G  <- abs_change_u_map_L2G  / max(l2_gamma_norm(b0$u_map0,  C), 1e-12)
  
  rel_change_p_mean_L2G <- abs_change_p_mean_L2G / max(l2_gamma_norm(b0$p_mean0, C), 1e-12)
  rel_change_p_map_L2G  <- abs_change_p_map_L2G  / max(l2_gamma_norm(b0$p_map0,  C), 1e-12)
  
  rel_change_a_mean_L2G <- abs_change_a_mean_L2G / max(l2_gamma_norm(a_mean0, C), 1e-12)
  rel_change_a_map_L2G  <- abs_change_a_map_L2G  / max(l2_gamma_norm(a_map0,  C), 1e-12)
  
  # RMSE changes (Euclidean) relative to baseline estimators
  rmse_change_u_mean <- rmse(bd$u_mean, b0$u_mean0)
  rmse_change_u_map  <- rmse(bd$u_map,  b0$u_map0)
  rmse_change_p_mean <- rmse(bd$p_mean, b0$p_mean0)
  rmse_change_p_map  <- rmse(bd$p_map,  b0$p_map0)
  
  # Amplification metrics A_* : abs_change / (delta * ||sigma_true||_{L2G})
  A_u_mean <- abs_change_u_mean_L2G / denom_A
  A_u_map  <- abs_change_u_map_L2G  / denom_A
  A_p_mean <- abs_change_p_mean_L2G / denom_A
  A_p_map  <- abs_change_p_map_L2G  / denom_A
  A_a_mean <- abs_change_a_mean_L2G / denom_A
  A_a_map  <- abs_change_a_map_L2G  / denom_A
  
  data.frame(
    delta = delta,
    trial = tid,
    
    data_abs_eucl = data_abs_eucl,
    data_rel_eucl = data_rel_eucl,
    data_abs_L2G  = data_abs_L2G,
    data_rel_L2G  = data_rel_L2G,
    
    abs_change_u_mean_L2G = abs_change_u_mean_L2G,
    abs_change_u_map_L2G  = abs_change_u_map_L2G,
    abs_change_p_mean_L2G = abs_change_p_mean_L2G,
    abs_change_p_map_L2G  = abs_change_p_map_L2G,
    abs_change_a_mean_L2G = abs_change_a_mean_L2G,
    abs_change_a_map_L2G  = abs_change_a_map_L2G,
    
    rel_change_u_mean_L2G = rel_change_u_mean_L2G,
    rel_change_u_map_L2G  = rel_change_u_map_L2G,
    rel_change_p_mean_L2G = rel_change_p_mean_L2G,
    rel_change_p_map_L2G  = rel_change_p_map_L2G,
    rel_change_a_mean_L2G = rel_change_a_mean_L2G,
    rel_change_a_map_L2G  = rel_change_a_map_L2G,
    
    rmse_change_u_mean = rmse_change_u_mean,
    rmse_change_u_map  = rmse_change_u_map,
    rmse_change_p_mean = rmse_change_p_mean,
    rmse_change_p_map  = rmse_change_p_map,
    
    A_u_mean = A_u_mean,
    A_u_map  = A_u_map,
    A_p_mean = A_p_mean,
    A_p_map  = A_p_map,
    A_a_mean = A_a_mean,
    A_a_map  = A_a_map
  )
}

# Build change_df for all deltas and all trials
change_rows <- list()
for (delta in DELTAS_NEW) {
  for (tid in BASE_TRIAL_IDS) {
    rr <- one_change_row(delta, tid)
    if (!is.null(rr)) change_rows[[length(change_rows) + 1]] <- rr
  }
}
change_df <- if (length(change_rows) > 0) do.call(rbind, change_rows) else data.frame()

safe_write_csv(change_df, file.path(REPORT_DIR, "wellposedness_change_metrics_trials.csv"))

# Summary "delta, metric, mean, sd"
if (nrow(change_df) > 0) {
  num_cols <- setdiff(names(change_df), c("delta", "trial"))
  summary_change <- do.call(rbind, lapply(split(change_df, change_df$delta), function(df) {
    data.frame(
      delta = unique(df$delta),
      metric = num_cols,
      mean = sapply(df[, num_cols, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
      sd   = sapply(df[, num_cols, drop = FALSE], function(z) sd(z, na.rm = TRUE))
    )
  }))
  safe_write_csv(summary_change, file.path(REPORT_DIR, "wellposedness_change_metrics_summary_by_delta.csv"))
}

cat("\nDONE.\n",
    "New delta folders written under:\n  ", PERT_BASEDIR, "\n",
    "using suffix: ", RUN_SUFFIX, "\n\n",
    "Report folder:\n  ", REPORT_DIR, "\n\n",
    "Key files:\n",
    "  - wellposedness_change_metrics_trials.csv\n",
    "  - wellposedness_change_metrics_summary_by_delta.csv\n", sep = "")



# ============================================================
# 10) PLOT STAGE: boxplots + summary plots saved under REPORT_DIR
# ============================================================

# ----------------------------
# 10.0) Load change_df 
# ----------------------------
if (!exists("change_df") || !is.data.frame(change_df) || nrow(change_df) == 0) {
  csv_trials <- file.path(REPORT_DIR, "wellposedness_change_metrics_trials.csv")
  stopifnot(file.exists(csv_trials))
  change_df <- read.csv(csv_trials, stringsAsFactors = FALSE)
}

# Make sure delta is numeric and a factor for plotting order
change_df$delta <- as.numeric(change_df$delta)
change_df$delta_f <- factor(change_df$delta, levels = sort(unique(change_df$delta)))

# ----------------------------
# 10.1)  (21pt bold ticks/labels)
# ----------------------------
theme_paper_21 <- function() {
  ggplot2::theme_bw(base_size = 21) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 21, face = "bold"),
      axis.text  = ggplot2::element_text(size = 21, face = "bold"),
      strip.text = ggplot2::element_text(size = 21, face = "bold"),
      plot.title = ggplot2::element_text(size = 21, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(size = 21, face = "bold"),
      legend.text  = ggplot2::element_text(size = 21),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ----------------------------
# 10.2) Helpers: long-format without extra packages
# ----------------------------
to_long_df <- function(df, metrics, id_cols = c("delta", "delta_f", "trial")) {
  metrics <- metrics[metrics %in% names(df)]
  if (length(metrics) == 0) return(data.frame())
  mat <- as.matrix(df[, metrics, drop = FALSE])
  out <- data.frame(
    df[, id_cols, drop = FALSE],
    metric = rep(metrics, each = nrow(df)),
    value  = as.vector(mat),
    stringsAsFactors = FALSE
  )
  out$metric <- factor(out$metric, levels = metrics)
  out
}

save_png_pdf <- function(plot_obj, out_base, width = 14, height = 9, dpi = 300) {
  ggplot2::ggsave(paste0(out_base, ".png"), plot_obj, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(out_base, ".pdf"), plot_obj, width = width, height = height)
}

plot_box_facets <- function(df_long, out_base, ncol = 2) {
  if (!is.data.frame(df_long) || nrow(df_long) == 0) return(invisible(NULL))
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = delta_f, y = value)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 2) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = "delta", y = "") +
    theme_paper_21()
  save_png_pdf(p, out_base, width = 16, height = 10)
  invisible(p)
}

plot_mean_sd_facets <- function(df, metrics, out_base, ncol = 2) {
  metrics <- metrics[metrics %in% names(df)]
  if (length(metrics) == 0) return(invisible(NULL))
  # summarize mean/sd by delta
  agg_list <- lapply(metrics, function(m) {
    dd <- df[, c("delta", "delta_f", m)]
    names(dd)[3] <- "value"
    dd <- dd[is.finite(dd$value), , drop = FALSE]
    if (nrow(dd) == 0) return(NULL)
    mean_by <- tapply(dd$value, dd$delta, mean, na.rm = TRUE)
    sd_by   <- tapply(dd$value, dd$delta, sd,   na.rm = TRUE)
    out <- data.frame(
      delta = as.numeric(names(mean_by)),
      mean  = as.numeric(mean_by),
      sd    = as.numeric(sd_by),
      metric = m,
      stringsAsFactors = FALSE
    )
    out
  })
  agg <- do.call(rbind, Filter(Negate(is.null), agg_list))
  if (is.null(agg) || nrow(agg) == 0) return(invisible(NULL))
  
  agg$metric <- factor(agg$metric, levels = metrics)
  agg$delta_f <- factor(agg$delta, levels = sort(unique(df$delta)))
  
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = delta, y = mean)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), width = 0.02) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = "delta", y = "") +
    theme_paper_21()
  save_png_pdf(p, out_base, width = 16, height = 10)
  invisible(p)
}

# ----------------------------
# 10.3) Choose metric groups 
# ----------------------------
metrics_data <- c("data_abs_eucl", "data_rel_eucl", "data_abs_L2G", "data_rel_L2G")

metrics_abs_change <- c(
  "abs_change_u_mean_L2G", "abs_change_u_map_L2G",
  "abs_change_p_mean_L2G", "abs_change_p_map_L2G",
  "abs_change_a_mean_L2G", "abs_change_a_map_L2G"
)

metrics_rel_change <- c(
  "rel_change_u_mean_L2G", "rel_change_u_map_L2G",
  "rel_change_p_mean_L2G", "rel_change_p_map_L2G",
  "rel_change_a_mean_L2G", "rel_change_a_map_L2G"
)

metrics_rmse_change <- c("rmse_change_u_mean", "rmse_change_u_map", "rmse_change_p_mean", "rmse_change_p_map")

metrics_A <- c("A_u_mean", "A_u_map", "A_p_mean", "A_p_map", "A_a_mean", "A_a_map")

# ----------------------------
# 10.4) Output folder
# ----------------------------
PLOTS_DIR <- file.path(REPORT_DIR, "plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 10.5) Faceted boxplots (main deliverables)
# ----------------------------
df_data_long <- to_long_df(change_df, metrics_data)
plot_box_facets(df_data_long, file.path(PLOTS_DIR, "box_data_metrics"), ncol = 2)

df_abs_long <- to_long_df(change_df, metrics_abs_change)
plot_box_facets(df_abs_long, file.path(PLOTS_DIR, "box_abs_change_L2G"), ncol = 2)

df_rel_long <- to_long_df(change_df, metrics_rel_change)
plot_box_facets(df_rel_long, file.path(PLOTS_DIR, "box_rel_change_L2G"), ncol = 2)

df_rmse_long <- to_long_df(change_df, metrics_rmse_change)
plot_box_facets(df_rmse_long, file.path(PLOTS_DIR, "box_rmse_change"), ncol = 2)

df_A_long <- to_long_df(change_df, metrics_A)
plot_box_facets(df_A_long, file.path(PLOTS_DIR, "box_amplification_A"), ncol = 2)

# ----------------------------
# 10.6) Optional: mean ± sd trend plots vs delta
# ----------------------------
plot_mean_sd_facets(change_df, metrics_data,       file.path(PLOTS_DIR, "trend_mean_sd_data"), ncol = 2)
plot_mean_sd_facets(change_df, metrics_abs_change, file.path(PLOTS_DIR, "trend_mean_sd_abs_change_L2G"), ncol = 2)
plot_mean_sd_facets(change_df, metrics_rel_change, file.path(PLOTS_DIR, "trend_mean_sd_rel_change_L2G"), ncol = 2)
plot_mean_sd_facets(change_df, metrics_rmse_change,file.path(PLOTS_DIR, "trend_mean_sd_rmse_change"), ncol = 2)
plot_mean_sd_facets(change_df, metrics_A,          file.path(PLOTS_DIR, "trend_mean_sd_A"), ncol = 2)

cat("\nPLOTS DONE.\nSaved under:\n  ", PLOTS_DIR, "\n", sep = "")

