suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ----------------------------
# 0) Set OUTDIR
# ----------------------------
OUTDIR <- "/Users/wenwenli/inverse_mcmc_metricgraph_20260121_144015"
stopifnot(dir.exists(OUTDIR))
cat("Using OUTDIR =", OUTDIR, "\n")

# ----------------------------
# 1) Locate saved bundles
# ----------------------------
GLOBAL_PATH <- file.path(OUTDIR, "global", "global_bundle.rds")
if (!file.exists(GLOBAL_PATH)) {
  stop("Cannot find global_bundle.rds at: ", GLOBAL_PATH)
}
global_bundle <- readRDS(GLOBAL_PATH)

if (!("C" %in% names(global_bundle))) {
  stop("global_bundle does not contain C (mass matrix). Cannot compute L2(Gamma) errors.")
}
C_numeric <- as(global_bundle$C, "dgCMatrix")

# ----scan trial bundles
trial_dirs_all <- list.dirs(OUTDIR, full.names = TRUE, recursive = FALSE)
trial_dirs <- trial_dirs_all[grepl("^trial_[0-9]{3}$", basename(trial_dirs_all))]

trial_paths <- file.path(trial_dirs, "trial_bundle.rds")
trial_paths <- trial_paths[file.exists(trial_paths)]

if (length(trial_paths) == 0) {
  stop("No trial_bundle.rds found under OUTDIR/trial_###/.")
}
cat("Found ", length(trial_paths), " trial_bundle.rds files.\n")

trial_id_from_path <- function(path) {
  d <- basename(dirname(path))  
  as.integer(sub("^trial_([0-9]{3})$", "\\1", d))
}

trial_ids <- sapply(trial_paths, trial_id_from_path)
ord <- order(trial_ids)
trial_paths <- trial_paths[ord]
trial_ids   <- trial_ids[ord]

# ----------------------------
# 2) Helper functions
# ----------------------------
rmse <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((est[ok] - truth[ok])^2))
}

nrmse_range <- function(est, truth, eps = 1e-12) {
  truth <- as.numeric(truth)
  ok <- is.finite(truth)
  rng <- diff(range(truth[ok], na.rm = TRUE))
  rmse(est, truth) / max(rng, eps)
}

exp_clip <- function(u, clip = 20) {
  u <- as.numeric(u)
  exp(pmin(pmax(u, -clip), clip))
}

l2_gamma_norm <- function(v, Cmat) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  v2 <- v
  v2[!ok] <- 0
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

mean_sd <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
fmt_mean_sd <- function(m, s, digits = 4) {
  if (!is.finite(m) || !is.finite(s)) return("")
  sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), m, s)
}

# ----------------------------
# 2.1 Plot style (21pt + bold)
# ----------------------------
theme_paper_21_bold <- function() {
  theme_bw(base_size = 21) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_text(size = 21, face = "bold"),
      axis.title.y = element_text(size = 21, face = "bold"),
      axis.text.x  = element_text(size = 21, face = "bold"),
      axis.text.y  = element_text(size = 21, face = "bold"),
      axis.ticks = element_line(linewidth = 0.6),
      axis.line  = element_line(linewidth = 0.6),
      panel.grid.major = element_line(linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

# ----------------------------
# 2.2 Boxplot savers
# ----------------------------
save_boxplot_2 <- function(df_long, outpath_png, outpath_pdf, ylab_parse) {
  # enforce ordering mean then mode for each variable type (if present)
  if (all(c("hat(u)[plain(mean)]","hat(u)[plain(mode)]") %in% unique(df_long$estimator))) {
    df_long$estimator <- factor(df_long$estimator,
                                levels = c("hat(u)[plain(mean)]","hat(u)[plain(mode)]"))
  }
  if (all(c("hat(p)[plain(mean)]","hat(p)[plain(mode)]") %in% unique(df_long$estimator))) {
    df_long$estimator <- factor(df_long$estimator,
                                levels = c("hat(p)[plain(mean)]","hat(p)[plain(mode)]"))
  }
  if (all(c("hat(a)[plain(mean)]","hat(a)[plain(mode)]") %in% unique(df_long$estimator))) {
    df_long$estimator <- factor(df_long$estimator,
                                levels = c("hat(a)[plain(mean)]","hat(a)[plain(mode)]"))
  }
  
  p <- ggplot(df_long, aes(x = estimator, y = value)) +
    geom_boxplot(outlier.size = 1.0, linewidth = 0.6) +
    scale_x_discrete(labels = scales::label_parse()) +
    labs(x = "", y = parse(text = ylab_parse)[[1]]) +
    theme_paper_21_bold()
  
  ggsave(outpath_png, p, width = 6.5, height = 4.2, dpi = 220)
  ggsave(outpath_pdf, p, width = 6.5, height = 4.2)
  invisible(p)
}

save_boxplot_1 <- function(df_one, outpath_png, outpath_pdf, ylab_parse) {
  p <- ggplot(df_one, aes(x = " ", y = value)) +
    geom_boxplot(outlier.size = 1.0, linewidth = 0.6) +
    labs(x = "", y = parse(text = ylab_parse)[[1]]) +
    theme_paper_21_bold()
  
  ggsave(outpath_png, p, width = 6.5, height = 4.2, dpi = 220)
  ggsave(outpath_pdf, p, width = 6.5, height = 4.2)
  invisible(p)
}

# ----------------------------
# 3) Read each trial bundle and compute metrics
# ----------------------------
read_one_trial_metrics <- function(path, C_numeric) {
  b <- readRDS(path)
  
  required <- c("u_true","p_true","u_post_mean","u_map_raw","p_post_mean","p_map_raw","obs_indices")
  miss <- setdiff(required, names(b))
  if (length(miss) > 0) {
    warning("Skipping (missing fields): ", path, " | missing: ", paste(miss, collapse=", "))
    return(NULL)
  }
  
  tid <- if ("trial" %in% names(b) && is.finite(b$trial)) as.integer(b$trial) else trial_id_from_path(path)
  
  u_true <- as.numeric(b$u_true)
  p_true <- as.numeric(b$p_true)
  
  u_mean <- as.numeric(b$u_post_mean) # \hat u_mean
  u_mode <- as.numeric(b$u_map_raw)   # \hat u_mode (marginal mode)
  
  p_mean <- as.numeric(b$p_post_mean) # \hat p_mean = p(\hat u_mean)
  p_mode <- as.numeric(b$p_map_raw)   # \hat p_mode = p(\hat u_mode)
  
  obs_idx <- as.integer(b$obs_indices)
  obs_idx <- obs_idx[is.finite(obs_idx) & obs_idx >= 1 & obs_idx <= length(p_true)]
  
  # a = exp(u)
  a_true <- exp_clip(u_true)
  a_mean <- exp_clip(u_mean)
  a_mode <- exp_clip(u_mode)
  
  out <- data.frame(
    trial = tid,
    trial_path = path,
    
    # diagnostics 
    acc_rate_all = if ("acc_rate_all" %in% names(b)) as.numeric(b$acc_rate_all) else NA_real_,
    beta_final   = if ("beta_final"   %in% names(b)) as.numeric(b$beta_final)   else NA_real_,
    
    # saved posterior spread summaries (if present)
    sd_u_mean = if ("sd_u_mean" %in% names(b)) as.numeric(b$sd_u_mean) else NA_real_,
    sd_u_max  = if ("sd_u_max"  %in% names(b)) as.numeric(b$sd_u_max)  else NA_real_,
    
    # RMSE(u)
    rmse_u_mean = rmse(u_mean, u_true),
    rmse_u_mode = rmse(u_mode, u_true),
    
    # RMSE(a)
    rmse_a_mean = rmse(a_mean, a_true),
    rmse_a_mode = rmse(a_mode, a_true),
    
    # NRMSE_range(u)
    nrmse_u_range_mean = nrmse_range(u_mean, u_true),
    nrmse_u_range_mode = nrmse_range(u_mode, u_true),
    
    # RMSE(p) all + obs
    rmse_p_mean_all = rmse(p_mean, p_true),
    rmse_p_mode_all = rmse(p_mode, p_true),
    
    rmse_p_mean_obs = if (length(obs_idx) > 0) rmse(p_mean[obs_idx], p_true[obs_idx]) else NA_real_,
    rmse_p_mode_obs = if (length(obs_idx) > 0) rmse(p_mode[obs_idx], p_true[obs_idx]) else NA_real_,
    
    # RelErr_{L2(Gamma)}(a)
    relL2G_a_mean = rel_l2_gamma(a_mean, a_true, C_numeric),
    relL2G_a_mode = rel_l2_gamma(a_mode, a_true, C_numeric),
    
    # RelErr_{L2(Gamma)}(p)
    relL2G_p_mean_all = rel_l2_gamma(p_mean, p_true, C_numeric),
    relL2G_p_mode_all = rel_l2_gamma(p_mode, p_true, C_numeric),
    
    # record truth scales
    u_range_true = diff(range(u_true, na.rm = TRUE)),
    p_range_true = diff(range(p_true, na.rm = TRUE))
  )
  
  out
}

metrics_list <- lapply(trial_paths, read_one_trial_metrics, C_numeric = C_numeric)
metrics_list <- Filter(Negate(is.null), metrics_list)
metrics_df <- do.call(rbind, metrics_list)

metrics_df <- metrics_df[order(metrics_df$trial), ]
cat("\nComputed metrics for ", nrow(metrics_df), " trials.\n")

# ----------------------------
# 4) Output dir 
# ----------------------------
REPORT_DIR <- file.path(OUTDIR, paste0("report_stage_from_saved_UPDATED_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Saving report artifacts under:\n  ", REPORT_DIR, "\n\n")

write.csv(metrics_df, file.path(REPORT_DIR, "metrics_trials_from_saved_UPDATED.csv"), row.names = FALSE)

# Summary (mean ± sd across trials)
num_cols <- setdiff(names(metrics_df), c("trial","trial_path"))
summary_df <- data.frame(
  metric = num_cols,
  mean   = sapply(metrics_df[, num_cols, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
  sd     = sapply(metrics_df[, num_cols, drop = FALSE], function(z) sd(z, na.rm = TRUE))
)
write.csv(summary_df, file.path(REPORT_DIR, "metrics_summary_from_saved_UPDATED.csv"), row.names = FALSE)

cat("========== Head(metrics_df) ==========\n")
print(head(metrics_df))
cat("\n========== Summary(mean ± sd) ==========\n")
print(summary_df)

# ----------------------------
# 4.1 Table 
# ----------------------------
get_ms <- function(df, col) {
  if (!col %in% names(df)) return(c(mean = NA_real_, sd = NA_real_))
  mean_sd(df[[col]])
}

ms <- list(
  rmse_u_mean = get_ms(metrics_df, "rmse_u_mean"),
  rmse_u_mode = get_ms(metrics_df, "rmse_u_mode"),
  
  rmse_a_mean = get_ms(metrics_df, "rmse_a_mean"),
  rmse_a_mode = get_ms(metrics_df, "rmse_a_mode"),
  
  rmse_p_all_mean = get_ms(metrics_df, "rmse_p_mean_all"),
  rmse_p_all_mode = get_ms(metrics_df, "rmse_p_mode_all"),
  
  sd_u_mean = get_ms(metrics_df, "sd_u_mean"),
  sd_u_max  = get_ms(metrics_df, "sd_u_max"),
  
  nrmse_u_range_mean = get_ms(metrics_df, "nrmse_u_range_mean"),
  nrmse_u_range_mode = get_ms(metrics_df, "nrmse_u_range_mode"),
  
  rel_a_L2_mean = get_ms(metrics_df, "relL2G_a_mean"),
  rel_a_L2_mode = get_ms(metrics_df, "relL2G_a_mode"),
  
  rel_p_L2_mean = get_ms(metrics_df, "relL2G_p_mean_all"),
  rel_p_L2_mode = get_ms(metrics_df, "relL2G_p_mode_all")
)

paper_table <- data.frame(
  Metric = c(
    "RMSE(u)",
    "RMSE(a)",
    "RMSE(p) (all nodes)",
    "mean(sd_{u|y}) [saved]",
    "max(sd_{u|y})  [saved]",
    "NRMSE_range(u) = RMSE(u)/range(u_true)",
    "RelErr_{L2(Gamma)}(a)",
    "RelErr_{L2(Gamma)}(p)"
  ),
  `u_hat_mean / p_hat_mean / a_hat_mean` = c(
    fmt_mean_sd(ms$rmse_u_mean["mean"], ms$rmse_u_mean["sd"], 4),
    fmt_mean_sd(ms$rmse_a_mean["mean"], ms$rmse_a_mean["sd"], 4),
    fmt_mean_sd(ms$rmse_p_all_mean["mean"], ms$rmse_p_all_mean["sd"], 4),
    fmt_mean_sd(ms$sd_u_mean["mean"], ms$sd_u_mean["sd"], 4),
    fmt_mean_sd(ms$sd_u_max["mean"],  ms$sd_u_max["sd"],  4),
    fmt_mean_sd(ms$nrmse_u_range_mean["mean"], ms$nrmse_u_range_mean["sd"], 4),
    fmt_mean_sd(ms$rel_a_L2_mean["mean"], ms$rel_a_L2_mean["sd"], 4),
    fmt_mean_sd(ms$rel_p_L2_mean["mean"], ms$rel_p_L2_mean["sd"], 4)
  ),
  `u_hat_mode / p_hat_mode / a_hat_mode` = c(
    fmt_mean_sd(ms$rmse_u_mode["mean"], ms$rmse_u_mode["sd"], 4),
    fmt_mean_sd(ms$rmse_a_mode["mean"], ms$rmse_a_mode["sd"], 4),
    fmt_mean_sd(ms$rmse_p_all_mode["mean"], ms$rmse_p_all_mode["sd"], 4),
    "",
    "",
    fmt_mean_sd(ms$nrmse_u_range_mode["mean"], ms$nrmse_u_range_mode["sd"], 4),
    fmt_mean_sd(ms$rel_a_L2_mode["mean"], ms$rel_a_L2_mode["sd"], 4),
    fmt_mean_sd(ms$rel_p_L2_mode["mean"], ms$rel_p_L2_mode["sd"], 4)
  ),
  check.names = FALSE
)
write.csv(paper_table, file.path(REPORT_DIR, "paper_summary_table_from_saved_UPDATED.csv"), row.names = FALSE)
cat("\n========== Paper summary table (mean ± sd) ==========\n")
print(paper_table)

# ----------------------------
# 5) Boxplots 
# ----------------------------

# 5.1 RMSE(u)
save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mean)]", value = metrics_df$rmse_u_mean),
    data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mode)]", value = metrics_df$rmse_u_mode)
  ),
  file.path(REPORT_DIR, "fig_rmse_u_boxplot_UPDATED.png"),
  file.path(REPORT_DIR, "fig_rmse_u_boxplot_UPDATED.pdf"),
  ylab_parse = "plain(RMSE)(hat(u), u[0])"
)

# 5.2 RMSE(a)
save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mean)]", value = metrics_df$rmse_a_mean),
    data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mode)]", value = metrics_df$rmse_a_mode)
  ),
  file.path(REPORT_DIR, "fig_rmse_a_boxplot_UPDATED.png"),
  file.path(REPORT_DIR, "fig_rmse_a_boxplot_UPDATED.pdf"),
  ylab_parse = "plain(RMSE)(hat(a), a[0])"
)

# 5.3 RMSE(p) all nodes
save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mean)]", value = metrics_df$rmse_p_mean_all),
    data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mode)]", value = metrics_df$rmse_p_mode_all)
  ),
  file.path(REPORT_DIR, "fig_rmse_p_all_boxplot_UPDATED.png"),
  file.path(REPORT_DIR, "fig_rmse_p_all_boxplot_UPDATED.pdf"),
  ylab_parse = "plain(RMSE)(hat(p), p[0])"
)

# 5.4 RMSE(p) obs nodes 
if (any(is.finite(metrics_df$rmse_p_mean_obs)) && any(is.finite(metrics_df$rmse_p_mode_obs))) {
  save_boxplot_2(
    rbind(
      data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mean)]", value = metrics_df$rmse_p_mean_obs),
      data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mode)]", value = metrics_df$rmse_p_mode_obs)
    ),
    file.path(REPORT_DIR, "fig_rmse_p_obs_boxplot_UPDATED.png"),
    file.path(REPORT_DIR, "fig_rmse_p_obs_boxplot_UPDATED.pdf"),
    ylab_parse = "plain(RMSE)(hat(p), p[0])"
  )
}

# 5.5 NRMSE_range(u)
save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mean)]", value = metrics_df$nrmse_u_range_mean),
    data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mode)]", value = metrics_df$nrmse_u_range_mode)
  ),
  file.path(REPORT_DIR, "fig_nrmse_u_range_boxplot_UPDATED.png"),
  file.path(REPORT_DIR, "fig_nrmse_u_range_boxplot_UPDATED.pdf"),
  ylab_parse = "plain(NRMSE)[plain(range)](hat(u), u[0])"
)

# 5.6 RelErr_{L2(Gamma)}(a)
save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mean)]", value = metrics_df$relL2G_a_mean),
    data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mode)]", value = metrics_df$relL2G_a_mode)
  ),
  file.path(REPORT_DIR, "fig_relL2G_a_boxplot_UPDATED.png"),
  file.path(REPORT_DIR, "fig_relL2G_a_boxplot_UPDATED.pdf"),
  ylab_parse = "plain(RelErr)[L^2(Gamma)](hat(a), a[0])"
)

# 5.7 RelErr_{L2(Gamma)}(p)
save_boxplot_2(
  rbind(
    data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mean)]", value = metrics_df$relL2G_p_mean_all),
    data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mode)]", value = metrics_df$relL2G_p_mode_all)
  ),
  file.path(REPORT_DIR, "fig_relL2G_p_boxplot_UPDATED.png"),
  file.path(REPORT_DIR, "fig_relL2G_p_boxplot_UPDATED.pdf"),
  ylab_parse = "plain(RelErr)[L^2(Gamma)](hat(p), p[0])"
)

# 5.8 sd summaries of u 
if (any(is.finite(metrics_df$sd_u_mean))) {
  save_boxplot_1(
    data.frame(value = metrics_df$sd_u_mean),
    file.path(REPORT_DIR, "fig_sd_u_mean_boxplot_UPDATED.png"),
    file.path(REPORT_DIR, "fig_sd_u_mean_boxplot_UPDATED.pdf"),
    ylab_parse = "plain(mean)(widehat(sd)[u*'|'*y])"
  )
}
if (any(is.finite(metrics_df$sd_u_max))) {
  save_boxplot_1(
    data.frame(value = metrics_df$sd_u_max),
    file.path(REPORT_DIR, "fig_sd_u_max_boxplot_UPDATED.png"),
    file.path(REPORT_DIR, "fig_sd_u_max_boxplot_UPDATED.pdf"),
    ylab_parse = "plain(max)(widehat(sd)[u*'|'*y])"
  )
}

# ----------------------------
# 6) Save bundle
# ----------------------------
report_bundle <- list(
  timestamp = Sys.time(),
  outdir = OUTDIR,
  report_dir = REPORT_DIR,
  global_path = GLOBAL_PATH,
  trial_paths = trial_paths,
  metrics_df = metrics_df,
  summary_df = summary_df,
  paper_table = paper_table
)
saveRDS(report_bundle, file.path(REPORT_DIR, "report_bundle_from_saved_UPDATED.rds"))

cat("\nDONE. Report artifacts saved under:\n  ", REPORT_DIR, "\n")
cat("Key files:\n",
    "  - metrics_trials_from_saved_UPDATED.csv\n",
    "  - metrics_summary_from_saved_UPDATED.csv\n",
    "  - paper_summary_table_from_saved_UPDATED.csv\n",
    "  - fig_*_UPDATED.pdf and fig_*_UPDATED.png\n",
    "  - report_bundle_from_saved_UPDATED.rds\n", sep = "")

