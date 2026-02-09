suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE)

# ----------------------------
# 0) User: absolute OUTDIR
# ----------------------------
OUTDIR <- "/Users/wenwenli/inverse_mcmc_fractional_20260118_105241"
stopifnot(dir.exists(OUTDIR))

# ----------------------------
# 1) Helpers: unique save path
# ----------------------------
make_unique_path <- function(path) {
  if (!file.exists(path)) return(path)
  ext <- tools::file_ext(path)
  base <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", path) else path
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (nzchar(ext)) paste0(base, "_", stamp, ".", ext) else paste0(base, "_", stamp)
}

safe_write_csv <- function(df, path) {
  path2 <- make_unique_path(path)
  write.csv(df, path2, row.names = FALSE)
  message("Wrote: ", path2)
  invisible(path2)
}

safe_write_lines <- function(lines, path) {
  path2 <- make_unique_path(path)
  writeLines(lines, con = path2)
  message("Wrote: ", path2)
  invisible(path2)
}

safe_ggsave_pair <- function(path_png, plot, width = 6.5, height = 4.2, dpi = 220) {
  # normalize extension
  if (!grepl("\\.png$", path_png, ignore.case = TRUE)) {
    path_png <- paste0(path_png, ".png")
  }
  
  path_png2 <- make_unique_path(path_png)
  
  path_pdf2 <- sub("\\.png$", ".pdf", path_png2, ignore.case = TRUE)
  if (file.exists(path_pdf2)) path_pdf2 <- make_unique_path(path_pdf2)
  
  ggsave(path_png2, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved: ", path_png2)
  
  ggsave(path_pdf2, plot = plot, width = width, height = height)
  message("Saved: ", path_pdf2)
  
  invisible(list(png = path_png2, pdf = path_pdf2))
}

# ----------------------------
# 2) Load global bundle for mass matrix C
# ----------------------------
gb_path <- file.path(OUTDIR, "global", "global_bundle.rds")
stopifnot(file.exists(gb_path))
gb <- readRDS(gb_path)

Cmat <- gb$C
if (!inherits(Cmat, "dgCMatrix")) Cmat <- as(Cmat, "dgCMatrix")

n_obs <- gb$n_obs
n_iter <- gb$n_iter
burn_in_user <- gb$burn_in_user

# ----------------------------
# 3) Metrics
# ----------------------------
rmse <- function(est, truth) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((est[ok] - truth[ok])^2))
}

# NRMSE_range(u) = RMSE(u)/range(u_true)
nrmse_u_range <- function(est, truth, eps_range = 1e-12) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  r <- diff(range(truth[ok], na.rm = TRUE))
  rmse(est[ok], truth[ok]) / max(r, eps_range)
}

# RelErr_{L2(Gamma)} using FEM mass matrix C:
#   ||e||_C / ||truth||_C, where ||v||_C^2 = v^T C v
relerr_L2Gamma <- function(est, truth, Cmat, eps = 1e-12) {
  est <- as.numeric(est); truth <- as.numeric(truth)
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) return(NA_real_)
  e <- est[ok] - truth[ok]
  t <- truth[ok]
  Csub <- Cmat[ok, ok, drop = FALSE]
  num <- sqrt(as.numeric(crossprod(e, Csub %*% e)))
  den <- sqrt(as.numeric(crossprod(t, Csub %*% t)))
  num / max(den, eps)
}

exp_clip <- function(u, clip = 20) {
  u <- as.numeric(u)
  exp(pmin(pmax(u, -clip), clip))
}

# ----------------------------
# 4) Read trial bundles
# ----------------------------
trial_dirs <- list.dirs(OUTDIR, recursive = FALSE, full.names = TRUE)
trial_dirs <- trial_dirs[grepl("^trial_[0-9]{3}$", basename(trial_dirs))]
tb_paths <- file.path(trial_dirs, "trial_bundle.rds")
tb_paths <- tb_paths[file.exists(tb_paths)]
if (length(tb_paths) == 0) stop("No trial_bundle.rds found under OUTDIR.")

read_one_trial_metrics <- function(tb_path) {
  tb <- readRDS(tb_path)
  
  u_true <- tb$u_true
  p_true <- tb$p_true
  
  # mean vs mode (marginal posterior mode)
  u_mean <- tb$u_post_mean
  u_mode <- tb$u_map_raw
  p_mean <- tb$p_post_mean
  p_mode <- tb$p_map_raw
  
  # a = exp(u)
  a_true <- exp_clip(u_true)
  a_mean <- exp_clip(u_mean)
  a_mode <- exp_clip(u_mode)
  
  data.frame(
    trial = tb$trial,
    
    # RMSE
    rmse_u_mean = rmse(u_mean, u_true),
    rmse_u_mode = rmse(u_mode, u_true),
    rmse_p_mean = rmse(p_mean, p_true),
    rmse_p_mode = rmse(p_mode, p_true),
    rmse_a_mean = rmse(a_mean, a_true),
    rmse_a_mode = rmse(a_mode, a_true),
    
    # scale-free errors
    nrmse_u_range_mean = nrmse_u_range(u_mean, u_true),
    nrmse_u_range_mode = nrmse_u_range(u_mode, u_true),
    
    relerr_p_L2G_mean  = relerr_L2Gamma(p_mean, p_true, Cmat),
    relerr_p_L2G_mode  = relerr_L2Gamma(p_mode, p_true, Cmat),
    
    relerr_a_L2G_mean  = relerr_L2Gamma(a_mean, a_true, Cmat),
    relerr_a_L2G_mode  = relerr_L2Gamma(a_mode, a_true, Cmat),
    
    # optional run fields 
    acc_rate_all = tb$acc_rate_all,
    beta_final   = tb$beta_final,
    burn_used    = tb$burn_used,
    i_star       = tb$i_star
  )
}

metrics_df <- do.call(rbind, lapply(tb_paths, read_one_trial_metrics))
metrics_df <- metrics_df[order(metrics_df$trial), ]
n_trials_done <- nrow(metrics_df)

# ----------------------------
# 5) Summary (mean ± sd across trials)
# ----------------------------
num_cols <- setdiff(names(metrics_df), "trial")
summary_df <- data.frame(
  metric = num_cols,
  mean   = sapply(metrics_df[, num_cols, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
  sd     = sapply(metrics_df[, num_cols, drop = FALSE], function(z) sd(z, na.rm = TRUE))
)

# ----------------------------
# 6) Table
# ----------------------------
fmt_pm <- function(m, s, digits = 3) {
  if (!is.finite(m) || !is.finite(s)) return("NA")
  sprintf(paste0("%.", digits, "f $\\pm$ %.", digits, "f"), m, s)
}
get_ms <- function(colname) {
  z <- metrics_df[[colname]]
  c(mean = mean(z, na.rm = TRUE), sd = sd(z, na.rm = TRUE))
}

u_rmse_mean <- get_ms("rmse_u_mean"); u_rmse_mode <- get_ms("rmse_u_mode")
p_rmse_mean <- get_ms("rmse_p_mean"); p_rmse_mode <- get_ms("rmse_p_mode")
a_rmse_mean <- get_ms("rmse_a_mean"); a_rmse_mode <- get_ms("rmse_a_mode")

u_rng_mean  <- get_ms("nrmse_u_range_mean"); u_rng_mode <- get_ms("nrmse_u_range_mode")
p_rel_mean  <- get_ms("relerr_p_L2G_mean");  p_rel_mode <- get_ms("relerr_p_L2G_mode")
a_rel_mean  <- get_ms("relerr_a_L2G_mean");  a_rel_mode <- get_ms("relerr_a_L2G_mode")

latex_lines <- c(
  "% Auto-generated referee metrics table (no rerun MCMC)",
  "\\begin{table}[t]",
  "\\centering",
  sprintf(paste0(
    "\\caption{Quantitative errors over %d repetitions ($n_{\\mathrm{obs}}=%d$). ",
    "We compare the posterior mean estimator (subscript \\texttt{mean}) and the pointwise marginal posterior mode (subscript \\texttt{mode}). ",
    "For $u$, we report $\\mathrm{NRMSE}_{\\mathrm{range}}(\\hat u,u_0)=\\mathrm{RMSE}(\\hat u,u_0)/(\\max u_0-\\min u_0)$. ",
    "For $p$ and $a=\\exp(u)$, we report $\\mathrm{RelErr}_{L^2(\\Gamma)}(\\hat v,v_0)=\\|\\hat v-v_0\\|_C/\\|v_0\\|_C$ using the FEM mass matrix.}"
  ), n_trials_done, n_obs),
  "\\label{tab:referee-metrics}",
  "\\begin{tabular}{lcc}",
  "\\hline",
  "Metric & $\\hat{\\cdot}_{\\rm mean}$ & $\\hat{\\cdot}_{\\rm mode}$ \\\\",
  "\\hline",
  sprintf("$\\mathrm{RMSE}(\\hat u,u_0)$ & %s & %s \\\\",
          fmt_pm(u_rmse_mean["mean"], u_rmse_mean["sd"]),
          fmt_pm(u_rmse_mode["mean"], u_rmse_mode["sd"])),
  sprintf("$\\mathrm{NRMSE}_{\\rm range}(\\hat u,u_0)$ & %s & %s \\\\",
          fmt_pm(u_rng_mean["mean"], u_rng_mean["sd"]),
          fmt_pm(u_rng_mode["mean"], u_rng_mode["sd"])),
  sprintf("$\\mathrm{RMSE}(\\hat p,p_0)$ & %s & %s \\\\",
          fmt_pm(p_rmse_mean["mean"], p_rmse_mean["sd"]),
          fmt_pm(p_rmse_mode["mean"], p_rmse_mode["sd"])),
  sprintf("$\\mathrm{RelErr}_{L^2(\\Gamma)}(\\hat p,p_0)$ & %s & %s \\\\",
          fmt_pm(p_rel_mean["mean"], p_rel_mean["sd"]),
          fmt_pm(p_rel_mode["mean"], p_rel_mode["sd"])),
  "\\hline",
  sprintf("$\\mathrm{RMSE}(\\hat a,a_0)$ & %s & %s \\\\",
          fmt_pm(a_rmse_mean["mean"], a_rmse_mean["sd"]),
          fmt_pm(a_rmse_mode["mean"], a_rmse_mode["sd"])),
  sprintf("$\\mathrm{RelErr}_{L^2(\\Gamma)}(\\hat a,a_0)$ & %s & %s \\\\",
          fmt_pm(a_rel_mean["mean"], a_rel_mean["sd"]),
          fmt_pm(a_rel_mode["mean"], a_rel_mode["sd"])),
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)

# ----------------------------
# 7) Plots: boxplots 
# ----------------------------
report_dir <- file.path(OUTDIR, "report_stage_referee_metrics_relerr_NOTATION_21BOLD")
dir.create(report_dir, showWarnings = FALSE)

theme_paper_21_bold <- function() {
  theme_bw(base_size = 21) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_text(size = 21, face = "bold"),
      axis.title.y = element_text(size = 21, face = "bold"),
      axis.text.x  = element_text(size = 21, face = "bold"),
      axis.text.y  = element_text(size = 21, face = "bold"),
      panel.grid.major = element_line(linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

# estimator labels as plotmath strings 
save_boxplot <- function(df, ylab_parse, fname_png, width = 6.5, height = 4.2) {
  p <- ggplot(df, aes(x = estimator, y = value)) +
    geom_boxplot(outlier.size = 1.0, linewidth = 0.6) +
    scale_x_discrete(labels = scales::label_parse()) +
    labs(x = "", y = parse(text = ylab_parse)[[1]]) +
    theme_paper_21_bold()
  
  # save BOTH png and pdf
  safe_ggsave_pair(file.path(report_dir, fname_png), p, width = width, height = height, dpi = 220)
  invisible(p)
}

# ---- RMSE(u)
df_rmse_u <- rbind(
  data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mean)]", value = metrics_df$rmse_u_mean),
  data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mode)]", value = metrics_df$rmse_u_mode)
)
save_boxplot(df_rmse_u, "plain(RMSE)(hat(u), u[0])", "box_rmse_u.png")

# ---- NRMSE_range(u)
df_nrmse_u <- rbind(
  data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mean)]", value = metrics_df$nrmse_u_range_mean),
  data.frame(trial = metrics_df$trial, estimator = "hat(u)[plain(mode)]", value = metrics_df$nrmse_u_range_mode)
)
save_boxplot(df_nrmse_u, "plain(NRMSE)[plain(range)](hat(u), u[0])", "box_nrmse_u_range.png")

# ---- RMSE(p)
df_rmse_p <- rbind(
  data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mean)]", value = metrics_df$rmse_p_mean),
  data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mode)]", value = metrics_df$rmse_p_mode)
)
save_boxplot(df_rmse_p, "plain(RMSE)(hat(p), p[0])", "box_rmse_p.png")

# ---- RelErr_{L2(Gamma)}(p)
df_rel_p <- rbind(
  data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mean)]", value = metrics_df$relerr_p_L2G_mean),
  data.frame(trial = metrics_df$trial, estimator = "hat(p)[plain(mode)]", value = metrics_df$relerr_p_L2G_mode)
)
save_boxplot(df_rel_p, "plain(RelErr)[L^2(Gamma)](hat(p), p[0])", "box_relerr_p_L2Gamma.png")

# ---- RMSE(a)
df_rmse_a <- rbind(
  data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mean)]", value = metrics_df$rmse_a_mean),
  data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mode)]", value = metrics_df$rmse_a_mode)
)
save_boxplot(df_rmse_a, "plain(RMSE)(hat(a), a[0])", "box_rmse_a.png")

# ---- RelErr_{L2(Gamma)}(a)
df_rel_a <- rbind(
  data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mean)]", value = metrics_df$relerr_a_L2G_mean),
  data.frame(trial = metrics_df$trial, estimator = "hat(a)[plain(mode)]", value = metrics_df$relerr_a_L2G_mode)
)
save_boxplot(df_rel_a, "plain(RelErr)[L^2(Gamma)](hat(a), a[0])", "box_relerr_a_L2Gamma.png")

# ----------------------------
# 8) Save tables 
# ----------------------------
safe_write_csv(metrics_df, file.path(report_dir, "metrics_trials_referee_NOTATION_21BOLD.csv"))
safe_write_csv(summary_df, file.path(report_dir, "metrics_summary_referee_NOTATION_21BOLD.csv"))
safe_write_lines(latex_lines, file.path(report_dir, "table_referee_metrics_NOTATION_21BOLD.tex"))

note_lines <- c(
  sprintf("OUTDIR: %s", OUTDIR),
  sprintf("Completed trials: %d", n_trials_done),
  sprintf("n_obs (data locations): %d", n_obs),
  sprintf("n_iter: %s, burn_in_user: %s", as.character(n_iter), as.character(burn_in_user)),
  "",
  "Key outputs written under:",
  report_dir,
  "",
  "Plot outputs:",
  "  - For each metric plot, BOTH .png and .pdf are saved (same basename).",
  "",
  "Plot styling:",
  "  - axis titles and tick labels: 21pt, bold",
  "  - y-axis label: metric notation only (no extra text)",
  "  - x-axis: parsed estimator notation (hat{u}_mean vs hat{u}_mode etc.)"
)
safe_write_lines(note_lines, file.path(report_dir, "README_referee_metrics_NOTATION_21BOLD.txt"))

message("Done. Referee-ready metrics written under: ", report_dir)

