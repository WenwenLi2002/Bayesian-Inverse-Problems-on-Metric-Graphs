# ------------------------------------------------
# 0. Load required packages
# ------------------------------------------------
library(rSPDE)
library(MetricGraph)
library(Matrix)
library(ggplot2)
library(viridis)
library(scales)

# ------------------------------------------------
# 1. Create a metric graph and build a mesh
# ------------------------------------------------
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance = list(edge_edge = 1e-3,
                   vertex_vertex = 1e-3,
                   edge_vertex   = 1e-3)
)
graph$build_mesh(h = 0.2)

# Extract vertex coordinates
meshV <- as.matrix(graph$mesh$V)
if (ncol(meshV) < 2) stop("Mesh vertex matrix has fewer than 2 columns!")
colnames(meshV) <- c("x", "y")
graph$mesh$V <- meshV
x <- meshV[, "x"]
y <- meshV[, "y"]

# ------------------------------------------------
# 2. Compute finite element matrices and cache
# ------------------------------------------------
graph$compute_fem()
G    <- graph$mesh$G         # stiffness matrix
C    <- graph$mesh$C         # mass matrix
fbar <- x^2 - y^2            # forcing function
f    <- as.vector(C %*% fbar)
Cmat <- as(C, "dgCMatrix") # cache sparse mass matrix
kappa <- 1                   # PDE parameter

# ------------------------------------------------
# 3. Fractional PDE solver via eigen-decomposition
# ------------------------------------------------
# Solves L_u^(-frac_order) * f, where L_u = kappa^2 * C + G * diag(e^u)
solve_fractional_pde <- function(u, frac_order) {
  u   <- as.numeric(u)
  Du  <- Diagonal(length(u), exp(u))
  Gu  <- (G %*% Du + Du %*% G) * 0.5
  L   <- kappa^2 * Cmat + Gu
  eig <- eigen(L, symmetric = TRUE)
  vals <- eig$values^(-frac_order)
  V    <- eig$vectors
  # Compute p = V * diag(vals) * V^T * f
  return(as.numeric(V %*% (vals * crossprod(V, f))))
}

# ------------------------------------------------
# 4. Sample a "true" parameter and compute true solution
# ------------------------------------------------
frac_order <- 1.5
mu_0       <- as.numeric(
  sample_spde(range = 3, sigma = 0.2, alpha = 1,
              graph = graph, type = "mesh")
)
p0         <- solve_fractional_pde(mu_0, frac_order)

# ------------------------------------------------
# 5. Generate noisy observations
# ------------------------------------------------
alpha_noise <- 0.03
beta_noise  <- 0.2
set.seed(123)
n_obs       <- min(507, nrow(meshV))
obs_indices <- sample.int(nrow(meshV), n_obs)
epsilon     <- rnorm(n_obs)
observations <- p0[obs_indices] +
  (alpha_noise * abs(p0[obs_indices]) + beta_noise) * epsilon


# ------------------------------------------------
# 6. Negative log-likelihood for MCMC
# ------------------------------------------------
neg_log_likelihood <- function(u, T = 1) {
  p_est <- solve_fractional_pde(u, frac_order)
  p_obs <- p_est[obs_indices]
  sd_v  <- alpha_noise * abs(p_obs) + beta_noise
  ll    <- sum(log(sd_v)) +
    0.5 * sum(((p_obs - observations) / sd_v)^2)
  return(ll / T)
}

# ------------------------------------------------
# 7. pCN MCMC sampler
# ------------------------------------------------
pCN_MCMC <- function(n_iter, beta, u_init,
                     T0 = 5, cooling = 0.95, adapt = 500) {
  u_curr <- as.numeric(u_init)
  n_dim  <- length(u_curr)
  samples<- matrix(NA, nrow = n_iter, ncol = n_dim)
  acc    <- logical(n_iter)
  
  for (i in seq_len(n_iter)) {
    T_curr <- max(1, T0 * cooling^(i / adapt))
    phi_c  <- neg_log_likelihood(u_curr, T_curr)
    xi     <- head(
      as.numeric(sample_spde(range = 3, sigma = 0.2, alpha = 1,
                             graph = graph, type = "mesh")),
      n_dim
    )
    u_prop <- sqrt(1 - beta^2) * u_curr + beta * xi
    phi_p  <- neg_log_likelihood(u_prop, T_curr)
    if (log(runif(1)) < -(phi_p - phi_c)) {
      u_curr <- u_prop
      phi_c  <- phi_p
      acc[i] <- TRUE
    }
    samples[i, ] <- u_curr
    
    if (i %% adapt == 0) {
      rate <- mean(acc[(i - adapt + 1):i])
      beta <- beta * ifelse(rate < 0.36, 0.9,
                            ifelse(rate > 0.44, 1.2, 1))
      cat(sprintf("Iter %d | rate %.3f | beta %.3f | T %.2f\n",
                  i, rate, beta, T_curr))
    }
  }
  
  return(list(samples = samples, acceptance = acc))
}

# ------------------------------------------------
# 8. Run MCMC and summarize posterior
# ------------------------------------------------
set.seed(456)
u_init <- head(
  as.numeric(sample_spde(range = 3, sigma = 0.2, alpha = 1,
                         graph = graph, type = "mesh")),
  length(mu_0)
)
n_iter <- 1e5; beta <- 0.3
mcmc_res  <- pCN_MCMC(n_iter, beta, u_init = u_init)
post_samps<- mcmc_res$samples

# Discard burn-in and compute statistics
burn      <- 7000
post      <- post_samps[(burn + 1):n_iter, ]
post_mean <- colMeans(post)
get_mode  <- function(x) density(x)$x[which.max(density(x)$y)]
map_est   <- apply(post, 2, get_mode)

# ------------------------------------------------
# 9. Visualization: μ fields with shared color scale
# ------------------------------------------------
mu_limits  <- range(c(post_mean, map_est, mu_0))
shared_col <- scale_color_viridis_c(
  option = "D", limits = mu_limits,
  oob = scales::squish
)

# Plot posterior mean, MAP, and true μ side by side
graph$plot_function(data = post_mean, continuous = TRUE,
                    scale_color = shared_col, vertex_size = 0,
                    main = "Posterior Mean of μ")
graph$plot_function(data = map_est,   continuous = TRUE,
                    scale_color = shared_col, vertex_size = 0,
                    main = "MAP Estimate of μ")
graph$plot_function(data = mu_0,       continuous = TRUE,
                    scale_color = shared_col, vertex_size = 0,
                    main = "True μ")

# ------------------------------------------------
# 10. Error metrics: RMSE and MAE
# ------------------------------------------------
p_est_mean <- solve_fractional_pde(post_mean, frac_order)
rmse_p     <- sqrt(mean((p_est_mean - p0)^2))
rmse_mu    <- sqrt(mean((post_mean - mu_0)^2))
mae_mu     <- mean(abs(map_est - mu_0))
cat(sprintf("RMSE(p): %.4f | RMSE(μ): %.4f | MAE(μ): %.4f\n",
            rmse_p, rmse_mu, mae_mu))

# ------------------------------------------------
# 11. Plot pressure fields: true, posterior mean, MAP
# ------------------------------------------------
graph$plot_function(
  data = solve_fractional_pde(mu_0,       frac_order),
  vertex_size = 0, main = "True p-field"
)
graph$plot_function(
  data = solve_fractional_pde(post_mean,  frac_order),
  vertex_size = 0, main = "Posterior Mean p-field"
)
graph$plot_function(
  data = solve_fractional_pde(map_est,    frac_order),
  vertex_size = 0, main = "MAP p-field"
)

# ------------------------------------------------
# 12. Posterior variance fields for μ and p
# ------------------------------------------------
var_mu   <- apply(post,    2, var)
zlim_mu  <- range(var_mu)
graph$plot_function(
  data = var_mu, continuous = TRUE,
  scale_color = scale_color_viridis_c(option = "D",
                                      limits = zlim_mu,
                                      oob = scales::squish),
  vertex_size = 0, main = "Posterior Variance of μ"
)

# Variance of pressure samples
p_samples <- t(apply(post, 1,
                     function(u) solve_fractional_pde(u, frac_order)))
var_p    <- apply(p_samples, 2, var)
zlim_p   <- range(var_p)
graph$plot_function(
  data = var_p, continuous = TRUE,
  scale_color = scale_color_viridis_c(option = "D",
                                      limits = zlim_p,
                                      oob = scales::squish),
  vertex_size = 0, main = "Posterior Variance of p"
)


# ------------------------------------------------
# 12. Posterior std fields for μ and p
# ------------------------------------------------
# Plot the posterior standard deviation of μ over the metric graph
std_mu <- sqrt(var_mu)  
zlim_mu_std <- range(std_mu) 
graph$plot_function(
  data        = std_mu,     
  continuous  = TRUE,       # treat data as a continuous field
  scale_color = scale_color_viridis_c(option = "D",  
                                      limits = zlim_mu_std,  
                                      oob = scales::squish),  
  vertex_size = 0,          # hide vertex markers
  main        = "Posterior Standard Deviation of u"  
) 

# Plot the posterior standard deviation of p over the metric graph
std_p <- sqrt(var_p)  
zlim_p_std <- range(std_p)
graph$plot_function(
  data        = std_p,     
  continuous  = TRUE,      
  scale_color = scale_color_viridis_c(option = "D",  
                                      limits = zlim_p_std,  
                                      oob = scales::squish),  
  vertex_size = 0,         # hide vertex markers
  main        = "Posterior Standard Deviation of p"  
)