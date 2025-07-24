# ------------------------------------------------
# 0. Load required packages
# ------------------------------------------------
library(rSPDE)
library(MetricGraph)
library(Matrix)
library(ggplot2)


# ------------------------------------------------
# 1. Create a metric graph and build a mesh
# ------------------------------------------------
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance = list(edge_edge = 1e-3, vertex_vertex = 1e-3, edge_vertex = 1e-3)
)
graph$build_mesh(h = 0.05)  # A smaller h improves mesh accuracy

# Extract mesh vertex coordinates
meshV <- as.matrix(graph$mesh$V)
if (ncol(meshV) < 2) {
  stop("Mesh node matrix has fewer than 2 columns, unable to extract x and y coordinates!")
} else {
  colnames(meshV) <- c("x", "y")
  graph$mesh$V <- meshV
}
x <- meshV[, "x"]
y <- meshV[, "y"]

# ------------------------------------------------
# 2. Compute the finite element matrices and define the right-hand side f(s) = x^2 - y^2
# ------------------------------------------------
graph$compute_fem()
G <- graph$mesh$G  # Stiffness matrix
C <- graph$mesh$C  # Mass matrix

fbar <- x^2 - y^2
f <- as.vector(C %*% fbar)  # Approximate integration using the mass matrix
kappa <- 1  # PDE parameter

# ------------------------------------------------
# 3. Define the function to solve the elliptic PDE: solve_pde
# ------------------------------------------------
solve_pde <- function(u) {
  n_mesh <- nrow(meshV)
  u <- as.numeric(u)
  # Note: The coefficient in the PDE is exp(u)
  Du <- Matrix::Diagonal(n_mesh, exp(u))
  Gu <- (G %*% Du + Du %*% G) / 2         # Modified stiffness matrix
  C_numeric <- as(C, "dgCMatrix")
  L <- kappa^2 * C_numeric + Gu           # Assemble system matrix
  p_solution <- solve(L, f)               # Solve linear system for p
  return(p_solution)
}

# ------------------------------------------------
# 4. Draw the "true" parameter mu_0 from the GP prior and compute the true solution p0
# ------------------------------------------------
mu_0 <- as.numeric(sample_spde(range = 3, sigma = 0.2, alpha = 1,
                               graph = graph, type = "mesh"))
p0 <- solve_pde(mu_0)

# ------------------------------------------------
# 5. Generate observational data using a mixed noise model (relative + absolute noise)
# ------------------------------------------------
alpha_noise <- 0.03   # Control the proportion of relative noise
beta_noise  <- 0.2    # Control the scale of absolute noise, ensuring reasonable noise when p is small

set.seed(123)
obs_indices <- sample(1:nrow(meshV), 1998)
epsilon <- rnorm(length(obs_indices), mean = 0, sd = 1) #Full observations
observations <- p0[obs_indices] + (alpha_noise * abs(p0[obs_indices]) + beta_noise) * epsilon

# ------------------------------------------------
# 6. Define the negative log-likelihood function with temperature annealing
# ------------------------------------------------
neg_log_likelihood <- function(u, T = 1) {
  # Solve the PDE and extract the p values at the observation points
  p_est <- solve_pde(u)
  p_est_obs <- p_est[obs_indices]
  
  # Compute the standard deviation vector based on the mixed noise model
  sd_vec <- alpha_noise * abs(p_est_obs) + beta_noise
  
  # Note: The constant term log(sqrt(2*pi)) is omitted, which does not affect MCMC
  likelihood_term <- sum(log(sd_vec)) + 0.5 * sum(((p_est_obs - observations) / sd_vec)^2)
  
  # Enhanced regularization to constrain the amplitude of u (alleviating nonlinear amplification issues). 
  lambda_reg <- 0.0 #Here we set the regularization term to 0.
  regularization_term <- 0.5 * lambda_reg * sum(u^2)
  
  # Return the negative log-likelihood with annealing incorporated
  return((likelihood_term + regularization_term) / T)
}

# ------------------------------------------------
# 7. Modified pCN MCMC algorithm with temperature annealing and adaptive beta
# ------------------------------------------------
pCN_MCMC <- function(n_iter, beta, u_init, T0 = 5, cooling_factor = 0.95, adapt_interval = 500) {
  u_current <- as.numeric(u_init)
  n_dim <- length(u_current)
  samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
  acceptance <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    # Temperature annealing: calculate the temperature based on the current iteration (gradually decreasing to 1)
    T_current <- max(1, T0 * cooling_factor^(i / adapt_interval))
    
    # Compute the current negative log-likelihood (with annealing factor)
    phi_current <- neg_log_likelihood(u_current, T = T_current)
    
    # Generate a pCN proposal
    xi <- as.numeric(sample_spde(range = 3, sigma = 0.2, alpha = 1,
                                 graph = graph, type = "mesh"))
    if (length(xi) != n_dim) {
      xi <- xi[1:n_dim]
    }
    u_proposed <- sqrt(1 - beta^2) * u_current + beta * xi
    phi_proposed <- neg_log_likelihood(u_proposed, T = T_current)
    
    # Compute acceptance probability (the annealed T is already incorporated in the likelihood)
    log_alpha <- -(phi_proposed - phi_current)
    if (log(runif(1)) < log_alpha) {
      u_current <- u_proposed
      phi_current <- phi_proposed
      acceptance[i] <- 1
    } else {
      acceptance[i] <- 0
    }
    
    samples[i, ] <- u_current
    
    # Adaptive adjustment of beta every adapt_interval iterations
    if (i %% adapt_interval == 0) {
      current_accept_rate <- mean(acceptance[(i - adapt_interval + 1):i])
      if (current_accept_rate < 0.4 * 0.9) {
        beta <- max(beta * 0.9, 0.01)
      } else if (current_accept_rate > 0.4 * 1.1) {
        beta <- beta * 1.2
      }
      cat("Iteration:", i, "Acceptance rate:", current_accept_rate, 
          "Adjusted beta:", beta, "Current Temperature:", T_current, "\n")
    }
  }
  return(list(samples = samples, acceptance = acceptance))
}

# ------------------------------------------------
# 8. Initial guess: sample from the prior, and run the MCMC algorithm
# ------------------------------------------------
u_init <- as.numeric(sample_spde(range = 3, sigma = 0.2, alpha = 1,
                                 graph = graph, type = "mesh"))
if (length(u_init) != length(mu_0)) {
  u_init <- u_init[1:length(mu_0)]
}
n_iter <- 100000
beta <- 0.3
set.seed(456)
mcmc_result <- pCN_MCMC(n_iter, beta, u_init)
posterior_samples <- mcmc_result$samples

# ------------------------------------------------
# 9. Analyze the posterior samples and apply burn-in
# ------------------------------------------------
burn_in <- 7000
posterior_samples_post_burnin <- posterior_samples[(burn_in + 1):n_iter, ]

# ------------------------------------------------
# 10. Plot Posterior Mean, MAP, and True μ₀ – all on the same color scale
# ------------------------------------------------

# (a) Compute posterior mean and MAP
posterior_mean <- apply(posterior_samples_post_burnin, 2, mean)

get_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
map_est <- apply(posterior_samples_post_burnin, 2, get_mode)

# (b) Determine the global colour limits over all three μ fields
all_mu <- c(posterior_mean, map_est, mu_0)
zlim_mu <- range(all_mu)   # [min, max] of all μ values

# (c) Build a shared viridis color scale
library(ggplot2)
shared_mu_scale <- scale_color_viridis_c(
  option = "D",         # viridis “D” palette
  limits = zlim_mu,     # same min/max across plots
  oob    = scales::squish
)

# (d) Plot each field with the shared scale

# (1) Posterior Mean of μ
graph$plot_function(
  data        = posterior_mean,
  continuous  = TRUE,
  scale_color = shared_mu_scale,
  vertex_size = 0,
  main        = "Posterior Mean of μ"
)
try(
  points(meshV[, "x"], meshV[, "y"],
         col = "blue", pch = 19, cex = 0.5),
  silent = TRUE
)

# (2) MAP Estimate of μ
graph$plot_function(
  data        = map_est,
  continuous  = TRUE,
  scale_color = shared_mu_scale,
  vertex_size = 0,
  main        = "MAP Estimate of μ"
)
try(
  points(meshV[, "x"], meshV[, "y"],
         col = "blue", pch = 19, cex = 0.5),
  silent = TRUE
)

# (3) True μ₀
graph$plot_function(
  data        = mu_0,
  continuous  = TRUE,
  scale_color = shared_mu_scale,
  vertex_size = 0,
  main        = "True μ₀"
)
try(
  points(meshV[, "x"], meshV[, "y"],
         col = "red", pch = 19, cex = 0.5),
  silent = TRUE
)

# ------------------------------------------------
# 11. Compute the RMSE of mu (root-mean-square error between the posterior mean and the true mu_0)
# ------------------------------------------------
rmse_mu <- sqrt(mean((posterior_mean - mu_0)^2))
cat("RMSE for mu over all mesh points:", rmse_mu, "\n")

# ------------------------------------------------
# 12. Compute and report the mean absolute error of mu (mean absolute difference between the MAP and the true mu_0)
# ------------------------------------------------
map_est <- apply(posterior_samples_post_burnin, 2, get_mode)
mae_mu <- mean(abs(map_est - mu_0))
cat("Mean absolute error (MAP vs mu_0) over all mesh points:", mae_mu, "\n")


# ------------------------------------------------
# 13. Plot the solution p for True μ₀, Posterior Mean, and Smoothed MAP μ
# ------------------------------------------------

# Compute the three pressure fields
p_true       <- solve_pde(mu_0)            # True parameter field
p_post_mean  <- solve_pde(posterior_mean)  # Posterior mean field

# Smooth the MAP estimate μ over the graph via the unnormalized Laplacian
E            <- graph$mesh$E               
n            <- length(map_est)
A            <- sparseMatrix(
  i    = c(E[,1], E[,2]),
  j    = c(E[,2], E[,1]),
  x    = 1,
  dims = c(n, n)
)
D            <- Diagonal(x = rowSums(A))   # Degree matrix
L            <- D - A                      # Graph Laplacian
gamma        <- 0.1                        # Smoothing strength
u_smooth     <- solve(D + gamma * L, map_est)
p_map_smooth <- solve_pde(u_smooth)        # Smoothed MAP field

# Create a shared viridis colour scale across all three fields
all_p        <- c(p_true, p_post_mean, p_map_smooth)
zlim         <- range(all_p, na.rm = TRUE)
shared_scale <- scale_color_viridis_c(
  option = "D",
  limits = zlim,
  oob    = scales::squish
)

# Plot True μ₀ → p_true
plot_true <- graph$plot_function(
  data        = p_true,
  continuous  = TRUE,
  scale_color = shared_scale,
  vertex_size = 0,
  main        = "Solution p — True μ₀"
)
plot_true

# Plot Posterior Mean μ → p_post_mean
plot_post <- graph$plot_function(
  data        = p_post_mean,
  continuous  = TRUE,
  scale_color = shared_scale,
  vertex_size = 0,
  main        = "Solution p — Posterior Mean μ"
)
plot_post

# Plot Smoothed MAP μ → p_map_smooth
plot_map <- graph$plot_function(
  data        = p_map_smooth,
  continuous  = TRUE,
  scale_color = shared_scale,
  vertex_size = 0,
  main        = "Solution p — Smoothed MAP μ"
)
plot_map


# ------------------------------------------------
# 14. Compute posterior variance at each mesh node and plot it
# ------------------------------------------------

# Calculate the variance of μ at each node from the MCMC samples after burn-in
var_mu <- apply(posterior_samples_post_burnin, 2, var)

#  Define a shared viridis colour scale for the variance plot
#  Here we use the variance range as the limits
zlim_var <- range(var_mu, na.rm = TRUE)
shared_var_scale <- scale_color_viridis_c(
  option = "D",        # use the “D” viridis palette
  limits = zlim_var,   # set colour limits to [min(var_mu), max(var_mu)]
  oob    = scales::squish  # squash values outside the limits
)

# Plot the posterior variance field over the metric graph
graph$plot_function(
  data        = var_mu,          # the vector of variances at each node
  continuous  = TRUE,            # treat data as a continuous field
  scale_color = shared_var_scale,# apply the shared viridis scale
  vertex_size = 0,               # hide vertex markers
  main        = "Posterior Variance of μ"
)

# ------------------------------------------------
# 15. Compute posterior variance of the solution field p at each node and plot
# ------------------------------------------------

# For each posterior μ sample, solve the PDE to get p, then collect into a matrix
n_post <- nrow(posterior_samples_post_burnin)
n_nodes <- nrow(meshV)
p_samples <- matrix(NA, nrow = n_post, ncol = n_nodes)
for (i in seq_len(n_post)) {
  # Solve PDE for the i-th posterior sample of μ
  p_samples[i, ] <- solve_pde(posterior_samples_post_burnin[i, ])
}

#  Compute variance of p across the posterior samples at each node
var_p <- apply(p_samples, 2, var)

#  Define a separate viridis colour scale for the p-variance plot
zlim_pvar <- range(var_p, na.rm = TRUE)
pvar_scale <- scale_color_viridis_c(
  option = "D",        # use viridis “D” palette
  limits = zlim_pvar,  # set colour limits to [min(var_p), max(var_p)]
  oob    = scales::squish
)

#  Plot the posterior variance of p over the metric graph
graph$plot_function(
  data        = var_p,     # vector of pressure variances at each node
  continuous  = TRUE,      # treat data as a continuous field
  scale_color = pvar_scale,# apply the new viridis scale
  vertex_size = 0,         # hide vertex markers
  main        = "Posterior Variance of p"
)

# ------------------------------------------------
# 16. Compute posterior std of u and p at each node and plot
# ------------------------------------------------
# Check if p_samples exists and is not empty
if (exists("p_samples") && !all(is.na(p_samples))) {
  
  # 1. Remove NA values before computing standard deviation
  p_samples_clean <- p_samples[complete.cases(p_samples), ]  # Remove rows with NA values
  
  if (nrow(p_samples_clean) > 0) {  
    
    # 2. Compute the standard deviation of p across the posterior samples at each node
    sd_p <- apply(p_samples_clean, 2, sd)  # Compute standard deviation
    
    # 3. Define a separate viridis color scale for the p-standard deviation plot
    zlim_psd <- range(sd_p, na.rm = TRUE)
    psd_scale <- scale_color_viridis_c(
      option = "D",        
      limits = zlim_psd,   # set color limits to [min(sd_p), max(sd_p)]
      oob    = scales::squish
    )
    
    # 4. Plot the posterior standard deviation of p over the metric graph
    graph$plot_function(
      data        = sd_p,     # vector of pressure standard deviations at each node
      continuous  = TRUE,     
      scale_color = psd_scale, 
      vertex_size = 0,        
      main        = "Posterior Standard Deviation of p"
    )
    
  } else {
    cat("No valid data left in p_samples after removing NAs.")
  }
} else {
  cat("p_samples is missing or empty.")
}

# Check if posterior_samples_post_burnin exists and is not empty
if (exists("posterior_samples_post_burnin") && !all(is.na(posterior_samples_post_burnin))) {
  
  # 1. Remove NA values before computing standard deviation
  posterior_samples_clean <- posterior_samples_post_burnin[complete.cases(posterior_samples_post_burnin), ]  # Remove rows with NA values
  
  if (nrow(posterior_samples_clean) > 0) {  
    
    # 2. Compute the standard deviation of posterior samples across the posterior samples at each node
    sd_posterior <- apply(posterior_samples_clean, 2, sd)  # Compute standard deviation
    
    # 3. Define a separate viridis color scale for the posterior standard deviation plot
    zlim_posterior_sd <- range(sd_posterior, na.rm = TRUE)
    posterior_sd_scale <- scale_color_viridis_c(
      option = "D",       
      limits = zlim_posterior_sd,   # set color limits to [min(sd_posterior), max(sd_posterior)]
      oob    = scales::squish
    )
    
    # 4. Plot the posterior standard deviation of posterior_samples_post_burnin over the metric graph
    graph$plot_function(
      data        = sd_posterior,     # vector of standard deviations at each node
      continuous  = TRUE,             
      scale_color = posterior_sd_scale, 
      vertex_size = 0,                
      main        = "Posterior Standard Deviation of Posterior Samples"
    ) 
    
  } else {
    cat("No valid data left in posterior_samples_post_burnin after removing NAs.")
  }
} else {
  cat("posterior_samples_post_burnin is missing or empty.")
}

