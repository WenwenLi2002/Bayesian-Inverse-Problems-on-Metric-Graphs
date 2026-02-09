
# Bayesian Inverse Problems on Metric Graphs

This repository contains R code for solving **Bayesian inverse problems** on **metric graphs** using **finite‑element methods (FEM)**.  
A preconditioned Crank–Nicolson (**pCN**) MCMC sampler with temperature annealing is used to draw posterior samples of the unknown parameter field.

---

## Overview

1. Build a metric graph and FEM mesh.  
2. Solve elliptic or fractional PDEs on the graph.  
3. Generate noisy observations with a mixed (relative + absolute) noise model.  
4. Sample from the posterior distribution using a pCN MCMC algorithm.  
5. Post‑process samples to compute posterior mean, MAP, variances, RMSE/MAE, and graphical summaries.

---

## Requirements

To install the required dependencies, run the following R code:

```r
install.packages(c(
  "rSPDE",      # stochastic PDE priors
  "MetricGraph",# graph + mesh utilities
  "Matrix",     # sparse matrix algebra
  "ggplot2",    # plotting
  "viridis",    # colour scales
  "scales"      # helper scales
))
```

---

## Usage

### 1. Create the Metric Graph and Mesh

We use the `metric_graph` class to create the graph and mesh, which forms the foundation for solving the PDEs.

```r
graph <- metric_graph$new(
  perform_merges = TRUE,
  tolerance      = list(edge_edge   = 1e-3,
                        vertex_vertex = 1e-3,
                        edge_vertex   = 1e-3)
)

graph$build_mesh(h = 0.05)   # finer mesh ⇒ higher accuracy
```

### 2. Define the PDE Solver

The `solve_pde` function solves the elliptic PDE using the FEM method, given the graph and finite element matrices.

```r
solve_pde <- function(u) {
  ## Assemble the system matrix L(u) and solve L(u) p = f
  ## Returns the pressure / state vector p.
}
```

### 3. Generate Observational Data

Observational data is generated with a mixed noise model, combining both relative and absolute noise.

```r
observations <- p0[obs_indices] +
  (alpha_noise * abs(p0[obs_indices]) + beta_noise) * epsilon
```

### 4. MCMC Sampling

The `pCN_MCMC` function implements the pCN MCMC algorithm with temperature annealing to sample from the posterior distribution.

```r
mcmc_result <- pCN_MCMC(
  n_iter = 100000,
  beta   = 0.3,
  u_init = u_init
)
posterior_samples <- mcmc_result$samples
```

### 5. Analyze Posterior Samples

After the burn-in phase, posterior samples are analyzed, and the posterior mean and MAP estimates are computed.

```r
burn_in <- 7000
post_samples <- posterior_samples[(burn_in + 1):n_iter, ]

posterior_mean <- colMeans(post_samples)
get_mode       <- function(x) density(x)$x[which.max(density(x)$y)]
map_est        <- apply(post_samples, 2, get_mode)
```

---

## Results

The script outputs the following:

| **Quantity**                | **Description**                                            |
|-----------------------------|------------------------------------------------------------|
| **Posterior mean / MAP of μ** | Point estimates of the parameter field                    |
| **RMSE, MAE**               | Errors between estimates and ground truth                  |
| **Variance & Std-Dev fields** | Nodal uncertainty for μ and pressure p                    |
| **Pressure visualisations**  | True, posterior-mean, and MAP-based p-fields               |

### Example Plots

- **Posterior Mean**: The posterior mean of the unknown parameter \( u \), estimated from MCMC samples.
- **MAP Estimate**: The maximum a posteriori estimate of \( u \), computed using the mode of the posterior distribution.
- **Pressure Field**: The solutions \( p \) for true \( u_0 \), posterior mean \( u \), and smoothed MAP \( u \).

---

## Data

The reference data for this project can be found in the following [OneDrive link](https://uchicagoedu-my.sharepoint.com/:u:/g/personal/wenwenli_uchicago_edu/EcJ_tJHCalNGnN-ZdKqYFAwBbd_o9RTG_ear0MGDPypknw?e=Yl2ubI). 

---

## Updates (see folder `complement_reference_code`)

We have updated the reference code to include **outer-loop drivers** (multi-trial runs) to demonstrate the robustness of the inference-and-reconstruction pipeline across different problem settings; see **`Final_standard_with_outloop`** and **`Final_fractional_with_outloop`**.  
We also consider an additional **oscillatory forcing term** as a further robustness check; see **`Final_standard_oscillary_rhs`**.

The corresponding reference outputs are packaged as **`Standard_data.zip`** and **`Fractional_data.zip`**. Each archive contains saved per-trial bundles (e.g., `trial_bundle.rds`) with ground truth fields, observations, posterior summaries, and diagnostics, so the reported tables and box plots can be reproduced directly from these bundles.  
In addition, the **well-posedness / data-perturbation analysis** (together with its saved trial bundles) can be found in **`Final_standard_wp`** and **`Final_fractional_wp`**.

---

## License

Released under the **MIT License** – see [`LICENSE`](LICENSE) for details.
