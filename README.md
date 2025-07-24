# Bayesian Inverse Problems on Metric Graphs

This repository contains R code for solving **Bayesian inverse problems** on **metric graphs** with **finite‑element methods (FEM)**.  
A preconditioned Crank–Nicolson (**pCN**) MCMC sampler with temperature annealing is used to draw posterior samples of the unknown parameter field.

---

## Overview

1. Build a metric graph and FEM mesh.  
2. Solve elliptic or fractional PDEs on the graph.  
3. Generate noisy observations with a mixed (relative + absolute) noise model.  
4. Sample from the posterior distribution with a pCN MCMC algorithm.  
5. Post‑process samples to obtain posterior mean, MAP, variances, RMSE / MAE, and graphical summaries.

---
## Usage

### 1. Create the Metric Graph and Mesh

We use the `metric_graph` class to create the graph and mesh, which forms the foundation for solving the PDEs.

### 2. Define the PDE Solver

The `solve_pde` function solves the elliptic PDE using the FEM method, given the graph and finite element matrices.

### 3. Generate Observational Data

Observational data is generated with a mixed noise model, combining both relative and absolute noise.

### 4. MCMC Sampling

The `pCN_MCMC` function implements the pCN MCMC algorithm with temperature annealing to sample from the posterior distribution.

### 5. Analyze Posterior Samples

After the burn-in phase, posterior samples are analyzed, and the posterior mean and MAP estimates are computed.

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

- **Posterior Mean**: The posterior mean of the unknown parameter \( \mu \), estimated from MCMC samples.
- **MAP Estimate**: The maximum a posteriori estimate of \( \mu \), computed using the mode of the posterior distribution.
- **Pressure Field**: The solutions \( p \) for true \( \mu_0 \), posterior mean \( \mu \), and smoothed MAP \( \mu \).

---

## Requirements

```r
install.packages(c(
  "rSPDE",      # stochastic PDE priors
  "MetricGraph",# graph + mesh utilities
  "Matrix",     # sparse matrix algebra
  "ggplot2",    # plotting
  "viridis",    # colour scales
  "scales"      # helper scales
))


