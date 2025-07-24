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
