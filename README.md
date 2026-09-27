# VIPRS-x — Variational Inference for Polygenic Risk Scores across Ancestries

Variational EM for Bayesian polygenic risk scores, implemented from scratch in R, extended
from a single population to joint modelling across two ancestries. 

**The problem.** Polygenic risk scores predict a phenotype from hundreds of thousands of
genetic markers. Bayesian formulations give calibrated uncertainty but are usually fit by
MCMC, which does not scale. VIPRS (Zabad & Li, 2023) replaces MCMC with variational inference,
approximating the posterior over effect sizes by a tractable family and maximising an ELBO.

It fits one population at a time. Most PRS methods are trained on European-ancestry cohorts
and transfer poorly elsewhere, so the question here is whether **jointly** modelling effect
sizes across ancestries — letting related populations share information while keeping
population-specific effects — predicts better than fitting each population separately.

## Model

**Baseline (VIPRS).** For one population, `y = Xβ + ε` with a spike-and-slab prior on each
effect size, `β_j ~ π N(0, σ²_β) + (1 − π) δ₀`. The variational family factorises as

    q(β, s) = Π_j N(β_j; μ_j, σ²_j) · Bern(s_j; γ_j)

and the E-step updates (μ_j, σ²_j, γ_j) marker by marker while the M-step updates the
hyperparameters σ²_β, π and σ²_ε in closed form. Convergence is monitored by the ELBO.

**Extension (VIPRS-x).** For P populations, `y_p = X_p β_p + ε_p`. Effect sizes are stacked
into a P-vector per marker and given a **multivariate Gaussian mixture** prior with K
components,

    q(β, s) = Π_j Π_k [ N_p(β_j; μ_jk, Σ_jk) ]^{s_jk} [ γ_jk ]^{s_jk}

so that μ_jk is a P-vector and Σ_jk a P × P matrix carrying the cross-population covariance of
effects. The E-step updates μ_jk, Σ_jk and γ_jk; the M-step updates Σ_k, π_k and σ²_{ε,p}.
Here P = 2 (East Asian, European) and K = 2 mixture components.

Everything — E-step, M-step, ELBO, the full variational EM loop — is written in base R with
`MASS` for the generalized inverse. No modelling packages.

## Data

1000 Genomes Project, chromosome 22, M = 100 SNPs, two populations: East Asian (n = 352
train) and European (n = 273 train), with held-out test sets. Phenotypes simulated as
`y_p = X_p β_p + ε_p`. Per population the repository carries the genotype matrices, phenotypes,
the LD matrix R, and the marginal effect estimates, under `data/`.

## Results

Pearson correlation between predicted and true phenotype:

| Method | Population | Train | Test |
|---|---|---:|---:|
| VIPRS (single) | East Asian | 0.7315 | 0.7315 |
| VIPRS (single) | European | 0.6263 | 0.6263 |
| VIPRS-x (joint) | East Asian | 0.6111 | 0.5660 |
| VIPRS-x (joint) | European | 0.4887 | 0.3812 |

Squared error, ‖ŷ − y‖₂:

| Method | Population | Train | Test |
|---|---|---:|---:|
| VIPRS | East Asian | 17.7760 | 10.7644 |
| VIPRS | European | 15.7962 | 9.8137 |
| VIPRS-x | East Asian | 18.7325 | 11.3123 |
| VIPRS-x | European | 16.4901 | 10.1970 |

**Joint modelling did not help here.** VIPRS-x is worse than the single-population baseline on
correlation in both populations and on both splits, while the squared errors are close enough
to be uninformative. The gap is largest in the European population, where the predicted
phenotypes are compressed toward zero — the signature of underestimated effect sizes.

The report attributes this to σ²_{ε,p} not being updated well in the joint M-step. See
[`NOTES.md`](NOTES.md) for an additional candidate explanation found while preparing this
repository, concerning the initialisation of the prior covariance Σ_k.

Convergence: VIPRS reached the 10⁻⁴ ELBO tolerance after 7 iterations for East Asian and 64
for European; VIPRS-x after roughly 30–40 iterations.

Full write-up, including the derivations and a discussion of hyperparameter-tuning strategies
(grid search, Bayesian optimisation, Bayesian model averaging):
[`report/final-project.pdf`](report/final-project.pdf)

## Layout

```
R/original/   the code as submitted (see NOTES.md before running)
data/         1000 Genomes-derived inputs, by population
figures/      Pearson correlation plots, train and test, both methods
report/       project write-up
```

## Running it

Requires R with `MASS`, `psych`, `ggplot2`, `gridExtra`.

The scripts in `R/original/` use absolute paths from the machine they were written on and need
editing before they will run — see [`NOTES.md`](NOTES.md), which also documents several issues
found in the submitted code. **Read it before using these results.**

## Reference

Zabad, S., Gravel, S. & Li, Y. (2023). Fast and accurate Bayesian polygenic risk modeling with
variational inference. *American Journal of Human Genetics*.
