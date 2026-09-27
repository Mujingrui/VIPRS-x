# Notes on the submitted code

The code in `R/original/` is exactly as submitted for the course — nothing has been changed.
This file records what is in it, what needs fixing before anyone can run it, and one
statistical issue worth investigating before the results are relied on.

## What is where

| File | Contents |
|---|---|
| `FinalProject_VIPRS.Rmd` | Single-population VIPRS: `E_step_single`, `M_step_single`, `ELBO_single`, `VIPRS_single`, data loading, fitting both populations, Pearson plots |
| `FinalProject_VIPRS-x.Rmd` | Multi-population VIPRS-x: `E_step`, `M_step`, `ELBO`, the EM loop run inline, Pearson plots |
| `Final_Project.R` | Data loading, the VIPRS-x functions again, a `VIPRS_x()` wrapper, the EM loop, plots, and scratch work |
| `Final_Project2.R` | The VIPRS_single functions again, plus fitting and plotting |
| `Untitled.R` | A loose copy of the E-step body, not wrapped in a function — scratch |

The two `.Rmd` files are the authoritative versions: they are what produced the results in the
report. `E_step`/`M_step`/`ELBO` are duplicated between `FinalProject_VIPRS-x.Rmd` and
`Final_Project.R`; `VIPRS_single` and friends are duplicated between `FinalProject_VIPRS.Rmd`
and `Final_Project2.R`.

## Before anything will run

**Absolute paths.** Every data read points at
`/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/...`. With the data now under
`data/` in this repository, those become e.g. `data/eas/LD_eas.csv`. Note also that
`Final_Project.R` reads `LD.csv`, `X_train.csv` etc. while the `.Rmd` files read
`LD_eas.csv`, `X_train_eas.csv` — only the latter names exist.

**`Final_Project2.R` does not parse.** Line 212 contains `round(_PPC_eur_single, 4)`;
`_PPC_eur_single` is not a valid R name. Should be `test_PPC_eur_single`.

**`VIPRS_x()` in `Final_Project.R` cannot execute.** Four problems:

- `Sigma_k_current = Simga_k_path[[iter-1]]` — typo for `Sigma_k_path`
- the call to `E_step()` omits `beta_marginal`, which is a required argument
- `Sigma_k = Sigma_k[[iter]]` in the ELBO call should be `Sigma_k_path[[iter]]`
- `numiter` is only assigned inside the convergence `break`, so the return fails if the loop
  runs to `maxiter`; and `Sigma_jk_path[[1:numiter]]` is invalid — `[[ ]]` takes a single index

This is presumably why `FinalProject_VIPRS-x.Rmd` runs the EM loop inline rather than calling
the wrapper. Either fix the wrapper or delete it.

**Load order.** `Final_Project2.R` references `LD_eas`, `scaled.X_train_eas`, `beta_marginal`,
`X_data_test` and `VIPRS_plot2` before they exist in that file; it only works if
`Final_Project.R` has been run first in the same session.

## Statistical issues to check

**1. The prior covariance Σ_k is not a valid covariance matrix.**

```r
Sigma_k_initial[,,1] <- matrix(c(0.08^2, 0.8, 0.8, 0.08^2), 2, 2)   # diag 0.0064, off-diag 0.8
Sigma_k_initial[,,2] <- matrix(c(1e-5^2, 0.8, 0.8, 1e-5^2), 2, 2)   # diag 1e-10, off-diag 0.8
```

Both matrices have a negative eigenvalue (for the first, 0.0064 ± 0.8), so neither is positive
semi-definite and the implied cross-population correlation far exceeds 1. Because the code
uses `ginv()` rather than `solve()`, R never errors — it returns a generalized inverse and the
algorithm proceeds on an invalid prior.

If the off-diagonal was intended as a *correlation* of 0.8, the covariance would be
0.8 × 0.08 × 0.08 ≈ 0.0051. This is a plausible alternative explanation for VIPRS-x
underperforming, alongside the σ²_{ε,p} issue identified in the report. **Worth re-running
with a valid Σ_k before treating the comparison as settled.**

**2. A dropped term in the single-population E-step.**

```r
uj_update[j] <- log(pii/(1-pii)) + 0.5*log(tau_beta/tau_betaj_update[j])
+0.5*tau_betaj_update[j]*(mu_betaj_update[j]^2)
```

The second line begins a new complete expression, so R evaluates it and discards the result.
The third term of u_j never enters the update, which affects γ_j and therefore the VIPRS
baseline results. The fix is to end the first line with `+`.

**3. `temp33` squares a variance.** In `M_step`:

```r
temp33[p,k,j] <- mu_jk[p,k,j]^2 + Sigma_jk[[j]][p,p,k]^2
```

`Sigma_jk[[j]][p,p,k]` is already a variance; the second moment would use it unsquared. Worth
checking against the derivation in the report.

**4. Convergence is not enforced in VIPRS-x.** The ELBO tolerance check in the EM loop is
commented out and the loop runs a fixed 50 iterations. The report states convergence at
30–40 iterations, so 50 is probably past it, but the check should be live.

**5. Prediction iterations are hard-coded.** `gammaj_out[,7]` for East Asian and
`gammaj_out[,64]` for European are the convergence points read off by hand. They should come
from the returned iteration count. Relatedly, `count <- i` sits after `break` in
`VIPRS_single`, so `count` is always 0 — which is why the indices had to be found manually.

## Smaller things

- `Untitled.R` is scratch and can be deleted.
- γ is clipped to [0.01, 0.99] in most places and [0.001, 0.999] in `Untitled.R`.
- `beta_marginal` is recomputed from the data in the `.Rmd` files even though
  `beta_marginal_eas.csv` and `beta_marginal_eur.csv` exist in `data/`.
- Axis limits in `create_plot` are hard-coded per method (`±0.2` for VIPRS, `±0.001` for
  VIPRS-x), which is the compression of the VIPRS-x predictions showing up in the plotting
  code.

## Data provenance

The CSVs under `data/` were provided with the course project and derive from 1000 Genomes
Project chromosome 22 genotypes; phenotypes are simulated. Confirm that redistributing them
is permitted before making this repository public.
