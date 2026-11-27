# ci_LoA
Simulations for the coverage of various confidence intervals for the LoA limits

The repo contains R scripts for the simulations in the manuscript Confidence Intervals 
for Bland-Altman Limits of Agreementon, by Oke Gerke and Erik Thorlund Parner. 

## Simulation scripts: organization

Monte Carlo simulations that compare methods for prediction intervals (PI) /
limits of agreement (LoA) and for confidence intervals of the median. They
vary sample size and data-generating distribution (`normal`, `lognormal`,
`chisq4`, `mixture`). 

### Workflow

1. `Simulate <dist> non-bootstrap.R` (4 scripts)
2. `Simulate <dist> bootstrap batch <k>.R` (one script per batch)
3. `Combine simulations.R`

### Methods compared
- `confintr` (nonparametric)
- `blandr` / `t.test` (normal theory)
- Non-central-t exact CI (PI/LoA only)
- Harrell-Davis (`Hmisc::hdquantile`)
- Sasaki-Wang (`exquantile_R.R`)
- Bootstrap (`boot` / `boot.ci`: `norm`, `perc`, `bca`; bootstrap scripts only)

### Simulation scripts
Each script computes PI/LoA and median results together for every replicate.
It sets its own `true_lower`, `true_upper`, `true_median` and the
`outcome <- r...()` call. The sample-size grid `n_values` is
`seq(10, 1000, by = 10)` plus 5, 6, 145 and 146.

| Script type | Behaviour |
|---|---|
| `Simulate <dist> non-bootstrap.R` | Loops through all batches up to `target_total = 25000` (25 batches). Skips batches whose file already exists. A single run completes everything remaining. |
| `Simulate <dist> bootstrap batch <k>.R` | One batch per run, set by `batch <- <k>`. Stops if the batch file exists. Intended to run on separate computers because `boot::boot()` is slow. |

Batch scheme (identical across all scripts):
- Batch `k` covers replicates `((k-1)*1000 + 1):(k*1000)`.
- The seed is `1` for batch 1 and `(k-1)*1000` for `k > 1`.
- Each batch is saved to its own `<dist>_nonboot_sims_batch_<k>.RData` or
  `<dist>_boot_sims_batch_<k>.RData` file (never cumulative).
- `boot_replicates` is 1100 and must stay above `max(n_values)` for `bca` to be valid.

### Helper scripts
- `Combine simulations.R`: for each distribution, row-binds the batch files
  (up to `max_batch`, the furthest completed bootstrap batch). It then
  `full_join`s the non-bootstrap and bootstrap results on `(n, replicate)` and
  saves `Data/<dist>_combined_sims.RData`.
  Batch files must be in `Data/` first.
- `exquantile_R.R`: Sasaki-Wang quantile CI implementation.

### Conventions
- Result columns are named `<method>_<lower/upper>_<coverage/width>`, so figure code works for any distribution.
- `batch_seed()` and the batch size (1000) are identical across all distributions.
- To add a distribution, copy the closest existing script pair and edit the true values and the `r...()` call.
`Data/` (inputs and combined results), `Results/`, `Papers/`, `Communication/`, `Develop exact CI(PI)/`, `Test and development/`, `Old versions/`.

