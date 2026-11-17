library(tidyverse)
library(confintr)   # ci_quantile
library(blandr)     # blandr.statistics (PI)
library(tictoc)
library(Hmisc)      # hdquantile

options(warn = -1)

source("exquantile_R.R")

# Settings ----------------------------------------------------------------
batch_size   <- 1000    # replicates saved per batch file
target_total <- 50000   # total replicates to accumulate across all batches
n_batches    <- target_total / batch_size

# Batch k covers replicates ((k - 1) * batch_size + 1):(k * batch_size).
# Batch 1 uses seed 1; batch k > 1 uses seed (k - 1) * batch_size.
batch_seed <- function(batch) if (batch == 1L) 1L else (batch - 1L) * batch_size

# True values for N(0,1)
true_lower  <- qnorm(0.025)
true_upper  <- qnorm(0.975)
true_median <- 0

# Combined n grid (union of PI boundary 145/146 and median boundary 5/6)
n_values <- sort(unique(c(seq(10, 1000, by = 10), 5, 6, 145, 146)))
nsetting <- length(n_values)

# Helper functions --------------------------------------------------------
ci_LoA_exact <- function(x, p = 0.95, cilevel = 0.95) {
  N    <- length(x)
  mean <- mean(x)
  sd   <- sd(x)

  kloa        <- qnorm(1 - (1 - p) / 2)
  df          <- N - 1
  np          <- kloa * sqrt(N)
  smallgamma  <- (1 - cilevel) / 2
  biggamma    <- (1 + cilevel) / 2
  kouterciloa <- qt(biggamma,   df, np) / sqrt(N)
  kinnerciloa <- qt(smallgamma, df, np) / sqrt(N)

  data.frame(
    cilloalow = mean - kouterciloa * sd,
    cilloaupp = mean - kinnerciloa * sd,
    ciuloalow = mean + kinnerciloa * sd,
    ciuloaupp = mean + kouterciloa * sd
  )
}

# Exact normal-quantile CI (Carkeet 2015 / Owen 1968), generalized to any
# target probability q via the noncentral t distribution. At q = 0.5 the
# noncentrality parameter is zero, so this reduces to the ordinary
# one-sample t interval around the mean.
ci_quantile_exact <- function(x, q, cilevel = 0.95) {
  N    <- length(x)
  mean <- mean(x)
  sd   <- sd(x)

  df  <- N - 1
  ncp <- qnorm(q) * sqrt(N)
  smallgamma <- (1 - cilevel) / 2
  biggamma   <- (1 + cilevel) / 2

  data.frame(
    lower = mean + qt(smallgamma, df, ncp) * sd / sqrt(N),
    upper = mean + qt(biggamma,   df, ncp) * sd / sqrt(N)
  )
}

ci_hd <- function(x, q, conf = 0.95) {
  hd <- hdquantile(x, probs = q, se = TRUE)
  se <- attr(hd, "se")
  z  <- qnorm(1 - (1 - conf) / 2)
  c(estimate = unname(hd), lower = unname(hd) - z * unname(se), upper = unname(hd) + z * unname(se))
}

# Run one batch of `batch_size` replicates across all n values ------------
run_batch <- function(batch) {
  batch_file <- paste0("normal_nonboot_sims_batch_", batch, ".RData")
  if (file.exists(batch_file)) {
    cat(sprintf("Batch %d / %d already completed (%s); skipping.\n",
                batch, n_batches, batch_file))
    return(invisible(NULL))
  }

  rep_start <- (batch - 1) * batch_size + 1
  rep_end   <- batch * batch_size
  cat(sprintf("\n=== Batch %d / %d: replicates %d-%d, seed = %d ===\n",
              batch, n_batches, rep_start, rep_end, batch_seed(batch)))
  set.seed(batch_seed(batch))

  replicates <- batch_size

  # Simulations -------------------------------------------------------------
  sims_list <- vector("list", nsetting)

  tic()
  for (i in seq_len(nsetting)) {
  n <- n_values[i]
  cat(sprintf("n = %d (%d / %d)\n", n, i, nsetting))

  # PI lower (q = 0.025) ---------------------------------------------------
  np_lower_estimate           <- rep(NA, replicates)
  np_lower_coverage           <- rep(NA, replicates)
  np_lower_finite             <- rep(NA, replicates)
  np_lower_width              <- rep(NA, replicates)
  normal_lower_estimate       <- rep(NA, replicates)
  normal_lower_coverage       <- rep(NA, replicates)
  normal_lower_width          <- rep(NA, replicates)
  normal_exact_lower_coverage <- rep(NA, replicates)
  normal_exact_lower_width    <- rep(NA, replicates)
  hd_lower_estimate           <- rep(NA, replicates)
  hd_lower_coverage           <- rep(NA, replicates)
  hd_lower_width              <- rep(NA, replicates)
  sw_lower_estimate           <- rep(NA, replicates)
  sw_lower_coverage           <- rep(NA, replicates)
  sw_lower_width              <- rep(NA, replicates)

  # PI upper (q = 0.975) ---------------------------------------------------
  np_upper_estimate           <- rep(NA, replicates)
  np_upper_coverage           <- rep(NA, replicates)
  np_upper_finite             <- rep(NA, replicates)
  np_upper_width              <- rep(NA, replicates)
  normal_upper_estimate       <- rep(NA, replicates)
  normal_upper_coverage       <- rep(NA, replicates)
  normal_upper_width          <- rep(NA, replicates)
  normal_exact_upper_coverage <- rep(NA, replicates)
  normal_exact_upper_width    <- rep(NA, replicates)
  hd_upper_estimate           <- rep(NA, replicates)
  hd_upper_coverage           <- rep(NA, replicates)
  hd_upper_width              <- rep(NA, replicates)
  sw_upper_estimate           <- rep(NA, replicates)
  sw_upper_coverage           <- rep(NA, replicates)
  sw_upper_width              <- rep(NA, replicates)

  # Median (q = 0.5) -------------------------------------------------------
  np_estimate     <- rep(NA, replicates)
  np_coverage     <- rep(NA, replicates)
  np_finite       <- rep(NA, replicates)
  np_width        <- rep(NA, replicates)
  normal_estimate <- rep(NA, replicates)
  normal_coverage <- rep(NA, replicates)
  normal_width    <- rep(NA, replicates)
  normal_exact_coverage <- rep(NA, replicates)
  normal_exact_width    <- rep(NA, replicates)
  hd_estimate     <- rep(NA, replicates)
  hd_coverage     <- rep(NA, replicates)
  hd_width        <- rep(NA, replicates)
  sw_estimate     <- rep(NA, replicates)
  sw_coverage     <- rep(NA, replicates)
  sw_width        <- rep(NA, replicates)

  for (j in seq_len(replicates)) {
    outcome <- rnorm(n)

    # Non-parametric PI lower (confintr).
    pi_lower <- ci_quantile(outcome, q = 0.025)
    np_lower_estimate[j] <- pi_lower$estimate
    np_lower_coverage[j] <- as.numeric(pi_lower$interval[1] <= true_lower & true_lower <= pi_lower$interval[2])
    np_lower_finite[j]   <- as.numeric(pi_lower$interval[1] > -Inf)
    np_lower_width[j]    <- pi_lower$interval[2] - pi_lower$interval[1]

    # Normal PI lower (blandr).
    pi_obs <- blandr.statistics(outcome, rep(0, n))
    normal_lower_estimate[j] <- pi_obs$lowerLOA
    normal_lower_coverage[j] <- as.numeric(pi_obs$lowerLOA_lowerCI <= true_lower & true_lower <= pi_obs$lowerLOA_upperCI)
    normal_lower_width[j]    <- pi_obs$lowerLOA_upperCI - pi_obs$lowerLOA_lowerCI

    # Non-central t exact PI lower (Carkeet 2015 / Owen 1968).
    pi_exact <- ci_LoA_exact(outcome)
    normal_exact_lower_coverage[j] <- as.numeric(pi_exact$cilloalow <= true_lower & true_lower <= pi_exact$cilloaupp)
    normal_exact_lower_width[j]    <- pi_exact$cilloaupp - pi_exact$cilloalow

    # Harrell-Davis PI lower.
    hd_ci_lower <- ci_hd(outcome, q = 0.025)
    hd_lower_estimate[j] <- hd_ci_lower["estimate"]
    hd_lower_coverage[j] <- as.numeric(hd_ci_lower["lower"] <= true_lower & true_lower <= hd_ci_lower["upper"])
    hd_lower_width[j]    <- hd_ci_lower["upper"] - hd_ci_lower["lower"]

    # Sasaki-Wang PI lower.
    sw_ci_lower <- tryCatch(exquantile(outcome, q = 0.025), error = function(e) NULL)
    if (!is.null(sw_ci_lower)) {
      sw_lower_estimate[j] <- unname(sw_ci_lower$estimate)
      sw_lower_coverage[j] <- as.numeric(sw_ci_lower$conf_int["lower"] <= true_lower & true_lower <= sw_ci_lower$conf_int["upper"])
      sw_lower_width[j]    <- unname(sw_ci_lower$conf_int["upper"] - sw_ci_lower$conf_int["lower"])
    }

    # Non-parametric PI upper (confintr).
    pi_upper <- ci_quantile(outcome, q = 0.975)
    np_upper_estimate[j] <- pi_upper$estimate
    np_upper_coverage[j] <- as.numeric(pi_upper$interval[1] <= true_upper & true_upper <= pi_upper$interval[2])
    np_upper_finite[j]   <- as.numeric(pi_upper$interval[2] < Inf)
    np_upper_width[j]    <- pi_upper$interval[2] - pi_upper$interval[1]

    # Normal PI upper (blandr) — pi_obs already computed above.
    normal_upper_estimate[j] <- pi_obs$upperLOA
    normal_upper_coverage[j] <- as.numeric(pi_obs$upperLOA_lowerCI <= true_upper & true_upper <= pi_obs$upperLOA_upperCI)
    normal_upper_width[j]    <- pi_obs$upperLOA_upperCI - pi_obs$upperLOA_lowerCI

    # Non-central t exact PI upper — pi_exact already computed above.
    normal_exact_upper_coverage[j] <- as.numeric(pi_exact$ciuloalow <= true_upper & true_upper <= pi_exact$ciuloaupp)
    normal_exact_upper_width[j]    <- pi_exact$ciuloaupp - pi_exact$ciuloalow

    # Harrell-Davis PI upper.
    hd_ci_upper <- ci_hd(outcome, q = 0.975)
    hd_upper_estimate[j] <- hd_ci_upper["estimate"]
    hd_upper_coverage[j] <- as.numeric(hd_ci_upper["lower"] <= true_upper & true_upper <= hd_ci_upper["upper"])
    hd_upper_width[j]    <- hd_ci_upper["upper"] - hd_ci_upper["lower"]

    # Sasaki-Wang PI upper.
    sw_ci_upper <- tryCatch(exquantile(outcome, q = 0.975), error = function(e) NULL)
    if (!is.null(sw_ci_upper)) {
      sw_upper_estimate[j] <- unname(sw_ci_upper$estimate)
      sw_upper_coverage[j] <- as.numeric(sw_ci_upper$conf_int["lower"] <= true_upper & true_upper <= sw_ci_upper$conf_int["upper"])
      sw_upper_width[j]    <- unname(sw_ci_upper$conf_int["upper"] - sw_ci_upper$conf_int["lower"])
    }

    # Non-parametric median (confintr).
    median_ci <- ci_quantile(outcome, q = 0.50)
    np_estimate[j] <- median_ci$estimate
    np_coverage[j] <- as.numeric(median_ci$interval[1] <= true_median & true_median <= median_ci$interval[2])
    np_finite[j]   <- as.numeric(-Inf < median_ci$interval[1] & median_ci$interval[2] < Inf)
    np_width[j]    <- median_ci$interval[2] - median_ci$interval[1]

    # Normal median (t-test).
    ttest_res <- t.test(outcome)
    normal_estimate[j] <- ttest_res$estimate
    normal_coverage[j] <- as.numeric(ttest_res$conf.int[1] <= true_median & true_median <= ttest_res$conf.int[2])
    normal_width[j]    <- ttest_res$conf.int[2] - ttest_res$conf.int[1]

    # Non-central t exact median (q = 0.5; equals the one-sample t interval).
    median_exact <- ci_quantile_exact(outcome, q = 0.50)
    normal_exact_coverage[j] <- as.numeric(median_exact$lower <= true_median & true_median <= median_exact$upper)
    normal_exact_width[j]    <- median_exact$upper - median_exact$lower

    # Harrell-Davis median.
    hd_ci <- ci_hd(outcome, q = 0.50)
    hd_estimate[j] <- hd_ci["estimate"]
    hd_coverage[j] <- as.numeric(hd_ci["lower"] <= true_median & true_median <= hd_ci["upper"])
    hd_width[j]    <- hd_ci["upper"] - hd_ci["lower"]

    # Sasaki-Wang median.
    sw_ci <- tryCatch(suppressWarnings(exquantile(outcome, q = 0.50)), error = function(e) NULL)
    if (!is.null(sw_ci)) {
      sw_estimate[j] <- unname(sw_ci$estimate)
      sw_coverage[j] <- as.numeric(sw_ci$conf_int["lower"] <= true_median & true_median <= sw_ci$conf_int["upper"])
      sw_width[j]    <- unname(sw_ci$conf_int["upper"] - sw_ci$conf_int["lower"])
    }
  }

  sims_list[[i]] <- data.frame(
    n                           = n,
    replicate                   = rep_start - 1 + seq_len(replicates),
    np_lower_estimate           = np_lower_estimate,
    np_lower_coverage           = np_lower_coverage,
    np_lower_finite             = np_lower_finite,
    np_lower_width              = np_lower_width,
    normal_lower_estimate       = normal_lower_estimate,
    normal_lower_coverage       = normal_lower_coverage,
    normal_lower_width          = normal_lower_width,
    normal_exact_lower_coverage = normal_exact_lower_coverage,
    normal_exact_lower_width    = normal_exact_lower_width,
    hd_lower_estimate           = hd_lower_estimate,
    hd_lower_coverage           = hd_lower_coverage,
    hd_lower_width              = hd_lower_width,
    sw_lower_estimate           = sw_lower_estimate,
    sw_lower_coverage           = sw_lower_coverage,
    sw_lower_width              = sw_lower_width,
    np_upper_estimate           = np_upper_estimate,
    np_upper_coverage           = np_upper_coverage,
    np_upper_finite             = np_upper_finite,
    np_upper_width              = np_upper_width,
    normal_upper_estimate       = normal_upper_estimate,
    normal_upper_coverage       = normal_upper_coverage,
    normal_upper_width          = normal_upper_width,
    normal_exact_upper_coverage = normal_exact_upper_coverage,
    normal_exact_upper_width    = normal_exact_upper_width,
    hd_upper_estimate           = hd_upper_estimate,
    hd_upper_coverage           = hd_upper_coverage,
    hd_upper_width              = hd_upper_width,
    sw_upper_estimate           = sw_upper_estimate,
    sw_upper_coverage           = sw_upper_coverage,
    sw_upper_width              = sw_upper_width,
    np_estimate                 = np_estimate,
    np_coverage                 = np_coverage,
    np_finite                   = np_finite,
    np_width                    = np_width,
    normal_estimate             = normal_estimate,
    normal_coverage             = normal_coverage,
    normal_width                = normal_width,
    normal_exact_coverage       = normal_exact_coverage,
    normal_exact_width          = normal_exact_width,
    hd_estimate                 = hd_estimate,
    hd_coverage                 = hd_coverage,
    hd_width                    = hd_width,
    sw_estimate                 = sw_estimate,
    sw_coverage                 = sw_coverage,
    sw_width                    = sw_width
  )
  }
  sim_time <- toc()

  sims <- bind_rows(sims_list)
  save(sims, sim_time, file = batch_file)
  cat("Saved:", batch_file, "\n")
}

for (batch in seq_len(n_batches)) {
  run_batch(batch)
}
