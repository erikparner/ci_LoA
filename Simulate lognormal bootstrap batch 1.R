library(tidyverse)
library(boot)    # boot, boot.ci
library(tictoc)

options(warn = -1)


# Settings ----------------------------------------------------------------
# Batch number for this run: 1, 2, 3, ... Each batch is independent and is
# meant to be run as its own script invocation (e.g. on a separate
# computer); change this value and re-run to produce a different batch.
batch <- 1

batch_size      <- 1000    # replicates per batch
boot_replicates <- 1100    # must exceed max(n) for bca to be valid
boot_types      <- c("norm", "perc", "bca")

# Batch k covers replicates ((k - 1) * batch_size + 1):(k * batch_size).
# Batch 1 uses seed 1; batch k > 1 uses seed (k - 1) * batch_size, matching
# the seed scheme used by "Simulate lognormal non-bootstrap.R".
batch_seed <- function(batch) if (batch == 1L) 1L else (batch - 1L) * batch_size

# True values for the standard log-normal distribution (meanlog=0, sdlog=1)
true_lower  <- qlnorm(0.025)
true_upper  <- qlnorm(0.975)
true_median <- qlnorm(0.5)    # = exp(0) = 1

# Combined n grid (union of PI boundary 145/146 and median boundary 5/6)
n_values <- sort(unique(c(seq(10, 1000, by = 10), 5, 6, 145, 146)))
nsetting <- length(n_values)

# Output file for this batch -----------------------------------------------
batch_file <- paste0("lognormal_boot_sims_batch_", batch, ".RData")
if (file.exists(batch_file)) {
  stop("Batch ", batch, " already completed (", batch_file, "); ",
       "change `batch` above to run a different batch.")
}

rep_start <- (batch - 1) * batch_size + 1
rep_end   <- batch * batch_size
cat(sprintf("=== Batch %d: replicates %d-%d, seed = %d ===\n",
            batch, rep_start, rep_end, batch_seed(batch)))
set.seed(batch_seed(batch))

replicates <- batch_size

# Simulations -------------------------------------------------------------
sims_list <- vector("list", nsetting)

tic()
for (i in seq_len(nsetting)) {
  n <- n_values[i]
  cat(sprintf("n = %d (%d / %d)\n", n, i, nsetting))

  # PI lower (q = 0.025) ---------------------------------------------------
  boot_norm_lower_coverage <- rep(NA, replicates)
  boot_norm_lower_width    <- rep(NA, replicates)
  boot_perc_lower_coverage <- rep(NA, replicates)
  boot_perc_lower_width    <- rep(NA, replicates)
  boot_bca_lower_coverage  <- rep(NA, replicates)
  boot_bca_lower_width     <- rep(NA, replicates)

  # PI upper (q = 0.975) ---------------------------------------------------
  boot_norm_upper_coverage <- rep(NA, replicates)
  boot_norm_upper_width    <- rep(NA, replicates)
  boot_perc_upper_coverage <- rep(NA, replicates)
  boot_perc_upper_width    <- rep(NA, replicates)
  boot_bca_upper_coverage  <- rep(NA, replicates)
  boot_bca_upper_width     <- rep(NA, replicates)

  # Median (q = 0.5) -------------------------------------------------------
  boot_norm_coverage <- rep(NA, replicates)
  boot_norm_width    <- rep(NA, replicates)
  boot_perc_coverage <- rep(NA, replicates)
  boot_perc_width    <- rep(NA, replicates)
  boot_bca_coverage  <- rep(NA, replicates)
  boot_bca_width     <- rep(NA, replicates)

  for (j in seq_len(replicates)) {
    outcome <- rlnorm(n)

    # Single resample for PI lower (q=0.025), median (q=0.5), PI upper (q=0.975).
    boot_out    <- boot(outcome,
                        statistic = function(x, id) quantile(x[id], probs = c(0.025, 0.5, 0.975), names = FALSE),
                        R = boot_replicates, parallel = "multicore")
    boot_ci_lower  <- boot.ci(boot_out, conf = 0.95, type = boot_types, index = 1)
    boot_ci_med    <- boot.ci(boot_out, conf = 0.95, type = boot_types, index = 2)
    boot_ci_upper  <- boot.ci(boot_out, conf = 0.95, type = boot_types, index = 3)

    boot_norm_lower_coverage[j] <- as.numeric(boot_ci_lower$normal[2] <= true_lower & true_lower <= boot_ci_lower$normal[3])
    boot_norm_lower_width[j]    <- boot_ci_lower$normal[3] - boot_ci_lower$normal[2]
    boot_norm_upper_coverage[j] <- as.numeric(boot_ci_upper$normal[2] <= true_upper & true_upper <= boot_ci_upper$normal[3])
    boot_norm_upper_width[j]    <- boot_ci_upper$normal[3] - boot_ci_upper$normal[2]

    boot_perc_lower_coverage[j] <- as.numeric(boot_ci_lower$percent[1, 4] <= true_lower & true_lower <= boot_ci_lower$percent[1, 5])
    boot_perc_lower_width[j]    <- boot_ci_lower$percent[1, 5] - boot_ci_lower$percent[1, 4]
    boot_perc_upper_coverage[j] <- as.numeric(boot_ci_upper$percent[1, 4] <= true_upper & true_upper <= boot_ci_upper$percent[1, 5])
    boot_perc_upper_width[j]    <- boot_ci_upper$percent[1, 5] - boot_ci_upper$percent[1, 4]

    boot_bca_lower_coverage[j]  <- as.numeric(boot_ci_lower$bca[1, 4] <= true_lower & true_lower <= boot_ci_lower$bca[1, 5])
    boot_bca_lower_width[j]     <- boot_ci_lower$bca[1, 5] - boot_ci_lower$bca[1, 4]
    boot_bca_upper_coverage[j]  <- as.numeric(boot_ci_upper$bca[1, 4] <= true_upper & true_upper <= boot_ci_upper$bca[1, 5])
    boot_bca_upper_width[j]     <- boot_ci_upper$bca[1, 5] - boot_ci_upper$bca[1, 4]

    boot_norm_coverage[j] <- as.numeric(boot_ci_med$normal[2] <= true_median & true_median <= boot_ci_med$normal[3])
    boot_norm_width[j]    <- boot_ci_med$normal[3] - boot_ci_med$normal[2]

    boot_perc_coverage[j] <- as.numeric(boot_ci_med$percent[1, 4] <= true_median & true_median <= boot_ci_med$percent[1, 5])
    boot_perc_width[j]    <- boot_ci_med$percent[1, 5] - boot_ci_med$percent[1, 4]

    boot_bca_coverage[j]  <- as.numeric(boot_ci_med$bca[1, 4] <= true_median & true_median <= boot_ci_med$bca[1, 5])
    boot_bca_width[j]     <- boot_ci_med$bca[1, 5] - boot_ci_med$bca[1, 4]
  }

  sims_list[[i]] <- data.frame(
    n                        = n,
    replicate                = rep_start - 1 + seq_len(replicates),
    boot_norm_lower_coverage = boot_norm_lower_coverage,
    boot_norm_lower_width    = boot_norm_lower_width,
    boot_perc_lower_coverage = boot_perc_lower_coverage,
    boot_perc_lower_width    = boot_perc_lower_width,
    boot_bca_lower_coverage  = boot_bca_lower_coverage,
    boot_bca_lower_width     = boot_bca_lower_width,
    boot_norm_upper_coverage = boot_norm_upper_coverage,
    boot_norm_upper_width    = boot_norm_upper_width,
    boot_perc_upper_coverage = boot_perc_upper_coverage,
    boot_perc_upper_width    = boot_perc_upper_width,
    boot_bca_upper_coverage  = boot_bca_upper_coverage,
    boot_bca_upper_width     = boot_bca_upper_width,
    boot_norm_coverage       = boot_norm_coverage,
    boot_norm_width          = boot_norm_width,
    boot_perc_coverage       = boot_perc_coverage,
    boot_perc_width          = boot_perc_width,
    boot_bca_coverage        = boot_bca_coverage,
    boot_bca_width           = boot_bca_width
  )
}
sim_time <- toc()

sims <- bind_rows(sims_list)
save(sims, sim_time, file = batch_file)
cat("Saved:", batch_file, "\n")
