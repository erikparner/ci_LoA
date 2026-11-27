library(tidyverse)

# -------------------------------------------------------------------------
# For each distribution:
#   1. Load every "<dist>_nonboot_sims_batch_<k>.RData" file found in the
#      `data_dir` subfolder (with k <= max_batch) and row-bind them into
#      one non-bootstrap data frame (replicate numbers are globally unique
#      across batches by design, see "Simulate <dist> non-bootstrap.R").
#   2. Load every "<dist>_boot_sims_batch_<k>.RData" file found in
#      `data_dir` (with k <= max_batch; produced by "Simulate <dist>
#      bootstrap batch.R", often run on separate computers) and row-bind
#      them the same way.
#   3. Join on (n, replicate) and save as "<data_dir>/<dist>_combined_sims.RData".
# Only the batches/replicates that currently exist (and are <= max_batch)
# are combined; re-run this script after more batch files have been
# produced to refresh the combined files.
# -------------------------------------------------------------------------

# Settings ------------------------------------------------------------------
# Subfolder (relative to the working directory) holding the batch files and
# where the combined files are saved.
data_dir <- "Data"

# Highest batch number to include when combining. Bootstrap batches are far
# more expensive than non-bootstrap batches and so typically lag behind;
# set this to the highest batch number completed on the bootstrap side (or
# Inf to include every batch file found) so the two sides line up.
max_batch <- 25

# Load and row-bind every "<prefix>_batch_<k>.RData" file matching `pattern`
# under `data_dir`, with k <= max_batch, ordered by batch number. Returns
# NULL (with a message) if none are found.
load_batches <- function(pattern, label, data_dir, max_batch = Inf) {
  batch_files <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
  if (length(batch_files) == 0) {
    cat(sprintf("  No %s batch files found; skipping.\n", label))
    return(NULL)
  }
  batch_numbers <- as.integer(sub(".*_batch_(\\d+)\\.RData$", "\\1", batch_files))

  keep         <- batch_numbers <= max_batch
  batch_files   <- batch_files[keep]
  batch_numbers <- batch_numbers[keep]
  if (length(batch_files) == 0) {
    cat(sprintf("  No %s batch files with batch <= %s; skipping.\n", label, max_batch))
    return(NULL)
  }
  batch_files <- batch_files[order(batch_numbers)]

  sims_list <- vector("list", length(batch_files))
  for (k in seq_along(batch_files)) {
    load(batch_files[k])   # loads `sims`, `sim_time`
    sims_list[[k]] <- sims
  }
  sims_out <- bind_rows(sims_list)
  cat(sprintf("  %s: %d batch file(s) (<= batch %s), %d replicates, %d rows\n",
              label, length(batch_files), max_batch,
              length(unique(sims_out$replicate)), nrow(sims_out)))
  sims_out
}

distributions <- c("normal", "lognormal", "chisq4", "mixture")

for (dist in distributions) {
  cat("\n---", dist, "---\n")

  sims_nonboot <- load_batches(
    pattern   = paste0("^", dist, "_nonboot_sims_batch_\\d+\\.RData$"),
    label     = "Non-bootstrap",
    data_dir  = data_dir,
    max_batch = max_batch
  )
  if (is.null(sims_nonboot)) next

  sims_boot <- load_batches(
    pattern   = paste0("^", dist, "_boot_sims_batch_\\d+\\.RData$"),
    label     = "Bootstrap",
    data_dir  = data_dir,
    max_batch = max_batch
  )
  if (is.null(sims_boot)) next

  sims_combined <- full_join(sims_nonboot, sims_boot, by = c("n", "replicate"))

  cat(sprintf("  Combined: %d rows, %d columns\n", nrow(sims_combined), ncol(sims_combined)))

  out_file <- file.path(data_dir, paste0(dist, "_combined_sims.RData"))
  save(sims_combined, file = out_file)
  cat("  Saved:", out_file, "\n")
}

cat("\nDone.\n")


