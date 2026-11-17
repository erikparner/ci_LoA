# Faithful R port of the Sasaki & Wang (2022) Stata command exquantile
#
# Based directly on the supplied exquantile.ado source.
#
# Main function:
#   exquantile(y, x = NULL, q = 0.99, k = 0, xval = 999999,
#              id = NULL, time = NULL)
#
# Conditional version:
#   y is the dependent variable
#   x is the conditioning variable
#   id identifies panel units / repeated cross-sectional units
#   time is accepted for interface compatibility; it is not otherwise used
#   by the supplied Mata code.
#
# The Stata program requires xtset when x is supplied. Thus id and time
# are required here for the conditional version. The algorithm itself uses
# id to obtain one y observation per unique id, namely the observation whose
# x is nearest to xval.

exquantile <- function(y, x = NULL, q = 0.99, k = 0,
                       xval = 999999, id = NULL, time = NULL) {

  # ---- Basic validation corresponding to Stata's numeric varlist/options ----
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(q) != 1L || !is.numeric(q) || !is.finite(q))
    stop("q must be a finite scalar.")
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k))
    stop("k must be a finite scalar.")
  if (length(xval) != 1L || !is.numeric(xval) || !is.finite(xval))
    stop("xval must be a finite scalar.")

  conditional <- !is.null(x)

  if (conditional) {
    if (!is.numeric(x)) stop("x must be numeric.")
    if (length(x) != length(y))
      stop("y and x must have the same length.")
    if (is.null(id))
      stop("id must be supplied for the conditional version.")
    if (length(id) != length(y))
      stop("id must have the same length as y.")
    if (is.null(time)) {
      # The supplied Mata code only uses time to obtain xtset's time variable;
      # it does not subsequently use it. Keep the R interface permissive.
      time <- seq_along(y)
    }
    if (length(time) != length(y))
      stop("time must have the same length as y.")
  }

  # Stata marksample touse; for an R vector implementation, complete cases
  # are the natural equivalent for the variables actually used.
  if (conditional) {
    keep <- stats::complete.cases(y, x, id, time)
    y0 <- y[keep]
    x0 <- x[keep]
    id0 <- id[keep]
  } else {
    keep <- is.finite(y)
    y0 <- y[keep]
  }

  if (conditional) {
    # ---------------------------------------------------------------
    # Faithful translation of:
    # uniqueids = uniqrows(id)
    # for each id:
    #   curr_xy = select((x,y),id:==uniqueids[idx])
    #   minindex(abs(curr_xy[.,1] :- xval),1,minidx=.,.)
    #   nnY[idx,1] = curr_xy[minidx[1,1],2]
    # ---------------------------------------------------------------
    #
    # Stata/Mata uniqrows sorts unique values. That ordering does not affect
    # the subsequent estimator, because nnY is sorted before estimation.
    unique_ids <- unique(id0)
    unique_ids <- sort(unique_ids)

    nnY <- numeric(length(unique_ids))

    for (j in seq_along(unique_ids)) {
      ii <- which(id0 == unique_ids[j])
      d <- abs(x0[ii] - xval)

      # Mata minindex(...,1,...) returns the first minimum in case of a tie.
      jj <- ii[which.min(d)]
      nnY[j] <- y0[jj]
    }

    y_used <- nnY
  } else {
    y_used <- y0
  }

  n <- length(y_used)

  if (n < 3L)
    stop("Too few observations.")

  # The supplied Stata code does not guard against invalid k values.
  # These checks prevent R from silently producing nonsensical indexing.
  if (k <= 0) {
    k_used <- floor(0.05 * n^0.99)
  } else {
    k_used <- k
  }

  # Stata's Mata indices are integer indices. Use the integer value supplied
  # by the user, matching the intended Stata usage.
  if (length(k_used) != 1L || !is.finite(k_used))
    stop("Invalid k.")
  k_used <- as.integer(k_used)

  if (k_used < 1L || k_used >= n)
    stop("k must satisfy 1 <= k < n.")

  # Median normalization:
  # normy = y :- sort(y,1)[floor(n/2),1]
  ys <- sort(y_used)
  med_index <- floor(n / 2)
  med_value <- ys[med_index]
  normy <- y_used - med_value
  normys <- sort(normy)

  if (q <= 0 || q >= 1)
    stop("q must be strictly between 0 and 1.")

  if (0.05 < q && q < 0.95) {
    warning("q is not an extreme quantile (the Stata command prints a warning).")
  }

  # ---------------------------------------------------------------
  # Faithful translation of the estimator
  # ---------------------------------------------------------------

  # Mata:
  # hill_right = mean(log(sort(normy)[(n-k+1..n)] /
  #                       sort(normy)[n-k+1]))
  upper_idx <- (n - k_used + 1L):n
  upper_threshold_idx <- n - k_used + 1L

  hill_right <- mean(
    log(normys[upper_idx] / normys[upper_threshold_idx])
  )

  # Mata:
  # hill_quantile_right_normalized =
  #   sort(normy)[n-k] * (k/(n*(1-q)))^hill_right
  #
  # Note carefully: n-k is the R index n-k (not n-k+1).
  right_base_idx <- n - k_used

  hill_quantile_right_normalized <-
    normys[right_base_idx] *
    (k_used / (n * (1 - q)))^hill_right

  # Add back the same median used for normalization.
  hill_quantile_right <-
    hill_quantile_right_normalized + med_value

  hill <- hill_right
  hill_quantile <- hill_quantile_right

  # Left-tail calculation. The supplied Stata code uses hill_right in the
  # exponent in BOTH lines 120/121, rather than hill_left. This is preserved
  # deliberately because the objective here is a faithful port.
  left_idx <- 1L:k_used
  left_threshold_idx <- k_used

  hill_left <- mean(
    log(normys[left_idx] / normys[left_threshold_idx])
  )

  left_base_idx <- k_used + 1L

  hill_quantile_left_normalized <-
    normys[left_base_idx] *
    (k_used / (n * q))^hill_right

  hill_quantile_left <-
    hill_quantile_left_normalized + med_value

  if (q < 0.5) {
    hill <- hill_left
    hill_quantile <- hill_quantile_left
  }

  # ---------------------------------------------------------------
  # Faithful translation of the variance calculation
  # ---------------------------------------------------------------

  hill_quantile_right_variance <-
    hill_right^2 *
    hill_quantile_right_normalized^2 *
    (log(k_used / n / (1 - q)))^2 /
    k_used

  hill_quantile_variance <- hill_quantile_right_variance

  hill_quantile_left_variance <-
    hill_left^2 *
    hill_quantile_left_normalized^2 *
    (log(k_used / n / q))^2 /
    k_used

  if (q < 0.5) {
    hill_quantile_variance <- hill_quantile_left_variance
  }

  # Stata returns b = hill_quantile and V = hill_quantile_variance.
  # Add useful R-side summaries without altering the underlying calculation.
  se <- sqrt(hill_quantile_variance)

  # Stata's ereturn display uses the usual 95% normal approximation.
  ci <- c(
    lower = hill_quantile - stats::qnorm(0.975) * se,
    upper = hill_quantile + stats::qnorm(0.975) * se
  )

  structure(
    list(
      estimate = unname(hill_quantile),
      variance = unname(hill_quantile_variance),
      std_error = unname(se),
      conf_int = ci,
      q = q,
      k = k_used,
      n = n,
      conditional = conditional,
      xval = if (conditional) xval else NULL,
      hill = unname(hill),
      hill_right = unname(hill_right),
      hill_left = unname(hill_left),
      median = unname(med_value),
      call = match.call()
    ),
    class = "exquantile"
  )
}


print.exquantile <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...) {
  cat("\n")
  if (isTRUE(x$conditional)) {
    cat(sprintf("Conditional %.6f-th Extreme Quantile Given X=%.6f\n",
                x$q, x$xval))
  } else {
    cat(sprintf("Unconditional %.6f-th Extreme Quantile\n", x$q))
  }
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("           k: %d\n", x$k))
  cat(sprintf("           q: %.6f\n", x$q))
  if (isTRUE(x$conditional))
    cat(sprintf("           x: %.6f\n", x$xval))
  cat("\n")
  cat(sprintf("Quantile: %.*f\n", digits, x$estimate))
  cat(sprintf("Std. Err.: %.*f\n", digits, x$std_error))
  cat(sprintf("95%% CI: [%.*f, %.*f]\n",
              digits, x$conf_int[["lower"]], digits, x$conf_int[["upper"]]))
  invisible(x)
}


coef.exquantile <- function(object, ...) {
  object$estimate
}

vcov.exquantile <- function(object, ...) {
  matrix(object$variance, 1L, 1L,
         dimnames = list("Quantile", "Quantile"))
}


# -------------------------------------------------------------------------
# Optional convenience function for reproducing the Stata-style calculation
# without the print method:
# -------------------------------------------------------------------------
summary.exquantile <- function(object, ...) {
  out <- c(
    Quantile = object$estimate,
    `Std. Error` = object$std_error,
    `95% CI lower` = object$conf_int[["lower"]],
    `95% CI upper` = object$conf_int[["upper"]],
    n = object$n,
    k = object$k,
    q = object$q
  )
  if (isTRUE(object$conditional))
    out <- c(out, xval = object$xval)
  out
}
