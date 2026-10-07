# ---------------------------------------------------------------------------
# Chiang II
#
# Chiang II：analytic sampling variance of ex, propagated from binomial death
# counts. It answers "how precise is this estimate given the number of deaths
# observed", and it is the right tool for a descriptive life table built from
# registered deaths.
#
# It is NOT the right tool for a life table built from modelled or smoothed
# death counts. Posterior means from a Bayesian model are not observed counts;
# feeding them in here produces an interval that describes a sampling process
# that never happened, and it will typically be far too narrow because the
# smoothing has already removed the variability. For modelled input, build the
# life table across posterior draws and take quantiles of ex instead.
# ---------------------------------------------------------------------------


#' Chiang standard errors and confidence intervals for life expectancy
#'
#' The variance of `ex` accumulates the contribution of each age group's
#' binomial sampling error in `qx`, weighted by how much that age group moves
#' life expectancy:
#'
#' \deqn{Var(e_x) = \frac{1}{l_x^2} \sum_{i \ge x} l_i^2 [n_i - a_i + e_{i+1}]^2 Var(q_i)}
#'
#' with \eqn{Var(q_i) = q_i^2 (1 - q_i) / D_i} for closed age groups. Here
#' \eqn{a_i} is in years, as everywhere in this package; Chiang's original
#' \eqn{(1 - a_i) n_i} is the same quantity with \eqn{a_i} as a fraction of the
#' interval. The
#' open-ended group is handled separately: its term reduces to
#' \eqn{l_n^2 e_n^2 / D_n}, the variance of a mean survival time estimated from
#' \eqn{D_n} deaths.
#'
#' That open-group term assumes \eqn{e_n = 1/m_n}, which holds under the
#' `"constant_hazard"` closeout but not under `"kannisto"` or `"gompertz"`.
#' The function checks this directly, so a supplied `ax` that satisfies it
#' passes silently. Otherwise it warns rather than refuses, since the resulting
#' inconsistency is small relative to the interval width, but the cleanest
#' combination is `closeout = "constant_hazard"`.
#'
#' Age groups with no deaths contribute zero variance rather than `NaN`. That is
#' a convenience. Mortality in such a group is genuinely unknown, and treating it
#' as known understates the interval. This matters most in small-area work, where
#' empty cells are common and are usually the reason for reaching for a smoothing
#' model in the first place.
#'
#' @param x A `letools_lt` object.
#' @param deaths Observed death counts, matching `x` in shape: either a vector
#'   of length `n_age` or an `n_age x n_draw` matrix.
#' @param conf Confidence level.
#' @param truncate_at_zero Clamp the lower limit at zero, since life
#'   expectancy cannot be negative.
#' @return A list of `n_age x n_draw` matrices: `se`, `lower`, `upper`, plus the
#'   variance components `var_qx` and `var_Tx`.
#'
#' @references
#' Chiang CL (1984). *The Life Table and its Applications*. Krieger.
#'
#' Eayres D, Williams ES (2004). Evaluation of methodologies for small area
#' life expectancy estimation. *J Epidemiol Community Health* 58(3):243-249.
#' @export
lt_chiang_ci <- function(x, deaths, conf = 0.95, truncate_at_zero = TRUE) {

  if (!inherits(x, "letools_lt")) {
    stop("`x` must be a life table from lt().", call. = FALSE)
  }
  if (conf <= 0 || conf >= 1) {
    stop("`conf` must be strictly between 0 and 1.", call. = FALSE)
  }

  n <- nrow(x$ex)
  m <- ncol(x$ex)

  D <- if (is.matrix(deaths)) deaths else matrix(deaths, ncol = 1L)
  if (ncol(D) == 1L && m > 1L) D <- D[, rep(1L, m), drop = FALSE]

  if (!identical(dim(D), c(n, m))) {
    stop(
      sprintf(
        "`deaths` is %d x %d but the life table is %d x %d.",
        nrow(D), ncol(D), n, m
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(D)) || any(D < 0)) {
    stop("`deaths` must be finite and non-negative.", call. = FALSE)
  }

  # Check the assumption itself rather than the closeout label: a supplied ax
  # (a published table, say) can satisfy e_open = 1/mx just as well.
  if (any(abs(x$ax[n, ] * x$mx[n, ] - 1) > 1e-8)) {
    warning(
      "Chiang's open-group variance assumes e_open = 1/mx, which this life ",
      'table does not satisfy (closeout = "', attr(x, "closeout"), '"). The ',
      "interval remains usable but is internally inconsistent at the oldest ",
      "age group.",
      call. = FALSE
    )
  }

  nx <- matrix(attr(x, "nx"), nrow = n, ncol = m)

  # --- variance of qx ---
  var_qx <- matrix(0, nrow = n, ncol = m)

  closed <- seq_len(n - 1L)
  ok <- D[closed, , drop = FALSE] > 0
  vq <- x$qx[closed, , drop = FALSE]^2 *
    (1 - x$qx[closed, , drop = FALSE]) / D[closed, , drop = FALSE]
  vq[!ok] <- 0
  var_qx[closed, ] <- vq

  # open group: variance of a mean survival time estimated from D_n deaths
  ok_open <- D[n, ] > 0 & x$mx[n, ] > 0
  vq_open <- numeric(m)
  vq_open[ok_open] <- 4 / (D[n, ok_open] * x$mx[n, ok_open]^2)
  var_qx[n, ] <- vq_open

  # --- weight each age group's error by its leverage on ex ---
  ex_lead <- rbind(x$ex[-1, , drop = FALSE], matrix(0, nrow = 1L, ncol = m))

  w <- matrix(0, nrow = n, ncol = m)
  w[closed, ] <- var_qx[closed, , drop = FALSE] *
    x$lx[closed, , drop = FALSE]^2 *
    (nx[closed, , drop = FALSE] - x$ax[closed, , drop = FALSE] +
       ex_lead[closed, , drop = FALSE])^2
  w[n, ] <- (x$lx[n, ] / 2)^2 * var_qx[n, ]

  var_Tx <- apply_rev_cumsum(w)

  se <- matrix(0, nrow = n, ncol = m)
  pos <- x$lx > 0
  se[pos] <- sqrt(var_Tx[pos] / x$lx[pos]^2)

  z <- stats::qnorm(1 - (1 - conf) / 2)
  lower <- x$ex - z * se
  upper <- x$ex + z * se

  if (isTRUE(truncate_at_zero)) lower[lower < 0] <- 0

  list(se = se, lower = lower, upper = upper,
       var_qx = var_qx, var_Tx = var_Tx, conf = conf)
}
