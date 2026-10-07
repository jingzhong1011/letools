# ---------------------------------------------------------------------------
# ax
#
# ax：average time lived within an age interval by those who die in it.
# Assuming deaths fall uniformly (ax = n/2) is badly wrong at both ends of the
# age range: infant deaths cluster in the first days of life, and at old ages
# mortality rises steeply enough within a five-year band that deaths are
# front-loaded. The two ends are handled separately, and the middle can
# optionally be refined by iteration.
# ---------------------------------------------------------------------------


#' Normalise a sex argument to 1 (male) or 2 (female)
#'
#' @param sex `1`, `2`, `"male"`, `"female"`, `"m"` or `"f"`.
#' @return `1L` or `2L`.
#' @export
normalise_sex <- function(sex) {
  if (length(sex) != 1L) {
    stop("`sex` must be a single value.", call. = FALSE)
  }
  if (is.character(sex)) {
    sex <- switch(
      tolower(sex),
      male = 1L, m = 1L, female = 2L, f = 2L,
      stop('`sex` must be 1, 2, "male" or "female". Got "', sex, '".',
           call. = FALSE)
    )
  }
  sex <- as.integer(sex)
  if (!sex %in% c(1L, 2L)) {
    stop("`sex` must be 1 (male) or 2 (female). Got ", sex, ".", call. = FALSE)
  }
  sex
}


#' Coale-Demeny separation factors for infancy and early childhood
#'
#' @param m0 Infant mortality rate; a numeric vector, one value per draw.
#' @param sex `1`/`"male"` or `2`/`"female"`.
#' @return A two-row matrix, rows `a0` and `a1`, one column per draw.
#'
#' @references
#' Coale AJ, Demeny P, Vaughan B (1983). *Regional Model Life Tables and Stable
#' Populations*, 2nd ed. Academic Press.
#'
#' Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
#' Modeling Population Processes*, p. 48. Blackwell.
#' @export
ax_coale_demeny <- function(m0, sex) {
  sex <- normalise_sex(sex)
  m0  <- as.numeric(m0)

  if (any(!is.finite(m0)) || any(m0 < 0)) {
    stop("`m0` must be finite and non-negative.", call. = FALSE)
  }

  high <- m0 >= 0.107

  if (sex == 1L) {
    a0 <- ifelse(high, 0.330, 0.045 + 2.684 * m0)
    a1 <- ifelse(high, 1.352, 1.651 - 2.816 * m0)
  } else {
    a0 <- ifelse(high, 0.350, 0.053 + 2.800 * m0)
    a1 <- ifelse(high, 1.361, 1.522 - 1.518 * m0)
  }

  rbind(a0 = a0, a1 = a1)
}


#' Old-age separation factors by single-year extrapolation
#'
#' The two supported old-age models work identically and differ only in what is
#' taken to be linear in age:
#'
#' \describe{
#'   \item{`"kannisto"`}{`logit(qx)` is linear in age, a logistic hazard.
#'     Mortality decelerates at the highest ages, which is what the empirical
#'     record above 95 shows.}
#'   \item{`"gompertz"`}{`log(mu)` is linear in age. Mortality keeps rising
#'     exponentially without limit, which overstates it at the highest ages.}
#' }
#'
#' The procedure is the same for both: take a constant hazard within each
#' observed old-age band, fit the chosen model by least squares to one
#' single-year probability per band at the band midpoint, extrapolate year by
#' year past the open-group boundary, then aggregate the implied deaths back to
#' the original bands. Each band's separation factor is the death-weighted mean
#' time lived within it; for the open-ended band, that quantity is life
#' expectancy at the boundary.
#'
#' Note what the fit actually sees: three closed bands give three points. One
#' point per band is deliberate. Replicating each band's value across its single
#' years would add variance to the predictor without adding information to the
#' response, attenuating the slope by 25/28 for three five-year bands and
#' overstating life expectancy.
#'
#' @param lx Survivorship at the old-age band boundaries, ages x draws. The
#'   last row is the open-group boundary.
#' @param age Ages for the rows of `lx`; must be evenly spaced.
#' @param method `"kannisto"` or `"gompertz"`.
#' @param max_age Age at which the single-year extrapolation stops.
#' @param b_min Lower bound on the fitted slope. A flat or negative slope implies
#'   old-age mortality that does not rise with age, which is a small-numbers
#'   artefact and makes the closeout diverge.
#' @return A list with `ax` (one row per band, one column per draw; the last row
#'   is life expectancy at the open-group boundary), the fitted coefficients `a`
#'   and `b`, and the `method` used.
#'
#' @references
#' Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at Ages 80
#' to 120*. Odense University Press.
#' @export
ax_old_age <- function(lx,
                       age,
                       method = c("kannisto", "gompertz"),
                       max_age = 130,
                       b_min = 0.02) {

  method <- match.arg(method)
  lx <- as.matrix(lx)

  n_band <- nrow(lx)
  m      <- ncol(lx)

  if (length(age) != n_band) {
    stop("`age` must have one value per row of `lx`.", call. = FALSE)
  }
  if (n_band < 4L) {
    stop(
      "Old-age extrapolation needs at least four band boundaries (three ",
      "usable intervals). Got ", n_band, ".",
      call. = FALSE
    )
  }

  width_all <- diff(age)
  if (length(unique(width_all)) != 1L) {
    stop(
      "Old-age bands must be evenly spaced; got widths ",
      paste(width_all, collapse = ", "), ".",
      call. = FALSE
    )
  }
  width <- width_all[1]

  if (max_age <= age[n_band]) {
    stop(
      sprintf("`max_age` (%s) must exceed the open-group boundary (%s).",
              max_age, age[n_band]),
      call. = FALSE
    )
  }
  if (any(lx <= 0)) {
    stop(
      "Survivorship must be strictly positive at old ages. A zero usually ",
      "means the radix was too small relative to the mortality level.",
      call. = FALSE
    )
  }

  n_closed <- n_band - 1L

  # --- constant hazard within each observed band ---
  mu_band <- (log(lx[-n_band, , drop = FALSE]) -
                log(lx[-1, , drop = FALSE])) / width

  if (any(!is.finite(mu_band)) || any(mu_band <= 0)) {
    stop(
      "Non-positive force of mortality at old ages. Check for survivorship ",
      "that is flat or increasing with age.",
      call. = FALSE
    )
  }

  band_of <- rep(seq_len(n_closed), each = width)
  offset  <- rep(seq_len(width) - 1L, times = n_closed)
  age_obs <- rep(age[-n_band], each = width) + offset

  lx_obs <- lx[band_of, , drop = FALSE] *
    exp(-mu_band[band_of, , drop = FALSE] * offset)
  qx_obs <- 1 - exp(-mu_band[band_of, , drop = FALSE])
  dx_obs <- lx_obs * qx_obs

  # --- fit the chosen model, one observation per band ---
  qx_fit <- 1 - exp(-mu_band)
  xv <- age[-n_band] + width / 2

  link <- function(q) {
    switch(method,
           kannisto = log(q / (1 - q)),
           gompertz = log(-log(1 - q)))
  }

  y <- link(qx_fit)

  if (any(!is.finite(y))) {
    stop(
      "Non-finite values when transforming old-age probabilities; mortality ",
      "may be zero or one in some draws.",
      call. = FALSE
    )
  }

  xbar <- mean(xv)
  xc   <- xv - xbar

  b <- as.numeric(crossprod(xc, y) / sum(xc^2))
  b <- pmax(b, b_min)

  a <- colMeans(y) - b * xbar

  # --- extrapolate past the open-group boundary ---
  age_ext <- seq.int(age[n_band], max_age - 1L)
  n_ext <- length(age_ext)

  eta <- outer(age_ext + 0.5, b) + rep(a, each = n_ext)
  qx_ext <- switch(
    method,
    kannisto = 1 / (1 + exp(-eta)),
    gompertz = 1 - exp(-exp(eta))
  )
  qx_ext <- pmin(pmax(qx_ext, 1e-12), 1 - 1e-12)

  lx_ext <- matrix(0, nrow = n_ext, ncol = m)
  lx_ext[1, ] <- lx[n_band, ]
  if (n_ext > 1L) {
    surv <- apply_cumprod(1 - qx_ext[-n_ext, , drop = FALSE])
    lx_ext[-1, ] <- sweep(surv, 2L, lx[n_band, ], "*")
  }
  dx_ext <- lx_ext * qx_ext

  # --- aggregate deaths back to the original bands ---
  age_all <- c(age_obs, age_ext)
  dx_all <- rbind(dx_obs, dx_ext)

  band_index <- findInterval(age_all, age)
  G <- matrix(0, nrow = length(age_all), ncol = n_band)
  G[cbind(seq_along(age_all), band_index)] <- 1

  yl  <- age_all - age[band_index] + 0.5
  num <- crossprod(G, dx_all * yl)
  den <- crossprod(G, dx_all)

  ax <- matrix(width / 2, nrow = n_band, ncol = m)
  usable <- den > 0
  ax[usable] <- num[usable] / den[usable]

  list(ax = ax, a = a, b = b, method = method)
}


#' Keyfitz iterative refinement of separation factors
#'
#' Where deaths are not evenly distributed across neighbouring age groups, the
#' average time lived within an interval shifts away from its midpoint. This
#' applies the correction
#' \deqn{a_x = n/2 + (n/24)(d_{x+n} - d_{x-n}) / d_x}
#' in its usual form.
#'
#' Results are clamped to `[0, n]`. Values outside that range almost always come
#' from an age group with too few deaths to say anything.
#'
#' @param dx Deaths by age, ages x draws.
#' @param nx Interval widths, one per age group.
#' @param rows Rows to refine. All other rows are returned unchanged.
#' @param ax Current separation factors.
#' @return A matrix the same shape as `ax`.
#'
#' @references
#' Keyfitz N (1966). A life table that agrees with the data. *Journal of the
#' American Statistical Association* 61(314):305-312.
#' @export
ax_keyfitz <- function(dx, nx, rows, ax) {
  out <- ax
  if (!length(rows)) return(out)

  n <- nrow(dx)
  for (k in rows) {
    if (k <= 1L || k >= n) next

    d_prev <- dx[k - 1L, ]
    d_here <- dx[k, ]
    d_next <- dx[k + 1L, ]

    new <- nx[k] / 2 + (nx[k] / 24) * (d_next - d_prev) / d_here
    new[!is.finite(new) | d_here <= 0] <- nx[k] / 2

    out[k, ] <- pmin(pmax(new, 0), nx[k])
  }
  out
}


#' Column-wise cumulative product that always returns a matrix
#' @noRd
apply_cumprod <- function(x) {
  out <- apply(x, 2L, cumprod)
  if (!is.matrix(out)) out <- matrix(out, nrow = nrow(x), ncol = ncol(x))
  out
}


#' Column-wise cumulative sum from the bottom up
#' @noRd
apply_rev_cumsum <- function(x) {
  n <- nrow(x)
  out <- apply(x[n:1, , drop = FALSE], 2L, cumsum)
  if (!is.matrix(out)) out <- matrix(out, nrow = n, ncol = ncol(x))
  out[n:1, , drop = FALSE]
}
