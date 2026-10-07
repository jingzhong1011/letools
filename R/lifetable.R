# ---------------------------------------------------------------------------
# lt
#
# lt：abridged life table built from age-specific mortality rates.
# Everything is an n_age x n_draw matrix. A single life table is just the
# one-column case, so descriptive analysis and posterior samples run through
# the same code with no branching.
# ---------------------------------------------------------------------------


#' Validate an abridged age structure
#'
#' @param age Left-hand boundaries of the age groups, strictly ascending,
#'   starting at 0.
#' @return Invisibly, the interval widths, with `NA` for the open group.
#' @export
validate_age <- function(age) {
  if (length(age) < 3L) {
    stop("An abridged life table needs at least three age groups.", call. = FALSE)
  }
  if (any(!is.finite(age)) || any(age < 0)) {
    stop("`age` must be finite and non-negative.", call. = FALSE)
  }
  if (is.unsorted(age, strictly = TRUE)) {
    stop("`age` must be strictly increasing.", call. = FALSE)
  }
  if (age[1] != 0) {
    stop(
      "The first age group must start at 0; life expectancy at birth is not ",
      "defined otherwise. Got ", age[1], ".",
      call. = FALSE
    )
  }
  invisible(c(diff(age), NA_real_))
}


#' Build abridged life tables from age-specific mortality rates
#'
#' @param mx Mortality rates: a numeric vector of length `n_age`, or an
#'   `n_age x n_draw` matrix, ordered youngest to oldest with the open group
#'   last.
#' @param sex `1`/`"male"` or `2`/`"female"`, used for the Coale-Demeny
#'   separation factors at ages 0 and 1-4. The 1-4 factor is applied only when
#'   the second age group is exactly 1-4; otherwise that group starts at the
#'   interval midpoint like the rest.
#' @param age Left-hand boundaries of the age groups. Defaults to the standard
#'   abridged structure `0, 1, 5, ..., 85`.
#' @param closeout How to handle old age. `"kannisto"` fits a logistic hazard,
#'   `"gompertz"` an exponential one, each extrapolated year by year past the
#'   open-group boundary; `"constant_hazard"` simply assumes the observed rate
#'   persists, giving `e_open = 1/mx`. Kannisto is the better-supported model
#'   above 95 and is the default; `"constant_hazard"` is the only choice fully
#'   consistent with [lt_chiang_ci()].
#' @param ax Optional separation factors supplied directly, as a vector of
#'   length `n_age` or an `n_age x n_draw` matrix. When given, every rule for
#'   deriving `ax` is bypassed: `closeout`, `n_fit_ages` and `keyfitz_iter` are
#'   ignored, and the last value is taken as life expectancy in the open
#'   group.
#'
#'   Use this to reproduce a published life table exactly. With the published
#'   `mx` and `ax` both supplied, all that is left is arithmetic. If you compare
#'   `ex` against a published table without fixing `ax`, any gap you see is a
#'   difference of assumptions rather than an error.
#' @param n_fit_ages Number of oldest age groups used to fit the old-age model,
#'   counting the open group. Ignored for `"constant_hazard"`.
#' @param keyfitz_iter Rounds of Keyfitz refinement applied to the middle age groups.
#'   `0` leaves them at the interval midpoint; a few rounds are enough to
#'   converge. The young-age and old-age separation factors are never
#'   overwritten by the iteration.
#' @param mx_floor Optional lower bounds on `mx`: a single value applied to the
#'   open group only, or a vector of length `n_age` with `NA` where no bound
#'   applies. Floors stabilise the closeout when projected old-age rates fall
#'   very low, but they also cap life expectancy from above. The number of cells
#'   actually bound is recorded in the `floored` attribute; a large count means
#'   the floor is driving the result.
#' @param radix Life table radix.
#' @param max_age Age at which the old-age extrapolation stops.
#' @param b_min Lower bound on the fitted old-age slope.
#' @return An object of class `letools_lt`: a list of `n_age x n_draw` matrices
#'   `mx`, `ax`, `qx`, `px`, `lx`, `dx`, `Lx`, `Tx`, `ex`, carrying `age`, `nx`,
#'   `sex`, `closeout`, `old_age_fit` and `floored` attributes.
#'
#' @examples
#' mx <- c(0.006, 0.0004, 0.0002, 0.0002, 0.0004, 0.0006, 0.0007, 0.0009,
#'         0.0012, 0.0018, 0.0028, 0.0044, 0.0068, 0.0105, 0.0165, 0.0270,
#'         0.0450, 0.0760, 0.1500)
#' x <- lt(mx, sex = 1)
#' x
#' x$ex[1, ]
#' @export
lt <- function(mx,
               sex,
               age = c(0, 1, seq(5, 85, 5)),
               closeout = c("kannisto", "gompertz", "constant_hazard"),
               ax = NULL,
               n_fit_ages = 4L,
               keyfitz_iter = 0L,
               mx_floor = NULL,
               radix = 1e5,
               max_age = 130,
               b_min = 0.02) {

  closeout <- match.arg(closeout)
  sex <- normalise_sex(sex)
  nx <- validate_age(age)

  mx <- if (is.matrix(mx)) mx else matrix(mx, ncol = 1L)
  n <- nrow(mx)
  m <- ncol(mx)

  if (n != length(age)) {
    stop(
      sprintf("`mx` has %d age groups but `age` has %d.", n, length(age)),
      call. = FALSE
    )
  }
  if (any(!is.finite(mx)) || any(mx < 0)) {
    stop("`mx` must be finite and non-negative.", call. = FALSE)
  }
  if (any(mx[n, ] <= 0)) {
    stop(
      "Mortality in the open-ended age group must be strictly positive; the ",
      "group has no upper bound, so a zero rate implies infinite life.",
      call. = FALSE
    )
  }
  if (!is.numeric(keyfitz_iter) || length(keyfitz_iter) != 1L || keyfitz_iter < 0) {
    stop("`keyfitz_iter` must be a single non-negative number.", call. = FALSE)
  }
  keyfitz_iter <- as.integer(keyfitz_iter)

  ax_supplied <- !is.null(ax)
  if (ax_supplied) {
    if (!is.matrix(ax)) {
      if (length(ax) != n) {
        stop(
          sprintf("`ax` vector must have length %d to match `age`, got %d.",
                  n, length(ax)),
          call. = FALSE
        )
      }
      ax <- matrix(ax, nrow = n, ncol = m)
    }

    if (!identical(dim(ax), c(n, m))) {
      stop(
        sprintf("`ax` is %d x %d but `mx` is %d x %d.",
                nrow(ax), ncol(ax), n, m),
        call. = FALSE
      )
    }
    if (any(!is.finite(ax)) || any(ax < 0)) {
      stop("`ax` must be finite and non-negative.", call. = FALSE)
    }
    too_long <- which(apply(ax[-n, , drop = FALSE] > nx[-n], 1L, any))
    if (length(too_long)) {
      stop(
        "`ax` exceeds the interval width at age ",
        paste(age[too_long], collapse = ", "),
        "; those dying in an interval cannot live longer than it lasts.",
        call. = FALSE
      )
    }
  }

  fits_old <- !ax_supplied && closeout != "constant_hazard"
  if (fits_old) {
    if (n_fit_ages < 4L) {
      stop(
        "`n_fit_ages` must be at least 4 to fit an old-age model.",
        call. = FALSE
      )
    }
    if (n - n_fit_ages < 2L) {
      stop(
        sprintf(
          paste0("`n_fit_ages` = %d leaves too few younger age groups; the ",
                 "fit would reach into childhood mortality, which follows ",
                 "neither model."),
          n_fit_ages
        ),
        call. = FALSE
      )
    }
  }

  # --- optional floors on mx ---
  floored <- 0L
  if (!is.null(mx_floor)) {
    if (length(mx_floor) == 1L) {
      floor_vec <- c(rep(NA_real_, n - 1L), mx_floor)
    } else if (length(mx_floor) == n) {
      floor_vec <- as.numeric(mx_floor)
    } else {
      stop(
        sprintf("`mx_floor` must have length 1 or %d, got %d.",
                n, length(mx_floor)),
        call. = FALSE
      )
    }

    for (i in which(!is.na(floor_vec))) {
      bound <- mx[i, ] < floor_vec[i]
      floored <- floored + sum(bound)
      mx[i, bound] <- floor_vec[i]
    }
  }

  nx_mat <- matrix(nx, nrow = n, ncol = m)

  # --- separation factors, first pass ---
  if (!ax_supplied) {
    ax <- nx_mat / 2
    cd <- ax_coale_demeny(mx[1, ], sex)
    ax[1, ] <- cd["a0", ]
    # Coale-Demeny's a1 is 4a1, for ages 1-4 only; with any other second
    # interval (single years, say) it would exceed the interval width
    if (all(age[2:3] == c(1, 5))) ax[2, ] <- cd["a1", ]
    ax[n, ] <- 1 / mx[n, ]
  }

  built <- build_from_ax(mx, ax, nx_mat, radix)

  # --- old-age closeout ---
  old_rows <- if (fits_old) seq.int(n - n_fit_ages + 1L, n) else n
  old_fit  <- NULL

  if (fits_old) {
    old_fit <- ax_old_age(
      built$lx[old_rows, , drop = FALSE],
      age           = age[old_rows],
      method        = closeout,
      max_age       = max_age,
      b_min         = b_min
    )
    ax[old_rows, ] <- old_fit$ax
    built <- build_from_ax(mx, ax, nx_mat, radix)
  }

  # --- optional Keyfitz refinement of the middle ---
  protected <- unique(c(1L, 2L, old_rows))
  middle <- setdiff(seq_len(n - 1L), protected)

  if (!ax_supplied && keyfitz_iter > 0L && length(middle)) {
    for (i in seq_len(keyfitz_iter)) {
      ax <- ax_keyfitz(built$dx, nx, middle, ax)
      built <- build_from_ax(mx, ax, nx_mat, radix)
    }
  }

  # --- Lx, Tx, ex ---
  Lx <- rbind(
    nx_mat[-n, , drop = FALSE] * built$lx[-1, , drop = FALSE] +
      ax[-n, , drop = FALSE] * built$dx[-n, , drop = FALSE],
    built$lx[n, ] * ax[n, ]
  )

  Tx <- apply_rev_cumsum(Lx)
  ex <- Tx / built$lx

  structure(
    list(mx = mx, ax = ax, qx = built$qx, px = built$px,
         lx = built$lx, dx = built$dx, Lx = Lx, Tx = Tx, ex = ex),
    class = "letools_lt",
    age = age,
    nx = nx,
    sex = sex,
    closeout = if (ax_supplied) "supplied" else closeout,
    keyfitz_iter = if (ax_supplied) 0L else keyfitz_iter,
    old_age_fit = old_fit,
    floored = floored,
    n_cells = n * m
  )
}


#' Rebuild qx, px, lx and dx from a set of separation factors
#' @noRd
build_from_ax <- function(mx, ax, nx_mat, radix) {
  n <- nrow(mx)
  m <- ncol(mx)

  qx <- (nx_mat * mx) / (1 + (nx_mat - ax) * mx)
  qx[-n, ] <- pmin(pmax(qx[-n, , drop = FALSE], 0), 1 - 1e-10)
  qx[n, ]  <- 1

  px <- 1 - qx
  lx <- rbind(rep(radix, m), radix * apply_cumprod(px[-n, , drop = FALSE]))
  dx <- lx * qx

  list(qx = qx, px = px, lx = lx, dx = dx)
}


#' @export
print.letools_lt <- function(x, ...) {
  age <- attr(x, "age")
  e0 <- x$ex[1, ]

  cat(sprintf(
    "<letools_lt>  %d age groups (%g-%g+), %d draw%s, sex = %s\n",
    length(age), age[1], age[length(age)], ncol(x$ex),
    if (ncol(x$ex) == 1L) "" else "s",
    if (attr(x, "sex") == 1L) "male" else "female"
  ))
  cat(sprintf("  closeout: %s", attr(x, "closeout")))
  if (attr(x, "keyfitz_iter") > 0L) {
    cat(sprintf(", keyfitz_iter: %d", attr(x, "keyfitz_iter")))
  }
  cat("\n")

  if (ncol(x$ex) == 1L) {
    cat(sprintf("  e0: %.2f\n", e0))
  } else {
    q <- stats::quantile(e0, c(0.025, 0.5, 0.975))
    cat(sprintf("  e0: %.2f (95%% CrI %.2f-%.2f)\n", q[2], q[1], q[3]))
  }

  fl <- attr(x, "floored")
  if (fl > 0L) {
    cat(sprintf(
      "  note: mx floored in %d of %d cells (%.1f%%)\n",
      fl, attr(x, "n_cells"), 100 * fl / attr(x, "n_cells")
    ))
  }
  invisible(x)
}


#' Convert a life table to a long data frame
#'
#' @param x A `letools_lt` object.
#' @param row.names,optional Ignored; present for method consistency.
#' @param ... Ignored.
#' @return A data frame with one row per age group per draw.
#' @export
as.data.frame.letools_lt <- function(x, row.names = NULL, optional = FALSE, ...) {
  age <- attr(x, "age")
  n <- length(age)
  m <- ncol(x$ex)

  out <- data.frame(
    age  = rep(age, times = m),
    draw = rep(seq_len(m), each = n),
    stringsAsFactors = FALSE
  )
  for (nm in c("mx", "ax", "qx", "lx", "dx", "Lx", "Tx", "ex")) {
    out[[nm]] <- as.vector(x[[nm]])
  }
  out
}
