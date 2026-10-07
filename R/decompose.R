# ---------------------------------------------------------------------------
# decomp
#
# decomp_*：additive contribution of each age group, or cause, to a difference
# in life expectancy at birth. Everything here takes life tables as produced
# by lt(), so each component is an n_age x n_draw matrix and draws come along
# for free.
#
# Direction convention, used consistently throughout: `from` is the reference
# (baseline, worse-off group, earlier year) and `to` is the comparison
# (scenario, better-off group, later year). Contributions therefore sum to
#     e0(to) - e0(from)
# and a positive contribution means the age group or cause moved life
# expectancy upward from `from` to `to`.
# ---------------------------------------------------------------------------


#' Coerce a life table component to an n_age x n_draw matrix
#' @noRd
as_draw_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  matrix(x, ncol = 1L)
}


#' Check that a matrix argument has the expected n_age x n_draw shape
#'
#' R would otherwise recycle a short matrix silently, or fail later with an
#' error that does not name the argument.
#' @noRd
check_shape <- function(x, expected, arg) {
  d <- dim(as_draw_matrix(x))
  if (!identical(as.integer(d), as.integer(expected))) {
    stop(
      sprintf("%s is %d x %d but should be %d x %d.",
              arg, d[1], d[2], expected[1], expected[2]),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Check that two life tables can be compared
#' @noRd
check_conformable <- function(lt_from, lt_to) {
  needed <- c("lx", "Lx", "ex")
  missing_from <- setdiff(needed, names(lt_from))
  missing_to <- setdiff(needed, names(lt_to))

  if (length(missing_from) || length(missing_to)) {
    stop(
      "Life tables must contain lx, Lx and ex. Missing: ",
      paste(unique(c(missing_from, missing_to)), collapse = ", "),
      call. = FALSE
    )
  }

  d_from <- dim(as_draw_matrix(lt_from$lx))
  d_to <- dim(as_draw_matrix(lt_to$lx))

  if (!identical(d_from, d_to)) {
    stop(
      sprintf(
        "Life tables have different shapes: `from` is %d x %d, `to` is %d x %d.",
        d_from[1], d_from[2], d_to[1], d_to[2]
      ),
      call. = FALSE
    )
  }

  if (d_from[1] < 2L) {
    stop("A life table needs at least two age groups to decompose.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Arriaga decomposition of a life expectancy difference by age
#'
#' The decomposition is exact, so contributions sum to the total difference up
#' to floating-point error. `check = TRUE` verifies this and is cheap, so it is
#' on by default. Turn it off only inside a hot loop.
#'
#' @param lt_from Reference life table, as returned by [lt()].
#' @param lt_to Comparison life table with the same age structure and number of
#'   draws.
#' @param check Verify that contributions sum to the total difference.
#' @param tol Absolute tolerance for that check, in years.
#' @return An `n_age x n_draw` matrix of contributions in years.
#'
#' @references
#' Arriaga EE (1984). Measuring and explaining the change in life expectancies.
#' *Demography* 21(1):83-96.
#'
#' @examples
#' \dontrun{
#' lt_2000 <- lt(mx_2000, sex = 1)
#' lt_2021 <- lt(mx_2021, sex = 1)
#' contrib <- decomp_arriaga(lt_2000, lt_2021)
#' colSums(contrib)  # equals e0(2021) - e0(2000)
#' }
#' @export
decomp_arriaga <- function(lt_from, lt_to, check = TRUE, tol = 1e-8) {
  check_conformable(lt_from, lt_to)

  lx_from <- as_draw_matrix(lt_from$lx)
  Lx_from <- as_draw_matrix(lt_from$Lx)
  ex_from <- as_draw_matrix(lt_from$ex)
  lx_to <- as_draw_matrix(lt_to$lx)
  Lx_to <- as_draw_matrix(lt_to$Lx)
  ex_to <- as_draw_matrix(lt_to$ex)

  n <- nrow(lx_from)
  radix <- lx_from[1, ]

  contrib <- matrix(0, nrow = n, ncol = ncol(lx_from))
  i <- seq_len(n - 1L)

  # person-years gained within each closed age group
  direct <- lx_from[i, , drop = FALSE] *
    (Lx_to[i, , drop = FALSE] / lx_to[i, , drop = FALSE] -
       Lx_from[i, , drop = FALSE] / lx_from[i, , drop = FALSE])

  # person-years gained beyond the group because more survivors reach it
  indirect <- (lx_from[i, , drop = FALSE] *
                 lx_to[i + 1L, , drop = FALSE] / lx_to[i, , drop = FALSE] -
                 lx_from[i + 1L, , drop = FALSE]) *
    ex_to[i + 1L, , drop = FALSE]

  contrib[i, ] <- sweep(direct + indirect, 2L, radix, "/")

  # the open-ended group has no ages beyond it, so no indirect effect
  contrib[n, ] <- lx_from[n, ] *
    (Lx_to[n, ] / lx_to[n, ] - Lx_from[n, ] / lx_from[n, ]) / radix

  if (isTRUE(check)) {
    total <- ex_to[1, ] - ex_from[1, ]
    err <- max(abs(colSums(contrib) - total))
    if (!is.finite(err) || err > tol) {
      stop(
        sprintf(
          paste0("Arriaga decomposition did not close: largest discrepancy ",
                 "%.3e years exceeds tol = %.3e. This usually means lx, Lx ",
                 "and ex came from different life tables, or that the open-",
                 "ended age group was handled inconsistently upstream."),
          err, tol
        ),
        call. = FALSE
      )
    }
  }

  contrib
}


#' Split age-specific contributions across causes of death
#'
#' Partitions each age group's contribution in proportion to how much each
#' cause contributed to the change in the all-cause mortality rate at that age.
#'
#' Where the all-cause rate barely moves, the denominator approaches zero and
#' the split is not identified; those cells are returned as zero rather than as
#' a very large number. `min_diff` sets the threshold, expressed relative to the
#' larger of the two all-cause rates, so it does not need retuning when rates
#' are on a different scale.
#'
#' @param contrib Age contributions from [decomp_arriaga()] or
#'   [decomp_pollard()].
#' @param mx_from,mx_to All-cause mortality rate matrices for the two life
#'   tables, same shape as `contrib`.
#' @param mx_cause_from,mx_cause_to Named lists of cause-specific mortality rate
#'   matrices. Both must use the same cause names.
#' @param min_diff Relative threshold below which an age group's all-cause
#'   change is treated as too small to attribute.
#' @return A named list of `n_age x n_draw` contribution matrices, one per
#'   cause.
#' @export
decomp_by_cause <- function(contrib,
                            mx_from, mx_to,
                            mx_cause_from, mx_cause_to,
                            min_diff = 1e-10) {

  causes <- names(mx_cause_from)
  if (is.null(causes) || !length(causes)) {
    stop("`mx_cause_from` must be a named list of rate matrices.", call. = FALSE)
  }
  if (!setequal(causes, names(mx_cause_to))) {
    stop(
      "Cause names differ between `mx_cause_from` and `mx_cause_to`: ",
      paste(symmetric_diff(causes, names(mx_cause_to)), collapse = ", "),
      call. = FALSE
    )
  }

  mx_from <- as_draw_matrix(mx_from)
  mx_to   <- as_draw_matrix(mx_to)
  contrib <- as_draw_matrix(contrib)

  expected <- dim(contrib)
  check_shape(mx_from, expected, "`mx_from`")
  check_shape(mx_to, expected, "`mx_to`")
  for (cause in causes) {
    check_shape(mx_cause_from[[cause]], expected,
                sprintf('`mx_cause_from[["%s"]]`', cause))
    check_shape(mx_cause_to[[cause]], expected,
                sprintf('`mx_cause_to[["%s"]]`', cause))
  }

  diff_all <- mx_from - mx_to
  scale    <- pmax(abs(mx_from), abs(mx_to))
  usable   <- abs(diff_all) > min_diff * pmax(scale, .Machine$double.eps)

  out <- lapply(causes, function(cause) {
    diff_c <- as_draw_matrix(mx_cause_from[[cause]]) -
      as_draw_matrix(mx_cause_to[[cause]])

    share <- matrix(0, nrow = nrow(diff_all), ncol = ncol(diff_all))
    share[usable] <- diff_c[usable] / diff_all[usable]

    contrib * share
  })

  names(out) <- causes
  out
}


#' Symmetric difference of two character vectors
#' @noRd
symmetric_diff <- function(a, b) {
  unique(c(setdiff(a, b), setdiff(b, a)))
}


#' Pollard decomposition of a life expectancy difference by age
#'
#' An alternative to Arriaga that expresses the difference as a weighted
#' integral of the difference in age-specific mortality rates. Because the
#' contribution is linear in the rate difference, splitting it across causes is
#' exact rather than proportional. Prefer Pollard when the cause split is the
#' quantity you actually care about.
#'
#' The cost is that Pollard is exact only in continuous time. The discrete
#' residual is first order in the age interval width: on the published Taiwan
#' 2023 male-female comparison it overshoots by 9 per cent on five-year bands
#' and 1.8 per cent on single years, where Arriaga closes to floating-point
#' precision at either width. If the totals matter as much as the cause split,
#' work in single years or take the totals from Arriaga.
#'
#' @param lt_from,lt_to Life tables as returned by [lt()].
#' @param mx_from,mx_to All-cause mortality rate matrices for those life tables.
#' @param check Verify that contributions sum to the total difference.
#' @param rel_tol Relative tolerance for that check, as a fraction of the total
#'   difference. Discretisation alone runs to 9 per cent on five-year bands and
#'   1.8 per cent on single years, while rates that do not belong to these life
#'   tables start at about 33 per cent, so the default separates the two.
#' @return An `n_age x n_draw` matrix of contributions in years.
#'
#' @references
#' Pollard JH (1988). On the decomposition of changes in expectation of life
#' and differentials in life expectancy. *Demography* 25(2):265-276.
#' @export
decomp_pollard <- function(lt_from, lt_to,
                           mx_from, mx_to,
                           check = TRUE, rel_tol = 0.2) {
  check_conformable(lt_from, lt_to)

  Lx_from <- as_draw_matrix(lt_from$Lx)
  ex_from <- as_draw_matrix(lt_from$ex)
  lx_from <- as_draw_matrix(lt_from$lx)
  Lx_to <- as_draw_matrix(lt_to$Lx)
  ex_to <- as_draw_matrix(lt_to$ex)
  lx_to <- as_draw_matrix(lt_to$lx)

  check_shape(mx_from, dim(Lx_from), "`mx_from`")
  check_shape(mx_to, dim(Lx_from), "`mx_to`")

  # person-years lived are expressed per survivor of each table's own radix so
  # that the weights are on the same scale as the rate differences, whatever
  # radix either table was built with
  Lx_from_scaled <- sweep(Lx_from, 2L, lx_from[1, ], "/")
  Lx_to_scaled <- sweep(Lx_to, 2L, lx_to[1, ], "/")

  weights <- (Lx_from_scaled * ex_to + Lx_to_scaled * ex_from) / 2

  contrib <- (as_draw_matrix(mx_from) - as_draw_matrix(mx_to)) * weights

  if (isTRUE(check)) {
    total <- ex_to[1, ] - ex_from[1, ]
    scale <- pmax(abs(total), .Machine$double.eps)
    err <- max(abs(colSums(contrib) - total) / scale)
    if (!is.finite(err) || err > rel_tol) {
      warning(
        sprintf(
          paste0("Pollard decomposition closed to %.1f%% of the total ",
                 "difference, above rel_tol = %.1f%%. Discretisation alone ",
                 "accounts for a few per cent per year of age interval ",
                 "width; more than that suggests the rate matrices do not ",
                 "match the life tables."),
          100 * err, 100 * rel_tol
        ),
        call. = FALSE
      )
    }
  }

  contrib
}
