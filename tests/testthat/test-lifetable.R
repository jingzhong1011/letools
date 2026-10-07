# ---------------------------------------------------------------------------
# Most of these check internal identities: relationships that must hold for any
# life table whatever the input, like lx declining by exactly dx or Tx being the
# reverse cumulative sum of Lx. They catch structural mistakes -- an off-by-one,
# a mishandled open group, cross-contamination between draws -- without needing
# reference data.
#
# They cannot catch a wrong assumption applied consistently. A life table built
# with the wrong separation factors is still internally consistent; it is just
# wrong. Only the golden test at the bottom, against a published life table,
# settles that.
# ---------------------------------------------------------------------------

ages19 <- c(0, 1, seq(5, 85, 5))

fake_mx <- function(level = 1, n_draw = 1L, seed = NULL, age = ages19) {
  if (!is.null(seed)) set.seed(seed)
  base <- 0.00005 * exp(0.085 * age)
  base[1] <- 0.006
  base[2] <- 0.0004

  mx <- matrix(base * level, nrow = length(base), ncol = n_draw)
  if (n_draw > 1L) {
    mx <- mx * matrix(stats::runif(length(mx), 0.9, 1.1), nrow = length(base))
  }
  mx
}

all_closeouts <- c("kannisto", "gompertz", "constant_hazard")


# --- internal identities ---------------------------------------------------

test_that("survivorship declines by exactly the death counts", {
  for (co in all_closeouts) {
    x <- lt(fake_mx(), sex = 1, closeout = co)
    n <- nrow(x$lx)
    expect_equal(x$lx[-1, ], (x$lx - x$dx)[-n, ], tolerance = 1e-10,
                 info = co)
  }
})


test_that("deaths sum to the radix", {
  for (co in all_closeouts) {
    x <- lt(fake_mx(), sex = 2, closeout = co, radix = 1e5)
    expect_equal(sum(x$dx), 1e5, tolerance = 1e-8, info = co)
  }
})


test_that("Tx is the reverse cumulative sum of Lx, and ex is Tx over lx", {
  x <- lt(fake_mx(), sex = 1)
  n <- nrow(x$Lx)

  expect_equal(x$Tx[n, ], x$Lx[n, ])
  expect_equal(
    x$Tx[-n, , drop = FALSE],
    x$Tx[-1, , drop = FALSE] + x$Lx[-n, , drop = FALSE]
  )
  expect_equal(x$ex, x$Tx / x$lx)
})


test_that("the open group is self-consistent under every closeout", {
  for (co in all_closeouts) {
    x <- lt(fake_mx(), sex = 1, closeout = co)
    n <- nrow(x$ex)
    expect_equal(x$qx[n, ], 1, info = co)
    expect_equal(x$Lx[n, ], x$lx[n, ] * x$ex[n, ], info = co)
    expect_equal(x$ax[n, ], x$ex[n, ], info = co)
  }
})


test_that("qx stays in the unit interval and lx strictly declines", {
  x <- lt(fake_mx(level = 5), sex = 1)   # deliberately harsh mortality
  expect_true(all(x$qx >= 0 & x$qx <= 1))
  expect_true(all(diff(x$lx[, 1]) < 0))
})


# --- invariances -----------------------------------------------------------

test_that("the radix cancels out of ex, qx and ax", {
  a <- lt(fake_mx(), sex = 1, radix = 1e5)
  b <- lt(fake_mx(), sex = 1, radix = 1)

  expect_equal(a$ex, b$ex, tolerance = 1e-9)
  expect_equal(a$qx, b$qx, tolerance = 1e-9)
  expect_equal(a$ax, b$ax, tolerance = 1e-9)
  expect_equal(a$lx / 1e5, b$lx, tolerance = 1e-9)
})


test_that("draws do not leak into each other", {
  mx <- fake_mx(n_draw = 8, seed = 11)

  for (co in all_closeouts) {
    together <- lt(mx, sex = 1, closeout = co)
    separate <- vapply(
      seq_len(ncol(mx)),
      function(j) lt(mx[, j, drop = FALSE], sex = 1, closeout = co)$ex[1, ],
      numeric(1)
    )
    expect_equal(together$ex[1, ], separate, tolerance = 1e-9, info = co)
  }
})


test_that("uniformly lower mortality raises life expectancy", {
  for (co in all_closeouts) {
    hi <- lt(fake_mx(level = 1.0), sex = 1, closeout = co)$ex[1, ]
    lo <- lt(fake_mx(level = 0.7), sex = 1, closeout = co)$ex[1, ]
    expect_gt(lo, hi)
  }
})


# --- closeout --------------------------------------------------------------

test_that("constant-hazard closeout gives exactly 1/mx in the open group", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  n  <- nrow(mx)
  expect_equal(x$ex[n, ], 1 / mx[n, ], tolerance = 1e-10)
})


test_that("Kannisto gives a longer open-group life than Gompertz", {
  # the logistic hazard decelerates at the highest ages where the exponential
  # one keeps climbing, so survival past the boundary is longer under Kannisto
  mx <- fake_mx()
  n  <- nrow(mx)
  k <- lt(mx, sex = 1, closeout = "kannisto")$ex[n, ]
  g <- lt(mx, sex = 1, closeout = "gompertz")$ex[n, ]
  expect_gt(k, g)
})


test_that("the closeout is recorded and reaches the fitted coefficients", {
  x <- lt(fake_mx(), sex = 1, closeout = "kannisto")
  fit <- attr(x, "old_age_fit")

  expect_equal(attr(x, "closeout"), "kannisto")
  expect_equal(fit$method, "kannisto")
  expect_true(is.finite(fit$b))
  expect_gt(fit$b, 0)

  expect_null(attr(lt(fake_mx(), sex = 1, closeout = "constant_hazard"),
                   "old_age_fit"))
})


test_that("extrapolating further changes the open group only slightly", {
  mx <- fake_mx()
  n  <- nrow(mx)
  short <- lt(mx, sex = 1, max_age = 110)$ex[n, ]
  long  <- lt(mx, sex = 1, max_age = 150)$ex[n, ]
  expect_equal(short, long, tolerance = 0.05)
})


# --- Keyfitz iteration -----------------------------------------------------

test_that("iteration leaves young and old separation factors untouched", {
  mx <- fake_mx()
  plain <- lt(mx, sex = 1, keyfitz_iter = 0L)
  iter  <- lt(mx, sex = 1, keyfitz_iter = 4L)

  protected <- c(1L, 2L, 16L, 17L, 18L, 19L)
  expect_equal(iter$ax[protected, ], plain$ax[protected, ], tolerance = 1e-9)
})


test_that("iteration moves the middle separation factors off the midpoint", {
  mx <- fake_mx()
  iter <- lt(mx, sex = 1, keyfitz_iter = 4L)

  middle <- 3:15
  expect_false(isTRUE(all.equal(iter$ax[middle, ], rep(2.5, length(middle)),
                                check.attributes = FALSE)))
  expect_true(all(iter$ax[middle, ] >= 0 & iter$ax[middle, ] <= 5))
})


test_that("iteration converges rather than drifting", {
  mx <- fake_mx()
  a <- lt(mx, sex = 1, keyfitz_iter = 4L)$ex[1, ]
  b <- lt(mx, sex = 1, keyfitz_iter = 8L)$ex[1, ]
  expect_equal(a, b, tolerance = 1e-6)
})


test_that("iteration is recorded and rejects nonsense", {
  expect_equal(attr(lt(fake_mx(), sex = 1, keyfitz_iter = 3L), "keyfitz_iter"), 3L)
  expect_error(lt(fake_mx(), sex = 1, keyfitz_iter = -1), "non-negative")
})


# --- mortality floors ------------------------------------------------------

test_that("a floor that never binds changes nothing and is recorded as zero", {
  mx <- fake_mx()
  free    <- lt(mx, sex = 1)
  floored <- lt(mx, sex = 1, mx_floor = 1e-8)

  expect_equal(free$ex, floored$ex, tolerance = 1e-12)
  expect_equal(attr(floored, "floored"), 0L)
})


test_that("a floor on the open group bites only under constant_hazard", {
  # Under kannisto or gompertz the open group's own rate never enters the
  # calculation: qx is 1 there by construction, and the closeout is driven by
  # the survivorship pattern across the fitting bands. Flooring mx in the open
  # group alone therefore changes the count but not the answer. This is worth
  # knowing before setting a floor and assuming it did something.
  mx <- fake_mx(level = 0.3, n_draw = 4, seed = 7)

  free_k    <- lt(mx, sex = 1, closeout = "kannisto")
  floored_k <- lt(mx, sex = 1, closeout = "kannisto", mx_floor = 0.09)

  expect_equal(attr(floored_k, "floored"), 4L)
  expect_equal(floored_k$ex, free_k$ex, tolerance = 1e-12)

  free_c    <- lt(mx, sex = 1, closeout = "constant_hazard")
  floored_c <- lt(mx, sex = 1, closeout = "constant_hazard", mx_floor = 0.09)

  expect_equal(attr(floored_c, "floored"), 4L)
  expect_true(all(floored_c$ex[1, ] < free_c$ex[1, ]))
})


test_that("a floor on the fitting bands caps life expectancy under kannisto", {
  mx <- fake_mx(level = 0.3, n_draw = 4, seed = 7)

  # bound the 80-84 band, which does feed the old-age fit
  bands <- c(rep(NA_real_, 17), 0.025, NA_real_)

  free    <- lt(mx, sex = 1, closeout = "kannisto")
  floored <- lt(mx, sex = 1, closeout = "kannisto", mx_floor = bands)

  expect_equal(attr(floored, "floored"), 4L)
  expect_true(all(floored$ex[1, ] < free$ex[1, ]))
})


test_that("a per-age floor vector must match the age structure", {
  expect_error(lt(fake_mx(), sex = 1, mx_floor = c(0.01, 0.02)),
               "length 1 or 19")
})


# --- input validation ------------------------------------------------------

test_that("malformed age structures are rejected with a useful message", {
  expect_error(validate_age(c(0, 5, 1)), "strictly increasing")
  expect_error(validate_age(c(1, 5, 10)), "must start at 0")
  expect_error(validate_age(c(0, 5)), "at least three")
})


test_that("bad mortality input is rejected before any arithmetic", {
  bad <- fake_mx(); bad[5, 1] <- -0.01
  expect_error(lt(bad, sex = 1), "finite and non-negative")

  zero_open <- fake_mx(); zero_open[19, 1] <- 0
  expect_error(lt(zero_open, sex = 1), "strictly positive")

  expect_error(lt(fake_mx(), sex = 1, age = c(0, 1, 5)), "19 age groups but")
})


test_that("a fit window that reaches into childhood is refused", {
  expect_error(lt(fake_mx(), sex = 1, n_fit_ages = 18L), "too few younger")
  expect_error(lt(fake_mx(), sex = 1, n_fit_ages = 3L), "at least 4")
})


test_that("constant_hazard needs no fit window", {
  expect_silent(
    lt(fake_mx(), sex = 1, closeout = "constant_hazard", n_fit_ages = 2L)
  )
})


# --- methods ---------------------------------------------------------------

test_that("the long data frame has one row per age group per draw", {
  df <- as.data.frame(lt(fake_mx(n_draw = 3, seed = 2), sex = 1))
  expect_equal(nrow(df), 19L * 3L)
  expect_true(all(c("age", "draw", "mx", "ax", "qx", "lx", "ex") %in% names(df)))
})


test_that("printing reports the closeout and flags floored cells", {
  expect_output(print(lt(fake_mx(), sex = 1)), "kannisto")
  expect_output(
    print(lt(fake_mx(level = 0.3), sex = 1, mx_floor = 0.09)),
    "floored"
  )
})


# --- supplied separation factors -------------------------------------------

test_that("supplied ax bypasses every rule for deriving it", {
  mx <- fake_mx()
  own <- lt(mx, sex = 1)

  same <- lt(mx, sex = 1, ax = own$ax)
  expect_equal(same$ex, own$ex, tolerance = 1e-12)
  expect_equal(same$ax, own$ax)
  expect_equal(attr(same, "closeout"), "supplied")

  # closeout and keyfitz_iter must have no effect once ax is given
  a <- lt(mx, sex = 1, ax = own$ax, closeout = "gompertz", keyfitz_iter = 4L)
  expect_equal(a$ex, same$ex, tolerance = 1e-12)
})


test_that("supplied ax still produces a self-consistent table", {
  mx <- fake_mx()
  ax <- c(0.1, 1.5, rep(2.5, 16), 6.5)
  x  <- lt(mx, sex = 1, ax = ax)
  n  <- nrow(x$lx)

  expect_equal(sum(x$dx), 1e5, tolerance = 1e-8)
  expect_equal(x$ex[n, ], 6.5)
  expect_equal(x$Lx[n, ], x$lx[n, ] * 6.5)
  expect_equal(x$ex, x$Tx / x$lx)
})


test_that("Coale-Demeny 4a1 is applied only to a 1-4 age group", {
  # standard abridged ages: second group is 1-4, so it gets 4a1
  ab <- lt(fake_mx(), sex = 1)
  expect_equal(ab$ax[2, ], ax_coale_demeny(fake_mx()[1, 1], 1)["a1", ],
               ignore_attr = TRUE)

  # single years: 4a1 would exceed the one-year interval, so the group starts
  # at the midpoint like any other
  ref <- taiwan_lt("male")
  x <- lt(ref$mx, sex = 1, age = ref$age)
  nx <- attr(x, "nx")
  n <- nrow(x$ax)

  expect_equal(x$ax[2, ], 0.5)
  expect_true(all(x$ax[-n, ] >= 0 & x$ax[-n, ] <= nx[-n]))
})


test_that("supplied ax is validated", {
  mx <- fake_mx()
  expect_error(lt(mx, sex = 1, ax = c(1, 2, 3)), "must have length.*to match `age`")
  expect_error(lt(mx, sex = 1, ax = rep(-1, 19)), "finite and non-negative")

  # 1.5 years lived by those dying in a one-year interval is impossible
  too_long <- c(1.5, 1.5, rep(2.5, 16), 6.5)
  expect_error(lt(mx, sex = 1, ax = too_long), "exceeds the interval width at age 0")

  # the open group has no width, so any non-negative value is allowed there
  expect_silent(lt(mx, sex = 1, ax = c(0.1, 1.5, rep(2.5, 16), 40)))
})


# ---------------------------------------------------------------------------
# Golden test: Taiwan 簡易生命表, 民國112年 (2023), Ministry of the Interior.
#
# This does NOT check that letools reproduces the published life expectancy,
# and it should not. e0 is a function of mx AND of a set of assumptions, and
# the Ministry's assumptions are not this package's. Their a0 is 0.157 for
# males where Coale-Demeny West would give 0.058 -- nearly three times apart.
# Failing to match e0 under different assumptions says nothing about the code.
#
# What it does check is the arithmetic, with the assumptions pinned down.
# Supplying the published mx and ax together leaves nothing to disagree about:
# qx follows from mx and ax, lx from qx, Lx from lx and dx, Tx by accumulation,
# ex by division. Every published column must come back.
#
# One honest caveat: the fixture's ax is derived from the published Lx via
# Lx = l(x+1) + ax*dx, so reproducing Lx is partly circular. The columns that
# are not circular are lx (from the qx chain), Tx (accumulation) and ex
# (division), and those are the ones that matter.
#
# The fixture is single-year, ages 0 to 85+, for all three published series.
# data-raw/build_taiwan_lt.py regenerates inst/extdata/taiwan_lt_2023.csv from
# the source workbook, and helper-fixtures.R reads it through taiwan_lt().
# ---------------------------------------------------------------------------

test_that("lt() reproduces the published table when given its separation factors", {
  for (series in c("total", "male", "female")) {
    ref <- read_reference_lt(series)

    expect_equal(nrow(ref), 86L, info = series)

    # `sex` is inert here: supplying ax bypasses the Coale-Demeny rules, which
    # are the only place it would be used. The published "total" series has no
    # sex to pass.
    x <- lt(ref$mx, sex = 1, age = ref$age, ax = ref$ax)

    for (col in c("qx", "lx", "dx", "Lx", "Tx", "ex")) {
      expect_equal(as.vector(x[[col]]), ref[[col]],
                   tolerance = 1e-10, info = paste(series, col))
    }
  }
})


test_that("the published table is not reproducible under this package's own assumptions", {
  # The mirror image of the test above, kept as documentation rather than as a
  # correctness check: with the same mx but letools' own separation factors,
  # e0 lands close to the published figure without matching it. That gap is the
  # difference between two sets of assumptions, and it is expected.
  ref <- read_reference_lt("male")
  own <- lt(ref$mx, sex = 1, age = ref$age, closeout = "kannisto")

  expect_equal(own$ex[1, 1], ref$ex[1], tolerance = 0.05)
  expect_false(isTRUE(all.equal(own$ex[1, 1], ref$ex[1], tolerance = 1e-8)))
})
