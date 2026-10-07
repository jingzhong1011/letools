# ---------------------------------------------------------------------------
# Separation factors tested on their own, away from lt(). Two of these are
# exact: a Gompertz survivorship curve must return its own slope, and a
# constant hazard must drive the fitted slope down to the floor.
# ---------------------------------------------------------------------------

old_ages <- c(70, 75, 80, 85)

# survivorship generated from a known Gompertz hazard mu(x) = A exp(Bx)
gompertz_lx <- function(age, A = 3e-5, B = 0.093, radix = 1e5) {
  matrix(radix * exp(-(A / B) * (exp(B * age) - 1)), ncol = 1)
}


# --- Coale-Demeny ----------------------------------------------------------

test_that("Coale-Demeny switches to fixed values above the threshold", {
  hi <- ax_coale_demeny(0.15, sex = 1)
  expect_equal(unname(hi["a0", ]), 0.330)
  expect_equal(unname(hi["a1", ]), 1.352)

  lo <- ax_coale_demeny(0.05, sex = 1)
  expect_equal(unname(lo["a0", ]), 0.045 + 2.684 * 0.05)
})


test_that("Coale-Demeny is sex-specific and vectorised over draws", {
  m <- ax_coale_demeny(c(0.02, 0.05, 0.20), sex = "male")
  f <- ax_coale_demeny(c(0.02, 0.05, 0.20), sex = "female")

  expect_equal(dim(m), c(2L, 3L))
  expect_false(isTRUE(all.equal(m, f)))
})


test_that("lower infant mortality concentrates deaths earlier", {
  low  <- ax_coale_demeny(0.01, sex = 1)["a0", ]
  high <- ax_coale_demeny(0.09, sex = 1)["a0", ]
  expect_lt(low, high)
})


test_that("the sex argument accepts the usual spellings and rejects the rest", {
  expect_equal(normalise_sex("Female"), 2L)
  expect_equal(normalise_sex(1), 1L)
  expect_error(normalise_sex(0), "must be 1")
  expect_error(normalise_sex("other"), "must be 1, 2")
})


# --- old-age closeout ------------------------------------------------------

test_that("the Gompertz slope is recovered from synthetic survivorship", {
  B_true <- 0.093
  lx  <- gompertz_lx(old_ages, B = B_true)
  fit <- ax_old_age(lx, age = old_ages, method = "gompertz")
  expect_equal(fit$b, B_true, tolerance = 1e-6)
})


test_that("the slope floor engages when old-age mortality is flat", {
  lx <- matrix(1e5 * exp(-0.05 * (old_ages - 70)), ncol = 1)
  for (meth in c("kannisto", "gompertz")) {
    fit <- ax_old_age(lx, age = old_ages, method = meth, b_min = 0.02)
    expect_equal(fit$b, 0.02, info = meth)
  }
})


test_that("separation factors stay inside their intervals", {
  lx <- gompertz_lx(old_ages)
  for (meth in c("kannisto", "gompertz")) {
    fit <- ax_old_age(lx, age = old_ages, method = meth)
    closed <- fit$ax[1:3, ]
    expect_true(all(closed > 0 & closed < 5), info = meth)
    expect_gt(fit$ax[4, ], 0)
  }
})


test_that("steeper old-age mortality front-loads deaths within each band", {
  shallow <- ax_old_age(gompertz_lx(old_ages, B = 0.06), age = old_ages)
  steep   <- ax_old_age(gompertz_lx(old_ages, B = 0.14), age = old_ages)

  expect_true(all(steep$ax[1:3, ] < shallow$ax[1:3, ]))
  expect_lt(steep$ax[4, ], shallow$ax[4, ])
})


test_that("Kannisto gives a longer open-group life than Gompertz", {
  lx <- gompertz_lx(old_ages)
  k <- ax_old_age(lx, age = old_ages, method = "kannisto")
  g <- ax_old_age(lx, age = old_ages, method = "gompertz")
  expect_gt(k$ax[4, ], g$ax[4, ])
})


test_that("the closeout is vectorised over draws without leakage", {
  lx <- cbind(gompertz_lx(old_ages, B = 0.08),
              gompertz_lx(old_ages, B = 0.11))

  both <- ax_old_age(lx, age = old_ages, method = "kannisto")
  one  <- ax_old_age(lx[, 1, drop = FALSE], age = old_ages, method = "kannisto")
  two  <- ax_old_age(lx[, 2, drop = FALSE], age = old_ages, method = "kannisto")

  expect_equal(dim(both$ax), c(4L, 2L))
  expect_equal(both$ax[, 1, drop = FALSE], one$ax, tolerance = 1e-10)
  expect_equal(both$ax[, 2, drop = FALSE], two$ax, tolerance = 1e-10)
})


test_that("malformed old-age input is rejected", {
  lx <- gompertz_lx(old_ages)

  expect_error(ax_old_age(lx[1:3, , drop = FALSE], age = old_ages[1:3]),
               "at least four")
  expect_error(ax_old_age(lx, age = c(70, 75, 80, 90)), "evenly spaced")
  expect_error(ax_old_age(lx, age = old_ages, max_age = 80), "must exceed")
  expect_error(ax_old_age(lx, age = old_ages[1:3]), "one value per row")

  flat <- matrix(1e5, nrow = 4, ncol = 1)
  expect_error(ax_old_age(flat, age = old_ages), "Non-positive force")
})


# --- Keyfitz ---------------------------------------------------------------

test_that("a flat death distribution leaves separation factors at the midpoint", {
  dx <- matrix(100, nrow = 5, ncol = 1)
  ax <- matrix(2.5, nrow = 5, ncol = 1)
  out <- ax_keyfitz(dx, nx = rep(5, 5), rows = 2:4, ax = ax)
  expect_equal(out, ax)
})


test_that("deaths rising with age push the separation factor later", {
  dx <- matrix(c(50, 100, 200, 400, 800), ncol = 1)
  ax <- matrix(2.5, nrow = 5, ncol = 1)
  out <- ax_keyfitz(dx, nx = rep(5, 5), rows = 2:4, ax = ax)
  expect_true(all(out[2:4, ] > 2.5))
})


test_that("Keyfitz clamps to the interval and survives empty age groups", {
  dx <- matrix(c(1, 1e-12, 1e6, 1, 1), ncol = 1)
  ax <- matrix(2.5, nrow = 5, ncol = 1)
  out <- ax_keyfitz(dx, nx = rep(5, 5), rows = 2:4, ax = ax)

  expect_true(all(out >= 0 & out <= 5))
  expect_true(all(is.finite(out)))

  zero <- matrix(c(1, 0, 1, 1, 1), ncol = 1)
  expect_equal(ax_keyfitz(zero, rep(5, 5), 2L, ax)[2, ], 2.5)
})


test_that("rows outside the requested set are untouched", {
  dx <- matrix(c(50, 100, 200, 400, 800), ncol = 1)
  ax <- matrix(c(0.1, 0.4, 2.5, 2.5, 2.5), ncol = 1)
  out <- ax_keyfitz(dx, nx = rep(5, 5), rows = 3L, ax = ax)

  expect_equal(out[c(1, 2, 4, 5), ], ax[c(1, 2, 4, 5), ])
  expect_false(isTRUE(all.equal(out[3, ], ax[3, ])))
})
