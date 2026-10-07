# ---------------------------------------------------------------------------
# Chiang's variance has two properties that pin it down exactly, so these are
# not approximate comparisons: the standard error must scale as the inverse
# square root of the death count, and at the open-ended age group it must
# reduce to e_open / sqrt(D). Either one breaking means the weighting or the
# accumulation is wrong.
# ---------------------------------------------------------------------------

ages19 <- c(0, 1, seq(5, 85, 5))

fake_mx <- function(level = 1, n_draw = 1L, age = ages19) {
  base <- 0.00005 * exp(0.085 * age)
  base[1] <- 0.006
  base[2] <- 0.0004
  matrix(base * level, nrow = length(base), ncol = n_draw)
}

# Deaths consistent with the rates, for a population of the given size. Not
# rounded: Chiang's variance does not require integers, and rounding would break
# the exact scaling relationship the first test relies on.
fake_deaths <- function(mx, pop = 5e4) mx * pop


test_that("standard errors scale as one over the square root of deaths", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")

  small <- lt_chiang_ci(x, fake_deaths(mx, pop = 5e4))
  large <- lt_chiang_ci(x, fake_deaths(mx, pop = 5e4 * 4))

  # four times the deaths, half the standard error, at every age
  expect_equal(large$se, small$se / 2, tolerance = 1e-6)
})


test_that("the open group reduces to e_open over the square root of deaths", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  D  <- fake_deaths(mx)

  ci <- lt_chiang_ci(x, D)
  n  <- nrow(mx)

  expect_equal(ci$se[n, ], x$ex[n, ] / sqrt(D[n, ]), tolerance = 1e-8)
})


test_that("closed-group variance matches a numerical delta method", {
  # The two tests above hold whatever weight each age group gets, so they
  # cannot catch a wrong leverage term. This one can: differentiate e0 with
  # respect to each closed-group qx numerically, holding ax fixed in years,
  # and the sum of squared gradients times Var(qx) must equal Chiang's.
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  D  <- fake_deaths(mx)
  n  <- nrow(mx)
  nx <- attr(x, "nx")
  ax <- x$ax[, 1]

  e0_of_q <- function(q) {
    l <- cumprod(c(1, 1 - q))
    d <- l[-n] * q
    sum(nx[-n] * l[-1] + ax[-n] * d) + l[n] / mx[n, 1]
  }

  q <- x$qx[-n, 1]
  h <- 1e-7
  grad <- vapply(seq_along(q), function(i) {
    up <- q; up[i] <- up[i] + h
    dn <- q; dn[i] <- dn[i] - h
    (e0_of_q(up) - e0_of_q(dn)) / (2 * h)
  }, numeric(1))

  delta <- sum(grad^2 * q^2 * (1 - q) / D[-n, 1])

  ci <- lt_chiang_ci(x, D)
  open_part <- (x$lx[n, 1] / 2)^2 * ci$var_qx[n, 1]
  chiang <- (ci$var_Tx[1, 1] - open_part) / x$lx[1, 1]^2

  expect_equal(chiang, delta, tolerance = 1e-6)
})


test_that("the interval is symmetric and has the requested width", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  ci <- lt_chiang_ci(x, fake_deaths(mx), conf = 0.95)

  expect_equal((ci$lower + ci$upper) / 2, x$ex, tolerance = 1e-10)
  expect_equal(ci$upper - ci$lower,
               2 * stats::qnorm(0.975) * ci$se, tolerance = 1e-10)
})


test_that("a wider confidence level gives a wider interval", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  D  <- fake_deaths(mx)

  narrow <- lt_chiang_ci(x, D, conf = 0.80)
  wide   <- lt_chiang_ci(x, D, conf = 0.99)

  expect_true(all(wide$upper >= narrow$upper))
  expect_true(all(wide$lower <= narrow$lower))
  expect_equal(wide$se, narrow$se)   # se does not depend on conf
})


test_that("variance accumulates from the oldest age upward", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  ci <- lt_chiang_ci(x, fake_deaths(mx))
  n  <- nrow(mx)

  # var_Tx is a reverse cumulative sum, so it never increases with age
  expect_true(all(diff(ci$var_Tx[, 1]) <= 0))
  expect_gt(ci$var_Tx[1, 1], ci$var_Tx[n, 1])
})


test_that("age groups with no deaths contribute nothing instead of NaN", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")

  D <- fake_deaths(mx)
  D[4:6, ] <- 0

  ci <- lt_chiang_ci(x, D)

  expect_true(all(is.finite(ci$se)))
  expect_equal(ci$var_qx[4:6, ], rep(0, 3), ignore_attr = TRUE)
})


test_that("zero deaths everywhere gives a zero-width interval", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  ci <- lt_chiang_ci(x, matrix(0, nrow = 19, ncol = 1))

  expect_equal(ci$se, matrix(0, nrow = 19, ncol = 1))
  expect_equal(ci$lower, ci$upper)
})


test_that("the lower limit is clamped at zero only when asked", {
  mx <- fake_mx(level = 20)              # extreme mortality, tiny e0
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")
  D  <- matrix(2, nrow = 19, ncol = 1)   # almost no information

  clamped <- lt_chiang_ci(x, D, truncate_at_zero = TRUE)
  raw     <- lt_chiang_ci(x, D, truncate_at_zero = FALSE)

  expect_true(all(clamped$lower >= 0))
  expect_true(all(raw$lower <= clamped$lower))
})


test_that("a single vector of deaths is recycled across draws", {
  mx <- fake_mx(n_draw = 3)
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")

  ci <- lt_chiang_ci(x, fake_deaths(mx[, 1]))
  expect_equal(dim(ci$se), c(19L, 3L))
  expect_equal(ci$se[, 1], ci$se[, 2])
})


test_that("mismatched input is rejected", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "constant_hazard")

  expect_error(lt_chiang_ci(x, matrix(1, nrow = 5, ncol = 1)), "but the life table is")
  expect_error(lt_chiang_ci(x, fake_deaths(mx), conf = 1.5), "strictly between")
  expect_error(lt_chiang_ci(list(), 1), "must be a life table")

  neg <- fake_deaths(mx); neg[3, 1] <- -1
  expect_error(lt_chiang_ci(x, neg), "finite and non-negative")
})


test_that("a Gompertz closeout is flagged as inconsistent with the open-group variance", {
  mx <- fake_mx()
  x  <- lt(mx, sex = 1, closeout = "gompertz")
  expect_warning(lt_chiang_ci(x, fake_deaths(mx)), "e_open = 1/mx")
})


test_that("a supplied ax is judged by the assumption, not by its label", {
  # The published Taiwan table closes out with ax = 1/mx, so it satisfies
  # Chiang's assumption even though lt() labels its closeout "supplied".
  ref <- taiwan_lt("male")
  x   <- lt(ref$mx, sex = 1, age = ref$age, ax = ref$ax)
  expect_equal(attr(x, "closeout"), "supplied")
  expect_silent(lt_chiang_ci(x, ref$dx))

  # a supplied open-group ax that breaks the assumption still warns
  ax <- ref$ax
  ax[length(ax)] <- ax[length(ax)] * 1.2
  y <- lt(ref$mx, sex = 1, age = ref$age, ax = ax)
  expect_warning(lt_chiang_ci(y, ref$dx), "e_open = 1/mx")
})
