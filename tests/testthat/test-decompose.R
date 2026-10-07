# ---------------------------------------------------------------------------
# The property being tested here is the one worth having: Arriaga is an exact
# decomposition, so the contributions must sum to the difference in life
# expectancy at birth. If any change to the life table code breaks that -- a
# mishandled open-ended group, an off-by-one in lx, a radix that stops being
# 100,000 -- these tests fail immediately instead of producing plausible but
# wrong figures.
# ---------------------------------------------------------------------------

# A minimal life table builder, used only to generate test fixtures. It is
# deliberately independent of the package's own lt() so that a bug in lt()
# cannot hide by being present on both sides of the comparison.
make_lt <- function(mx, ax = NULL, nx = c(1, 4, rep(5, 16), NA)) {
  mx <- as.matrix(mx)
  n  <- nrow(mx)
  m  <- ncol(mx)

  if (is.null(ax)) ax <- c(0.33, 1.35, rep(2.5, n - 3), NA)
  ax_mat <- matrix(ax, nrow = n, ncol = m)
  nx_mat <- matrix(nx, nrow = n, ncol = m)

  qx <- (nx_mat * mx) / (1 + (nx_mat - ax_mat) * mx)
  qx[n, ] <- 1
  px <- 1 - qx

  lx <- rbind(rep(1e5, m), 1e5 * apply(px[-n, , drop = FALSE], 2, cumprod))
  dx <- lx * qx

  Lx_closed <- nx_mat[-n, , drop = FALSE] * lx[-1, , drop = FALSE] +
    ax_mat[-n, , drop = FALSE] * dx[-n, , drop = FALSE]
  Lx_open <- lx[n, ] / mx[n, ]
  Lx <- rbind(Lx_closed, Lx_open)

  Tx <- apply(Lx[n:1, , drop = FALSE], 2, cumsum)[n:1, , drop = FALSE]

  # unname everything: rbind() picks up row names from symbol arguments, which
  # would otherwise propagate into colSums() output and break comparisons
  lapply(list(lx = lx, Lx = Lx, ex = Tx / lx, dx = dx, qx = qx, mx = mx), unname)
}

# Gompertz-ish schedule over the standard 19 abridged age groups, with an
# infant excess. `level` scales overall mortality so two comparable tables can
# be generated from one function.
fake_mx <- function(level = 1, n_draw = 1L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ages <- c(0, 1, seq(5, 85, 5))
  base <- 0.00005 * exp(0.085 * ages)
  base[1] <- 0.006
  base[2] <- 0.0004

  mx <- outer(base * level, rep(1, n_draw))
  if (n_draw > 1L) mx <- mx * matrix(runif(length(mx), 0.9, 1.1), nrow = length(base))
  mx
}


test_that("Arriaga contributions sum to the difference in e0", {
  lt_a <- make_lt(fake_mx(level = 1.00))
  lt_b <- make_lt(fake_mx(level = 0.75))

  contrib <- decomp_arriaga(lt_a, lt_b, check = FALSE)
  total   <- lt_b$ex[1, ] - lt_a$ex[1, ]

  expect_equal(colSums(contrib), total, tolerance = 1e-9)
})


test_that("Arriaga closes for every draw when given a posterior sample", {
  lt_a <- make_lt(fake_mx(level = 1.00, n_draw = 50, seed = 1))
  lt_b <- make_lt(fake_mx(level = 0.75, n_draw = 50, seed = 2))

  contrib <- decomp_arriaga(lt_a, lt_b, check = FALSE)
  total   <- lt_b$ex[1, ] - lt_a$ex[1, ]

  expect_equal(dim(contrib), c(19L, 50L))
  expect_equal(colSums(contrib), total, tolerance = 1e-9)
})


test_that("a single life table is just the one-column case", {
  lt_a <- make_lt(fake_mx(level = 1.00))
  lt_b <- make_lt(fake_mx(level = 0.75))

  wide_a <- lapply(lt_a, function(x) cbind(x, x, x))
  wide_b <- lapply(lt_b, function(x) cbind(x, x, x))

  one  <- decomp_arriaga(lt_a, lt_b)
  many <- decomp_arriaga(wide_a, wide_b)

  expect_equal(many[, 2, drop = FALSE], one, ignore_attr = TRUE)
})


test_that("decomposing a table against itself gives exactly zero", {
  lt_a <- make_lt(fake_mx())
  expect_equal(
    decomp_arriaga(lt_a, lt_a),
    matrix(0, nrow = 19, ncol = 1)
  )
})


test_that("the direction convention is the documented one", {
  lt_worse  <- make_lt(fake_mx(level = 1.00))
  lt_better <- make_lt(fake_mx(level = 0.75))

  # from worse to better, e0 rises, so contributions are positive on net
  expect_gt(sum(decomp_arriaga(lt_worse, lt_better)), 0)
  # reversing the arguments flips the sign of the total
  expect_lt(sum(decomp_arriaga(lt_better, lt_worse)), 0)
})


test_that("the closure check catches mismatched life tables", {
  lt_a <- make_lt(fake_mx(level = 1.00))
  lt_b <- make_lt(fake_mx(level = 0.75))

  # ex taken from a third table it does not belong to
  lt_broken <- lt_b
  lt_broken$ex <- make_lt(fake_mx(level = 0.50))$ex

  expect_error(decomp_arriaga(lt_a, lt_broken), "did not close")
})


test_that("mismatched shapes are rejected before any arithmetic", {
  lt_a <- make_lt(fake_mx(n_draw = 5, seed = 3))
  lt_b <- make_lt(fake_mx(n_draw = 3, seed = 4))
  expect_error(decomp_arriaga(lt_a, lt_b), "different shapes")
})


test_that("cause contributions sum back to the age contributions", {
  mx_a <- fake_mx(level = 1.00)
  mx_b <- fake_mx(level = 0.75)
  lt_a <- make_lt(mx_a)
  lt_b <- make_lt(mx_b)

  # an arbitrary but exhaustive split of all-cause mortality into three causes
  shares <- list(cvd = 0.5, cancer = 0.3, other = 0.2)
  cause_a <- lapply(shares, function(s) mx_a * s)
  cause_b <- lapply(shares, function(s) mx_b * s)

  contrib <- decomp_arriaga(lt_a, lt_b)
  by_cause <- decomp_by_cause(contrib, mx_a, mx_b, cause_a, cause_b)

  expect_named(by_cause, c("cvd", "cancer", "other"))
  expect_equal(Reduce(`+`, by_cause), contrib, tolerance = 1e-9)
})


test_that("causes that did not change contribute nothing", {
  mx_a <- fake_mx(level = 1.00)
  mx_b <- mx_a
  mx_b[10:19, ] <- mx_b[10:19, ] * 0.7  # only older ages improve

  lt_a <- make_lt(mx_a)
  lt_b <- make_lt(mx_b)

  cause_a <- list(moving = mx_a * 0.6, static = mx_a * 0.4)
  cause_b <- list(moving = mx_b - mx_a * 0.4, static = mx_a * 0.4)

  by_cause <- decomp_by_cause(
    decomp_arriaga(lt_a, lt_b), mx_a, mx_b, cause_a, cause_b
  )

  expect_equal(sum(by_cause$static), 0, tolerance = 1e-9)
})


test_that("mismatched cause names are rejected", {
  mx_a <- fake_mx(1.00); mx_b <- fake_mx(0.75)
  contrib <- decomp_arriaga(make_lt(mx_a), make_lt(mx_b))

  expect_error(
    decomp_by_cause(contrib, mx_a, mx_b,
                    list(a = mx_a), list(b = mx_b)),
    "Cause names differ"
  )
})


test_that("cause split rejects matrices that do not match the contributions", {
  # Without the check, a short cause matrix is recycled by R and the split
  # comes back with plausible-looking but meaningless numbers.
  mx_a <- fake_mx(1.00); mx_b <- fake_mx(0.75)
  contrib <- decomp_arriaga(make_lt(mx_a), make_lt(mx_b))
  short <- mx_a[1:10, , drop = FALSE]

  expect_error(
    decomp_by_cause(contrib, short, mx_b, list(a = mx_a), list(a = mx_b)),
    "`mx_from` is 10 x 1 but should be 19 x 1"
  )
  expect_error(
    decomp_by_cause(contrib, mx_a, mx_b,
                    list(a = mx_a * 0.5, b = short),
                    list(a = mx_b * 0.5, b = mx_b * 0.5)),
    'mx_cause_from\\[\\["b"\\]\\]'
  )
})


# --- cause split, the paths the proportional tests above do not reach ---

test_that("causes moving in opposite directions are attributed separately", {
  mx_from <- fake_mx(level = 1.00)
  mx_to <- fake_mx(level = 0.75)

  # one cause improves faster than all-cause, the other worsens
  cause_from <- list(better = mx_from * 0.5, worse = mx_from * 0.5)
  cause_to <- list(better = mx_to * 0.5 - mx_from * 0.15,
                   worse = mx_to * 0.5 + mx_from * 0.15)

  contrib <- decomp_arriaga(make_lt(mx_from), make_lt(mx_to))
  by_cause <- decomp_by_cause(contrib, mx_from, mx_to, cause_from, cause_to)

  expect_gt(sum(by_cause$better), sum(by_cause$worse))
  expect_equal(Reduce(`+`, by_cause), contrib, tolerance = 1e-9)
})


test_that("age groups with no all-cause change give zero, not NaN", {
  mx_from <- fake_mx(level = 1.00)
  mx_to <- mx_from
  mx_to[10:19, ] <- mx_to[10:19, ] * 0.7   # ages under 45 do not move at all

  contrib <- decomp_arriaga(make_lt(mx_from), make_lt(mx_to))
  by_cause <- decomp_by_cause(
    contrib, mx_from, mx_to, list(only = mx_from), list(only = mx_to)
  )

  expect_true(all(is.finite(by_cause$only)))
  expect_equal(by_cause$only[1:9, ], contrib[1:9, ] * 0, tolerance = 1e-12)
})


test_that("the cause split carries through a posterior sample", {
  mx_from <- fake_mx(level = 1.00, n_draw = 4, seed = 11)
  mx_to <- fake_mx(level = 0.75, n_draw = 4, seed = 12)

  contrib <- decomp_arriaga(make_lt(mx_from), make_lt(mx_to))
  by_cause <- decomp_by_cause(
    contrib, mx_from, mx_to,
    list(a = mx_from * 0.6, b = mx_from * 0.4),
    list(a = mx_to * 0.6, b = mx_to * 0.4)
  )

  expect_equal(dim(by_cause$a), dim(contrib))
  expect_equal(Reduce(`+`, by_cause), contrib, tolerance = 1e-9)
})


# --- Pollard ---

test_that("Pollard gives exactly zero against an identical table", {
  mx <- fake_mx()
  x <- make_lt(mx)
  expect_equal(decomp_pollard(x, x, mx, mx), matrix(0, nrow = nrow(mx), ncol = 1))
})


test_that("Pollard follows the same direction convention as Arriaga", {
  mx_from <- fake_mx(level = 1.00)
  mx_to <- fake_mx(level = 0.75)
  lt_from <- make_lt(mx_from)
  lt_to <- make_lt(mx_to)

  forward <- decomp_pollard(lt_from, lt_to, mx_from, mx_to, check = FALSE)
  reverse <- decomp_pollard(lt_to, lt_from, mx_to, mx_from, check = FALSE)

  expect_gt(sum(forward), 0)                       # `to` has lower mortality
  expect_equal(sum(forward), -sum(reverse), tolerance = 1e-9)
})


test_that("Pollard converges on the true difference as the intervals narrow", {
  # Pollard is exact only in continuous time. The discrete residual is first
  # order in the interval width, so moving from five-year bands to single years
  # must shrink it by roughly a factor of five. A residual that does not shrink
  # means the weights are wrong, not that the approximation is coarse.
  residual <- function(abridged) {
    ref_from <- taiwan_lt("male", abridged = abridged)
    ref_to <- taiwan_lt("female", abridged = abridged)

    lt_from <- lt(ref_from$mx, sex = 1, age = ref_from$age, ax = ref_from$ax)
    lt_to <- lt(ref_to$mx, sex = 2, age = ref_to$age, ax = ref_to$ax)

    contrib <- decomp_pollard(lt_from, lt_to, ref_from$mx, ref_to$mx,
                              check = FALSE)
    abs(sum(contrib) - (lt_to$ex[1, 1] - lt_from$ex[1, 1]))
  }

  wide <- residual(TRUE)
  narrow <- residual(FALSE)

  expect_gt(wide / narrow, 3)
  expect_lt(narrow, 0.2)
})


test_that("the closure check passes discretisation but catches a mismatched mx", {
  ref_from <- taiwan_lt("male", abridged = TRUE)
  ref_to <- taiwan_lt("female", abridged = TRUE)
  lt_from <- lt(ref_from$mx, sex = 1, ax = ref_from$ax)
  lt_to <- lt(ref_to$mx, sex = 2, ax = ref_to$ax)

  # five-year bands leave 9 per cent, which discretisation accounts for
  expect_silent(decomp_pollard(lt_from, lt_to, ref_from$mx, ref_to$mx))

  # rates off by a tenth is the mildest mismatch worth catching, at 33 per cent
  expect_warning(
    decomp_pollard(lt_from, lt_to, ref_from$mx * 1.1, ref_to$mx),
    "rel_tol"
  )
  expect_silent(
    decomp_pollard(lt_from, lt_to, ref_from$mx * 1.1, ref_to$mx, check = FALSE)
  )
})


test_that("splitting Pollard across causes is exact, not proportional", {
  # This is the property that makes Pollard the right choice for cause of
  # death work. Its contribution is linear in the rate difference, so the
  # proportional split done by decomp_by_cause coincides exactly with running
  # Pollard on the cause-specific rates directly. Arriaga has no such identity.
  mx_from <- fake_mx(level = 1.00)
  mx_to <- fake_mx(level = 0.75)
  lt_from <- make_lt(mx_from)
  lt_to <- make_lt(mx_to)

  # a deliberately uneven split, so a proportional shortcut would show up
  cause_from <- list(early = mx_from * 0.8, late = mx_from * 0.2)
  cause_to <- list(early = mx_to * 0.3, late = mx_to * 0.7)

  contrib <- decomp_pollard(lt_from, lt_to, mx_from, mx_to, check = FALSE)
  by_cause <- decomp_by_cause(contrib, mx_from, mx_to, cause_from, cause_to)

  for (cause in names(cause_from)) {
    direct <- decomp_pollard(lt_from, lt_to,
                             cause_from[[cause]], cause_to[[cause]],
                             check = FALSE)
    expect_equal(by_cause[[cause]], direct, tolerance = 1e-12, info = cause)
  }
})


test_that("Pollard does not depend on either table's radix", {
  # Each table's person-years must be scaled by its own radix. Scaling both by
  # the `from` radix shrinks the `to` half of every weight by the ratio of the
  # two, which here is 1e5.
  mx_from <- fake_mx(level = 1.00)
  mx_to <- fake_mx(level = 0.75)

  same <- decomp_pollard(lt(mx_from, sex = 1, radix = 1e5),
                         lt(mx_to, sex = 1, radix = 1e5),
                         mx_from, mx_to, check = FALSE)
  mixed <- decomp_pollard(lt(mx_from, sex = 1, radix = 1e5),
                          lt(mx_to, sex = 1, radix = 1),
                          mx_from, mx_to, check = FALSE)

  expect_equal(mixed, same, tolerance = 1e-9)
})


test_that("Pollard rejects rate matrices that do not match the life tables", {
  mx_from <- fake_mx(level = 1.00)
  mx_to <- fake_mx(level = 0.75)
  lt_from <- make_lt(mx_from)
  lt_to <- make_lt(mx_to)

  expect_error(decomp_pollard(lt_from, lt_to, mx_from[1:18, , drop = FALSE], mx_to),
               "`mx_from` is 18 x 1 but should be 19 x 1")
  expect_error(decomp_pollard(lt_from, lt_to, mx_from, cbind(mx_to, mx_to)),
               "`mx_to` is 19 x 2")
})


test_that("Pollard handles a posterior sample", {
  mx_from <- fake_mx(level = 1.00, n_draw = 6, seed = 21)
  mx_to <- fake_mx(level = 0.80, n_draw = 6, seed = 22)
  lt_from <- make_lt(mx_from)
  lt_to <- make_lt(mx_to)

  contrib <- decomp_pollard(lt_from, lt_to, mx_from, mx_to, check = FALSE)

  expect_equal(dim(contrib), c(19L, 6L))
  expect_true(all(colSums(contrib) > 0))
})
