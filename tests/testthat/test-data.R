# ---------------------------------------------------------------------------
# The bundled table is only useful if it is internally consistent, so these
# check the identities it claims rather than hard-coding values. The abridged
# form gets the same treatment: collapsing must preserve total person-years and
# total deaths exactly, because Lx and dx are additive across ages.
# ---------------------------------------------------------------------------

test_that("every published series loads with the expected shape", {
  for (series in c("total", "male", "female")) {
    ref <- taiwan_lt(series)
    expect_equal(nrow(ref), 86L, info = series)
    expect_equal(ref$age, 0:85, info = series)
    expect_named(ref, c("age", "mx", "ax", "qx", "lx", "dx", "Lx", "Tx", "ex"))
    expect_true(all(vapply(ref, function(x) all(is.finite(x)), logical(1))),
                info = series)
  }
})


test_that("the recovered mx and ax satisfy the identities they came from", {
  ref <- taiwan_lt("male")
  n <- nrow(ref)

  expect_equal(ref$mx, ref$dx / ref$Lx, tolerance = 1e-12)
  expect_equal(ref$Lx[-n], ref$lx[-1] + ref$ax[-n] * ref$dx[-n], tolerance = 1e-8)
  expect_equal(ref$ax[n], ref$ex[n], tolerance = 1e-10)
})


test_that("the published table is a valid life table", {
  ref <- taiwan_lt("total")
  n <- nrow(ref)

  expect_equal(ref$lx[1], 1e5)
  expect_equal(sum(ref$dx), 1e5, tolerance = 1e-6)
  expect_true(all(diff(ref$lx) < 0))
  expect_equal(ref$qx[n], 1)
  expect_equal(ref$Tx, rev(cumsum(rev(ref$Lx))), tolerance = 1e-6)
})


test_that("female life expectancy exceeds male at every age", {
  m <- taiwan_lt("male")
  f <- taiwan_lt("female")
  expect_true(all(f$ex > m$ex))
})


test_that("abridging preserves person-years and deaths exactly", {
  for (series in c("total", "male", "female")) {
    full <- taiwan_lt(series)
    ab   <- taiwan_lt(series, abridged = TRUE)

    expect_equal(nrow(ab), 19L, info = series)
    expect_equal(ab$age, c(0, 1, seq(5, 85, 5)), info = series)

    # Lx and dx are additive across ages, so nothing may be lost
    expect_equal(sum(ab$Lx), sum(full$Lx), tolerance = 1e-8, info = series)
    expect_equal(sum(ab$dx), sum(full$dx), tolerance = 1e-8, info = series)

    # e0 is a property of the whole table and must survive collapsing
    expect_equal(ab$ex[1], full$ex[1], tolerance = 1e-8, info = series)
  }
})


test_that("the abridged table satisfies the same identities", {
  ab <- taiwan_lt("male", abridged = TRUE)
  n  <- nrow(ab)
  nx <- c(diff(ab$age), NA)

  expect_equal(ab$mx, ab$dx / ab$Lx, tolerance = 1e-12)
  expect_equal(ab$Lx[-n],
               nx[-n] * ab$lx[-1] + ab$ax[-n] * ab$dx[-n], tolerance = 1e-8)
  expect_equal(ab$qx, ab$dx / ab$lx, tolerance = 1e-12)
  expect_true(all(ab$ax[3:(n - 1)] > 0 & ab$ax[3:(n - 1)] < 5))
})


test_that("lt() round-trips the abridged table when given its own ax", {
  ab <- taiwan_lt("male", abridged = TRUE)
  x  <- lt(ab$mx, sex = 1, ax = ab$ax)

  for (col in c("qx", "lx", "dx", "Lx", "Tx", "ex")) {
    expect_equal(as.vector(x[[col]]), ab[[col]], tolerance = 1e-8, info = col)
  }
})


test_that("an unknown series is rejected", {
  expect_error(taiwan_lt("everyone"), "'arg' should be one of")
})
