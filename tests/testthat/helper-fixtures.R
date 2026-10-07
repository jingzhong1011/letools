# ---------------------------------------------------------------------------
# The reference life table now ships with the package, so tests read it through
# the same accessor users do. That is deliberate: if taiwan_lt() breaks, the
# golden tests should fail rather than quietly fall back to a private copy.
# ---------------------------------------------------------------------------

read_reference_lt <- function(series) {
  taiwan_lt(series)
}
