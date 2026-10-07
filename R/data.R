#' Taiwan National Life Table, 2023
#'
#' The complete life table published by Taiwan's Ministry of the Interior for
#' 2023, covering single years of age from 0 to 85+, for the total population
#' and stratified by sex.
#'
#' The published table provides `qx`, `lx`, `dx`, `Lx`, `Tx`, and `ex`.
#' To allow exact reproduction of the published life expectancies using `lt()`,
#' the central death rate (`mx`) and average years lived by those dying in the interval (`ax`)
#' were derived algebraically from the published parameters using the identities
#' `mx = dx / Lx` and `ax = (Lx - n * l(x+n)) / dx`.
#'
#' Running `lt()` on the derived `mx` alone overshoots the published `e0` by
#' about 0.8 years, from the open-group closeout rather than from `a0`. The
#' Ministry uses `ax[85+] = 1/mx[85+]`, so `closeout = "constant_hazard"`
#' reproduces it; supplying `ax` reproduces every column exactly.
#'
#' @param series Character string indicating the target population: `"total"`, `"male"`, or `"female"`.
#' @param abridged Logical. If `TRUE`, collapses the single-year age groups into the
#'   standard abridged structure (0, 1, 5, 10, ..., 85+). `Lx` and `dx` are summed
#'   within each interval, and `mx` and `ax` are recalculated accordingly.
#' @return A data frame containing the life table columns `age`, `mx`, `ax`, `qx`,
#'   `lx`, `dx`, `Lx`, `Tx`, and `ex`, ordered by age.
#'
#' @source Department of Statistics, Ministry of the Interior, Taiwan.
#'   National Life Tables, 2023.
#'
#' @examples
#' ref <- taiwan_lt("male")
#' head(ref)
#'
#' # Reproduce the published table exactly by supplying the derived ax
#' x <- lt(ref$mx, sex = 1, age = ref$age, ax = ref$ax)
#' max(abs(as.vector(x$ex) - ref$ex))
#'
#' # Calculate life expectancy using the package's default ax assumptions
#' lt(ref$mx, sex = 1, age = ref$age)$ex[1, ]
#' ref$ex[1]
#' @export
#'
taiwan_lt <- function(series = c("total", "male", "female"), abridged = FALSE) {
  series <- match.arg(series)

  path <- system.file("extdata", "taiwan_lt_2023.csv", package = "letools")
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "Could not find the bundled Taiwan life table. If you are working from ",
      "a source checkout, load the package with devtools::load_all() rather ",
      "than sourcing the R/ files directly.",
      call. = FALSE
    )
  }

  ref <- utils::read.csv(path, stringsAsFactors = FALSE)
  ref <- ref[ref$sex == series, setdiff(names(ref), "sex")]
  ref <- ref[order(ref$age), ]
  rownames(ref) <- NULL

  if (!abridged) return(ref)

  abridge_lt(ref)
}


#' Collapse a single-year life table to the standard abridged structure
#'
#' `Lx` and `dx` are additive across ages, so they sum within each band; `lx` is
#' taken at the band boundary. `mx` and `ax` then follow from the same
#' identities used to build the single-year table, which keeps the abridged
#' version internally consistent with the published one.
#'
#' @param ref A single-year life table as returned by [taiwan_lt()].
#' @param edges Left-hand boundaries of the abridged bands.
#' @return A data frame with the same columns, one row per band.
#' @noRd
abridge_lt <- function(ref, edges = c(0, 1, seq(5, 85, 5))) {
  n_band <- length(edges)
  top <- max(ref$age)

  upper <- c(edges[-1], top + 1L)
  nx <- c(diff(edges), NA_real_)

  lx <- Lx <- dx <- numeric(n_band)
  for (i in seq_len(n_band)) {
    inside <- ref$age >= edges[i] & ref$age < upper[i]
    lx[i] <- ref$lx[ref$age == edges[i]]
    Lx[i] <- sum(ref$Lx[inside])
    dx[i] <- sum(ref$dx[inside])
  }

  ax <- numeric(n_band)
  ax[-n_band] <- (Lx[-n_band] - nx[-n_band] * lx[-1]) / dx[-n_band]
  ax[n_band] <- Lx[n_band] / lx[n_band]

  Tx <- rev(cumsum(rev(Lx)))

  data.frame(
    age = edges,
    mx  = dx / Lx,
    ax  = ax,
    qx  = dx / lx,
    lx  = lx,
    dx  = dx,
    Lx  = Lx,
    Tx  = Tx,
    ex  = Tx / lx
  )
}
