# letools

Abridged life tables and exact decomposition of life expectancy differences, built for population health metrics where input data can be a posterior sample (`n_age x n_draw`) rather than just a single set of rates.

Every function operates column-wise over an `n_age x n_draw` matrix. Descriptive analyses and Bayesian models run through identical code, propagating uncertainty seamlessly without changing calling syntax.

## Install

```r
# install.packages("remotes")
remotes::install_github("jingzhong1011/letools")
```

## Quick Start

A published life table ships with the package for immediate use:

```r
library(letools)

ref <- taiwan_lt("male", abridged = TRUE) # Taiwan, 2023, Ministry of the Interior

# Default assumptions vs. official published table
lt(ref$mx, sex = 1)
lt(ref$mx, sex = 1, ax = ref$ax)          # Reproduces official results exactly

# Compare open-age closeout methods
vapply(c("kannisto", "gompertz", "constant_hazard"),
       function(m) lt(ref$mx, sex = 1, closeout = m)$ex[1, ], numeric(1))
```

## Usage

### Single Tables & Posterior Draws

```r
# 1. Single life table
x <- lt(mx, sex = 1)
x$ex[1, ] # e0

# 2. Posterior sample (mx is n_age x n_draw)
x <- lt(mx_draws, sex = 1)
quantile(x$ex[1, ], c(0.025, 0.975))
```

### Age Decomposition
```r
m <- taiwan_lt("male",   abridged = TRUE)
f <- taiwan_lt("female", abridged = TRUE)

contrib <- decomp_arriaga(
  lt(m$mx, sex = 1, ax = m$ax),
  lt(f$mx, sex = 2, ax = f$ax)
)

sum(contrib)           # Total female advantage in e0 (years)
round(contrib[, 1], 3) # Age-specific contributions
```

## Methodology

### Old-Age Closeout

The choice of method for the open-ended age group moves e0 noticeably:

| `closeout` | Assumption | When to use it |
| --- | --- | --- |
| `"kannisto"` | Logit mortality deceleration | Default. Best supported empirically past age 95. |
| `"gompertz"` | Log hazard rises linearly | For comparing against existing Gompertz-based literature. |
| `"constant_hazard"` | $e_{\text{open}} = 1/m_{\text{open}}$ | Matches tables that close out with $1/m$ (e.g., Taiwan MOI); the assumption behind `lt_chiang_ci()`. |

`mx_floor` sets lower bounds on rates to stabilize extrapolations during low-mortality projections. The printed summary reports how many cells the floor bound, as a count and a percentage; a large share means the floor is driving the result.

### Uncertainty: Draws vs. Chiang CI

- **Posterior draws:** Run `lt()` directly across draws and take quantiles of ex. Uncertainty is already in the input.
- **Chiang standard errors:** For observed count data (e.g., vital registration), use `lt_chiang_ci(x, deaths)` to propagate binomial sampling error.

> Note: Do not pass smoothed or modeled counts into `lt_chiang_ci()`, as prior smoothing artificially narrows binomial sampling intervals.

### Decomposition: Arriaga vs. Pollard

- `decomp_arriaga()`: Exact additive age decomposition summing precisely to $\Delta e_0$. Cause-specific splits are proportional approximations.
- `decomp_pollard()`: Exact cause attribution due to linearity in rate differences, but discrete abridged intervals carry residual approximation error for total $\Delta e_0$.
- *Recommendation*: Use Arriaga for total age group contributions, and Pollard when exact cause-of-death attribution is the primary metric.

### Validation

The test suite enforces mathematical consistency and reproduces empirical baselines:

- Internal identities (decrement accounting, $T_x = \sum L_x$, floating-point decomposition closure).
- Chiang standard errors: they scale as $1/\sqrt{D}$, and the closed-age-group variance matches a numerical delta method.
- Exact reproduction of the Taiwan Ministry of the Interior 2023 life table, both the published single-year table (0–85+) and its abridged version.

## License

MIT