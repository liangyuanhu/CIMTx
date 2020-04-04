
trunc_fun <- function(x, trim_alpha = 0.05) {
  pmin(stats::quantile(x, (1-trim_alpha)), pmax(stats::quantile(x, trim_alpha), x))
}
