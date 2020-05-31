
#' Truncation
#'
#' This function implements the truncation feature when estimand is ATT. Please use our main function causal_multi_treat.R.
#'
#' @param x vector to be trimmed
#' @param trim_alpha alpha values for IPTW weight trimming, inherited from causal_multi_treat.R
#'
#' @return vector trimmed
#' @export
#' @examples
#' library(CIMTx)
#' trunc_fun(1:10)
trunc_fun <- function(x, trim_alpha = 0.05) {
  pmin(stats::quantile(x, (1-trim_alpha)), pmax(stats::quantile(x, trim_alpha), x))
}
