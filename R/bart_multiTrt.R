#' Bayesian Additive Regression Trees (BART)
#'
#' This function implements the BART method. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
#' @param discard discarding rules for BART method, inherited from causal_multi_treat.R
#' @param estimand causal estimands. Please select "ATT" or "ATE"
#' @param k For binary y, k is the number of prior standard deviations f(x) is away from +/-3. The bigger k is, the more conservative the fitting will be.
#' @param ntree The number of trees in the sum
#' @param ndpost The number of posterior draws returned
#' @param nskip Number of MCMC iterations to be treated as burn in
#' @param reference Reference group for ATT
#'
#' @return list with 2 elements for ATT effect. It contains
#' \item{ATT12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATT13:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' list with 3 elements for ATE effect. It contains
#' \item{ATE12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#' \item{ATE13:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATE23:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' @export
#' @examples
#'library(CIMTx)
#'set.seed(3242019)
#'idata = data_gen(n = 3, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'reference_trt <- 2
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATT", discard = "No", ndpost = 10, reference_trt = 2)
bart_multiTrt = function(y, x, trt, discard = FALSE, estimand="ATE", k=2, ntree=100, ndpost=parent.frame()$ndpost, nskip=1000, reference = parent.frame()$reference_trt) {
  x <- x[, -1]
  # Data structure
  #        Y(1) Y(2) Y(3)
  # trt=1   *    ?    ?
  # trt=2   ?    *    ?
  # trt=3   ?    ?    *

  #        Y(1) Y(2) Y(3)
  # trt=1  y11  y12  y13
  # trt=2  y21  y22  y23
  # trt=3  y31  y32  y33

  if (estimand=="ATE") {
    bart_est = bart_multiTrt_ate(y, x, trt, discard = FALSE, k=2, ntree=100, ndpost=parent.frame()$ndpost, nskip=1000)
  }

  if (estimand=="ATT") {
    bart_est = bart_multiTrt_att(y, x, trt, discard = FALSE, k=2, ntree=100, ndpost=parent.frame()$ndpost, nskip=1000, reference = parent.frame()$reference_trt)
  }

  return(bart_est)
}
