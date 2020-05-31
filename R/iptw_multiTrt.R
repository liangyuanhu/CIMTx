
#' Inverse probability of treatment weighting (IPTW)
#'
#' This function implements the IPTW method. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param trt numeric vector for the treatment indicator
#' @param psdat data frame containing the treatment indicator and covariates
#' @param estimand causal estimands, "ATT" or "ATE"
#' @param method methods for causal inference with multiple treatments, inherited from causal_multi_treat.R
#' @param trim_alpha alpha values for IPTW weight trimming, inherited from causal_multi_treat.R
#' @param SL.library methods specified with SL.library in Superlearner package, inherited from causal_multi_treat.R
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
#'set.seed(1)
#'idata = data_gen(n = 500, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'causal_multi_treat(y = y,trt = trt_ind,
#'method = "IPTW-Logistics", estimand = "ATT", reference_trt = 2)

iptw_multiTrt = function(y, trt, psdat, estimand = "ATE", method, trim_alpha = parent.frame()$trim_alpha, SL.library = parent.frame()$SL.library, reference = parent.frame()$reference_trt) {
  if (estimand == "ATE") {
    iptw_est = iptw_multiTrt_ate(y, trt, psdat, method, trim_alpha = parent.frame()$trim_alpha,SL.library= parent.frame()$SL.library)
  }
  if (estimand == "ATT") {
    iptw_est = iptw_multiTrt_att(y, trt, psdat,method, trim_alpha = parent.frame()$trim_alpha, SL.library= parent.frame()$SL.library, reference = parent.frame()$reference_trt)
  }
  return(iptw_est)
}
