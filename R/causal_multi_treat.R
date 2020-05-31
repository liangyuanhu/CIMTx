#' Estimation of causal effects of multiple treatments
#'
#' This function estimates the causal effects of multiple treatments with a binary outcome.
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
#' @param method methods for causal inference with multiple treatments. Please select one of the following methods:
#'\enumerate{
#'\item {Regression Adjustment: }{Logistics regression to impute missing outcomes}
#'\item {VM Matching: }{vector matching}
#'\item {BART: }{Bayesian Additive Regression Trees}
#'\item {TMLE: }{Targeted maximum likelihood}
#'\item {IPTW-Logistics: }{Inverse probability of treatment weighting (IPTW) with weights from logistics regression}
#'\item {IPTW-Logistics-Trim: }{IPTW with trimmed weights from logistics regression}
#'\item {IPTW-GBM: }{IPTW with weights from generalized boosted method}
#'\item {IPTW-GBM-Trim: }{IPTW with trimmed weights from generalized boosted method}
#'\item {IPTW-Superlearner: }{IPTW with weights from superlearner}
#'\item {IPTW-Superlearner-Trim: }{IPTW with trimmed weights from superlearner}
#' }
#' @param trim_alpha alpha values for IPTW weight trimming. The default is 0.05, which means we truncate upper 95\% and lower 5\% of the weights for further IPTW estimation. The default is a combination of SL.glm, SL.gam and SL.knn.
#' @param SL.library methods specified with SL.library in Superlearner package
#' @param discard discarding rules for BART method. Please select "No", "Lenient" or "Stringent". The default is "No".
#' @param estimand causal estimands. Please select "ATT" or "ATE"
#' @param reference_trt Reference group for ATT
#' @param ndpost number of independent simulation draws to create
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
#'
#' @export causal_multi_treat
#'
#' @examples
#' library(CIMTx)
#'set.seed(3242019)
#'idata = data_gen(n = 12, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'
#'# Regression Adjustment
#'causal_multi_treat(y = y, x = idata$trtdat,ndpost = 10,
#'trt = trt_ind, method ="Regression Adjustment", estimand = "ATT", reference_trt = 3)
#'causal_multi_treat(y = y, x = idata$trtdat,ndpost = 10,
#'trt = trt_ind, method ="Regression Adjustment",
#'estimand = "ATE")
#'
#'
#'
#'# BART with and without discarding
#'\dontrun{
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATE", discard = "No")
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATE", discard = "No")
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATT", discard = "Stringent")
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATT", discard = "Stringent")
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATT", discard = "Lenient")
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATT", discard = "Lenient")
#'
#'# VM Matching
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind,method = "VM Matching", estimand = "ATT")
#'
#'# IPTW-related methods
#'causal_multi_treat(y = y,trt = trt_ind,
#'method = "IPTW-Logistics", estimand = "ATT")
#'causal_multi_treat(y = y,trt = trt_ind,
#'method = "IPTW-Logistics", estimand = "ATE")
#'causal_multiple_treatment(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "IPTW-GBM", estimand = "ATE")
#'causal_multiple_treatment(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "IPTW-GBM-Trim", estimand = "ATE")
#'causal_multiple_treatment(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "IPTW-Superlearner", estimand = "ATE")
#'causal_multiple_treatment(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "IPTW-Superlearner-Trim", estimand = "ATE")
#'causal_multiple_treatment(y = y, x = idata$trtdat,
#'trt = trt_ind,method = "IPTW-Superlearner", estimand = "ATT")
#' causal_multiple_treatment(y = y, x = idata$trtdat,
#' trt = trt_ind, method = "IPTW-Superlearner-Trim", estimand = "ATT")
#'
#'}


causal_multi_treat <- function(y, x, trt, method, discard = "No", estimand, trim_alpha = 0.05, SL.library = c("SL.glm", "SL.gam", "SL.knn"), reference_trt = 1, ndpost = 1000){
  if (method == "Regression Adjustment" && estimand == "ATE") {
    result <- regadj_multiTrt(
      y,
      x,
      trt,
      ndpost = parent.frame()$ndpost,
      estimand = "ATE"
    )
  } else if (method == "Regression Adjustment" && estimand == "ATT"){
    result <- regadj_multiTrt(
      y,
      x,
      trt,
      ndpost = parent.frame()$ndpost,
      reference = parent.frame()$reference_trt,
      estimand = "ATT"
    )
  } else if (method == "VM Matching" && estimand == "ATT") {
    result <- vm_multiTrt_att(y = y, x, trt)
  } else if (method == "BART" && estimand == "ATE" && discard == "No"){
    result <- bart_multiTrt(
      y,
      x,
      trt,
      estimand = "ATE",
      discard = "No",
      ndpost = parent.frame()$ndpost
    )
  } else if (method == "BART" && estimand == "ATE" && discard == "Stringent"){
    result <- bart_multiTrt(
      y,
      x,
      trt,
      estimand = "ATE",
      ndpost = parent.frame()$ndpost,
      discard = "Stringent"
    )
  } else if (method == "BART" && estimand == "ATE" && discard == "Lenient"){
    result <- bart_multiTrt(
      y,
      x,
      trt,
      estimand = "ATE",
      ndpost = parent.frame()$ndpost,
      discard = "Lenient"
    )
  } else if (method == "BART" && estimand == "ATT" && discard == "No"){
    result <- bart_multiTrt(
      y,
      x,
      trt,
      estimand = "ATT",
      ndpost = parent.frame()$ndpost,
      reference = parent.frame()$reference_trt,
      discard = "No"
    )
  } else if (method == "BART" && estimand == "ATT" && discard == "Stringent"){
    result <- bart_multiTrt(
      y,
      x,
      trt,
      estimand = "ATT",
      ndpost = parent.frame()$ndpost,
      reference = parent.frame()$reference_trt,
      discard = "Stringent"
    )
  } else if (method == "BART" && estimand == "ATT" && discard == "Lenient"){
    result <- bart_multiTrt(
      y,
      x,
      trt,
      estimand = "ATT",
      ndpost = parent.frame()$ndpost,
      reference = parent.frame()$reference_trt,
      discard = "Lenient"
    )
  } else if (method == "IPTW-Logistics" && estimand == "ATT" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATT",
      reference = parent.frame()$reference_trt,
      method = method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Logistics-Trim" && estimand == "ATT" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATT",
      reference = parent.frame()$reference_trt,
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-GBM" && estimand == "ATT" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATT",
      reference = parent.frame()$reference_trt,
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-GBM-Trim" && estimand == "ATT" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATT",
      reference = parent.frame()$reference_trt,
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Logistics" && estimand == "ATE" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATE",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Logistics-Trim" && estimand == "ATE" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATE",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-GBM" && estimand == "ATE" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATE",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-GBM-Trim" && estimand == "ATE" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATE",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Superlearner" && estimand == "ATT" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATT",
      reference = parent.frame()$reference_trt,
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Superlearner-Trim" && estimand == "ATT" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATT",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Superlearner" && estimand == "ATE" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATE",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "IPTW-Superlearner-Trim" && estimand == "ATE" ){
    result <- iptw_multiTrt(
      y,
      trt,
      trim_alpha =trim_alpha,
      estimand = "ATE",
      method= method,
      SL.library = SL.library
    )
  } else if (method == "TMLE" && estimand == "ATE" ){
    result <- tmle(
      y,
      trt,
      x
      # estimand = "ATE",
      # method
    )
  }
  return(result)
}
