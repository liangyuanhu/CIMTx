#######################################################################
# Function to summarize posterior samples of RD, RR and OR
# Author: Chenyang Gu
# Date: 02/23/2019
#######################################################################


#' Summarize posterior samples
#'
#' This function summarize posterior samples of RD, RR and OR. Please use our main function causal_multi_treat.R.
#' @param RD_est vector of estimation for RD
#' @param RR_est vector of estimation for RR
#' @param OR_est vector of estimation for OR
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
#' library(CIMTx)
#' postSumm(RD_est = 1:10, RR_est = 11:20, OR_est = 1:10)
postSumm = function(RD_est, RR_est, OR_est) {
    # Risk difference (RD)
    RD_mean = mean(RD_est)
    RD_se = stats::sd(RD_est)
    RD_lower = stats::quantile(RD_est, probs=0.025, na.rm = T)
    RD_upper = stats::quantile(RD_est, probs=0.975, na.rm = T)

    # Relative risk (RR)
    RR_mean = mean(RR_est)
    RR_se = stats::sd(RR_est)
    RR_lower = stats::quantile(RR_est, probs=0.025, na.rm = T)
    RR_upper = stats::quantile(RR_est, probs=0.975, na.rm = T)

    # Odds ratio (OR)
    OR_mean = mean(OR_est)
    OR_se = stats::sd(OR_est)
    OR_lower = stats::quantile(OR_est, probs=0.025, na.rm = T)
    OR_upper = stats::quantile(OR_est, probs=0.975, na.rm = T)

    # summarize results
    RD = c(RD_mean, RD_se, RD_lower, RD_upper)
    RR = c(RR_mean, RR_se, RR_lower, RR_upper)
    OR = c(OR_mean, OR_se, OR_lower, OR_upper)

    res = rbind(RD, RR, OR)
    colnames(res) = c("EST","SE","LOWER","UPPER")
    return(res)
}























