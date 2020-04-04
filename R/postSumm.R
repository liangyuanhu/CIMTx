#######################################################################
# Function to summarize posterior samples of RD, RR and OR
# Author: Chenyang Gu
# Date: 02/23/2019
#######################################################################


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























