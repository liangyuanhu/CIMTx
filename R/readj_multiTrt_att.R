# Estimate average treatment effect on the treated (ATT)
# The default reference group is 1st group
#' Regression Adjustment when estimand is ATT
#'
#' This function implements the regression adjustment method when estimand is ATT. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
#' @param ndpost number of independent simulation draws to create
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
#' library(CIMTx)
#' set.seed(3242019)
#' idata = data_gen(n = 12, ratio =1,scenario = 1)
#' trt_ind <- as.numeric(idata$trtdat$trt_ind)
#' all_vars <- idata$trtdat[, -1]
#' y <- idata$Yobs
#' reference_trt <- 2
#' regadj_multiTrt_att(y = y, x = idata$trtdat,trt = trt_ind, reference = 2, ndpost = 100)
regadj_multiTrt_att = function(y, x, trt, ndpost=parent.frame()$ndpost, reference = parent.frame()$reference_trt) {

  n_trt <- length(unique(trt))
  for (i in 1:n_trt){
    assign(paste0("n",i), sum(trt==i))
  }
  n = length(trt)
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference]
  # n1 = sum(trt==1)
  # n2 = sum(trt==2)
  # n3 = sum(trt==3)
  # n = n1 + n2 + n3

  xt = cbind(trt,x)
  mod_data = cbind(y,xt)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~., data = mod_data, family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = as.data.frame(stats::model.matrix(reg_mod))

  for (i in 1:n_trt){
    assign(paste0("x_tilde", reference, i), x_tilde[x_tilde[["trt"]] == reference, ])
  }

  for (i in 1:length(trt_indicator_no_reference)){
    assign(paste0("x_tilde", reference, trt_indicator_no_reference[i]), eval(parse(text = paste0("x_tilde", reference, trt_indicator_no_reference[i]))) %>%
             dplyr::mutate(trt = trt_indicator_no_reference[i]))

  }

  # x_tilde11 = x_tilde12 = x_tilde13 = x_tilde[x_tilde$trt==1,]
  # x_tilde12$trt = 2
  # x_tilde13$trt = 3

  for (i in 1:n_trt){
    assign(paste0("y",reference, i,"_tilde"), NULL)
  }

  # predictive simulation using the stats::binomial distribution
  # predict potential outcomes
  # y11_tilde = array(NA, c(ndpost, n1))
  # y12_tilde = array(NA, c(ndpost, n1))
  # y13_tilde = array(NA, c(ndpost, n1))

  for (s in 1:ndpost) {
    for (i in 1:n_trt) {
      assign(paste0("p",reference, i,"_tilde"), arm::invlogit(as.matrix(eval(parse(text = paste0("x_tilde",reference, i)))) %*% sim_beta[s,]))
      assign(paste0("y",reference, i,"_tilde_", s),stats::rbinom(eval(parse(text = paste0("n",reference))), 1, eval(parse(text = paste0("p",reference, i,"_tilde")))))
      assign(paste0("y",reference, i,"_tilde"),rbind(eval(parse(text = paste0("y",reference, i,"_tilde"))),
                                                     eval(parse(text = paste0("y",reference, i,"_tilde_", s)))))
    }
  }


  # for (s in 1:ndpost) {
  #   # predict potential outcome Y(1)
  #   p11_tilde = arm::invlogit(as.matrix(x_tilde11) %*% sim_beta[s,])
  #   y11_tilde[s,] = stats::rbinom(n1, 1, p11_tilde)
  #
  #   # predict potential outcome Y(2)
  #   p12_tilde = arm::invlogit(as.matrix(x_tilde12) %*% sim_beta[s,])
  #   y12_tilde[s,] = stats::rbinom(n1, 1, p12_tilde)
  #
  #   # predict potential outcome Y(3)
  #   p13_tilde = arm::invlogit(as.matrix(x_tilde13) %*% sim_beta[s,])
  #   y13_tilde[s,] = stats::rbinom(n1, 1, p13_tilde)
  # }


  # Estimate causal effects
  for (j in 1:length(trt_indicator_no_reference)){
    assign(paste0("RD",reference, trt_indicator_no_reference[j], "_est"), NULL)
    assign(paste0("RR",reference, trt_indicator_no_reference[j], "_est"), NULL)
    assign(paste0("OR",reference, trt_indicator_no_reference[j], "_est"), NULL)
  }

  # RD12_est = RR12_est = OR12_est = NULL
  # RD13_est = RR13_est = OR13_est = NULL

  for (m in 1:ndpost) {
    # Estimate E(Y1), E(Y2), E(Y3)
    for (i in 1:n_trt){
      assign(paste0("y", i, "_pred_",m), mean(eval(parse(text = paste0("y",reference,i,"_tilde"))) %>% as.data.frame %>% dplyr::slice(m) %>% as.numeric()))
    }

    for (j in 1:length(trt_indicator_no_reference)){
      assign(paste0("RD",reference, trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) - eval(parse(text =(paste0("y",j, "_pred_",m)))))
      assign(paste0("RR",reference, trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) / eval(parse(text =(paste0("y",j, "_pred_",m)))))
      assign(paste0("OR",reference, trt_indicator_no_reference[j], "_est_m"), (eval(parse(text =(paste0("y",i, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",i, "_pred_",m)))))) / (eval(parse(text =(paste0("y",j, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",j, "_pred_",m)))))))
      assign(paste0("RD",reference, trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RD",reference, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RD",reference, trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("RR",reference, trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RR",reference, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference, trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("OR",reference, trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("OR",reference, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference, trt_indicator_no_reference[j], "_est_m"))))))
    }

    # print(m)
  }

  # for (m in 1:ndpost) {
  #
  #   # Estimate E(Y1|trt=1), E(Y2|trt=1), E(Y3|trt=1)
  #   y1_pred = mean(y11_tilde[m,])
  #   y2_pred = mean(y12_tilde[m,])
  #   y3_pred = mean(y13_tilde[m,])
  #
  #   # Calculate risk difference (RD)
  #   RD12_est[m] = y1_pred - y2_pred
  #   RD13_est[m] = y1_pred - y3_pred
  #
  #   # Calculate relative risk (RR)
  #   RR12_est[m] = y1_pred / y2_pred
  #   RR13_est[m] = y1_pred / y3_pred
  #
  #   # Calculate odds ratio (OR)
  #   OR12_est[m] = (y1_pred / (1 - y1_pred)) / (y2_pred / (1 - y2_pred))
  #   OR13_est[m] = (y1_pred / (1 - y1_pred)) / (y3_pred / (1 - y3_pred))
  # }

  # att12 = postSumm(RD12_est, RR12_est, OR12_est)
  # att13 = postSumm(RD13_est, RR13_est, OR13_est)
  #
  # list(ATT12 = round(att12, digits=3),
  #      ATT13 = round(att13, digits=3))
  result <- NULL

  for (j in 1:length(trt_indicator_no_reference)){
    assign(paste0("att",reference, trt_indicator_no_reference[j]), postSumm(eval(parse(text =(paste0("RD",reference, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference, trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference, trt_indicator_no_reference[j], "_est"))))))
    assign(paste0("ATT",reference, trt_indicator_no_reference[j]), list(round(eval(parse(text =(paste0("att",reference, trt_indicator_no_reference[j])))), digits = 3)))
    assign(paste0("ATT",reference, trt_indicator_no_reference[j]), stats::setNames(eval(parse(text =(paste0("ATT",reference, trt_indicator_no_reference[j])))), paste0("ATT",reference, trt_indicator_no_reference[j])))
    result <- c(result, (eval(parse(text =(paste0("ATT",reference, trt_indicator_no_reference[j]))))))
  }
  return(result)
}
