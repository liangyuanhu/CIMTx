# Estimate average treatment effect (ATE)
#' Regression Adjustment when estimand is ATE
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
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
#' @export
#' @examples
#' library(CIMTx)
#' set.seed(3242019)
#' idata = data_gen(n = 12, ratio =1,scenario = 1)
#' trt_ind <- as.numeric(idata$trtdat$trt_ind)
#' all_vars <- idata$trtdat[, -1]
#' y <- idata$Yobs
#'causal_multi_treat(y = y, x = idata$trtdat,ndpost = 10,
#'trt = trt_ind, method ="Regression Adjustment", estimand = "ATE")
regadj_multiTrt_ate = function(y, x, trt, ndpost=parent.frame()$ndpost) {
  n_trt <- length(unique(trt))
  for (i in 1:n_trt){
    assign(paste0("n",i), sum(trt==i))
  }
  n = length(trt)

  xt = cbind(trt,x)
  mod_data = cbind(y,xt)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~ ., data = mod_data, family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = stats::model.matrix(reg_mod)

  for (i in 1:n_trt){
    assign(paste0("x_tilde",i), as.data.frame(x_tilde) %>%
             dplyr::mutate(trt = i))
  }

  # x_tilde1 = x_tilde2 = x_tilde3 = as.data.frame(x_tilde)
  # x_tilde1$trt = 1
  # x_tilde2$trt = 2
  # x_tilde3$trt = 3

  # predictive simulation using the stats::binomial distribution
  # predict potential outcomes
  for (i in 1:n_trt){
    assign(paste0("y",i,"_tilde"), NULL)
  }

  # y1_tilde = array(NA, c(ndpost, n))
  # y2_tilde = array(NA, c(ndpost, n))
  # y3_tilde = array(NA, c(ndpost, n))

  for (s in 1:ndpost) {
    for (i in 1:n_trt) {
      assign(paste0("p",i,"_tilde"), arm::invlogit(as.matrix(eval(parse(text = paste0("x_tilde",i)))) %*% sim_beta[s,]))
      assign(paste0("y",i,"_tilde_", s),stats::rbinom(n, 1, eval(parse(text = paste0("p",i,"_tilde")))))
      assign(paste0("y",i,"_tilde"),rbind(eval(parse(text = paste0("y",i,"_tilde"))),
                                          eval(parse(text = paste0("y",i,"_tilde_", s)))))
    }
  }
  #   # predict potential outcome Y(1)
  #   p1_tilde = arm::invlogit(as.matrix(x_tilde1) %*% sim_beta[s,])
  #   y1_tilde[s,] = stats::rbinom(n, 1, p1_tilde)
  #
  #   # predict potential outcome Y(2)
  #   p2_tilde = arm::invlogit(as.matrix(x_tilde2) %*% sim_beta[s,])
  #   y2_tilde[s,] = stats::rbinom(n, 1, p2_tilde)
  #
  #   # predict potential outcome Y(3)
  #   p3_tilde = arm::invlogit(as.matrix(x_tilde3) %*% sim_beta[s,])
  #   y3_tilde[s,] = stats::rbinom(n, 1, p3_tilde)
  # }


  # Estimate causal effects
  for (i in 1:(n_trt-1)){
    for (j in (i+1):(n_trt)){
      assign(paste0("RD",i,j, "_est"), NULL)
      assign(paste0("RR",i,j, "_est"), NULL)
      assign(paste0("OR",i,j, "_est"), NULL)
    }
  }

  # RD12_est = RR12_est = OR12_est = NULL
  # RD13_est = RR13_est = OR13_est = NULL
  # RD23_est = RR23_est = OR23_est = NULL

  for (m in 1:ndpost) {
    # Estimate E(Y1), E(Y2), E(Y3)
    for (i in 1:n_trt){
      assign(paste0("y",i, "_pred_",m), mean(eval(parse(text = paste0("y",i,"_tilde"))) %>% as.data.frame %>% dplyr::slice(m) %>% as.numeric()))
    }

    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j, "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) - eval(parse(text =(paste0("y",j, "_pred_",m)))))
        assign(paste0("RR",i,j, "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) / eval(parse(text =(paste0("y",j, "_pred_",m)))))
        assign(paste0("OR",i,j, "_est_m"), (eval(parse(text =(paste0("y",i, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",i, "_pred_",m)))))) / (eval(parse(text =(paste0("y",j, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",j, "_pred_",m)))))))
        assign(paste0("RD",i,j, "_est"), c(eval(parse(text =(paste0("RD",i,j, "_est")))), eval(parse(text =(paste0("RD",i,j, "_est_m"))))))

        assign(paste0("RR",i,j, "_est"), c(eval(parse(text =(paste0("RR",i,j, "_est")))), eval(parse(text =(paste0("RR",i,j, "_est_m"))))))

        assign(paste0("OR",i,j, "_est"), c(eval(parse(text =(paste0("OR",i,j, "_est")))), eval(parse(text =(paste0("OR",i,j, "_est_m"))))))
      }
    }
  # print(m)
  }

  result <- NULL
  for (i in 1:(n_trt-1)){
    for (j in (i + 1):n_trt){
      assign(paste0("ate",i,j), postSumm(eval(parse(text =(paste0("RD",i,j, "_est")))), eval(parse(text =(paste0("RR",i,j, "_est")))), eval(parse(text =(paste0("OR",i,j, "_est"))))))
      assign(paste0("ATE",i,j), list(round(eval(parse(text =(paste0("ate",i,j)))), digits = 3)))
      assign(paste0("ATE",i,j), stats::setNames(eval(parse(text =(paste0("ATE",i,j)))), paste0("ATE",i,j)))
      result <- c(result, (eval(parse(text =(paste0("ATE",i,j))))))
    }
  }
  return(result)

    # y1_pred = mean(y1_tilde[m,])
    # y2_pred = mean(y2_tilde[m,])
    # y3_pred = mean(y3_tilde[m,])
    #
    # # Calculate risk difference (RD)
    # RD12_est[m] = y1_pred - y2_pred
    # RD13_est[m] = y1_pred - y3_pred
    # RD23_est[m] = y2_pred - y3_pred
    #
    # # Calculate relative risk (RR)
    # RR12_est[m] = y1_pred / y2_pred
    # RR13_est[m] = y1_pred / y3_pred
    # RR23_est[m] = y2_pred / y3_pred
    #
    # # Calculate  odds ratio (OR)
    # OR12_est[m] = (y1_pred / (1 - y1_pred)) / (y2_pred / (1 - y2_pred))
    # OR13_est[m] = (y1_pred / (1 - y1_pred)) / (y3_pred / (1 - y3_pred))
    # OR23_est[m] = (y2_pred / (1 - y2_pred)) / (y3_pred / (1 - y3_pred))
  # }

  # ate12 = postSumm(RD12_est, RR12_est, OR12_est)
  # ate13 = postSumm(RD13_est, RR13_est, OR13_est)
  # ate23 = postSumm(RD23_est, RR23_est, OR23_est)
  #
  # list(ATE12 = round(ate12, digits=3),
  #      ATE13 = round(ate13, digits=3),
  #      ATE23 = round(ate23, digits=3))
}
