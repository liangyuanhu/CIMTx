#================================#
#Author: Liangyuan Hu, Jiayi Ji #
#================================#
#' Inverse probability of treatment weighting for ATE estimation (IPTW)
#'
#'This function implements the IPTW method when estimand is ATE. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param trt_ind numeric vector for the treatment indicator
#' @param psdat data frame containing the treatment indicator and covariates
#' @param method methods for causal inference with multiple treatments, inherited from causal_multi_treat.R
#' @param trim_alpha alpha values for IPTW weight trimming, inherited from causal_multi_treat.R
#' @param SL.library methods specified with SL.library in Superlearner package, inherited from causal_multi_treat.R
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
#'idata = data_gen(n = 50, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'iptw_multiTrt_ate(y=y, trt = trt_ind,SL.library = c("SL.glm"),
#'trim_alpha = 0.05, method = "IPTW-Logistics-Trim")
#'causal_multi_treat(y = y,trt = trt_ind,
#'method = "IPTW-Logistics", estimand = "ATE")
#'
iptw_multiTrt_ate<- function (y, trt_ind, psdat, method,trim_alpha,SL.library){
  psdat = idata$trtdat
  # temp<-noquote(names(psdat[,2:11]))
  # strFormula  = sprintf("trt_ind~%s", paste(temp, sep = "",collapse="+"))

  # use the GBM models (twang package) to estimate the ps
  # psmod<-twang::mnps(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, # formula ("strFormula")  ,
  #                    data=psdat,  n.trees = 1000, interaction.depth = 3, shrinkage = 0.01,
  #             bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
  #             verbose = F, estimand = "ATE", stop.method = c("es.max"),  sampw = NULL,
  #             treatATT = NULL)
  # #boxplot(psmod)
  # wt_hat<- twang::get.weights(psmod, stop.method = "es.max",estimand = "ATE")

  ###use logistic regression model with main effects only to estimate ps
  # psmod2 <-  nnet::multinom(formula ("strFormula"), data = psdat,trace = FALSE)
  n_trt <- length(unique(trt_ind))
  psmod2 <-  nnet::multinom(trt_ind~., data = psdat,trace = FALSE)
  pred_ps <- stats::fitted(psmod2)

  for (i in 1:n_trt){
    assign(paste0("ate_wt_",i), 1/pred_ps[,i])
  }

  weightit_superlearner <- WeightIt::weightit(trt_ind~., data = psdat,
                     method = "super", estimand = "ATE",
                     SL.library = SL.library)

  weight_superlearner <- weightit_superlearner$weights
  weight_superlearner_trim <- trunc(weightit_superlearner$weights)
  wt_hat <- weight_superlearner+0.1
  #1a_ to compute ATEs using LR estimated weights
  if (method == "IPTW-Logistics") {
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_iptw"), sum(y[trt_ind == i] * eval(parse(text = paste0("ate_wt_",i)))[trt_ind == i]) / sum(eval(parse(text = paste0("ate_wt_",i)))[trt_ind == i]))
    }
    result_list_multinomial <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw"))) - eval(parse(text = paste0("mu_",j, "_hat_iptw"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw"))) / eval(parse(text = paste0("mu_",j, "_hat_iptw"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_iptw"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_iptw"))))) / (eval(parse(text = paste0("mu_",j, "_hat_iptw"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_iptw"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_multinomial <- c(result_list_multinomial, result_once_list)
      }
    }
    return(result_list_multinomial)
  } else if (method == "IPTW-Logistics-Trim") {
    for (i in 1:n_trt){
      assign(paste0("ate_wt_",i,"_trunc"), trunc_fun(eval(parse(text = paste0("ate_wt_",i))), trim_alpha))
    }
    # ate_wt_1_trunc<- trunc_fun(ate_wt_1, trim_alpha)
    # ate_wt_2_trunc<- trunc_fun(ate_wt_2, trim_alpha)
    # ate_wt_3_trunc<- trunc_fun(ate_wt_3, trim_alpha)
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_iptw_trim"), sum(y[trt_ind == i] * eval(parse(text = paste0("ate_wt_",i,"_trunc")))[trt_ind == i]) / sum(eval(parse(text = paste0("ate_wt_",i,"_trunc")))[trt_ind == i]))
    }
    result_list_multinomial_trim <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))) - eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))) / eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_iptw_trim"))))) / (eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_iptw_trim"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_multinomial_trim <- c(result_list_multinomial_trim, result_once_list)
      }
    }
    return(result_list_multinomial_trim)

  } else if (method == "IPTW-GBM") {
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hatgbm"), sum(y[trt_ind == i] * wt_hat[trt_ind == i]) / sum(wt_hat[trt_ind == i]))
    }
    result_list_gbm <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm"))) - eval(parse(text = paste0("mu_",j, "_hatgbm"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm"))) / eval(parse(text = paste0("mu_",j, "_hatgbm"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hatgbm"))) /(1 - eval(parse(text = paste0("mu_",i, "_hatgbm"))))) / (eval(parse(text = paste0("mu_",j, "_hatgbm"))) /(1 - eval(parse(text = paste0("mu_",j, "_hatgbm"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_gbm <- c(result_list_gbm, result_once_list)
      }
    }
    return(result_list_gbm)
  } else if (method == "IPTW-GBM-Trim") {
    wt_hat_trunc <- trunc(wt_hat)
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hatgbm_trim"), sum(y[trt_ind == i] * wt_hat_trunc[trt_ind == i]) / sum(wt_hat_trunc[trt_ind == i]))
    }
    result_list_gbm_trim <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))) - eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))) / eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))) /(1 - eval(parse(text = paste0("mu_",i, "_hatgbm_trim"))))) / (eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))) /(1 - eval(parse(text = paste0("mu_",j, "_hatgbm_trim"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_gbm_trim <- c(result_list_gbm_trim, result_once_list)
      }
    }
    return(result_list_gbm_trim)
  } else if (method == "IPTW-Superlearner") {
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_superlearner"), sum(y[trt_ind == i] * weight_superlearner[trt_ind == i]) / sum(weight_superlearner[trt_ind == i]))
    }
    result_list_superlearner <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner"))) - eval(parse(text = paste0("mu_",j, "_hat_superlearner"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner"))) / eval(parse(text = paste0("mu_",j, "_hat_superlearner"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_superlearner"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_superlearner"))))) / (eval(parse(text = paste0("mu_",j, "_hat_superlearner"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_superlearner"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_superlearner <- c(result_list_superlearner, result_once_list)
      }
    }
    return(result_list_superlearner)
  } else if (method == "IPTW-Superlearner-Trim") {
    for (i in 1:n_trt){
      assign(paste0("mu_",i, "_hat_superlearner_trim"), sum(y[trt_ind == i] * weight_superlearner_trim[trt_ind == i]) / sum(weight_superlearner_trim[trt_ind == i]))
    }
    result_list_superlearner_trim <- NULL
    for (i in 1:(n_trt-1)){
      result_once <- NULL
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))) - eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))))
        assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))) / eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))))
        assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat_superlearner_trim"))))) / (eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat_superlearner_trim"))))))
        result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
        colnames(result_once) <- "EST"
        rownames(result_once) <- c("RD", "RR", "OR")
        result_once_list <- list(result_once)
        names(result_once_list) <- paste0("ATE",i,j)
        result_list_superlearner_trim <- c(result_list_superlearner_trim, result_once_list)
      }
    }
    return(result_list_superlearner_trim)

  }
}
