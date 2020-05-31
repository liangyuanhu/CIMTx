#===============================#
#Author: Liangyuan Hu, Jiayi Ji #
#===============================#
#' Inverse probability of treatment weighting for ATT estimation (IPTW)
#'
#'This function implements the IPTW method when estimand is ATT. Please use our main function causal_multi_treat.R.
#' @param y numeric vector for the binary outcome
#' @param trt numeric vector for the treatment indicator
#' @param psdat data frame containing the treatment indicator and covariates
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
#'idata = data_gen(n = 50, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'reference_trt <- 2
#'causal_multi_treat(y = y,trt = trt_ind,
#'method = "IPTW-Logistics", estimand = "ATT", reference_trt = 2)

iptw_multiTrt_att <- function (y, trt, psdat,method, trim_alpha, reference = parent.frame()$reference_trt, SL.library  = parent.frame()$SL.library){
  # use the GBM models (twang package) to estimate the ps
  psdat = idata$trtdat
  n_trt <- length(unique(trt))
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference]

  # temp<-noquote(names(psdat[,2:11]))
  # strFormula  = sprintf("trt_ind~%s", paste(temp, sep = "",collapse="+"))
  # GBM
  # psmod<-twang::mnps(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, # formula ("strFormula")  ,
  #                    data=psdat,  n.trees = 1000, interaction.depth = 3, shrinkage = 0.01,
  #             bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
  #             verbose = F, estimand = "ATT", stop.method = c("es.max"),  sampw = NULL,
  #             treatATT = 1)
  #boxplot(psmod)
  # wt_hat<- twang::get.weights(psmod, stop.method = "es.max",estimand = "ATT")

  weightit_superlearner <- WeightIt::weightit(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = psdat,focal = reference,
                     method = "super", estimand = "ATT",
                     SL.library = SL.library)


  weight_superlearner <- weightit_superlearner$weights
  # weight_superlearner_trim <- trunc(weightit_superlearner$weights)
  wt_hat <- weight_superlearner+0.1
  ###use logistic regression model with main effects only to estimate ps
  # psmod2 <-  nnet::multinom(formula ("trt_ind ~."), data = psdat,trace = FALSE)
  psmod2 <-  nnet::multinom(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = psdat,trace = FALSE)
  pred_ps <- stats::fitted(psmod2)
  for (j in 1:length((trt_indicator_no_reference))){
    assign(paste0("att_wt_",reference, trt_indicator_no_reference[j]), pred_ps[,reference]/pred_ps[,trt_indicator_no_reference[j]])
  }


  # att_wt_12 <- pred_ps[,1]/pred_ps[,2]
  # att_wt_13 <- pred_ps[,1]/pred_ps[,3]

  if (method == "IPTW-Superlearner") {
    assign(paste0("mu_",reference,"_hat_iptw_superlearner"), mean(y[trt_ind == reference]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_superlearner"), sum(y[trt_ind == trt_indicator_no_reference[i]] * weight_superlearner[trt_ind == trt_indicator_no_reference[i]]) / sum(weight_superlearner[trt_ind == trt_indicator_no_reference[i]]))
    }
    result_list_superlearner <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner")))) / (1 - eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
      result_list_superlearner <- c(result_list_superlearner, result_once_list)
    }
    return(result_list_superlearner)
  }

  #1a_ to compute ATTs using LR estimated weights
  if (method == "IPTW-Logistics") {
    assign(paste0("mu_",reference,"_hat_iptw"), mean(y[trt_ind == reference]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw"), sum(y[trt_ind == trt_indicator_no_reference[i]] * eval(parse(text = paste0("att_wt_",reference, trt_indicator_no_reference[i])))[trt_ind == trt_indicator_no_reference[i]]) / sum(eval(parse(text = paste0("att_wt_",reference,trt_indicator_no_reference[i])))[trt_ind == trt_indicator_no_reference[i]]))
    }
    result_list_multinomial <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference, "_hat_iptw")))) / (1 - eval(parse(text =(paste0("mu_",reference, "_hat_iptw")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
      result_list_multinomial <- c(result_list_multinomial, result_once_list)
    }
    return(result_list_multinomial)
  }


  #1b_ to compute ATTs using GBM estimated weights
  if (method == "IPTW-GBM") {
    assign(paste0("mu_",reference,"_hat_iptw_gbm"), mean(y[trt_ind == reference]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_gbm"), sum(y[trt_ind == trt_indicator_no_reference[i]] * wt_hat[trt_ind == trt_indicator_no_reference[i]]) / sum(wt_hat[trt_ind == trt_indicator_no_reference[i]]))
    }
    result_list_gbm <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm")))) / (1 - eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
      result_list_gbm <- c(result_list_gbm, result_once_list)
    }
    return(result_list_gbm)
  }

  #2a_ to compute ATTs using trimmed LR weights

  if (method == "IPTW-Logistics-Trim") {
    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("att_wt_",reference,trt_indicator_no_reference[i],"_trunc"), trunc_fun(eval(parse(text = paste0("att_wt_",reference, trt_indicator_no_reference[i]))), trim_alpha))
    }
    assign(paste0("mu_",reference,"_hat_iptw_trim"), mean(y[trt_ind == reference]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_trim"), sum(y[trt_ind == trt_indicator_no_reference[i]] * eval(parse(text = paste0("att_wt_",reference, trt_indicator_no_reference[i],"_trunc")))[trt_ind == trt_indicator_no_reference[i]]) / sum(eval(parse(text = paste0("att_wt_",reference,trt_indicator_no_reference[i],"_trunc")))[trt_ind == trt_indicator_no_reference[i]]))
    }
    result_list_multinomial_trim <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_trim")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_trim")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference, "_hat_iptw_trim")))) / (1 - eval(parse(text =(paste0("mu_",reference, "_hat_iptw_trim")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_trim")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
      result_list_multinomial_trim <- c(result_list_multinomial_trim, result_once_list)
    }
    return(result_list_multinomial_trim)
  }


  #2b_ to compute ATTs using trimmed GBM weights

  # wt_hat_trunc[trt_ind==2]<- trunc_fun(wt_hat_trunc[trt_ind==2])
  # wt_hat_trunc[trt_ind==3]<- trunc_fun(wt_hat_trunc[trt_ind==3])
  if (method == "IPTW-GBM-Trim") {
    wt_hat_trunc <- wt_hat
    for (i in 1:length((trt_indicator_no_reference))){
      wt_hat_trunc[trt_ind==i] <- trunc_fun(wt_hat_trunc[trt_ind==i], trim_alpha)
    }
    assign(paste0("mu_",reference,"_hat_iptw_gbm_trim"), mean(y[trt_ind == reference]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_gbm_trim"), sum(y[trt_ind == trt_indicator_no_reference[i]] * wt_hat_trunc[trt_ind == trt_indicator_no_reference[i]]) / sum(wt_hat_trunc[trt_ind == trt_indicator_no_reference[i]]))
    }
    result_list_gbm_trim <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm_trim")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm_trim")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm_trim")))) / (1 - eval(parse(text =(paste0("mu_",reference, "_hat_iptw_gbm_trim")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_gbm_trim")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
      result_list_gbm_trim <- c(result_list_gbm_trim, result_once_list)
    }
    return(result_list_gbm_trim)
  }


  if (method == "IPTW-Superlearner-Trim") {
    weight_superlearner_trunc <- weight_superlearner
    for (i in 1:length((trt_indicator_no_reference))){
      weight_superlearner_trunc[trt_ind==i] <- trunc_fun(weight_superlearner_trunc[trt_ind==i], trim_alpha)
    }
    assign(paste0("mu_",reference,"_hat_iptw_superlearner_trim"), mean(y[trt_ind == reference]))

    for (i in 1:length((trt_indicator_no_reference))){
      assign(paste0("mu_",trt_indicator_no_reference[i], "_hat_iptw_superlearner_trim"), sum(y[trt_ind == trt_indicator_no_reference[i]] * weight_superlearner_trunc[trt_ind == trt_indicator_no_reference[i]]) / sum(weight_superlearner_trunc[trt_ind == trt_indicator_no_reference[i]]))
    }
    result_list_superlearner_trim <- NULL
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner_trim")))) - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j]), eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner_trim")))) / eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j]), (eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner_trim")))) / (1 - eval(parse(text =(paste0("mu_",reference, "_hat_iptw_superlearner_trim")))))) / (eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))) / (1 - eval(parse(text =(paste0("mu_",trt_indicator_no_reference[j], "_hat_iptw_superlearner_trim")))))))
      result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j]))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j]))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
      result_list_superlearner_trim <- c(result_list_superlearner_trim, result_once_list)
    }
    return(result_list_superlearner_trim)
  }
}

