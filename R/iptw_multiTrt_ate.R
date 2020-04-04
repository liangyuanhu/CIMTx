#================================#
#Author: Liangyuan Hu, Jiayi Ji #
#================================#
iptw_multiTrt_ate<- function (y, trt_ind, psdat, wt1, wt2, wt3, method,trim_alpha,SL.library){
  psdat = idata$trtdat
  temp<-noquote(names(psdat[,2:11]))
  # strFormula  = sprintf("trt_ind~%s", paste(temp, sep = "",collapse="+"))

  # use the GBM models (twang package) to estimate the ps
  psmod<-twang::mnps(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, # formula ("strFormula")  ,
                     data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
              bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
              verbose = F, estimand = "ATE", stop.method = c("es.max"),  sampw = NULL,
              treatATT = NULL)
  #boxplot(psmod)
  wt_hat<- twang::get.weights(psmod, stop.method = "es.max",estimand = "ATE")

  ###use logistic regression model with main effects only to estimate ps
  # psmod2 <-  nnet::multinom(formula ("strFormula"), data = psdat,trace = FALSE)
  psmod2 <-  nnet::multinom(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = psdat,trace = FALSE)
  pred_ps <- stats::fitted(psmod2)
  ate_wt_1<- 1/pred_ps[,1]
  ate_wt_2<- 1/pred_ps[,2]
  ate_wt_3<- 1/pred_ps[,3]

  weightit_superlearner <- WeightIt::weightit(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = psdat,
                     method = "super", estimand = "ATE",
                     SL.library = SL.library)

  weight_superlearner <- weightit_superlearner$weights
  weight_superlearner_trim <- trunc(weightit_superlearner$weights)

  #1a_ to compute ATEs using LR estimated weights
  RD12_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[1]]
  RD13_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[2]]
  RD23_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[3]]
  RR12_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[4]]
  RR13_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[5]]
  RR23_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[6]]
  OR12_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[7]]
  OR13_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[8]]
  OR23_iptw_multinomial = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1, wt2 = ate_wt_2, wt3 = ate_wt_3)[[9]]

  #1b_ to compute ATTs using GBM estimated weights
  RD12_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[1]]
  RD13_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[2]]
  RD23_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[3]]
  RR12_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[4]]
  RR13_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[5]]
  RR23_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[6]]
  OR12_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[7]]
  OR13_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[8]]
  OR23_iptw_gbm = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat, wt2 = wt_hat, wt3 = wt_hat)[[9]]


  RD12_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[1]]
  RD13_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[2]]
  RD23_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[3]]
  RR12_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[4]]
  RR13_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[5]]
  RR23_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[6]]
  OR12_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[7]]
  OR13_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[8]]
  OR23_iptw_superlearner = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner, wt2 = weight_superlearner, wt3 = weight_superlearner)[[9]]


  #2a_ to compute ATEs using trimmed LR weights
  ate_wt_1_trunc<- trunc_fun(ate_wt_1, trim_alpha)
  ate_wt_2_trunc<- trunc_fun(ate_wt_2, trim_alpha)
  ate_wt_3_trunc<- trunc_fun(ate_wt_3, trim_alpha)
  RD12_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[1]]
  RD13_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[2]]
  RD23_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[3]]
  RR12_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[4]]
  RR13_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[5]]
  RR23_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[6]]
  OR12_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[7]]
  OR13_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[8]]
  OR23_iptw_multinomial_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = ate_wt_1_trunc, wt2 = ate_wt_2_trunc, wt3 = ate_wt_3_trunc)[[9]]

  #2b_ to compute ATEs using trimmed GBM weights
  wt_hat_trunc <- trunc(wt_hat)
  RD12_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[1]]
  RD13_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[2]]
  RD23_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[3]]
  RR12_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[4]]
  RR13_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[5]]
  RR23_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[6]]
  OR12_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[7]]
  OR13_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[8]]
  OR23_iptw_gbm_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = wt_hat_trunc, wt2 = wt_hat_trunc, wt3 = wt_hat_trunc)[[9]]

  RD12_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[1]]
  RD13_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[2]]
  RD23_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[3]]
  RR12_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[4]]
  RR13_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[5]]
  RR23_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[6]]
  OR12_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[7]]
  OR13_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[8]]
  OR23_iptw_superlearner_trim = ate_fun(y = y, trt_ind = trt_ind, wt1 = weight_superlearner_trim, wt2 = weight_superlearner_trim, wt3 = weight_superlearner_trim)[[9]]


  # ate_res = cbind(RD12_iptw_multinomial, RD13_iptw_multinomial, RD23_iptw_multinomial, RR12_iptw_multinomial, RR13_iptw_multinomial, RR23_iptw_multinomial, OR12_iptw_multinomial, OR13_iptw_multinomial, OR23_iptw_multinomial,
  #                 RD12_iptw_gbm, RD13_iptw_gbm, RD23_iptw_gbm, RR12_iptw_gbm, RR13_iptw_gbm, RR23_iptw_gbm, OR12_iptw_gbm, OR13_iptw_gbm, OR23_iptw_gbm,
  #                 RD12_iptw_multinomial_trim, RD13_iptw_multinomial_trim, RD23_iptw_multinomial_trim, RR12_iptw_multinomial_trim, RR13_iptw_multinomial_trim, RR23_iptw_multinomial_trim, OR12_iptw_multinomial_trim, OR13_iptw_multinomial_trim, OR23_iptw_multinomial_trim,
  #                 RD12_iptw_gbm_trim, RD13_iptw_gbm_trim, RD23_iptw_gbm_trim, RR12_iptw_gbm_trim, RR13_iptw_gbm_trim, RR23_iptw_gbm_trim, OR12_iptw_gbm_trim, OR13_iptw_gbm_trim, OR23_iptw_gbm_trim)

  if (method == "IPTW-Logistics") {
    res12 = rbind(RD = c(RD12_iptw_multinomial), RR = RR12_iptw_multinomial, OR = OR12_iptw_multinomial)
    res13 = rbind(RD = c(RD13_iptw_multinomial), RR = RR13_iptw_multinomial, OR = OR13_iptw_multinomial)
    res23 = rbind(RD = c(RD23_iptw_multinomial), RR = RR23_iptw_multinomial, OR = OR23_iptw_multinomial)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    colnames(res23) <- "EST"
    list(ATE12 = res12,
         ATE13 = res13,
         ATE23 = res23)
  } else if (method == "IPTW-Logistics-Trim") {
    res12 = rbind(RD = c(RD12_iptw_multinomial_trim), RR = RR12_iptw_multinomial_trim, OR = OR12_iptw_multinomial_trim)
    res13 = rbind(RD = c(RD13_iptw_multinomial_trim), RR = RR13_iptw_multinomial_trim, OR = OR13_iptw_multinomial_trim)
    res23 = rbind(RD = c(RD23_iptw_multinomial_trim), RR = RR23_iptw_multinomial_trim, OR = OR23_iptw_multinomial_trim)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    colnames(res23) <- "EST"
    list(ATE12 = res12,
         ATE13 = res13,
         ATE23 = res23)
  } else if (method == "IPTW-GBM") {
    res12 = rbind(RD = c(RD12_iptw_gbm), RR = RR12_iptw_gbm, OR = OR12_iptw_gbm)
    res13 = rbind(RD = c(RD13_iptw_gbm), RR = RR13_iptw_gbm, OR = OR13_iptw_gbm)
    res23 = rbind(RD = c(RD23_iptw_gbm), RR = RR23_iptw_gbm, OR = OR23_iptw_gbm)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    colnames(res23) <- "EST"
    list(ATE12 = res12,
         ATE13 = res13,
         ATE23 = res23)
  } else if (method == "IPTW-GBM-Trim") {
    res12 = rbind(RD = c(RD12_iptw_gbm_trim), RR = RR12_iptw_gbm_trim, OR = OR12_iptw_gbm_trim)
    res13 = rbind(RD = c(RD13_iptw_gbm_trim), RR = RR13_iptw_gbm_trim, OR = OR13_iptw_gbm_trim)
    res23 = rbind(RD = c(RD23_iptw_gbm_trim), RR = RR23_iptw_gbm_trim, OR = OR23_iptw_gbm_trim)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    colnames(res23) <- "EST"
    list(ATE12 = res12,
         ATE13 = res13,
         ATE23 = res23)
  } else if (method == "IPTW-Superlearner") {
    res12 = rbind(RD = c(RD12_iptw_superlearner), RR = RR12_iptw_superlearner, OR = OR12_iptw_superlearner)
    res13 = rbind(RD = c(RD13_iptw_superlearner), RR = RR13_iptw_superlearner, OR = OR13_iptw_superlearner)
    res23 = rbind(RD = c(RD23_iptw_superlearner), RR = RR23_iptw_superlearner, OR = OR23_iptw_superlearner)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    colnames(res23) <- "EST"
    list(ATE12 = res12,
         ATE13 = res13,
         ATE23 = res23)
  } else if (method == "IPTW-Superlearner-Trim") {
    res12 = rbind(RD = c(RD12_iptw_superlearner_trim), RR = RR12_iptw_superlearner_trim, OR = OR12_iptw_superlearner_trim)
    res13 = rbind(RD = c(RD13_iptw_superlearner_trim), RR = RR13_iptw_superlearner_trim, OR = OR13_iptw_superlearner_trim)
    res23 = rbind(RD = c(RD23_iptw_superlearner_trim), RR = RR23_iptw_superlearner_trim, OR = OR23_iptw_superlearner_trim)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    colnames(res23) <- "EST"
    list(ATE12 = res12,
         ATE13 = res13,
         ATE23 = res23)
  }
}
