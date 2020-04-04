#===============================#
#Author: Liangyuan Hu, Jiayi Ji #
#===============================#
iptw_multiTrt_att <- function (y, trt, psdat, wt12, wt13,method, trim_alpha, SL.library){
  # use the GBM models (twang package) to estimate the ps
  psdat = idata$trtdat
  temp<-noquote(names(psdat[,2:11]))
  # strFormula  = sprintf("trt_ind~%s", paste(temp, sep = "",collapse="+"))
  psmod<-twang::mnps(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, # formula ("strFormula")  ,
                     data=psdat,  n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
              bag.fraction = 1.0, perm.test.iters=0,  print.level = 2,  iterlim = 1000,
              verbose = F, estimand = "ATT", stop.method = c("es.max"),  sampw = NULL,
              treatATT = 1)
  #boxplot(psmod)
  wt_hat<- twang::get.weights(psmod, stop.method = "es.max",estimand = "ATT")

  weightit_superlearner <- WeightIt::weightit(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = psdat,focal = 1,
                     method = "super", estimand = "ATT",
                     SL.library = SL.library)


  weight_superlearner <- weightit_superlearner$weights
  weight_superlearner_trim <- trunc(weightit_superlearner$weights)

  ###use logistic regression model with main effects only to estimate ps
  # psmod2 <-  nnet::multinom(formula ("trt_ind ~."), data = psdat,trace = FALSE)
  psmod2 <-  nnet::multinom(trt_ind~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = psdat,trace = FALSE)

  pred_ps <- stats::fitted(psmod2)
  att_wt_12 <- pred_ps[,1]/pred_ps[,2]
  att_wt_13 <- pred_ps[,1]/pred_ps[,3]

  RD12_iptw_superlearner = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner, wt13 = weight_superlearner)[[1]]
  RD13_iptw_superlearner = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner, wt13 = weight_superlearner)[[2]]
  RR12_iptw_superlearner = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner, wt13 = weight_superlearner)[[3]]
  RR13_iptw_superlearner = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner, wt13 = weight_superlearner)[[4]]
  OR12_iptw_superlearner = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner, wt13 = weight_superlearner)[[5]]
  OR13_iptw_superlearner = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner, wt13 = weight_superlearner)[[6]]

  #1a_ to compute ATTs using LR estimated weights
  RD12_iptw_multinomial = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12, wt13 = att_wt_13)[[1]]
  RD13_iptw_multinomial = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12, wt13 = att_wt_13)[[2]]
  RR12_iptw_multinomial = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12, wt13 = att_wt_13)[[3]]
  RR13_iptw_multinomial = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12, wt13 = att_wt_13)[[4]]
  OR12_iptw_multinomial = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12, wt13 = att_wt_13)[[5]]
  OR13_iptw_multinomial = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12, wt13 = att_wt_13)[[6]]

  #1b_ to compute ATTs using GBM estimated weights
  RD12_iptw_gbm = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat, wt13 = wt_hat)[[1]]
  RD13_iptw_gbm = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat, wt13 = wt_hat)[[2]]
  RR12_iptw_gbm = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat, wt13 = wt_hat)[[3]]
  RR13_iptw_gbm = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat, wt13 = wt_hat)[[4]]
  OR12_iptw_gbm = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat, wt13 = wt_hat)[[5]]
  OR13_iptw_gbm = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat, wt13 = wt_hat)[[6]]

  #2a_ to compute ATTs using trimmed LR weights
  att_wt_12_trunc<- trunc_fun(att_wt_12, trim_alpha)
  att_wt_13_trunc<- trunc_fun(att_wt_13, trim_alpha)
  RD12_iptw_multinomial_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12_trunc, wt13 = att_wt_13_trunc)[[1]]
  RD13_iptw_multinomial_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12_trunc, wt13 = att_wt_13_trunc)[[2]]
  RR12_iptw_multinomial_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12_trunc, wt13 = att_wt_13_trunc)[[3]]
  RR13_iptw_multinomial_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12_trunc, wt13 = att_wt_13_trunc)[[4]]
  OR12_iptw_multinomial_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12_trunc, wt13 = att_wt_13_trunc)[[5]]
  OR13_iptw_multinomial_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = att_wt_12_trunc, wt13 = att_wt_13_trunc)[[6]]

  #2b_ to compute ATTs using trimmed GBM weights
  wt_hat_trunc <- wt_hat
  wt_hat_trunc[trt_ind==2]<- trunc_fun(wt_hat_trunc[trt_ind==2])
  wt_hat_trunc[trt_ind==3]<- trunc_fun(wt_hat_trunc[trt_ind==3])
  RD12_iptw_gbm_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat_trunc, wt13 = wt_hat_trunc)[[1]]
  RD13_iptw_gbm_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat_trunc, wt13 = wt_hat_trunc)[[2]]
  RR12_iptw_gbm_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat_trunc, wt13 = wt_hat_trunc)[[3]]
  RR13_iptw_gbm_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat_trunc, wt13 = wt_hat_trunc)[[4]]
  OR12_iptw_gbm_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat_trunc, wt13 = wt_hat_trunc)[[5]]
  OR13_iptw_gbm_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = wt_hat_trunc, wt13 = wt_hat_trunc)[[6]]

  RD12_iptw_superlearner_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner_trim, wt13 = weight_superlearner_trim)[[1]]
  RD13_iptw_superlearner_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner_trim, wt13 = weight_superlearner_trim)[[2]]
  RR12_iptw_superlearner_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner_trim, wt13 = weight_superlearner_trim)[[3]]
  RR13_iptw_superlearner_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner_trim, wt13 = weight_superlearner_trim)[[4]]
  OR12_iptw_superlearner_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner_trim, wt13 = weight_superlearner_trim)[[5]]
  OR13_iptw_superlearner_trim = att_fun(y = y, trt_ind = trt_ind, wt12 = weight_superlearner_trim, wt13 = weight_superlearner_trim)[[6]]

  # att_res = cbind(RD12_iptw_multinomial, RD13_iptw_multinomial, RR12_iptw_multinomial, RR13_iptw_multinomial, OR12_iptw_multinomial, OR13_iptw_multinomial,
  #                 RD12_iptw_gbm, RD13_iptw_gbm, RR12_iptw_gbm, RR13_iptw_gbm, OR12_iptw_gbm, OR13_iptw_gbm,
  #                 RD12_iptw_multinomial_trim, RD13_iptw_multinomial_trim, RR12_iptw_multinomial_trim, RR13_iptw_multinomial_trim, OR12_iptw_multinomial_trim, OR13_iptw_multinomial_trim,
  #                 RD12_iptw_gbm_trim, RD13_iptw_gbm_trim, RR12_iptw_gbm_trim, RR13_iptw_gbm_trim, OR12_iptw_gbm_trim, OR13_iptw_gbm_trim)
  # return(att_res)

  if (method == "IPTW-Logistics") {
    res12 = rbind(RD = c(RD12_iptw_multinomial), RR = RR12_iptw_multinomial, OR = OR12_iptw_multinomial)
    res13 = rbind(RD = c(RD13_iptw_multinomial), RR = RR13_iptw_multinomial, OR = OR13_iptw_multinomial)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    list(ATT12 = res12,
         ATT13 = res13)
  } else if (method == "IPTW-Logistics-Trim") {
    res12 = rbind(RD = c(RD12_iptw_multinomial_trim), RR = RR12_iptw_multinomial_trim, OR = OR12_iptw_multinomial_trim)
    res13 = rbind(RD = c(RD13_iptw_multinomial_trim), RR = RR13_iptw_multinomial_trim, OR = OR13_iptw_multinomial_trim)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    list(ATT12 = res12,
         ATT13 = res13)
  } else if (method == "IPTW-GBM") {
    res12 = rbind(RD = c(RD12_iptw_gbm), RR = RR12_iptw_gbm, OR = OR12_iptw_gbm)
    res13 = rbind(RD = c(RD13_iptw_gbm), RR = RR13_iptw_gbm, OR = OR13_iptw_gbm)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    list(ATT12 = res12,
         ATT13 = res13)
  } else if (method == "IPTW-GBM-Trim") {
    res12 = rbind(RD = c(RD12_iptw_gbm_trim), RR = RR12_iptw_gbm_trim, OR = OR12_iptw_gbm_trim)
    res13 = rbind(RD = c(RD13_iptw_gbm_trim), RR = RR13_iptw_gbm_trim, OR = OR13_iptw_gbm_trim)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    list(ATT12 = res12,
         ATT13 = res13)
  } else if (method == "IPTW-Superlearner") {
    res12 = rbind(RD = c(RD12_iptw_multinomial), RR = RR12_iptw_multinomial, OR = OR12_iptw_multinomial)
    res13 = rbind(RD = c(RD13_iptw_multinomial), RR = RR13_iptw_multinomial, OR = OR13_iptw_multinomial)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    list(ATT12 = res12,
         ATT13 = res13)
  } else if (method == "IPTW-Superlearner-Trim") {
    res12 = rbind(RD = c(RD12_iptw_multinomial_trim), RR = RR12_iptw_multinomial_trim, OR = OR12_iptw_multinomial_trim)
    res13 = rbind(RD = c(RD13_iptw_multinomial_trim), RR = RR13_iptw_multinomial_trim, OR = OR13_iptw_multinomial_trim)
    colnames(res12) <- "EST"
    colnames(res13) <- "EST"
    list(ATT12 = res12,
         ATT13 = res13)
  # res12_iptw_multinomial = rbind(RD = c(RD12_iptw_multinomial), RR = RR12_iptw_multinomial, OR = OR12_iptw_multinomial)
  # res13_iptw_multinomial = rbind(RD = c(RD13_iptw_multinomial), RR = RR13_iptw_multinomial, OR = OR13_iptw_multinomial)
  # colnames(res12_iptw_multinomial) <- "EST"
  # colnames(res13_iptw_multinomial) <- "EST"
  #
  # res12_iptw_multinomial_trim = rbind(RD = c(RD12_iptw_multinomial_trim), RR = RR12_iptw_multinomial_trim, OR = OR12_iptw_multinomial_trim)
  # res13_iptw_multinomial_trim = rbind(RD = c(RD13_iptw_multinomial_trim), RR = RR13_iptw_multinomial_trim, OR = OR13_iptw_multinomial_trim)
  # colnames(res12_iptw_multinomial_trim) <- "EST"
  # colnames(res13_iptw_multinomial_trim) <- "EST"
  #
  # res12_iptw_gbm = rbind(RD = c(RD12_iptw_gbm), RR = RR12_iptw_gbm, OR = OR12_iptw_gbm)
  # res13_iptw_gbm = rbind(RD = c(RD13_iptw_gbm), RR = RR13_iptw_gbm, OR = OR13_iptw_gbm)
  # colnames(res12_iptw_gbm) <- "EST"
  # colnames(res13_iptw_gbm) <- "EST"
  #
  # res12_iptw_gbm_trim = rbind(RD = c(RD12_iptw_gbm_trim), RR = RR12_iptw_gbm_trim, OR = OR12_iptw_gbm_trim)
  # res13_iptw_gbm_trim = rbind(RD = c(RD13_iptw_gbm_trim), RR = RR13_iptw_gbm_trim, OR = OR13_iptw_gbm_trim)
  # colnames(res12_iptw_gbm_trim) <- "EST"
  # colnames(res13_iptw_gbm_trim) <- "EST"
  #
  # res12_iptw_superlearner = rbind(RD = c(RD12_iptw_superlearner), RR = RR12_iptw_superlearner, OR = OR12_iptw_superlearner)
  # res13_iptw_superlearner = rbind(RD = c(RD13_iptw_superlearner), RR = RR13_iptw_superlearner, OR = OR13_iptw_superlearner)
  # colnames(res12_iptw_superlearner) <- "EST"
  # colnames(res13_iptw_superlearner) <- "EST"
  #
  # res12_iptw_superlearner_trim = rbind(RD = c(RD12_iptw_superlearner_trim), RR = RR12_iptw_superlearner_trim, OR = OR12_iptw_superlearner_trim)
  # res13_iptw_superlearner_trim = rbind(RD = c(RD13_iptw_superlearner_trim), RR = RR13_iptw_superlearner_trim, OR = OR13_iptw_superlearner_trim)
  # colnames(res12_iptw_superlearner_trim) <- "EST"
  # colnames(res13_iptw_superlearner_trim) <- "EST"
  #
  # list(ATT12_iptw_multinomial = res12_iptw_multinomial,
  #      ATT13_iptw_multinomial = res13_iptw_multinomial,
  #      ATT12_iptw_multinomial_trim = res12_iptw_multinomial_trim,
  #      ATT13_iptw_multinomial_trim = res13_iptw_multinomial_trim,
  #      ATT12_iptw_gbm  = res12_iptw_gbm,
  #      ATT13_iptw_gbm  = res13_iptw_gbm,
  #      ATT12_iptw_gbm_trim  = res12_iptw_gbm_trim,
  #      ATT13_iptw_gbm_trim  = res13_iptw_gbm_trim,
  #      ATT12_iptw_superlearner = res12_iptw_superlearner,
  #      ATT13_iptw_superlearner = res13_iptw_superlearner,
  #      ATT12_iptw_superlearner_trim = res12_iptw_superlearner_trim,
  #      ATT13_iptw_superlearner_trim = res13_iptw_superlearner_trim)


}
}
