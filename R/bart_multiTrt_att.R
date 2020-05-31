# Estimate average treatment effect on the treated (ATT)
# The default reference group is 1st group
# Author: Liangyuan Hu, Chenyang Gu
#' Bayesian Additive Regression Trees (BART) for ATT estimation
#'
#' This function implements the BART method when estimand is ATT. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
#' @param k For binary y, k is the number of prior standard deviations f(x) is away from +/-3. The bigger k is, the more conservative the fitting will be.
#' @param discard discarding rules for BART method, inherited from causal_multi_treat.R
#' @param ntree The number of trees in the sum
#' @param ndpost The number of posterior draws returned
#' @param nskip Number of MCMC iterations to be treated as burn in
#' @param reference Reference group for ATT
#'
#' @return list with 2 elements for ATT effect. It contains
#' \item{ATT12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATT13:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' list with 3 elements for ATT effect. It contains
#' \item{ATE12:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#' \item{ATE13:}{A dataframe containing the estimation,
#' standard error, lower and upper 95\% CI for RD/RR/OR}
#'\item{ATE23:}{A dataframe containing the estimation,
#'standard error, lower and upper 95\% CI for RD/RR/OR}
#' @export
#' @examples
#'library(CIMTx)
#'set.seed(3242019)
#'idata = data_gen(n = 5, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATT", discard = "No", ndpost = 10, reference_trt = 2)
bart_multiTrt_att = function(y, x, trt, k=2, discard = "No", ntree=100, ndpost=1000, nskip=1000, reference = parent.frame()$reference_trt) {
  n_trt <- length(unique(trt))
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference]

  for (i in 1:n_trt){
    assign(paste0("n",i), sum(trt==i))
  }
  # n1 = sum(trt==1)
  # n2 = sum(trt==2)
  # n3 = sum(trt==3)

  xt = cbind(trt,x)

  # Fit BART
  bart_mod = BART::pbart(x.train = xt, y.train = y, k = k, ntree = ntree, ndpost = ndpost, nskip = nskip)



  # Predict potential outcomes for trt=1
  # xp1 = xt[trt==1,]
  assign(paste0("xp",reference), xt[trt==reference,])
  for (i in trt_indicator[trt_indicator!=reference]){
    assign(paste0("xp",i), xt[trt==i,])
    assign(paste0("xp",i), eval(parse(text = paste0("xp",i))) %>%
             dplyr::mutate(trt = i))
  }

  for (j in 1:(n_trt)){
      assign(paste0("bart_pred",reference,j), BART::pwbart(eval(parse(text = paste0("xp",j))), bart_mod$treedraws, mu=mean(y)))
      assign(paste0("pred_prop",reference,j), stats::pnorm(eval(parse(text = paste0("bart_pred",reference,j)))))
    }


  # xp2 = xp1
  # xp3 = xp1
  # xp2[,1] = 2  # switch treatment label 1 to 2
  # xp3[,1] = 3  # switch treatment label 1 to 3
  #
  # bart_pred11 = BART::pwbart(xp1, bart_mod$treedraws)
  # bart_pred12 = BART::pwbart(xp2, bart_mod$treedraws)
  # bart_pred13 = BART::pwbart(xp3, bart_mod$treedraws)
  #
  # pred_prop11 = stats::pnorm(bart_pred11)
  # pred_prop12 = stats::pnorm(bart_pred12)
  # pred_prop13 = stats::pnorm(bart_pred13)

  if (discard == "No") {
    # pred_prop11 = stats::pnorm(bart_pred11)
    # pred_prop12 = stats::pnorm(bart_pred12)
    # pred_prop13 = stats::pnorm(bart_pred13)
    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference,j), stats::pnorm(eval(parse(text = paste0("bart_pred",reference,j)))))
    }
  } else if (discard == "Stringent"){

    for (j in 1:(n_trt)){
      assign(paste0("post.ind.sd",j), apply(eval(parse(text = paste0("pred_prop",reference,j))), 2, stats::sd))
    }

    # #posterior standard deviation of the predicted outcome among those treated with W=1
    # post.ind.sd1 = apply(pred_prop11, 2, stats::sd)
    #
    # #discard unit i with W_i =1 if posterior sd of his/her counterfactual outcomes exceeeds threshold
    # post.ind.sd2 = apply(pred_prop12, 2, stats::sd)
    # post.ind.sd3 = apply(pred_prop13, 2, stats::sd)
    eligible <- TRUE

    for (j in 1:length((trt_indicator_no_reference))) {
      assign(paste0("eligible", trt_indicator_no_reference[j]), eval(parse(text = paste0("post.ind.sd",trt_indicator_no_reference[j]))) <= 2 * eval(parse(text = paste0("post.ind.sd",reference))))
    }
    for (j in 1:length((trt_indicator_no_reference))){
      assign(eligible, eval(parse(text = paste0("eligible",i,trt_indicator_no_reference[j]))) &eligible)
      assign(paste0("n_",i, "_discard_att"), sum(eligible == F))
    }
    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference,j), eval(parse(text = paste0("pred_prop",reference,j))) %>%
               as.data.frame() %>%
               dplyr::select(which(eligible)) %>%
               as.matrix())
    }
  # }

    # #discard unit i if sd of counterfactual outcomes is more than 2 times sd of fitted outcomes
    # eligible = (post.ind.sd2 <= 2*post.ind.sd1) & (post.ind.sd3 <= 2*post.ind.sd1)
    # n_discard_att = sum(eligible == F)
    # #n_discard_att
    # #these units are within the common support region
    # pred_prop11 = pred_prop11[,eligible]
    # pred_prop12 = pred_prop12[,eligible]
    # pred_prop13 = pred_prop13[,eligible]
  } else if (discard == "Lenient"){
    for (j in 1:(n_trt)){
      assign(paste0("post.ind.sd",j), apply(eval(parse(text = paste0("pred_prop",reference,j))), 2, stats::sd))
    }
    threshold <- max(eval(parse(text = paste0("post.ind.sd",reference))))

    for (j in 1:length((trt_indicator_no_reference))) {
      assign(paste0("eligible", trt_indicator_no_reference[j]), eval(parse(text = paste0("post.ind.sd",trt_indicator_no_reference[j]))) <= threshold)
    }

    for (j in 1:length((trt_indicator_no_reference))){
      assign(eligible, eval(parse(text = paste0("eligible",i,trt_indicator_no_reference[j]))) &eligible)
      assign(paste0("n_",i, "_discard_att"), sum(eligible == F))
    }

    for (j in 1:(n_trt)){
      assign(paste0("pred_prop",reference,j), eval(parse(text = paste0("pred_prop",reference,j))) %>%
               as.data.frame() %>%
               dplyr::select(which(eligible)) %>%
               as.matrix())
    }

    # post.ind.sd1 = apply(pred_prop11, 2, stats::sd)
    # threshold = max(post.ind.sd1)
    # #discard unit i with W_i =1 if posterior sd of his/her counterfactual outcomes exceeeds threshold
    # post.ind.sd2 = apply(pred_prop12, 2, stats::sd)
    # post.ind.sd3 = apply(pred_prop13, 2, stats::sd)

    # eligible = (post.ind.sd2 <= threshold) & (post.ind.sd3 <= threshold)
    #discard unit i if sd of counterfactual outcomes is more than 2 times sd of fitted outcomes
    # n_discard_att = sum(eligible == F)
    #n_discard_att
    #these units are within the common support region
    # pred_prop11 = pred_prop11[,eligible]
    # pred_prop12 = pred_prop12[,eligible]
    # pred_prop13 = pred_prop13[,eligible]

  }
  # Estimate causal effects
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j], "_est"), NULL)
      assign(paste0("RR",reference,trt_indicator_no_reference[j], "_est"), NULL)
      assign(paste0("OR",reference,trt_indicator_no_reference[j], "_est"), NULL)
    }

  # RD12_est = RR12_est = OR12_est = NULL
  # RD13_est = RR13_est = OR13_est = NULL
  for (m in 1:ndpost) {
      for (j in 1:n_trt){
        assign(paste0("y",j,"_pred"), mean(stats::rbinom(eval(parse(text =(paste0("n",reference)))), 1, eval(parse(text =(paste0("pred_prop",reference,j)))))))
      }
    for (j in 1:length((trt_indicator_no_reference))){
      assign(paste0("RD",reference,trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",reference, "_pred")))) - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))))
      assign(paste0("RR",reference,trt_indicator_no_reference[j], "_est_m"), eval(parse(text =(paste0("y",reference, "_pred")))) / eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))))
      assign(paste0("OR",reference,trt_indicator_no_reference[j], "_est_m"), (eval(parse(text =(paste0("y",reference, "_pred")))) / (1 - eval(parse(text =(paste0("y",reference, "_pred")))))) / (eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))) / (1 - eval(parse(text =(paste0("y",trt_indicator_no_reference[j], "_pred")))))))
      assign(paste0("RD",reference,trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RD",reference,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RD",reference,trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("RR",reference,trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("RR",reference,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference,trt_indicator_no_reference[j], "_est_m"))))))

      assign(paste0("OR",reference,trt_indicator_no_reference[j], "_est"), c(eval(parse(text =(paste0("OR",reference,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference,trt_indicator_no_reference[j], "_est_m"))))))
    }

  }
  result <- NULL
  for (j in 1:length((trt_indicator_no_reference))){
    assign(paste0("att",reference,trt_indicator_no_reference[j]), postSumm(eval(parse(text =(paste0("RD",reference,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("RR",reference,trt_indicator_no_reference[j], "_est")))), eval(parse(text =(paste0("OR",reference,trt_indicator_no_reference[j], "_est"))))))
    assign(paste0("ATT",reference,trt_indicator_no_reference[j]), list(round(eval(parse(text =(paste0("att",reference,trt_indicator_no_reference[j])))), digits = 3)))
    assign(paste0("ATT",reference,trt_indicator_no_reference[j]), stats::setNames(eval(parse(text =(paste0("ATT",reference,trt_indicator_no_reference[j])))), paste0("ATT",reference,trt_indicator_no_reference[j])))
    result <- c(result, (eval(parse(text =(paste0("ATT",reference,trt_indicator_no_reference[j]))))))
  }

  return(result)

  # for (m in 1:ndpost) {
  #
  #   # Estimate E(Y1|trt=1), E(Y2|trt=1), E(Y3|trt=1)
  #   y1_pred = mean(stats::rbinom(n1, 1, pred_prop11))
  #   y2_pred = mean(stats::rbinom(n1, 1, pred_prop12))
  #   y3_pred = mean(stats::rbinom(n1, 1, pred_prop13))
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
}
