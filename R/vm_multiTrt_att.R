#===========================#
#Author: Michael Lopez      #
#===========================#
#' Vector matching Matching (VM matching)
#'
#' This function implements the VM matching method. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
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
#'idata = data_gen(n = 120, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind,method = "VM Matching", estimand = "ATT", reference = 1)
vm_multiTrt_att = function(y, x, trt, reference = parent.frame()$reference_trt) {
  n_trt <- length(unique(trt))
  if (n_trt > 3){
    stop("We do not recommend using VM for more than 3 treatments")
  }
  x <- x %>% dplyr::select(-trt_ind)
  trt <- factor(trt)
  # estimate generalized propensity scores using multinomial logistic regression
  xydata = cbind(trt, x)


  # WARNING: reorder the dataset by "Treatment 1", "Treatment 2", "Treatment 3"
  # it will affect the identification of matched indices in row 249
  xydata = xydata[order(xydata$trt),]


  # ps model 1
  ps_fit = nnet::multinom(trt ~ ., data = xydata, trace = FALSE)
  probs_logit1 = data.frame(stats::fitted(ps_fit))
  colnames_probs_logit1 <- NULL
  for (i in 1:n_trt){
    colnames_probs_logit1_once <- paste0("p",i)
    colnames_probs_logit1 <- c(colnames_probs_logit1, colnames_probs_logit1_once)
  }
  colnames(probs_logit1) = colnames_probs_logit1
  xydata = cbind(xydata, probs_logit1)

  # Determine eligibility
  # min_max_Ps <- xydata %>%
  #   dplyr::group_by(trt) %>%
  #   dplyr::summarise(min1 = min(p1), max1 = max(p1),
  #             min2 = min(p2), max2 = max(p2),
  #             min3 = min(p3), max3 = max(p3))
  min_max_Ps <- NULL
  for (i in 1: n_trt){
    xydata_summarise_once <- xydata %>%
      dplyr::group_by(trt) %>%
      dplyr::summarise(min(eval(parse(text = paste0("p",i)))), max(eval(parse(text = paste0("p",i)))))
    names(xydata_summarise_once)[c(2,3)] <- c(paste0("min",i), paste0("max",i))
    if (i == 1){
      min_max_Ps <- xydata_summarise_once
    } else {
      min_max_Ps <- min_max_Ps %>% dplyr::inner_join(xydata_summarise_once)
    }
  }

  #min_max_Ps
  Eligible <- TRUE
  for (i in 1:n_trt){
    Eligible_once <- xydata[[paste0("p", i)]] >= max(min_max_Ps[[paste0("min", i)]]) & xydata[[paste0("p", i)]] <= min(min_max_Ps[[paste0("max", i)]])
    Eligible <- Eligible & Eligible_once
  }
  xydata <- xydata %>%
    dplyr::mutate(Eligible = Eligible)

  xydata = cbind(y, xydata)
  xydata = dplyr::filter(xydata, Eligible)

  # Calculate new propensity scores for eligible subjects
  ps_fit_E = nnet::multinom(trt ~ ., data = xydata[,-1], trace = FALSE)
  probs_logit1_E = stats::fitted(ps_fit_E)

  colnames_probs_logit1_E <- NULL
  for (i in 1:n_trt){
    colnames_probs_logit1_E_once <- paste0("p",i)
    colnames_probs_logit1_E <- c(colnames_probs_logit1_E, colnames_probs_logit1_E_once)
  }
  colnames(probs_logit1_E) = colnames_probs_logit1_E
  # colnames(probs_logit1_E) = c("p1", "p2", "p3")
  xydata <- xydata %>%
    dplyr :: select(-"p1", -"p2", -"p3")
  xydata = cbind(xydata, probs_logit1_E)

  for (i in 1:n_trt){
    assign(paste0("n",i), sum(xydata$trt == i))
  }

  ### Vector Matching for ATT for outcome 1 (comp_resp_obs)
  # Stratify car::logit(r(ti, Xi)) using K-means clustering
  clustnum <- 5

  for (i in 1:n_trt){
    temp_once <-  stats::kmeans(car::logit(xydata[[paste0("p",i)]]), clustnum)
    xydata <- xydata %>%
      dplyr::mutate(temp_once$cluster)
    names(xydata)[length(names(xydata))] <- paste0("Quint",i)
  }

  colnames(xydata)[1:2] = c("Y_obs","treat")
  trt_indicator = 1:n_trt
  trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference]

  for (i in 1:length((trt_indicator_no_reference))){
    assign(paste0("temp", reference,trt_indicator_no_reference[i]), dplyr::filter(xydata, treat %in% c(reference,trt_indicator_no_reference[i])))
  }
  assign(paste0("temp", trt_indicator_no_reference[1],trt_indicator_no_reference[2]), dplyr::filter(xydata, treat %in% trt_indicator_no_reference))
  # matching
  for (i in 1:length((trt_indicator_no_reference))){
    trt_indicator_left <- trt_indicator[!(trt_indicator %in% c(reference, trt_indicator_no_reference[i]))]
    assign(paste0("match",reference,trt_indicator_no_reference[i]),   Matching::Matchby(Y = eval(parse(text = paste0("temp", reference,trt_indicator_no_reference[i])))[["Y_obs"]], Tr = eval(parse(text = paste0("temp", reference,trt_indicator_no_reference[i])))[["treat"]] == reference,
                                                                                        X = car::logit(eval(parse(text = paste0("temp", reference,trt_indicator_no_reference[i])))[[paste0("p",reference)]]), by = eval(parse(text = paste0("temp", reference,trt_indicator_no_reference[i])))[[paste0("Quint", trt_indicator_left)]],
                                                                                        caliper = 0.5*stats::sd(car::logit(eval(parse(text = paste0("temp", reference,trt_indicator_no_reference[i])))[[paste0("p",reference)]])),  replace = T, estimand = "ATT", print.level = 0))
  }

  # Identify the matched subgroups
  rownames(xydata) = 1:nrow(xydata)
  xydata$id = 1:nrow(xydata)
  eligible_matching <- TRUE
  for (i in 1:length((trt_indicator_no_reference))){
    eligible_matching_once <- xydata$id %in% eval(parse(text = paste0("match",reference,trt_indicator_no_reference[i])))[["index.treated"]]
    eligible_matching <- eligible_matching & eligible_matching_once
  }

  xydata$both_1 <- eligible_matching
  temp = xydata[xydata$both_1 == "TRUE", ]


  for (i in 1:length((trt_indicator_no_reference))){
    if ( i == 1){
      assign(paste0("m",reference,trt_indicator_no_reference[i]), cbind( eval(parse(text = paste0("match",reference,trt_indicator_no_reference[i])))[["index.treated"]], eval(parse(text = paste0("match",reference,trt_indicator_no_reference[i])))[["index.control"]]))
    } else {
      trt_indicator_left <- trt_indicator[!(trt_indicator %in% c(reference, trt_indicator_no_reference[i]))]
      assign(paste0("m",reference,trt_indicator_no_reference[i]), cbind( eval(parse(text = paste0("match",reference,trt_indicator_no_reference[i])))[["index.treated"]], eval(parse(text = paste0("match",reference,trt_indicator_no_reference[i])))[["index.control"]] + sum(xydata[["treat"]] == trt_indicator_left)))
    }
  }
  # WARNING: reorder the dataset by "Treatment 1", "Treatment 2", "Treatment 3" before running the following lines
  for (i in 1:length((trt_indicator_no_reference))){
    assign(paste0("m",reference,trt_indicator_no_reference[i]), eval(parse(text = paste0("m",reference,trt_indicator_no_reference[i])))[eval(parse(text = paste0("m",reference,trt_indicator_no_reference[i])))[,1] %in% rownames(temp),])
  }
  triplets <- NULL
  for (i in 1:length((trt_indicator_no_reference))){
    triplets_once <- eval(parse(text = paste0("m",reference,trt_indicator_no_reference[i])))[order(eval(parse(text = paste0("m",reference,trt_indicator_no_reference[i])))[,1]), ]
    triplets <- cbind(triplets,triplets_once)
  }
  triplets = as.matrix(triplets[,c(1, 2, 4)])
  n_trip = nrow(triplets)
  # Matching Estimator
  # For subjects receiving reference treatment
  for (i in 1:3){
    assign(paste0("Y",i,"_imp"), xydata$Y_obs[triplets[,i]])
    assign(paste0("y",i,"_hat"), mean(eval(parse(text = paste0("Y",i,"_imp")))))
  }

  for (i in 1:length((trt_indicator_no_reference))){
    assign(paste0("RD", reference, trt_indicator_no_reference[i],"_est"), eval(parse(text = paste0("y",reference,"_hat"))) - eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))))
    assign(paste0("RR", reference, trt_indicator_no_reference[i],"_est"), eval(parse(text = paste0("y",reference,"_hat"))) / eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))))
    assign(paste0("OR", reference, trt_indicator_no_reference[i],"_est"), (eval(parse(text = paste0("y",reference,"_hat"))) / ( 1- eval(parse(text = paste0("y",reference,"_hat"))))) / (eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))) / ( 1- eval(parse(text = paste0("y",trt_indicator_no_reference[i],"_hat"))))))
  }

  result <- NULL
  result_list <- NULL
  for (j in 1:length(trt_indicator_no_reference)){
    result_once <- rbind(eval(parse(text = paste0("RD",reference,trt_indicator_no_reference[j],"_est"))), eval(parse(text = paste0("RR",reference,trt_indicator_no_reference[j],"_est"))), eval(parse(text = paste0("OR",reference,trt_indicator_no_reference[j],"_est"))))
    colnames(result_once) <- "EST"
    rownames(result_once) <- c("RD", "RR", "OR")
    result_once_list <- list(result_once)
    names(result_once_list) <- paste0("ATT",reference,trt_indicator_no_reference[j])
    result_list <- c(result_list, result_once_list)
  }
return(result_list)

}


