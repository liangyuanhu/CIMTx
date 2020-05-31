
#===========================================================#
#Author: Jiayi Ji, Liangyuan Hu #
#Reference: Rose and Normand (2019) 75(1) 289-296 Biometrics#
#===========================================================#
#' Targeted maximum likelihood (TMLE)
#'
#' This function implements the TMLE method. Please use our main function causal_multi_treat.R.
#' @param y numeric vector for the binary outcome
#'
#' @param trt numeric vector for the treatment indicator
#' @param x data frame containing the treatment indicator and covariates
#' @param ... Other arguments
#'
#' @importFrom magrittr "%>%"
#' @import SuperLearner
#' @export
#' @examples
#'library(CIMTx)
#'set.seed(3242019)
#'idata = data_gen(n = 120, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#' x = idata$trtdat
#' tmle(y = y, trt = trt_ind, x = idata$trtdat, SL.library = c("SL.glm"))
#'causal_multi_treat(y = y, x = idata$trtdat, trt = trt_ind,
#'method = "TMLE", estimand = "ATE")
tmle <- function(y, trt, x,...){
  trt_ind <- trt
  n_trt <- length(unique(trt))
  Xmat <- x
  for (i in 1:n_trt){
    Xmat <- Xmat %>%
      dplyr::as_tibble() %>%
      # dplyr::mutate(!!sym(new_col_name) := dplyr::case_when(trt_ind == i ~ 1,
      #                                                       TRUE ~ 0))%>%
      dplyr::mutate(dplyr::case_when(trt_ind == i ~ 1,
                                                            TRUE ~ 0))%>%
        as.data.frame()
    names(Xmat)[length(Xmat)] <- paste0("trt",i)
  }
  # Xmat <- x %>%
  #   dplyr::as_tibble() %>%
  #   dplyr::mutate("trt1" = dplyr::case_when(trt_ind == 1 ~ 1,
  #                           TRUE ~ 0),
  #          "trt2" = dplyr::case_when(trt_ind == 2 ~ 1,
  #                           TRUE ~ 0),
  #          "trt3" = dplyr::case_when(trt_ind == 3 ~ 1,
  #                           TRUE ~ 0)) %>%
  #   as.data.frame()
  K <- n_trt
  # K <- 3
  n <- dim(Xmat)[1]
  t_mat <- NULL
  for (i in 1:n_trt){
    t_mat_once <- Xmat %>%
      dplyr::select(paste0("trt",i))
    t_mat <- dplyr::bind_cols(t_mat, t_mat_once)
  }
  # t_mat <- Xmat[,c("trt1", "trt2", "trt3")]
  W <- x[, -1]
  #---------------------------------------------#
  ###Create Counterfactual Treatment Scenarios###
  #---------------------------------------------#
  for (i in 1:n_trt){
    assign(paste0("trt", i, "_countfactual"), W)
  }

  for (j in 1:n_trt){
    new_col_name <- paste0("trt",j)
    for (i in 1:n_trt){
      if (i == j) {
        names_trt_countfactual <- names(eval(parse(text =(paste0("trt", i, "_countfactual")))))
        # assign(paste0("trt", i, "_countfactual"), eval(parse(text = paste0("trt", i, "_countfactual"))) %>% dplyr::mutate(!!sym(new_col_name) := 1))
        assign(paste0("trt", i, "_countfactual"), eval(parse(text = paste0("trt", i, "_countfactual"))) %>% dplyr::mutate(1))
        # names(eval(parse(text = paste0("trt", i, "_countfactual"))))[length(eval(parse(text = paste0("trt", i, "_countfactual"))))] <- paste0("trt", j)
        # stats::setNames(eval(parse(text =(paste0("trt", i, "_countfactual")))), c(names(W), paste0("trt", j))) %>% head
        assign(paste0("trt", i, "_countfactual"), stats::setNames(eval(parse(text =(paste0("trt", i, "_countfactual")))), c(names_trt_countfactual, paste0("trt", j))))
      } else {
        names_trt_countfactual <- names(eval(parse(text =(paste0("trt", i, "_countfactual")))))
        # assign(paste0("trt", i, "_countfactual"), eval(parse(text = paste0("trt", i, "_countfactual"))) %>% dplyr::mutate(!!sym(new_col_name) := 0))
        assign(paste0("trt", i, "_countfactual"), eval(parse(text = paste0("trt", i, "_countfactual"))) %>% dplyr::mutate(0))
        # names(eval(parse(text = paste0("trt", i, "_countfactual"))))[length(eval(parse(text = paste0("trt", i, "_countfactual"))))] <- paste0("trt", j)
        assign(paste0("trt", i, "_countfactual"), stats::setNames(eval(parse(text =(paste0("trt", i, "_countfactual")))), c(names_trt_countfactual, paste0("trt", j))))
      }
    }
  }
  # trt1_countfactual %>% head
  # trt2_countfactual %>% head
  # trt3_countfactual %>% head
  # trt1_countfactual <- W %>%  ## should use W since dummy variables of treatment indicators are included, cannot include categorical trt indicator again
  #   dplyr::mutate("trt1"= 1, "trt2" = 0, "trt3" = 0)  #had everyone received trt 1
  #
  # trt2_countfactual <- W %>%
  #   dplyr::mutate("trt1"= 0, "trt2" = 1, "trt3" = 0) #had everyone received trt 2
  #
  # trt3_countfactual <- W %>%
  #   dplyr::mutate("trt1"= 0, "trt2" = 0, "trt3" = 1) #had everyone received trt 3
  trt_countfactual_combined <- NULL
  for (i in 1:n_trt){
    trt_countfactual_combined <- as.data.frame(rbind(trt_countfactual_combined, eval(parse(text = paste0("trt", i, "_countfactual")))))
  }
  # trt_countfactual_combined <- as.data.frame(rbind(trt1_countfactual, trt2_countfactual, trt3_countfactual))

  #-------------------------------------------------------------------------------------------#
  ###Run Super Learner Once, Obtain Initial Predicted Values for All Counterfactual Settings###
  #-------------------------------------------------------------------------------------------#
  #Step 1: Estimating the outcome regression using super learner

  sl_fit <-
    SuperLearner::SuperLearner(
      Y = y,
      X = Xmat[,-1],
      newX = trt_countfactual_combined,
      SL.library = c("SL.glm"),
      family = stats::binomial(),
      method = "method.NNLS",
      verbose = F
    )

  q_0 <- rep(0, n * K)
  q_tvector <- cbind(q_0, sl_fit$SL.predict)
  q_tmat <-
    matrix(unlist(split(
      as.data.frame(q_tvector), rep(1:K, each = n)
    )), ncol = 2 * K)    # 1200 rows and 6 columns, columns 1-2 for treatment 1, columns 3-4 for trt 2, columns 5-6 for trt3

  #-----------------------------------------------------------------------------------#
  ###Run TMLE to Calculate Point Estimates of each T=t with Custom Cluster-Based SEs###
  #-----------------------------------------------------------------------------------#
  #Steps 2-5 are performed in this code chunk#

  trt_results <- matrix(NA, nrow = K, ncol = 4)
  trt_results_row_names <- NULL
  for (i in 1:n_trt){
    trt_results_row_names_once <- paste0("trt",i)
    trt_results_row_names <- c(trt_results_row_names, trt_results_row_names_once)
  }
  rownames(trt_results) <- trt_results_row_names
  # rownames(trt_results) <-
  #   c("trt1", "trt2", "trt3")
  colnames(trt_results) <- c("EYt", "SE", "CI1", "CI2")
  start <- 1
  end <- 2
  gbound <- 0.025
  for(t in 1:K) {
    #Step 2: Super learner fit for P(T_k=t|W) specified with g.SL.library=SL.library in tmle call#
    #Steps 3-4: Performed in tmle call, target  EYt parameter using A=NULL and Delta=Tmat[,t]#
    fit <- tmle::tmle(
      Y = y,
      A = NULL,
      Delta = t_mat[, t],
      W = W,
      Q = q_tmat[, c(start, end)],
      g.SL.library = c("SL.glm", "SL.gam",
                       "SL.knn"),
      family = "binomial",
      verbose = TRUE
    )
    #Step 5: The parameter estimates are stored in fit$estimates$EY1$psi#
    trt_results[t, 1] <- fit$estimates$EY1$psi
    trt_results[t, 2] <- fit$estimates$EY1$var.psi
    trt_results[t, 3] <- fit$estimates$EY1$CI[1]
    trt_results[t, 4] <- fit$estimates$EY1$CI[2]
    start <- start + 2
    end <- end + 2
  }

  trt_result_one_repetition <- trt_results %>%
    dplyr::as_tibble(rownames = "treatment")
    # dplyr::mutate(repetition = i)
  # trt_result_all_scenario_2_weak <- trt_result_all_scenario_2_weak %>%
  #   dplyr::bind_rows(trt_result_one_repetition)
  for (i in 1:n_trt){
    assign(paste0("mu_",i, "_hat"), trt_result_one_repetition %>% dplyr::select("EYt") %>% dplyr::slice(i) %>% dplyr::pull("EYt"))
  }
  result_list <- NULL
  for (i in 1:(n_trt-1)){
    result_once <- NULL
    for (j in (i + 1):n_trt){
      assign(paste0("RD",i,j),eval(parse(text = paste0("mu_",i, "_hat"))) - eval(parse(text = paste0("mu_",j, "_hat"))))
      assign(paste0("RR",i,j),eval(parse(text = paste0("mu_",i, "_hat"))) / eval(parse(text = paste0("mu_",j, "_hat"))))
      assign(paste0("OR",i,j), (eval(parse(text = paste0("mu_",i, "_hat"))) /(1 - eval(parse(text = paste0("mu_",i, "_hat"))))) / (eval(parse(text = paste0("mu_",j, "_hat"))) /(1 - eval(parse(text = paste0("mu_",j, "_hat"))))))
      result_once <- rbind(eval(parse(text = paste0("RD",i,j))), eval(parse(text = paste0("RR",i,j))), eval(parse(text = paste0("OR",i,j))))
      colnames(result_once) <- "EST"
      rownames(result_once) <- c("RD", "RR", "OR")
      result_once_list <- list(result_once)
      names(result_once_list) <- paste0("ATE",i,j)
      result_list <- c(result_list, result_once_list)
    }
  }
  return(result_list)

  # trt_result_one_repetition_result <- trt_result_one_repetition %>%
  #   dplyr::select("treatment",   "EYt" ) %>%
  #   tidyr::spread("treatment", "EYt") %>%
  #   dplyr::mutate(
  #     RD12 = trt1 - trt2,
  #     RD13 = trt1 - trt3,
  #     RD23 = trt2 - trt3,
  #     RR12 = trt1 / trt2,
  #     RR13 = trt1 / trt3,
  #     RR23 = trt2 / trt3,
  #     OR12 = (trt1 / (1 - trt1)) / (trt2 / (1 - trt2)),
  #     OR13 = (trt1 / (1 - trt1)) / (trt3 / (1 - trt3)),
  #     OR23 = (trt2 / (1 - trt2)) / (trt3 / (1 - trt3))
  #   ) %>%
  #   dplyr::select(-"trt1" , -"trt2" , -"trt3")
  #
  # res12 = rbind(RD = c(trt_result_one_repetition_result$RD12), RR = trt_result_one_repetition_result$RR12, OR = trt_result_one_repetition_result$OR12)
  # res13 = rbind(RD = c(trt_result_one_repetition_result$RD13), RR = trt_result_one_repetition_result$RR13, OR = trt_result_one_repetition_result$OR13)
  # res23 = rbind(RD = c(trt_result_one_repetition_result$RD23), RR = trt_result_one_repetition_result$RR23, OR = trt_result_one_repetition_result$OR23)
  #
  # colnames(res12) <- "EST"
  # colnames(res13) <- "EST"
  # colnames(res23) <- "EST"
  # list(ATE12 = res12,
  #      ATE13 = res13,
  #      ATE23 = res23)
}
