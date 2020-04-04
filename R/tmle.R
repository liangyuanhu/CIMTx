
#===========================================================#
#Author: Jiayi Ji, Liangyuan Hu #
#Reference: Rose and Normand (2019) 75(1) 289-296 Biometrics#
#===========================================================#
#' @importFrom magrittr "%>%"
#' @import SuperLearner
tmle <- function(y, trt, psdat,...){
  trt_ind <- trt
  Xmat <- psdat %>%
    dplyr::as_tibble() %>%
    dplyr::mutate("trt1" = dplyr::case_when(trt_ind == 1 ~ 1,
                            TRUE ~ 0),
           "trt2" = dplyr::case_when(trt_ind == 2 ~ 1,
                            TRUE ~ 0),
           "trt3" = dplyr::case_when(trt_ind == 3 ~ 1,
                            TRUE ~ 0)) %>%
    as.data.frame()
  K <- 3
  n <- dim(Xmat)[1]
  t_mat <- Xmat[,c("trt1", "trt2", "trt3")]
  W <- psdat[, -1]
  #---------------------------------------------#
  ###Create Counterfactual Treatment Scenarios###
  #---------------------------------------------#

  trt1_countfactual <- W %>%  ## should use W since dummy variables of treatment indicators are included, cannot include categorical trt indicator again
    dplyr::mutate("trt1"= 1, "trt2" = 0, "trt3" = 0)  #had everyone received trt 1

  trt2_countfactual <- W %>%
    dplyr::mutate("trt1"= 0, "trt2" = 1, "trt3" = 0) #had everyone received trt 2

  trt3_countfactual <- W %>%
    dplyr::mutate("trt1"= 0, "trt2" = 0, "trt3" = 1) #had everyone received trt 3

  trt_countfactual_combined <- as.data.frame(rbind(trt1_countfactual, trt2_countfactual, trt3_countfactual))

  #-------------------------------------------------------------------------------------------#
  ###Run Super Learner Once, Obtain Initial Predicted Values for All Counterfactual Settings###
  #-------------------------------------------------------------------------------------------#
  #Step 1: Estimating the outcome regression using super learner

  sl_fit <-
    SuperLearner::SuperLearner(
      Y = y,
      X = Xmat[,-1],
      newX = trt_countfactual_combined,
      SL.library = c("SL.glm", "SL.gam",
                     "SL.knn"),
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
  rownames(trt_results) <-
    c("trt1", "trt2", "trt3")
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
      family = "stats::binomial",
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

  trt_result_one_repetition_result <- trt_result_one_repetition %>%
    dplyr::select("treatment",   "EYt" ) %>%
    tidyr::spread("treatment", "EYt") %>%
    dplyr::mutate(
      RD12 = trt1 - trt2,
      RD13 = trt1 - trt3,
      RD23 = trt2 - trt3,
      RR12 = trt1 / trt2,
      RR13 = trt1 / trt3,
      RR23 = trt2 / trt3,
      OR12 = (trt1 / (1 - trt1)) / (trt2 / (1 - trt2)),
      OR13 = (trt1 / (1 - trt1)) / (trt3 / (1 - trt3)),
      OR23 = (trt2 / (1 - trt2)) / (trt3 / (1 - trt3))
    ) %>%
    dplyr::select(-"trt1" , -"trt2" , -"trt3")

  res12 = rbind(RD = c(trt_result_one_repetition_result$RD12), RR = trt_result_one_repetition_result$RR12, OR = trt_result_one_repetition_result$OR12)
  res13 = rbind(RD = c(trt_result_one_repetition_result$RD13), RR = trt_result_one_repetition_result$RR13, OR = trt_result_one_repetition_result$OR13)
  res23 = rbind(RD = c(trt_result_one_repetition_result$RD23), RR = trt_result_one_repetition_result$RR23, OR = trt_result_one_repetition_result$OR23)

  colnames(res12) <- "EST"
  colnames(res13) <- "EST"
  colnames(res23) <- "EST"
  list(ATE12 = res12,
       ATE13 = res13,
       ATE23 = res23)
}
