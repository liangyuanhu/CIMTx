#######################################################################
# Functions to estimate causal effects of multiple treatment using BART
# Assume 3 treatment
# Author: Liangyuan Hu, Chenyang Gu
#######################################################################
#' Bayesian Additive Regression Trees (BART) for ATE estimation
#'
#' This function implements the BART method when estimand is ATE. Please use our main function causal_multi_treat.R.
#'
#' @param y numeric vector for the binary outcome
#' @param x dataframe including the treatment indicator and the covariates
#' @param trt numeric vector for the treatment indicator
#' @param k For binary y, k is the number of prior standard deviations f(x) is away from +/-3. The bigger k is, the more conservative the fitting will be.
#' @param discard discarding rules for BART method, inherited from causal_multi_treat.R
#' @param ntree The number of trees in the sum
#' @param ndpost The number of posterior draws returned
#' @param nskip Number of MCMC iterations to be treated as burn in
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
#'set.seed(3242019)
#'idata = data_gen(n = 3, ratio =1,scenario = 1)
#'trt_ind <- as.numeric(idata$trtdat$trt_ind)
#'all_vars <- idata$trtdat[, -1] #exclude treatment indicator
#'y <- idata$Yobs
#'causal_multi_treat(y = y, x = idata$trtdat,
#'trt = trt_ind, method = "BART", estimand = "ATE", discard = "No", ndpost = 10)
bart_multiTrt_ate = function(y, x, trt, k=2, discard = "No", ntree=100, ndpost=parent.frame()$ndpost, nskip=1000) {
    n_trt <- length(unique(trt))
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
    for (i in 1:n_trt){
      assign(paste0("xp",i), xt[trt==i,])
      for (j in 1:(n_trt)){

        assign(paste0("xp",j), eval(parse(text = paste0("xp",i))) %>%
                 dplyr::mutate(trt = j))
        assign(paste0("bart_pred",i,j), BART::pwbart(eval(parse(text = paste0("xp",j))), bart_mod$treedraws,mu=mean(y)))
        assign(paste0("pred_prop",i,j), stats::pnorm(eval(parse(text = paste0("bart_pred",i,j)))))
      }
    }

    if (discard == "No") {
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
      assign(paste0("pred_prop",i,j), stats::pnorm(eval(parse(text = paste0("bart_pred",i,j)))))
        }
      }
    } else if (discard == "Stringent"){
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
      assign(paste0("post.ind.sd",i,j), apply(eval(parse(text = paste0("pred_prop",i,j))), 2, stats::sd))
        }
      }
      for (i in 1:n_trt){
        n_trt_no_i <- unique(trt)[unique(trt)!=i]
        assign(paste0("eligible",i), TRUE)
        for (j in 1:length((n_trt_no_i))) {
          assign(paste0("eligible",i,n_trt_no_i[j]), eval(parse(text = paste0("post.ind.sd",i,n_trt_no_i[j]))) <= 2 * eval(parse(text = paste0("post.ind.sd",i,i))))
        }
        assign(paste0("eligible",i), eval(parse(text = paste0("eligible",i,n_trt_no_i[j]))) & eval(parse(text = paste0("eligible",i))))
        assign(paste0("n_",i, "_discard"), sum(eval(parse(text = paste0("eligible",i)))== F))
      }
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("pred_prop",i,j), eval(parse(text = paste0("pred_prop",i,j))) %>%
                   as.data.frame() %>%
                   dplyr::select(which(eval(parse(text = paste0("eligible",i))))) %>%
                   as.matrix())
        }
      }
    } else if(discard == "Lenient"){
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("post.ind.sd",i,j), apply(eval(parse(text = paste0("pred_prop",i,j))), 2, stats::sd))
        }
        assign(paste0("threshold",i), max(eval(parse(text = paste0("pred_prop",i,i)))))
      }

      for (i in 1:n_trt){
        n_trt_no_i <- unique(trt)[unique(trt)!=i]
        assign(paste0("eligible",i), TRUE)
        for (j in 1:length((n_trt_no_i))) {
          assign(paste0("criteria",i,n_trt_no_i[j]), eval(parse(text = paste0("post.ind.sd",i,n_trt_no_i[j]))) <= eval(parse(text = paste0("threshold",i))))
        }
        assign(paste0("eligible",i), eval(parse(text = paste0("criteria",i,n_trt_no_i[j]))) & eval(parse(text = paste0("criteria",i))))
        assign(paste0("n_",i, "_discard"), sum(eval(parse(text = paste0("eligible",i)))== F))
      }
      for (i in 1:n_trt){
        for (j in 1:(n_trt)){
          assign(paste0("pred_prop",i,j), eval(parse(text = paste0("pred_prop",i,j))) %>%
                   as.data.frame() %>%
                   dplyr::select(which(eval(parse(text = paste0("eligible",i))))) %>%
                   as.matrix())
        }
      }
    }


    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        assign(paste0("RD",i,j, "_est"), NULL)
        assign(paste0("RR",i,j, "_est"), NULL)
        assign(paste0("OR",i,j, "_est"), NULL)
      }
    }

    for (m in 1:ndpost) {
      for (i in 1:n_trt){
        assign(paste0("y",i), NULL)
        for (j in 1:n_trt){
          assign(paste0("y",i,j), c(stats::rbinom(eval(parse(text =(paste0("n",j)))), 1, eval(parse(text =(paste0("pred_prop",i,j)))) %>% as.data.frame %>% dplyr::slice(m) %>% as.numeric())))
          assign(paste0("y",i), c(eval(parse(text =(paste0("y",i)))), eval(parse(text =(paste0("y",i,j))))))
        }
        assign(paste0("y",i, "_pred_",m), mean(eval(parse(text =(paste0("y",i))))))
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
    }
    # for (m in 1:ndpost) {
    #   for (i in 1:(n_trt-1)){
    #     for (j in (i + 1):n_trt){
    #       assign(paste0("RD",i,j, "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) - eval(parse(text =(paste0("y",j, "_pred_",m)))))
    #       assign(paste0("RR",i,j, "_est_m"), eval(parse(text =(paste0("y",i, "_pred_",m)))) / eval(parse(text =(paste0("y",j, "_pred_",m)))))
    #       assign(paste0("OR",i,j, "_est_m"), (eval(parse(text =(paste0("y",i, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",i, "_pred_",m)))))) / (eval(parse(text =(paste0("y",j, "_pred_",m)))) / (1 - eval(parse(text =(paste0("y",j, "_pred_",m)))))))
    #       assign(paste0("RD",i,j, "_est"), c(eval(parse(text =(paste0("RD",i,j, "_est")))), eval(parse(text =(paste0("RD",i,j, "_est_m"))))))
    #
    #       assign(paste0("RR",i,j, "_est"), c(eval(parse(text =(paste0("RR",i,j, "_est")))), eval(parse(text =(paste0("RR",i,j, "_est_m"))))))
    #
    #       assign(paste0("OR",i,j, "_est"), c(eval(parse(text =(paste0("OR",i,j, "_est")))), eval(parse(text =(paste0("OR",i,j, "_est_m"))))))
    #     }
    #   }
    # }
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
}
