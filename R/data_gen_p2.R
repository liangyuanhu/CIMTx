#==============================#
# Author: Liangyuan Hu        #
#=============================#
#' Data generation function for scenario 2
#' This function generates data to test different causal inference methods for scenario 2.  Please use our main function data_gen.R
#'
#' @param n total number of units for simulation
#' @param p number of predictors
#' @param overlap levels of covariate overlap: Please select: weak, strong, moderate
#' @param all_confounder TRUE or FALSE. overlap is lacking for a variable that is not predictive of the outcome (all_confounder equals to TRUE) or situations when it is lacking for a true confounder (all_confounder equals to FALSE)
#'
#' @return list with the 5 elements. Nested within each list, it contains
#' \item{n:}{Number of units for simulation}
#' \item{trt_ind:}{A data frame with number of rows equals to n and 11 columns}
#' \item{Y:}{Observed binary outcome for 3 treatments}
#' \item{Yobs:}{Observed binary outcome}
#' \item{Est:}{True ATE/ATT for RD/RR/OR}
#' @export
#' @examples
#' library(CIMTx)
#' set.seed(3242019)
#' data_gen_p2(n = 116, p =10, overlap = "weak", all_confounder = TRUE)
data_gen_p2 = function(n = 11600, p =10, overlap = "weak", all_confounder = TRUE ) {
  #generate treatment label W
  W = sample(c(1,2,3), size=n, replace=TRUE, prob=c(.034,.517,.449))
  #generate covariates distribution conditional on treatment
  #==================================#
  #===========weak overlap===========#
  #==================================#
  if (overlap == "weak") {
    Xcon<-matrix(NA, nrow=n, ncol=5); Xcat<-matrix(NA, nrow=n, ncol=5);
    for (i in 1:n){
      Xcon[i,] = matrix(stats::rnorm(5,0,1)*(W[i]==1) + stats::rnorm(5, 1, 1)*(W[i]==2) + stats::rnorm(5, 2, 1)*(W[i]==3),ncol=5);
      Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.3,.3,.4))*(W[i]==1)
                        + sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.6,.2,.2))*(W[i]==2)
                        + sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.8,.1,.1))*(W[i]==3),ncol=5)
    }
    x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5];
    x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
    trtdat<- data.frame(W,Xcon, Xcat)

    #*************************#
    #draw potential outcomes
    tau1 = .1; tau2 = .6; tau3= 1.0
    Yp1 = expit(tau1+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 +0.2*x10
                -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .27*x4^3 )
    Yp2 = expit(tau2+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 + 0.2*x10
                -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .27*x4^3 )
    Yp3 = expit(tau3+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 + 0.2*x10
                -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .27*x4^3 )

    #potential outcomes
    Y1 = Y2 = Y3= NULL
    for (i in 1:n) {
      Y1[i] = stats::rbinom(1, 1, Yp1[i])
      Y2[i] = stats::rbinom(1, 1, Yp2[i])
      Y3[i] = stats::rbinom(1, 1, Yp3[i])}
    # table(Y1)[2]/sum(table(Y1))
    # table(Y2)[2]/sum(table(Y2))
    # table(Y3)[2]/sum(table(Y3))
    # observed outcomes
    Y = cbind(Y1,Y2,Y3)
    YW = cbind(Y,W)
    Yobs = apply(YW, 1, function(x) x[1:3][x[4]])
  }
  #==================================#
  #=========strong overlap===========#
  #==================================#
  if (overlap == "strong") {
    Xcon<-matrix(NA, nrow=n, ncol=5); Xcat<-matrix(NA, nrow=n, ncol=5);
    for (i in 1:n) {
      Xcon[i,] = matrix(stats::rnorm(5, mean = 0.05*W[i], sd = 1- 0.05*W[i]),ncol=5);
      Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.3-0.001*W[i],.3+.001*W[i],.4)),ncol=5)
    }
    x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5];
    x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
    trtdat<- data.frame(W, Xcon, Xcat)

    #********************************#
    #draw potential outcomes
    tau1 = .1; tau2 = .6; tau3= 1.0
    Yp1 = expit(tau1+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 +0.2*x10
                -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .27*x4^3 )
    Yp2 = expit(tau2+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 + 0.2*x10
                -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .27*x4^3 )
    Yp3 = expit(tau3+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 + 0.2*x10
                -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .27*x4^3 )

    #potential outcomes
    Y1 = Y2 = Y3= NULL
    for (i in 1:n) {
      Y1[i] = stats::rbinom(1, 1, Yp1[i])
      Y2[i] = stats::rbinom(1, 1, Yp2[i])
      Y3[i] = stats::rbinom(1, 1, Yp3[i])}
    # observed outcomes
    Y = cbind(Y1,Y2,Y3)
    YW = cbind(Y,W)
    Yobs = apply(YW, 1, function(x) x[1:3][x[4]])
  }
  #==================================#
  #========moderate overlap==========#
  #==================================#
  if (overlap == "moderate"){
    Xcon<-matrix(NA, nrow=n, ncol=5);  Xcat<-matrix(NA, nrow=n, ncol=5)
    for (i in 1:n) {
      Xcon[i,] = matrix(stats::rnorm(5, mean = 0.05*W[i], sd = 1- 0.05*W[i]),ncol=5);
      Xcat[i,] = matrix(sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.3,.3,.4))*(W[i]==1)
                        + sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.6,.2,.2))*(W[i]==2)
                        + sample (c(1,2,3), size = 5,replace =TRUE, prob = c(.8,.1,.1))*(W[i]==3),ncol=5)
    }
    x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5  = Xcon[,5];
    x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
    trtdat<- data.frame(W, Xcon, Xcat)

    if (all_confounder == TRUE) {
      tau1 = .1; tau2 = .3; tau3= 0.5
      Yp1 = expit(tau1+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 -0.3*x10
                  -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .1*x4^3 -.2*x10^2)
      Yp2 = expit(tau2+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9 - 0.3*x10
                  -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .1*x4^3 -.2*x10^2)
      Yp3 = expit(tau3+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 + 0.2*x6 + 0.4*x7 - 0.6*x8 + 0.2*x9- 0.3*x10
                  -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2 + .4*x1*x7 - 0.5*x2*x8 + .2*x6*x9 + .1*x4^3 -.2*x10^2)

      #potential outcomes
      Y1 = Y2 = Y3= NULL
      for (i in 1:n) {
        Y1[i] = stats::rbinom(1, 1, Yp1[i])
        Y2[i] = stats::rbinom(1, 1, Yp2[i])
        Y3[i] = stats::rbinom(1, 1, Yp3[i])}
      # observed outcomes
      Y = cbind(Y1,Y2,Y3)
      YW = cbind(Y,W)
      Yobs = apply(YW, 1, function(x) x[1:3][x[4]])
    }


    if (all_confounder == F) {
      tau1 = .2; tau2 = 0.4; tau3= 0.5
      Yp1 = expit(tau1+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5  -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2  - 0.5*x2*x5  + .23*x4^3 )
      Yp2 = expit(tau2+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2  - 0.5*x2*x5  + .23*x4^3 )
      Yp3 = expit(tau3+ .2*x1 + .3*x2 - .4*x3 + .2*x4 - .2*x5 -.6*x1^2 - 0.7*x2^2 - 0.6*x3^2  - 0.5*x2*x5  + .23*x4^3 )

      #potential outcomes
      Y1 = Y2 = Y3= NULL
      for (i in 1:n) {
        Y1[i] = stats::rbinom(1, 1, Yp1[i])
        Y2[i] = stats::rbinom(1, 1, Yp2[i])
        Y3[i] = stats::rbinom(1, 1, Yp3[i])}

      # observed outcomes
      Y = cbind(Y1,Y2,Y3)
      YW = cbind(Y,W)
      Yobs = apply(YW, 1, function(x) x[1:3][x[4]])
    }
  }
  #need the following two lines if run GBM; ignore if run BART
  colnames(trtdat)<- c("trt_ind", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
  trtdat$trt_ind<- as.factor(trtdat$trt_ind)


  ###########################################################
  # true treatment effects
  # ATE(1,2), ATE(1,3), ATE(2,3)
  ATE12_RD = mean(Y[,1]) - mean(Y[,2])
  ATE13_RD = mean(Y[,1]) - mean(Y[,3])
  ATE23_RD = mean(Y[,2]) - mean(Y[,3])
  ATE12_RR = mean(Y[,1]) / mean(Y[,2])
  ATE13_RR = mean(Y[,1]) / mean(Y[,3])
  ATE23_RR = mean(Y[,2]) / mean(Y[,3])
  ATE12_OR = mean(Y[,1])/(1-mean(Y[,1])) / (mean(Y[,2])/(1-mean(Y[,2])))
  ATE13_OR = mean(Y[,1])/(1-mean(Y[,1])) / (mean(Y[,3])/(1-mean(Y[,3])))
  ATE23_OR = mean(Y[,2])/(1-mean(Y[,2])) / (mean(Y[,3])/(1-mean(Y[,3])))

  # ATT(1,2), ATT(1,3)
  ATT12_RD = mean(Y[W==1,1]) - mean(Y[W==1,2])
  ATT13_RD = mean(Y[W==1,1]) - mean(Y[W==1,3])
  ATT12_RR = mean(Y[W==1,1]) / mean(Y[W==1,2])
  ATT13_RR = mean(Y[W==1,1]) / mean(Y[W==1,3])
  ATT12_OR = mean(Y[W==1,1])/(1-mean(Y[W==1,1])) / (mean(Y[W==1,2])/(1-mean(Y[W==1,2])))
  ATT13_OR = mean(Y[W==1,1])/(1-mean(Y[W==1,1])) / (mean(Y[W==1,3])/(1-mean(Y[W==1,3])))

  Ests = list(ATE12_RD=ATE12_RD, ATE13_RD=ATE13_RD, ATE23_RD=ATE23_RD,ATE12_RR=ATE12_RR,ATE13_RR=ATE13_RR, ATE23_RR=ATE23_RR, ATE12_OR=ATE12_OR, ATE13_OR=ATE13_OR, ATE23_OR=ATE23_OR,
              ATT12_RD=ATT12_RD, ATT13_RD=ATT13_RD, ATT12_RR=ATT12_RR, ATT13_RR=ATT13_RR, ATT12_OR=ATT12_OR, ATT13_OR=ATT13_OR)

  return(list(trtdat = trtdat, Y=Y, Yobs=Yobs, n=n, Ests=Ests))
}
