#============================================================================#
# simulation design - part 1                                                 #
# sample size 1200, 4000, 11600 with ratio of units 1:1:1, 1:5:4 and 1:15:13 #
# number of confounders 10                                                   #
# Author: Liangyuan Hu                                                       #
# March 24, 2019                                                             #
#============================================================================#
#' Data generation function for scenario 1
#'
#' This function generates data to test different causal inference methods for scenario 1.  Please use our main function data_gen.R
#'
#' @param n total number of units for simulation
#' @param ratio ratio of units in the treatment groups
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
#' data_gen_p1(n = 116, ratio = 3,all_confounder=FALSE)
data_gen_p1 = function(n = 11600, ratio = 3,all_confounder=FALSE ) {
  p = 10
  #treatment assignment model
  if (ratio ==1) { alpha1=0.32; alpha2=0.28}
  if (ratio ==2) { alpha1 = -1.3; alpha2 = 0.7}
  if (ratio ==3) {alpha1 = -2.57; alpha2 = 0.52}

  #set.seed(3242019)
  X = matrix(stats::rnorm(p*n), nrow=n, ncol=p/2)
  x1 = X[,1]; x2 = X[,2]; x3 = X[,3]; x4 = X[,4]; x5= X[,5]
  C = matrix(sample(0:2, n*p,replace = T, prob = c(.3,.3,.4)), nrow=n,ncol=p/2)
  x6 = C[,1]; x7 = C[,2]; x8 = C[,3]; x9 = C[,4]; x10 = C[,5];

  ex1 = exp(alpha1+.2*x1+.4*x2+.3*x3+.4*x4+.1*x5+.2*x6 -.8*x7 - 1.1*x8+.5*x9+0.5*x10
            + .2*x1^2 + .3*x2^2 + .4*x4^2 + .8*x1*x2 + .5*x1*x6 + .4*x4*x8)
  ex2 = exp(alpha2+.5*x1+.8*x2+.6*x3+.2*x4+.25*x5+.4*x6 -1.2*x7 - 1.5*x8 -.3*x9+1.5*x10
            + .2*x1^2 + 0.7*x2^2 + .2*x4^2 + .25*x1*x2 + .3*x1*x6 + .6*x4*x8)

  Wp1 = ex1 / (1 + ex1 + ex2)
  Wp2 = ex2 / (1 + ex1 + ex2)
  Wp3 = 1 - Wp1 - Wp2

  W = NULL
  for (i in 1:n) W[i] = sample(c(1,2,3), size=1, replace=T, prob=c(Wp1[i],Wp2[i],Wp3[i]))
  # table(W)[2]/table(W)[1]
  # table(W)[3]/table(W)[1]


  trtdat<- data.frame(W, X, C)
  colnames(trtdat)<- c("trt_ind", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
  trtdat$trt_ind<- as.factor(trtdat$trt_ind)
  # require(nnet)
  # summary( multinom(W ~ ., data = trtdat))

  #parallel response surface model
  # res complication event rate = 33.3%; 30.1%, 33.6% 33.3% by trt group
  # prolonged los event rate =  8.22%
  # icu stay event rate = 67.9%
  # re admission event rate = 9.01%
  #set.seed(3242019)
  tau1 = 0.08; tau2 = 0.3; tau3= 0.5
  Yp1 = expit(tau1+0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
              +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 -0.1*x10^2)
  Yp2 = expit(tau2+0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
              +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3 -0.1*x10^2)
  Yp3 = expit(tau3+0.2*x1 + 0.3*x2 - 0.8*x3 - 0.5*x4 + 0.2*x5 - 1.2*x6 + 0.5*x7 + 0.5*x8 + 0.2*x9 - .2*x10
              +.5*x1^2 - 1.1*x2^2 - 1.5*x3^2 + .1*x1*x7 - 1.5*x2*x8 + .2*x6*x9 + .15*x4^3-0.1*x10^2)

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
  #table(Yobs)[2]/(sum(table(Yobs)))

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


# set.seed(3242019)
# mydata1 = data_gen(n = 4e5, ratio =1)
# save(mydata1, file="sim_p1_scen1.RData")
# mydata1 = data_gen(n = 4e5, ratio =2)
# save(mydata1, file="sim_p1_scen2.RData")
# mydata1 = data_gen(n = 4e5, ratio = 3)
# save(mydata1, file="sim_p1_scen3.RData")

