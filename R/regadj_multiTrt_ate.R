# Estimate average treatment effect (ATE)
regadj_multiTrt_ate = function(y, x, trt, ndpost=1000) {

  n1 = sum(trt==1)
  n2 = sum(trt==2)
  n3 = sum(trt==3)
  n = n1 + n2 + n3

  xt = cbind(trt,x)
  mod_data = cbind(y,xt)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~ ., data = mod_data, family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = stats::model.matrix(reg_mod)

  x_tilde1 = x_tilde2 = x_tilde3 = as.data.frame(x_tilde)
  x_tilde1$trt = 1
  x_tilde2$trt = 2
  x_tilde3$trt = 3

  # predictive simulation using the stats::binomial distribution
  # predict potential outcomes
  y1_tilde = array(NA, c(ndpost, n))
  y2_tilde = array(NA, c(ndpost, n))
  y3_tilde = array(NA, c(ndpost, n))

  for (s in 1:ndpost) {
    # predict potential outcome Y(1)
    p1_tilde = arm::invlogit(as.matrix(x_tilde1) %*% sim_beta[s,])
    y1_tilde[s,] = stats::rbinom(n, 1, p1_tilde)

    # predict potential outcome Y(2)
    p2_tilde = arm::invlogit(as.matrix(x_tilde2) %*% sim_beta[s,])
    y2_tilde[s,] = stats::rbinom(n, 1, p2_tilde)

    # predict potential outcome Y(3)
    p3_tilde = arm::invlogit(as.matrix(x_tilde3) %*% sim_beta[s,])
    y3_tilde[s,] = stats::rbinom(n, 1, p3_tilde)
  }


  # Estimate causal effects
  RD12_est = RR12_est = OR12_est = NULL
  RD13_est = RR13_est = OR13_est = NULL
  RD23_est = RR23_est = OR23_est = NULL

  for (m in 1:ndpost) {
    # Estimate E(Y1), E(Y2), E(Y3)
    y1_pred = mean(y1_tilde[m,])
    y2_pred = mean(y2_tilde[m,])
    y3_pred = mean(y3_tilde[m,])

    # Calculate risk difference (RD)
    RD12_est[m] = y1_pred - y2_pred
    RD13_est[m] = y1_pred - y3_pred
    RD23_est[m] = y2_pred - y3_pred

    # Calculate relative risk (RR)
    RR12_est[m] = y1_pred / y2_pred
    RR13_est[m] = y1_pred / y3_pred
    RR23_est[m] = y2_pred / y3_pred

    # Calculate  odds ratio (OR)
    OR12_est[m] = (y1_pred / (1 - y1_pred)) / (y2_pred / (1 - y2_pred))
    OR13_est[m] = (y1_pred / (1 - y1_pred)) / (y3_pred / (1 - y3_pred))
    OR23_est[m] = (y2_pred / (1 - y2_pred)) / (y3_pred / (1 - y3_pred))
  }

  ate12 = postSumm(RD12_est, RR12_est, OR12_est)
  ate13 = postSumm(RD13_est, RR13_est, OR13_est)
  ate23 = postSumm(RD23_est, RR23_est, OR23_est)

  list(ATE12 = round(ate12, digits=3),
       ATE13 = round(ate13, digits=3),
       ATE23 = round(ate23, digits=3))
}
