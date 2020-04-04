# Estimate average treatment effect on the treated (ATT)
# The default reference group is 1st group
regadj_multiTrt_att = function(y, x, trt, ndpost=1000) {

  n1 = sum(trt==1)
  n2 = sum(trt==2)
  n3 = sum(trt==3)
  n = n1 + n2 + n3

  xt = cbind(trt,x)
  mod_data = cbind(y,xt)

  # Fit Bayesian logistic regression
  reg_mod = arm::bayesglm(y ~., data = mod_data, family = stats::binomial(link="logit"), x = TRUE)

  mod_sims = arm::sim(reg_mod, n.sims = ndpost)
  sim_beta = as.matrix(stats::coef(mod_sims))
  x_tilde  = as.data.frame(stats::model.matrix(reg_mod))

  x_tilde11 = x_tilde12 = x_tilde13 = x_tilde[x_tilde$trt==1,]
  x_tilde12$trt = 2
  x_tilde13$trt = 3

  # predictive simulation using the stats::binomial distribution
  # predict potential outcomes
  y11_tilde = array(NA, c(ndpost, n1))
  y12_tilde = array(NA, c(ndpost, n1))
  y13_tilde = array(NA, c(ndpost, n1))

  for (s in 1:ndpost) {
    # predict potential outcome Y(1)
    p11_tilde = arm::invlogit(as.matrix(x_tilde11) %*% sim_beta[s,])
    y11_tilde[s,] = stats::rbinom(n1, 1, p11_tilde)

    # predict potential outcome Y(2)
    p12_tilde = arm::invlogit(as.matrix(x_tilde12) %*% sim_beta[s,])
    y12_tilde[s,] = stats::rbinom(n1, 1, p12_tilde)

    # predict potential outcome Y(3)
    p13_tilde = arm::invlogit(as.matrix(x_tilde13) %*% sim_beta[s,])
    y13_tilde[s,] = stats::rbinom(n1, 1, p13_tilde)
  }


  # Estimate causal effects
  RD12_est = RR12_est = OR12_est = NULL
  RD13_est = RR13_est = OR13_est = NULL

  for (m in 1:ndpost) {

    # Estimate E(Y1|trt=1), E(Y2|trt=1), E(Y3|trt=1)
    y1_pred = mean(y11_tilde[m,])
    y2_pred = mean(y12_tilde[m,])
    y3_pred = mean(y13_tilde[m,])

    # Calculate risk difference (RD)
    RD12_est[m] = y1_pred - y2_pred
    RD13_est[m] = y1_pred - y3_pred

    # Calculate relative risk (RR)
    RR12_est[m] = y1_pred / y2_pred
    RR13_est[m] = y1_pred / y3_pred

    # Calculate odds ratio (OR)
    OR12_est[m] = (y1_pred / (1 - y1_pred)) / (y2_pred / (1 - y2_pred))
    OR13_est[m] = (y1_pred / (1 - y1_pred)) / (y3_pred / (1 - y3_pred))
  }

  att12 = postSumm(RD12_est, RR12_est, OR12_est)
  att13 = postSumm(RD13_est, RR13_est, OR13_est)

  list(ATT12 = round(att12, digits=3),
       ATT13 = round(att13, digits=3))
}
