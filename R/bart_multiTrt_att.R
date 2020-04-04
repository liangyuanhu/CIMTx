# Estimate average treatment effect on the treated (ATT)
# The default reference group is 1st group
# Author: Liangyuan Hu, Chenyang Gu
bart_multiTrt_att = function(y, x, trt, k=2, discard = "No", ntree=100, ndpost=1000, nskip=1000) {

  n1 = sum(trt==1)
  n2 = sum(trt==2)
  n3 = sum(trt==3)

  xt = cbind(trt,x)

  # Fit BART
  bart_mod = BART::pbart(x.train = xt, y.train = y, k = k, ntree = ntree, ndpost = ndpost, nskip = nskip)


  # Predict potential outcomes for trt=1
  xp1 = xt[trt==1,]
  xp2 = xp1
  xp3 = xp1
  xp2[,1] = 2  # switch treatment label 1 to 2
  xp3[,1] = 3  # switch treatment label 1 to 3

  bart_pred11 = BART::pwbart(xp1, bart_mod$treedraws)
  bart_pred12 = BART::pwbart(xp2, bart_mod$treedraws)
  bart_pred13 = BART::pwbart(xp3, bart_mod$treedraws)

  pred_prop11 = stats::pnorm(bart_pred11)
  pred_prop12 = stats::pnorm(bart_pred12)
  pred_prop13 = stats::pnorm(bart_pred13)

  if (discard == "No") {
    pred_prop11 = stats::pnorm(bart_pred11)
    pred_prop12 = stats::pnorm(bart_pred12)
    pred_prop13 = stats::pnorm(bart_pred13)
  } else if (discard == "Stringent"){

    #posterior standard deviation of the predicted outcome among those treated with W=1
    post.ind.sd1 = apply(pred_prop11, 2, stats::sd)

    #discard unit i with W_i =1 if posterior sd of his/her counterfactual outcomes exceeeds threshold
    post.ind.sd2 = apply(pred_prop12, 2, stats::sd)
    post.ind.sd3 = apply(pred_prop13, 2, stats::sd)

    #discard unit i if sd of counterfactual outcomes is more than 2 times sd of fitted outcomes
    eligible = (post.ind.sd2 <= 2*post.ind.sd1) & (post.ind.sd3 <= 2*post.ind.sd1)
    n_discard_att = sum(eligible == F)
    #n_discard_att
    #these units are within the common support region
    pred_prop11 = pred_prop11[,eligible]
    pred_prop12 = pred_prop12[,eligible]
    pred_prop13 = pred_prop13[,eligible]
  } else if (discard == "Lenient"){
    post.ind.sd1 = apply(pred_prop11, 2, stats::sd)
    threshold = max(post.ind.sd1)
    #discard unit i with W_i =1 if posterior sd of his/her counterfactual outcomes exceeeds threshold
    post.ind.sd2 = apply(pred_prop12, 2, stats::sd)
    post.ind.sd3 = apply(pred_prop13, 2, stats::sd)
    eligible = (post.ind.sd2 <= threshold) & (post.ind.sd3 <= threshold)
    #discard unit i if sd of counterfactual outcomes is more than 2 times sd of fitted outcomes
    n_discard_att = sum(eligible == F)
    #n_discard_att
    #these units are within the common support region
    pred_prop11 = pred_prop11[,eligible]
    pred_prop12 = pred_prop12[,eligible]
    pred_prop13 = pred_prop13[,eligible]

  }
  # Estimate causal effects
  RD12_est = RR12_est = OR12_est = NULL
  RD13_est = RR13_est = OR13_est = NULL

  for (m in 1:ndpost) {

    # Estimate E(Y1|trt=1), E(Y2|trt=1), E(Y3|trt=1)
    y1_pred = mean(stats::rbinom(n1, 1, pred_prop11))
    y2_pred = mean(stats::rbinom(n1, 1, pred_prop12))
    y3_pred = mean(stats::rbinom(n1, 1, pred_prop13))

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
