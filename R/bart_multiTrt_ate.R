#######################################################################
# Functions to estimate causal effects of multiple treatment using BART
# Assume 3 treatment
# Author: Liangyuan Hu, Chenyang Gu
#######################################################################


# Estimate average treatment effect (ATE)
bart_multiTrt_ate = function(y, x, trt, k=2, discard = "No", ntree=100, ndpost=1000, nskip=1000) {

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


    # Predict potential outcomes for trt=2
    xp2 = xt[trt==2,]
    xp1 = xp2
    xp3 = xp2
    xp1[,1] = 1  # switch treatment label 2 to 1
    xp3[,1] = 3  # switch treatment label 2 to 3

    bart_pred21 = BART::pwbart(xp1, bart_mod$treedraws)
    bart_pred22 = BART::pwbart(xp2, bart_mod$treedraws)
    bart_pred23 = BART::pwbart(xp3, bart_mod$treedraws)

    pred_prop21 = stats::pnorm(bart_pred21)
    pred_prop22 = stats::pnorm(bart_pred22)
    pred_prop23 = stats::pnorm(bart_pred23)

    # Predict potential outcomes for trt=3
    xp3 = xt[trt==3,]
    xp1 = xp3
    xp2 = xp3
    xp1[,1] = 1  # switch treatment label 3 to 1
    xp2[,1] = 2  # switch treatment label 3 to 2

    bart_pred31 = BART::pwbart(xp1, bart_mod$treedraws)
    bart_pred32 = BART::pwbart(xp2, bart_mod$treedraws)
    bart_pred33 = BART::pwbart(xp3, bart_mod$treedraws)

    pred_prop31 = stats::pnorm(bart_pred31)
    pred_prop32 = stats::pnorm(bart_pred32)
    pred_prop33 = stats::pnorm(bart_pred33)

    if (discard == "No") {
      pred_prop11 = stats::pnorm(bart_pred11)
      pred_prop12 = stats::pnorm(bart_pred12)
      pred_prop13 = stats::pnorm(bart_pred13)

      pred_prop21 = stats::pnorm(bart_pred21)
      pred_prop22 = stats::pnorm(bart_pred22)
      pred_prop23 = stats::pnorm(bart_pred23)

      pred_prop31 = stats::pnorm(bart_pred31)
      pred_prop32 = stats::pnorm(bart_pred32)
      pred_prop33 = stats::pnorm(bart_pred33)
    } else if(discard == "Stringent"){
      #posterior standard deviation of the predicted outcome among those treated with W=1
      post.ind.sd11 = apply(pred_prop11, 2, stats::sd)
      post.ind.sd12 = apply(pred_prop12, 2, stats::sd)
      post.ind.sd13 = apply(pred_prop13, 2, stats::sd)


      post.ind.sd21 = apply(pred_prop21, 2, stats::sd)
      post.ind.sd22 = apply(pred_prop22, 2, stats::sd)
      post.ind.sd23 = apply(pred_prop23, 2, stats::sd)


      post.ind.sd31 = apply(pred_prop31, 2, stats::sd)
      post.ind.sd32 = apply(pred_prop32, 2, stats::sd)
      post.ind.sd33 = apply(pred_prop33, 2, stats::sd)


      eligible1 =  (post.ind.sd12 <= 2*post.ind.sd11) & (post.ind.sd13 <= 2*post.ind.sd11)
      eligible2 =  (post.ind.sd21 <= 2*post.ind.sd22) & (post.ind.sd23 <= 2*post.ind.sd22)
      eligible3 =  (post.ind.sd31 <= 2*post.ind.sd33) & (post.ind.sd32 <= 2*post.ind.sd33)
      n_1_discard <- sum(eligible1 == F)
      n_2_discard <- sum(eligible2 == F)
      n_3_discard <- sum(eligible3 == F)

      #these units are within the common support region
      pred_prop11 = pred_prop11[,eligible1]
      pred_prop12 = pred_prop12[,eligible1]
      pred_prop13 = pred_prop13[,eligible1]

      pred_prop21 = pred_prop21[,eligible2]
      pred_prop22 = pred_prop22[,eligible2]
      pred_prop23 = pred_prop23[,eligible2]

      pred_prop31 = pred_prop31[,eligible3]
      pred_prop32 = pred_prop32[,eligible3]
      pred_prop33 = pred_prop33[,eligible3]
    } else if(discard == "Lenient"){
      #posterior standard deviation of the predicted outcome among those treated with W=1
      post.ind.sd11 = apply(pred_prop11, 2, stats::sd)
      post.ind.sd12 = apply(pred_prop12, 2, stats::sd)
      post.ind.sd13 = apply(pred_prop13, 2, stats::sd)
      #discard unit i with W_i =1 if posterior sd of his/her counterfactual outcomes exceeeds threshold
      threshold1 = max(post.ind.sd11)

      post.ind.sd21 = apply(pred_prop21, 2, stats::sd)
      post.ind.sd22 = apply(pred_prop22, 2, stats::sd)
      post.ind.sd23 = apply(pred_prop23, 2, stats::sd)
      threshold2 = max(post.ind.sd22)

      post.ind.sd31 = apply(pred_prop31, 2, stats::sd)
      post.ind.sd32 = apply(pred_prop32, 2, stats::sd)
      post.ind.sd33 = apply(pred_prop33, 2, stats::sd)
      threshold3 = max(post.ind.sd33)

      eligible1 = post.ind.sd12 <= threshold1 & post.ind.sd13 <= threshold1
      eligible2 = post.ind.sd21 <= threshold2 & post.ind.sd23 <= threshold2
      eligible3 = post.ind.sd31 <= threshold3 & post.ind.sd32 <= threshold3
      n_1_discard <- sum(eligible1 == F)
      n_2_discard <- sum(eligible2 == F)
      n_3_discard <- sum(eligible3 == F)

      #these units are within the common support region
      pred_prop11 = pred_prop11[,eligible1]
      pred_prop12 = pred_prop12[,eligible1]
      pred_prop13 = pred_prop13[,eligible1]

      pred_prop21 = pred_prop21[,eligible2]
      pred_prop22 = pred_prop22[,eligible2]
      pred_prop23 = pred_prop23[,eligible2]

      pred_prop31 = pred_prop31[,eligible3]
      pred_prop32 = pred_prop32[,eligible3]
      pred_prop33 = pred_prop33[,eligible3]
    }

    # Estimate causal effects
    RD12_est = RR12_est = OR12_est = NULL
    RD13_est = RR13_est = OR13_est = NULL
    RD23_est = RR23_est = OR23_est = NULL

    for (m in 1:ndpost) {

        # Estimate E(Y1), E(Y2), E(Y3)
        y1 = c(stats::rbinom(n1, 1, pred_prop11[m,]), stats::rbinom(n2, 1, pred_prop21[m,]), stats::rbinom(n3, 1, pred_prop31[m,]))
        y2 = c(stats::rbinom(n1, 1, pred_prop12[m,]), stats::rbinom(n2, 1, pred_prop22[m,]), stats::rbinom(n3, 1, pred_prop32[m,]))
        y3 = c(stats::rbinom(n1, 1, pred_prop13[m,]), stats::rbinom(n2, 1, pred_prop23[m,]), stats::rbinom(n3, 1, pred_prop33[m,]))

        y1_pred = mean(y1)
        y2_pred = mean(y2)
        y3_pred = mean(y3)

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


































