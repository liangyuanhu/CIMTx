#===========================#
#Author: Michael Lopez      #
#===========================#
vm_multiTrt_att = function(y, x, trt) {
  x <- x[, -1]
  trt <- factor(trt, levels = c(1, 2, 3),
                    labels = c("Treatment 1", "Treatment 2", "Treatment 3"))
  # estimate generalized propensity scores using multinomial logistic regression
  xydata = cbind(trt, x)


  # WARNING: reorder the dataset by "Treatment 1", "Treatment 2", "Treatment 3"
  # it will affect the identification of matched indices in row 249
  xydata = xydata[order(xydata$trt),]


  # ps model 1
  ps_fit = nnet::multinom(trt ~ ., data = xydata, trace = FALSE)
  probs_logit1 = data.frame(stats::fitted(ps_fit))
  colnames(probs_logit1) = c("p1", "p2", "p3")
  xydata = cbind(xydata, probs_logit1)

  # Determine eligibility
  min_max_Ps <- xydata %>%
    dplyr::group_by(trt) %>%
    dplyr::summarise(min1 = min(p1), max1 = max(p1),
              min2 = min(p2), max2 = max(p2),
              min3 = min(p3), max3 = max(p3))
  #min_max_Ps

  xydata$"Eligible" <-
    xydata$"p1" >= max(min_max_Ps$min1) & xydata$"p1" <= min(min_max_Ps$max1) &
    xydata$"p2" >= max(min_max_Ps$min2) & xydata$"p2" <= min(min_max_Ps$max2) &
    xydata$"p3" >= max(min_max_Ps$min3) & xydata$"p3" <= min(min_max_Ps$max3)
  #table(xydata$Eligible)

  xydata = cbind(y, xydata)
  xydata = dplyr::filter(xydata, Eligible)

  # Calculate new propensity scores for eligible subjects
  ps_fit_E = nnet::multinom(trt ~ ., data = xydata[,-1], trace = FALSE)
  probs_logit1_E = stats::fitted(ps_fit_E)
  colnames(probs_logit1_E) = c("p1", "p2", "p3")
  xydata <- xydata %>%
    dplyr :: select(-"p1", -"p2", -"p3")
  xydata = cbind(xydata, probs_logit1_E)

  n1 = sum(xydata$trt == "Treatment 1")
  n2 = sum(xydata$trt == "Treatment 2")
  n3 = sum(xydata$trt == "Treatment 3")

  ### Vector Matching for ATT for outcome 1 (comp_resp_obs)
  # Stratify car::logit(r(ti, Xi)) using K-means clustering
  clustnum <- 5

  xydata$Quint1 <- 1
  temp1 <- stats::kmeans(car::logit(xydata$"p1"), clustnum)
  xydata$Quint1 <- temp1$cluster

  xydata$Quint2 <- 1
  temp2 <- stats::kmeans(car::logit(xydata$"p2"), clustnum)
  xydata$Quint2 <- temp2$cluster

  xydata$Quint3 <- 1
  temp3 <- stats::kmeans(car::logit(xydata$"p3"), clustnum)
  xydata$Quint3 <- temp3$cluster

  colnames(xydata)[1:2] = c("Y_obs","treat")
  temp12 <- dplyr::filter(xydata, treat != "Treatment 3")
  temp13 <- dplyr::filter(xydata, treat != "Treatment 2")
  temp23 <- dplyr::filter(xydata, treat != "Treatment 1")

  # matching
  # t1 = reference treatment
  match12 <- Matching::Matchby(Y = temp12$Y_obs, Tr = temp12$"treat" == "Treatment 1",
                     X = car::logit(temp12$"p1"), by = temp12$Quint3,
                     caliper = 0.5*stats::sd(car::logit(temp12$"p1")),  replace = T, estimand = "ATT", print.level = 0)

  match13 <- Matching::Matchby(Y = temp13$Y_obs, Tr = temp13$"treat" == "Treatment 1",
                     X = car::logit(temp13$"p1"), by = temp13$Quint2,
                     caliper = 0.5*stats::sd(car::logit(temp13$"p1")), replace = T, estimand = "ATT", print.level = 0)


  # Identify the matched subgroups
  rownames(xydata) = 1:nrow(xydata)
  xydata$id = 1:nrow(xydata)
  xydata$both_1 <- xydata$id %in% match12$index.treated & xydata$id %in% match13$index.treated
  temp = xydata[xydata$both_1 == "TRUE", ]
  m12 = cbind(match12$index.treated, match12$index.control)
  #===============================================================================================================#
  # WARNING: reorder the dataset by "Treatment 1", "Treatment 2", "Treatment 3" before running the following lines
  m13 = cbind(match13$index.treated, match13$index.control + sum(xydata$"treat" == "Treatment 2"))
  m12 = m12[m12[,1] %in% rownames(temp), ]
  m13 = m13[m13[,1] %in% rownames(temp), ]
  triplets = cbind(m12[order(m12[,1]), ], m13[order(m13[,1]), ])
  triplets = as.matrix(triplets[,c(1, 2, 4)])
  n_trip = nrow(triplets)

  # Matching Estimator
  # For subjects receiving reference treatment
  Y1_imp = xydata$Y_obs[triplets[,1]] # observed outcomes for Treatment 1
  Y2_imp = xydata$Y_obs[triplets[,2]] # imputed outcomes for Treatment 2
  Y3_imp = xydata$Y_obs[triplets[,3]] # imputed outcomes for Treatment 3

  y1_hat = mean(Y1_imp)
  y2_hat = mean(Y2_imp)
  y3_hat = mean(Y3_imp)

  # Calculate risk difference (RD)
  RD12_est = y1_hat - y2_hat
  RD13_est = y1_hat - y3_hat

  # Calculate relative risk (RR)
  RR12_est = y1_hat / y2_hat
  RR13_est = y1_hat / y3_hat

  # Calculate  odds ratio (OR)
  OR12_est = (y1_hat / (1 - y1_hat)) / (y2_hat / (1 - y2_hat))
  OR13_est = (y1_hat / (1 - y1_hat)) / (y3_hat / (1 - y3_hat))

  # list(RD = c(RD12_est, RD13_est),
  #      RR = c(RR12_est, RR13_est),
  #      OR = c(OR12_est, OR13_est))
  res12 = rbind(RD = c(RD12_est), RR = RR12_est, OR = OR12_est)
  res13 = rbind(RD = c(RD13_est), RR = RR13_est, OR = OR13_est)
  colnames(res12) <- "EST"
  colnames(res13) <- "EST"
  list(ATT12 = res12,
       ATT13 = res13)
}

# comp_resp_att = vm_multiTrt_att_updated(y=comp_resp_obs, x=all_vars1, trt=trt1)
