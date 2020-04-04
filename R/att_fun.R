
att_fun <- function(wt12, wt13,y, trt_ind) {
  mu_1_1_hat_iptw = mean(y[trt_ind == 1])
  mu_1_2_hat_iptw = sum(y[trt_ind == 2] * wt12[trt_ind == 2]) / sum(wt12[trt_ind == 2])
  mu_1_3_hat_iptw = sum(y[trt_ind == 3] * wt13[trt_ind == 3]) / sum(wt13[trt_ind ==  3])
  RD12 = mu_1_1_hat_iptw - mu_1_2_hat_iptw
  RD13  = mu_1_1_hat_iptw - mu_1_3_hat_iptw
  RR12 = mu_1_1_hat_iptw / mu_1_2_hat_iptw
  RR13 = mu_1_1_hat_iptw / mu_1_3_hat_iptw
  OR12 = (mu_1_1_hat_iptw / (1 - mu_1_1_hat_iptw)) / (mu_1_2_hat_iptw / (1 - mu_1_2_hat_iptw))
  OR13 = (mu_1_1_hat_iptw / (1 - mu_1_1_hat_iptw)) / (mu_1_3_hat_iptw / (1 - mu_1_3_hat_iptw))
  res = list(RD12, RD13, RR12, RR13, OR12, OR13)
  return (res)
}
