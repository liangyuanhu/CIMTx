
ate_fun <- function(wt1,
                    wt2,
                    wt3,
                    y,
                    trt_ind) {
  mu_1_hat_iptw = sum(y[trt_ind == 1] * wt1[trt_ind == 1]) / sum(wt1[trt_ind == 1])
  mu_2_hat_iptw = sum(y[trt_ind == 2] * wt2[trt_ind == 2]) / sum(wt2[trt_ind == 2])
  mu_3_hat_iptw = sum(y[trt_ind == 3] * wt3[trt_ind == 3]) / sum(wt3[trt_ind == 3])
  RD12 = mu_1_hat_iptw - mu_2_hat_iptw
  RD13 = mu_1_hat_iptw - mu_3_hat_iptw
  RD23 = mu_2_hat_iptw - mu_3_hat_iptw
  RR12 = mu_1_hat_iptw / mu_2_hat_iptw
  RR13 = mu_1_hat_iptw / mu_3_hat_iptw
  RR23 = mu_2_hat_iptw / mu_3_hat_iptw
  OR12 = (mu_1_hat_iptw / (1 - mu_1_hat_iptw)) / (mu_2_hat_iptw / (1 - mu_2_hat_iptw))
  OR13 = (mu_1_hat_iptw / (1 - mu_1_hat_iptw)) / (mu_3_hat_iptw / (1 - mu_3_hat_iptw))
  OR23 = (mu_2_hat_iptw / (1 - mu_2_hat_iptw)) / (mu_3_hat_iptw / (1 - mu_3_hat_iptw))
  res = list(RD12, RD13, RD23, RR12, RR13, RR23, OR12, OR13, OR23)
  return (res)
}
