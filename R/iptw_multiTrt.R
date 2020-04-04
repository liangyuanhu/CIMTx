
iptw_multiTrt = function(y, trt, psdat, estimand = "ATE", method, wt1, wt2, wt3,wt12, wt13, trim_alpha,SL.library) {
  if (estimand == "ATE") {
    iptw_est = iptw_multiTrt_ate(y, trt, psdat, wt1, wt2, wt3, method, trim_alpha,SL.library)
  }
  if (estimand == "ATT") {
    iptw_est = iptw_multiTrt_att(y, trt, psdat, wt12, wt13,method, trim_alpha,SL.library)
  }
  return(iptw_est)
}
