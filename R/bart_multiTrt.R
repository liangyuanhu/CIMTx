bart_multiTrt = function(y, x, trt, discard = FALSE, estimand="ATE", k=2, ntree=100, ndpost=1000, nskip=1000) {
  x <- x[, -1]
  # Data structure
  #        Y(1) Y(2) Y(3)
  # trt=1   *    ?    ?
  # trt=2   ?    *    ?
  # trt=3   ?    ?    *

  #        Y(1) Y(2) Y(3)
  # trt=1  y11  y12  y13
  # trt=2  y21  y22  y23
  # trt=3  y31  y32  y33

  if (estimand=="ATE") {
    bart_est = bart_multiTrt_ate(y, x, trt, discard = FALSE, k=2, ntree=100, ndpost=1000, nskip=1000)
  }

  if (estimand=="ATT") {
    bart_est = bart_multiTrt_att(y, x, trt, discard = FALSE, k=2, ntree=100, ndpost=1000, nskip=1000)
  }

  return(bart_est)
}
