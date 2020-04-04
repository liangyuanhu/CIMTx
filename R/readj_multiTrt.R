
regadj_multiTrt = function(y, x, trt, estimand="ATE", ndpost=1000) {
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
    regadj_est = regadj_multiTrt_ate(y, x, trt, ndpost=1000)
  }

  if (estimand=="ATT") {
    regadj_est = regadj_multiTrt_att(y, x, trt, ndpost=1000)
  }

  return(regadj_est)
}
