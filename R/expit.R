#' Inverse logit
#'
#' This function inverse the logit function.
#'
#' @param x a vector
#'
#' @return a vector
#' @export
#' @examples
#' library(CIMTx)
#' expit(1:5)
expit <- function(x) {exp(x)/(1+exp(x))}
