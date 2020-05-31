#' Data generation function
#'
#'This function generates data to test different causal inference methods.
#' @param n total number of units for simulation
#' @param scenario simulation scenario 1 or scenario 2
#' @param ratio ratio of units in the treatment groups
#' @param overlap levels of covariate overlap: Please select: weak, strong, moderate
#' @param all_confounder TRUE or FALSE. overlap is lacking for a variable that is not predictive of the outcome (all_confounder equals to TRUE) or situations when it is lacking for a true confounder (all_confounder equals to FALSE)
#'
#' @return  list with the 5 elements. Nested within each list, it contains
#' \item{n:}{Number of units for simulation}
#' \item{trt_ind:}{A data frame with number of rows equals to n and 11 columns}
#' \item{Y:}{Observed binary outcome for 3 treatments}
#' \item{Yobs:}{Observed binary outcome}
#' \item{Est:}{True ATE/ATT for RD/RR/OR}
#' @export data_gen
#'
#' @examples
#' library(CIMTx)
#' set.seed(3242019)
#' idata = data_gen(n = 120, ratio =1,scenario = 1)
data_gen <- function(n, scenario, ratio, overlap, all_confounder){
  if (scenario == 1) {
    data_gen_result <- data_gen_p1(n, ratio,all_confounder=FALSE)
  }
  if (scenario == 2) {
    data_gen_result <- data_gen_p1(n, overlap, all_confounder)
  }
  return(data_gen_result)
}
