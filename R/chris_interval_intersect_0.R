#' Calculate the interval intersect (Chris' code).
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the interval
#'  intersect used in the calculation of the intergral for the infection likelihood. 
#'  This is the length of the interval for which (I_i, R_i) (i in infecteds) and
#'  (0, I_j) (j in all) overlap.
#'  This is the interval \eqn{max{0, [ min(R_i, I_j) - max(I_i, 0) ]}} for i in infecteds and j in all.
#'  This is an equivalent calculation, but only in the case where the epidemic is assumed 
#'  to start at time 0. It can be thought of as calculating the amount of individual j's susceptible
#'  period (0, I_j) that is spent in individual i's infectious period (I_i, R_i).
#'
#' @param interval_i Is the matrix cbind(t_inf, t_rem)[infected, ]. This is a 2 column matrix of the infection
#'                   and removal times for the infected individuals.
#' @param interval_j Is the matrix cbind(0, t_inf). This is a 2 column matrix of 0s in the first column and
#'                   the infection times in the second.
#'
#' @keywords interval intersection infection GSE likelihood
#' @export
#'
#' @return Returns a matrix where each column is a vector for each j.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

# == Functions ==
# "pmax()" calculates the "parallel" maximum of two vectors ie. pmax(c(1,2,3), c(3,2,1)) = c(3,2,3)
# similarly for "pmin()"

chris_interval_intersect_0 = function(interval_i, interval_j){
  
  # Calculate the maximum of 0 and the infection time of infected individual i.
  # This is max(0, I_i).
  int_start <- sapply(interval_j[,1], function(x) pmax(x, interval_i[,1]))
  
  # Calculate the minimum of the removal time of infected individual i,
  # and the infection time of individual j.
  # This is min(R_i, I_j).
  int_end <- sapply(interval_j[,2], function(x) pmin(x, interval_i[,2]))
  
  # Calculate the maximum of 0 and the difference between these two values.
  # This is max{0, [ min(R_i, I_j) - max(I_i, 0) ]}
  pmax(int_end - int_start, 0)
}