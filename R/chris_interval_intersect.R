#' Calculate the interval intersect (Chris' code).
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic,
#'  and calculates the interval intersect used in the calculation of the intergral for the infection likelihood.
#'  This is the interval \eqn{[ (R_i ^ I_j) - (I_i ^ I_j) ]} for i in infecteds and j in all.
#'
#' @param interval_i Is the matrix cbind(t_inf, t_rem)[infected, ]. This is a 2 column matrix of the infection
#'                   and removal times for the infected individuals.
#' @param interval_j Is the matrix cbind(t_inf, t_rem). This is a 2 column matrix of the infection
#'                   and removal times for all individuals.
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

chris_interval_intersect = function(interval_i, interval_j){

  # Calculates the minimum of infection time I_j of j and the infection time I_i of
  # infected individual i.
  # This is min(I_j, I_i).
  int_start <- sapply(interval_j[,1], function(x) pmin(x, interval_i[,1]))

  # Calculates the minimum of infection time I_j of individual j and the
  # removal time R_i of infected individual i, for i in infecteds, j in all.
  # This is min(I_j, R_i).
  int_end <- sapply(interval_j[,1], function(x) pmin(x, interval_i[,2]))

  # Returns a matrix where each column is a vector for each j
  int_end - int_start
}
