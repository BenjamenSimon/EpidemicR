#' Calculate the interval intersect.
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the interval
#'  intersect used in the calculation of the intergral for the infection likelihood.
#'  This is the interval \eqn{[ min(R_i, I_j) - min(I_i, I_j) ]} for i in infecteds, j in all.
#'
#' @param inf_times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem_times A vector of the removal times of all infected individuals.
#' @param infected_inds A vector of the infected individuals.
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

interval_intersect = function(inf_times, rem_times, infected_inds){

  # Calculates the minimum of infection time I_j of j and the infection time I_i of
  # infected individual i, i in infecteds, j in all.
  # This is min(I_j, I_i).
  int_start <- sapply(inf_times, function(x) pmin(x, inf_times[infected_inds, ]))

  # Calculates the minimum of infection time I_j of individual j and the
  # removal time R_i of infected individual i, for i in infecteds, j in all.
  # This is min(I_j, R_i).
  int_end <- sapply(inf_times, function(x) pmin(x, rem_times[infected_inds, ]))

  # Returns a matrix where each column is a vector for each j
  int_end - int_start
}






