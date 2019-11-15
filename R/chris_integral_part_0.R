#' Calculate the integral part of the infection likelihood (Chris' code).
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the integral part
#'  of the infection likelihood, using the alternative double-sum parameterisation.
#'
#' @param t_inf_j A vector of the infection times of all individuals (Inf if not infected), ordered by ID.
#' @param events A 2 column matrix where the first column is the infection times, and the second is the paired removal times.
#' @param B The infection rate matrix.
#'
#' @keywords Integral infection GSE likelihood
#' @export
#'
#' @return Returns the value of the integral detailed above.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

chris_integral_part_0 = function(t_inf_j, events, B){
  
  # The indexes of those who were infected.
  i_infected = events[,1] < Inf
  
  # A matrix (calculated using the interval intersect function)
  # which denotes the amount of time that individual j (column)
  # spent in infected indivudual i's (row) infectious period.
  E = chris_interval_intersect_0(events[i_infected,], cbind(0, t_inf_j))
  
  # Pointwise multiply the periods/intersects by the appropriate infection rate.
  integral = E * B[i_infected,]
  
  # Sum up the matrix, equivalent to taking the double sum in i (in infected), j (in all).
  sum(integral)
}