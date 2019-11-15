#' Calculate the integral part of the infection likelihood.
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the integral part
#'  of the infection likelihood, using the alternative double-sum parameterisation.
#'
#' @param inf_times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem_times A vector of the removal times of all infected individuals.
#' @param B The infection rate matrix.
#' @param infected_inds A vector of the infected individuals.
#'
#' @keywords Integral infection GSE likelihood
#' @export
#'
#' @return Returns the value of the integral detailed above.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

integral_part = function(inf_times, rem_times, B, infected_inds){

  # A matrix (calculated using the interval intersect function)
  # which denotes the amount of time that individual j (column)
  # spent in infected indivudual i's (row) infectious period.
  E = interval_intersect(inf_times, rem_times, infected_inds)

  # Pointwise multiply the periods/intersects by the appropriate infection rate.
  integral = E * B[infected_inds,]

  # Sum up the matrix, equivalent to taking the double sum in i (in infected), j (in all).
  sum(integral)
}





