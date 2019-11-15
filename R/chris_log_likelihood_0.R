#' Calculate the log-likelihood for an epidemic, up to a constant of proportionality (Chris' code).
#'
#' This function is used tp calculate the likelihood of a General Stochastic Epidemic, and brings together all the previously 
#'  specified functions; prod_part, interval_intersect, and integral_part.
#'  It assumes that the epidemic began at time 0.
#'
#' @param t_inf A vector of the infection times of all individuals (Inf if not infected).
#' @param t_rem A vector of the removal times of all infected individuals.
#' @param B The infection rate matrix.
#'
#' @keywords likelihood log loglikelihood infection GSE 
#' @export
#'
#' @return Returns the value of the log likelihood, up to a constant of proportionality.
#'
#' @examples
#' This function is utilised in the inference functions.

chris_log_likelihood_0 = function(t_inf, t_rem, B){
  
  # Calculate the product part
  prod = chris_prod_part(t_inf, cbind(t_inf, t_rem), B)
  
  # Calculate the integral part
  integral = chris_integral_part_0(t_inf, cbind(t_inf, t_rem), B)
  
  # Calculate the log likelikehood
  prod - integral
}





