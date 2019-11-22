#' Calculate the log-likelihood for an epidemic, up to a constant of proportionality (Simple code).
#'
#' This function is used tp calculate the likelihood of a General Stochastic Epidemic, and brings together all the previously 
#'  specified functions; prod_part, interval_intersect, and integral_part.
#'
#' @param inf.times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem.times A vector of the removal times of all infected individuals.
#' @param beta.mat The infection rate matrix.
#'
#' @keywords likelihood log loglikelihood infection GSE 
#' @export
#'
#' @return Returns the value of the log likelihood, up to a constant of proportionality.
#'
#' @examples
#' This function is utilised in the inference functions.

simple_log_likelihood = function(inf.times, rem.times, beta.mat){
  
  # Calculate the product part
  prod = simple_prod_part(inf.times, rem.times, beta.mat)
  
  # Calculate the integral part
  integral = simple_integral_part(inf.times, rem.times, beta.mat)
  
  # Calculate the log likelikehood
  prod - integral
}


