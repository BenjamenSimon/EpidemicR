#' Calculate the integral part of the infection likelihood (Simple code).
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates the integral part
#'  of the infection likelihood, using the alternative double-sum parameterisation.
#'
#' @param I A vector of the infection times of all individuals (Inf if not infected), ordered by ID.
#' @param R A vector of the infection times of all individuals (Inf if not infected), ordered by ID.
#' @param beta.mat The infection rate matrix.
#'
#' @keywords Integral infection GSE likelihood
#' @export
#'
#' @return Returns the value of the integral detailed above.
#'
#' @examples
#' This function is utilised by the simple_log_likelihood function.

simple_integral_part = function(I = inf.times, R = rem.times, beta.mat = beta.mat){
  
  # Vector of infected individuals
  infecteds <- which(I < Inf)
  
  integral = 0
  
  n.I = sum(I < Inf)
  N = length(I)
  
  for(i in infecteds){
    for(j in 1:N){
      # For each i in the infected individuals
      # and each j in all the individuals
      # calculate the value of beta_ij multiplied by the difference of minimums
      integral <- integral + beta.mat[i,j] * ( min(c(R[i], I[j])) - min(c(I[i], I[j])) )
    }
  }
  
  return(integral)
}