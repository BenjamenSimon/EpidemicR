#' Calculate the product part of the infection likelihood (Chris' code)
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates 
#'  \eqn{\prod_{j != kappa}^{n_I} [ \sum_{ i in I_{n_{j-}} } [ beta_{i,j} ] ]}.
#' 
#' @param t_inf_j A vector of the infection times of all individuals (Inf if not infected), ordered by ID.
#' @param events A 2 column matrix where the first column is the infection times, and the second is the paired removal times.
#' @param B The infection rate matrix.
#'
#' @keywords Product infection GSE likelihood
#' @export
#'
#' @return Returns the value of the product detailed above.
#'
#' @examples
#' This function is utilised by the chris_log_likelihood function.

# == Components ==
# is_infected = 0/1 vector which says if each individual became infected or not
# waifw = "who acquired infection from whom" matrix. True false matrix where the value of element (i,j) is 1 if $i$ could have infected $j$ and 0 otherwise. Susceptibles just have vector of 0s.
# lambda_j = vector of the sums of the infectious pressure exererted on each individual who became infected
#          = in other words it is the sum of the beta_ij for i in (the set of infected individuals).
# I0 = the index of the initial infective

chris_prod_part <- function(t_inf_j, events, B) {
  
  # The indexes of those who were infected.
  is_infected <- t_inf_j < Inf
  
  # Calculate a matrix of `who acquired infection from whom`
  # It details whether the individual in row i could have infected the individuals in column j (given they were infected)
  waifw <- sapply(t_inf_j, function(t) events[,1] < t & t < events[,2])
  
  # Calculate the value of lambda_j for each infected individual
  lambdaj <- colSums(B[,is_infected] * waifw[, is_infected])
  
  # Which individual is the initial infected
  I0 <- which.min(t_inf_j[is_infected])
  
  # Calculate the final value
  sum(log(lambdaj[-I0]))
}