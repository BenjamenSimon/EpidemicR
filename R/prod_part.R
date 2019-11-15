#' Calculate the product part of the infection likelihood
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates
#'  \eqn{\prod_{j != kappa}^{n_I} [ \sum_{ i in I_{n_{j-}} } [ beta_{i,j} ] ]}.
#'
#' @param inf_times A vector of the infection times of all individuals (Inf if not infected).
#' @param rem_times A vector of the removal times of all infected individuals.
#' @param B The infection rate matrix.
#' @param infected_inds A vector of the infected individuals.
#'
#' @keywords Product infection GSE likelihood
#' @export
#'
#' @return Returns the value of the product detailed above.
#'
#' @examples
#' This function is utilised by the log_likelihood function.

prod_part <- function(inf_times, rem_times, B, infected_inds) {

  # Calculate a matrix of `who acquired infection from whom`
  # It details whether the individual in row i could have infected the individuals in column j (given they were infected)
  waifw <- sapply(inf_times, function(t) inf_times < t & t < rem_times)

  # Calculate the value of lambda_j for each infected individual
  lambdaj <- colSums(B[,infected_inds] * waifw[, infected_inds])

  # Which individual is the initial infected
  I0 <- which.min(inf_times[infected_inds])

  # Calculate the final value
  sum(log(lambdaj[-I0]))
}


