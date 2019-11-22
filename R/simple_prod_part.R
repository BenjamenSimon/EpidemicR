#' Calculate the product part of the infection likelihood (Simple code)
#'
#' This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates 
#'  \eqn{\prod_{j != kappa}^{n_I} [ \sum_{ i in I_{n_{j-}} } [ beta_{i,j} ] ]}.
#' 
#' @param inf.times A vector of the infection times of all individuals (Inf if not infected), ordered by ID.
#' @param rem.times A vector of the removal times of all individuals (Inf if not infected), ordered by ID.
#' @param beta.mat The infection rate matrix.
#'
#' @keywords Product infection GSE likelihood
#' @export
#'
#' @return Returns the value of the product detailed above.
#'
#' @examples
#' This function is utilised by the simple_log_likelihood function.


simple_prod_part <- function(inf.times, rem.times, beta.mat) {
  
  # Vector of infected individuals
  infecteds <- which(inf.times < Inf)
  
  # Initial infected
  I0 <- which.min(inf.times[infecteds])
  
  # Vector of infected individuals exc. initial infected
  new.infs <- infecteds[-I0]
  
  # Find the `who acquired infection from whom` values
  
  inf.before.j <- list()
  
  for(j in new.infs){
    
    # Which individuals were infected before the jth individual
    infed.before <- which(inf.times < inf.times[j])
    
    # Which individuals were removed before the jth individual was infected
    remed.before <- which(rem.times < inf.times[j])
    
    # Record the vector of individuals who were still infected 
    # just before the jth person became infected
    inf.before.j[[j]] <- infed.before[!(infed.before %in% remed.before)]
  }
  
  
  # Calculate the product
  
  sum_over_infs <- 0
  
  for(j in new.infs){
    
    sum_over_ibj <- 0
    
    for(i in inf.before.j[[j]]){
      # For each infected, sum up the beta_ij values associated with
      # the individuals who could have infected them (lambda_j)
      sum_over_ibj <- sum_over_ibj + beta.mat[i,j]
    }
    
    # Sum for all infected infected individuals (exc. initial infected)
    sum_over_infs <- sum_over_infs + log(sum_over_infs)
  }

  return(sum_over_infs)
}











