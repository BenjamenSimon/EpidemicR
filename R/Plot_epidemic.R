#' Plot a spatial epidemic.
#'
#' This function takes in an epidemic simulation generated by this package and outputs a colour coded plot.
#'  The circles are the remaining susceptibles, the red crosses are the infected individuals, and the blue crossed boxes are the
#'  initial infectives
#'
#' @param xy.coords The coordinates of all individuals in the population
#' @param sim The output of the epidemic simulation
#' @param init.inf A vector of the ID(s) of the initial invective(s). This just defaults to 1 as in the simulation.
#'
#' @keywords plot epidemic colour
#' @export
#'
#' @return A spatial plot of the epidemic
#'
#' @examples
#' # Generate a distance matrix
#'    C <- Dist_mat_unif(N=100, xlim = 20, ylim = 20)
#'    xy.coords <- C[[1]]
#'    distance_mat <- C[[2]]
#' 
#' # Generate an associated infection rate rate matrix
#'    rate_mat <- Beta_mat_form(distance_mat, c(0.004, 0.002), 10)
#'    
#' # Generate a simulated epidemic
#'    Hetero_sim <- GSE_sim(N = 100, beta.mat = rate_mat, gamma = 0.15)
#'    
#' # Plot the epidemic
#'    Plot_epidemic(xy.coords, Hetero_sim, init.inf = c(1))

Plot_epidemic <- function(xy.coords, sim, init.inf = 1){
  
  inf.id <- sim[which(sim[,2] != Inf), 1]
  
  col <- rep("black", nrow(sim))
  col[inf.id] <- "red"
  col[init.inf] <- "blue"
  
  pch <- rep(1, nrow(sim))
  pch[inf.id] <- 4
  pch[init.inf] <- 9
  
  plot(xy.coords[,1], xy.coords[,2], type  = "p", col = col, pch = pch, xlab = "X", ylab = "Y")
}
