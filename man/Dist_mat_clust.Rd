% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dist_mat_clust.R
\name{Dist_mat_clust}
\alias{Dist_mat_clust}
\title{Clustered population distribution.}
\usage{
Dist_mat_clust(N, xlim = 20, ylim = 20, centres = 5, spread = 1)
}
\arguments{
\item{N}{The total number of individuals in the population.}

\item{xlim}{The width of the plane on which individuals will be generated (defaults to 20 units wide).}

\item{ylim}{The height of the place on which individuals will be generated (defaults to 20 units high).}

\item{centres}{The number of centres (clusters). N must be divisable by the number of centres.}

\item{spread}{The variance of the isotropic Gaussian distribution around each cluster (defaults to identity).}
}
\value{
The function returns a list. The first object in the list is a Nx2 matrix of the coordinates of the individuals.
  The second object is an NxN distance matrix.
}
\description{
This function generates a distance matrix for a clustered population, across a specified rectangular plane.
}
\details{
The function designates a number of individuals specified to be clusters, around each it then generates an appropriate
 number of "children" using an isotrpoic distribution with variance as specificed (defaults to identity). 
 
 Attempts to ensure no points have negative coordinates by not allowing centres to have coordinate values below two.
}
\examples{
xy.coords <- Dist_mat_clust(N=100, xlim = 20, ylim = 20, centres = 5, spread = 1)[[1]]
distance_mat <- Dist_mat_clust(N=100, xlim = 20, ylim = 20, centres = 5, spread = 1)[[2]]
plot(xy.coords[,1], xy.coords[,2], type  = "p")

}
\keyword{Clustered}
\keyword{distance}
\keyword{distribution}
\keyword{matrix}
\keyword{population}
