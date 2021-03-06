% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chris_prod_part.R
\name{chris_prod_part}
\alias{chris_prod_part}
\title{Calculate the product part of the infection likelihood (Chris' code)}
\usage{
chris_prod_part(t_inf_j, events, B)
}
\arguments{
\item{t_inf_j}{A vector of the infection times of all individuals (Inf if not infected), ordered by ID.}

\item{events}{A 2 column matrix where the first column is the infection times, and the second is the paired removal times.}

\item{B}{The infection rate matrix.}
}
\value{
Returns the value of the product detailed above.
}
\description{
This function is used in the calculation of the likelihood for a General Stochastic Epidemic, and calculates 
 \eqn{\prod_{j != kappa}^{n_I} [ \sum_{ i in I_{n_{j-}} } [ beta_{i,j} ] ]}.
}
\examples{
This function is utilised by the chris_log_likelihood function.
}
\keyword{GSE}
\keyword{Product}
\keyword{infection}
\keyword{likelihood}
