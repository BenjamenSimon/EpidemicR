% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_reparam.R
\name{MCMC_reparam}
\alias{MCMC_reparam}
\title{Make inference on a General Stochastic Epidemic.}
\usage{
MCMC_reparam(N.its, N, inf.ids, rem.times, dist.mat, lambda.b1 = 0.001,
  nu.b1 = 1, lambda.g = 0.001, nu.g = 1, inc.beta1 = list(NA, T),
  inc.p = list(NA, T), inc.dist = list(NA, T),
  inc.inf.times = list(NA, T), inc.gamma = list(NA, T), d.upper,
  sigmap, sigmad, infupdate = 1)
}
\arguments{
\item{N.its}{The number of desired iterations of the MCMC algorithm.}

\item{N}{The total size of the population.}

\item{inf.ids}{A vectors of the IDs of the infected individuals.}

\item{rem.times}{A vector of the removal times, ordered by individual ID.}

\item{dist.mat}{An NxN distance matrix.}

\item{lambda.b1}{The rate parameter for beta1, assuming a Gamma prior.}

\item{nu.b1}{The shape parameter for beta1, assuming a Gamma prior.}

\item{lambda.g}{The rate parameter for gamma, assuming a Gamma prior.}

\item{nu.g}{The shape parameter for gamma, assuming a Gamma prior.}

\item{inc.beta1}{A list object of 2 levels, the true value of beta1, and T/F binary value that says 
whether to make inference on beta1 (T = make inference, see details for the different options).}

\item{inc.p}{See inc.beta1, but for p.}

\item{inc.dist}{See inc.beta1, but for the distance d.}

\item{inc.inf.times}{A list object of 2 levels; a vector of the true infection times, and T/F binary value that says
whether to make inference on the infection times (T = make inference, see details for the different options).}

\item{inc.gamma}{See inc.beta1, but for the removal rate gamma.}

\item{d.upper}{The upper bound on the distance d.}

\item{sigmap}{Tuning parameter. The standard deviation of the folded random walk for p.}

\item{sigmad}{Tuning parameter. The standard deviation of the folded random walk for d.}

\item{infupdate}{Tuning parameter. The number of infection times that should be updated at each iteration.}
}
\value{
This function returns a list object with elements; a matrix of results (which included all the accepted samples and
 the log-likelihood at the end of each iteration), the acceptance rate of the infection times, the acceptance rate of p, and
 the acceptance rate of d.
}
\description{
This function uses realistically available information from a GSE to perform Bayesian inference and attempt to recover the parameters,
 in the setting with a reparameterisation of beta2 = p * beta1, with p in (0,1).
}
\details{
The function has the functionality to allow the fixing of some (or all) parameters. For the function parameters that begin inc.
 a list can be fed to them of the form (true parameter value(s), TRUE/FALSE), where the TRUE dictates that inference should be
 made for that parameter. Possible options for this list are; list(parameters, F) in which case the function will fix that parameter
 value, list(NA, T) in which case the function will randomly generate a valid initial value for the parameter and then make inference,
 list(parameter, T) in which case the function will initialise the parameter at its true value and then make inference, or (NA, F) 
 which we do not suggest using as it will fix the parameter at a random value.
 
 Gamma is always assumed to be initialised at the value 0.15.
}
\examples{
infernce_r <- MCMC_reparam(N.its = 10000, N = 10, inf.ids, rem.times, dist.mat, 
lambda.b1 = 0.001, nu.b1 = 1, lambda.g = 0.001 , nu.g = 1,
inc.beta1 = list(0.004, T), inc.p = list(NA, T), inc.dist = list(NA, T),
inc.inf.times = list(inf.times, F), inc.gamma = list(NA, T),
d.upper = 1.5, sigmap = 0.1, sigmad = 2, infupdate = 1)
}
\keyword{Gibbs}
\keyword{Hastings}
\keyword{MCMC}
\keyword{MH}
\keyword{Metropolis}
\keyword{inference}
\keyword{reparameterisation}
