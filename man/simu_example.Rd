\name{simu_example}
\alias{simu_example}
\title{A simulated data example}
\description{This is version 0.2-1 updated in April, 2011. A list containing a simulated data example for package \code{dynamics}}

\value{
example is a list of ten components
\item{knots.tr} {a vector of bspline knots}
\item{beta.tr} {a  vector of corresponding basis coefficients}
\item{data.path}{a list of measurements}
\item{obstime}{a list of measurement times}
\item{numsub}{number of clusters}
\item{subject.id}{a vector of subject IDs}
\item{x0.tr}{a vector of initial conditions}
\item{sig.theta}{sd of the distribution of scale parameters}
\item{sig.x0}{sd of the distribution of initial conditions}
\item{sigeps.tr}{sd of the error distribution}
}



\references{
D. Paul, J. Peng, P. Burman, W. Sacks(2008). Semiparametric modelling of autonomous nonlinear dynamical systems with applications.
}

\usage{data(simu_example)}

\keyword{datasets}

