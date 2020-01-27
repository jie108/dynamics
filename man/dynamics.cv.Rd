\name{dynamics.cv}
\alias{dynamics.cv}
\title{A function to calculate the approximate CV score}

\description{This is version 0.2-1 updated in April, 2011. 
A function to calculate the approximate CV score for a given model (specified by method and knots.u) based on the final estimates of dynamics.fit (x0.final, theta.final, beta.final).
Note: two cv scores are returned: one based on the prediction error criterion (cvcur.alt); and another based on the loss function with penalty terms (cvcur.score); cvcur.alt is recommended.
}   
\usage{
dynamics.cv(x0.final,theta.final,beta.final,subid=NULL,data.path.u,obstime.u,N.u=1000,method="bspline",knots.u,lambda1.final,lambda2.final,Phi.pen=NULL,ini.fix=FALSE, maxloop=10,thre=0.0001)
}

\arguments{
    \item{x0.final}{numeric vector:final estimates of the initial conditions by dynamics.fit; vector of length N (total number of curves).}
    \item{theta.final}{numeric vector: final estimates of the subject-specific scale parameters by dynamics.fit; vector of length n (number of cluters);}
    \item{beta.final}{numeric vector: final estimates of the basis coefficients by dynamics.fit; vector of length M (number of basis functions)}
    \item{subid}{vector (numeric of string): cluster id for each curve, vector of length n (number of clusters); default(=NULL): no cluster}
    \item{data.path.u}{list (numeric): list of measurements for each curve, list length=N; each component is a numerical vector corresponding to measurements of one curve}
    \item{obstime.u}{list (numeric): list of measurement times for each curve, list length=N;each component is a numerical vector corresponding to measurement times of one curve}
    \item{N.u}{integer: number of steps for the Runge-Kunta method; default=1000}
    \item{method}{string: method for basis representation; two possible values: "bspline" (default), "natural" (truncated cubic polynomials)}
    \item{knots.u}{numeric vector: sequence of knots; length=M (number of basis funcitons)}
    \item{lambda1.final}{numeric scalar: final value for the penalty parameter of the initial conditions by dynamics.fit}
    \item{lambda2.final}{numeric scalar: final value for the penalty parameter of the scale parameters by dynamics.fit}
    \item{Phi.pen}{numeric matrix: penalty matrix of the basis coefficients. it is a M by M numeric matrix;defaul(=NULL): zero matrix}
    \item{ini.fix}{logic (TRUE or FALSE): indicating whether initial conditions are viewed as known or not; if known, then x0.ini will not be updated;default=FALSE (not known)}
    \item{maxloop}{integer: maximum number of iterations for calculating the drop-one-out initial condition $a^{-(i)}$; default=10}
    \item{thre}{numeric scalar: threshold to check the convergence of the drop-one-out initial condition  estimate; default=0.0001}
    }        
\details{
    \code{dynamics.cv} uses a second order Taylor expansion to approximate the leave-one-curve-out cv score}

\value{
    A list with five components
    \item{cvcur.alt}{numeric scalar: the approximate leave-one-curve-out cv score based on prediction error criterion}
    \item{cvcur.score}{numeric scalar: the approximate leave-one-curve-out cv score based on the loss function with penalty terms}
    \item{x0.drop}{numeric vector of length N: the drop-one-curve-out estimates of the initial conditions for each curve}
    \item{theta.drop}{numeric vector of length N: the drop-one-curve-out estimates of the scale parameters for each curve}
    \item{beta.drop}{numeric matrix of dimension M by N: the drop-one-curve-out estimates of the basis coefficients. Note: when method="natural", nrow(beta.NR)=length(knots.u)+3, since the first three terms of beta.NR corresponds to $(x,x^2,x^3)$}
    }

\references{
D. Paul, J. Peng, P. Burman, W. Sacks(2008). Semiparametric modelling of autonomous nonlinear dynamical systems with applications.
}

\author{J. Peng, D. Paul}

\keyword{methods}
