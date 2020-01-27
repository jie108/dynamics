\name{dynamics.fit}
\alias{dynamics.fit}
\title{A function to estimate the ``base gradient" function in an autonomous dynamical system via basis representation}
\description{This is version 0.2-1 updated in April, 2011. Note 1: In the dynamical system subject-specific scale parameters are allowed; 
Note 2: when method="natural", length(beta.ini)=length(knots.u)+3, since the first three terms of beta.ini corresponds to $(x,x^2,x^3)$ (note that in this case, the constant term is omitted and thus g(0)=0 is imposed); 
when method="bspline", length(beta.ini)=length(knots.u)}   
\usage{
dynamics.fit(x0.ini=NULL,theta.ini=NULL,beta.ini=NULL,subid=NULL,data.path.u,obstime.u,N.u=1000,method="bspline",knots.u,lambda1.ini=0,lambda2.ini=1,lambda3.v=NULL,Phi.pen=NULL,
adapt1.LM=FALSE, adapt2.LM=FALSE,adapt1.NR=TRUE, adapt2.NR=TRUE, ini.fix=FALSE, niter.LM=10,thre.LM=1e-3,conv.cri.LM="beta",niter.NR=10,thre.NR=1e-3,conv.cri.NR="grad",
tau.v.LM=NULL,tau.v.NR=NULL,NRadapt=FALSE,kstep=1)
}

\arguments{
    \item{x0.ini}{numeric vector:initial values of the initial conditions, vector of length N (total number of curves); default(=NULL): first measurements of each curve}
    \item{theta.ini}{numeric vector: initial values of the subject-specific scale parameters, vector of length n (number of cluters); default(=NULL): ZERO vector}
    \item{beta.ini}{numeric vector: initial values of the basis coefficients, vector of length M (number of basis functions); default(=NULL): ONE vector}
    \item{subid}{vector (numeric of string): cluster id for each curve, vector of length n (number of clusters); default(=NULL): no cluster}
    \item{data.path.u}{list (numeric): list of measurements for each curve, list length=N; each component is a numerical vector corresponding to measurements of one curve}
    \item{obstime.u}{list (numeric): list of measurement times for each curve, list length=N;each component is a numerical vector corresponding to measurement times of one curve}
    \item{N.u}{integer: number of steps for the Runge-Kunta method; default=1000}
    \item{method}{string: method for basis representation; two possible values: "bspline" (default), "natural" (truncated cubic polynomials)}
    \item{knots.u}{numeric vector: sequence of knots; length=M (number of basis funcitons)}
    \item{lambda1.ini}{numeric scalar: initial value for the penalty parameter of the initial conditions; default=0}
    \item{lambda2.ini}{numeric scalar: initial value for the penalty parameter of the scale parameters; default=1}
    \item{lambda3.v}{numeric vector: penalty parameters of the basis coefficients. it is a descreasing sequence of numbers; defaul(=NULL):0.1/(1:niter.LM)}
    \item{Phi.pen}{numeric matrix: penalty matrix of the basis coefficients. it is a M by M numeric matrix;defaul(=NULL): zero matrix}
    \item{adapt1.LM}{logic (TRUE or FALSE): indicating whether lambda1 is adaptively updated in the LM step; default=FALSE (no updating)}
    \item{adapt2.LM}{logic (TRUE or FALSE): indicating whether lambda2 is adaptively updated in the LM step; default=FALSE (no updating)}
    \item{adapt1.NR}{logic (TRUE or FALSE): indicating whether lambda1 is adaptively updated in the NR step; default=TRUE (updating)}
    \item{adapt2.NR}{logic (TRUE or FALSE): indicating whether lambda2 is adaptively updated in the NR step; default=TRUE (updating)}
    \item{ini.fix}{logic (TRUE or FALSE): indicating whether initial conditions are viewed as known or not; if known, then x0.ini will not be updated;default=FALSE (not known)}
    \item{niter.LM}{integer: maximum number of LM steps; default=10}
    \item{thre.LM}{numeric scalar: threshold for checking convergence of the LM step; default=0.001}
    \item{conv.cri.LM}{string: convergence criterion for the LM step; two possible values: "beta" (default)--- convergence based on basis coefficients; or "grad"--convergence based on maximum gradients; }
     \item{niter.NR}{integer: maximum number of NR steps; default=10}
    \item{thre.NR}{numeric scalar: threshold for checking convergence of the NR step; default=0.001}
    \item{conv.cri.NR}{string: convergence criterion for the NR step; two possible values: "beta"--- convergence based on basis coefficients; or "grad" (default)--convergence based on maximum gradients; } 
    \item{tau.v.LM}{numeric matrix: dampening factors in LM step for  the three updating circles for initial conditions, scale parameters and basis coefficients; a 3 by niter.LM numeric matrix; default(=NULL): ONE matrix (no dampening)}
    \item{tau.v.NR}{numeric matrix: dampening factors in NR step for  the three updating circles for initial conditions, scale parameters and basis coefficients; a 3 by niter.NR numeric matrix; default(=NULL): ONE matrix (no dampening)}
    \item{NRadapt}{logic (TRUE or FALSE): in the NR step, whether only updating the Hessian matrix of the basis coefficients  in every kstep steps; this is used to faciliate computation;default=FALSE;} 
    \item{kstep}{integer: if NRadapt=TRUE, only updating the Hessian matrix of the basis coefficients  in every kstep steps; default=1}
    }        
\details{
    \code{dynamics.fit} uses a numeric ODE solver (4th order Runge-Kunta method) and nonlinear optimization techniques (Levenberg-Marqardt and Newton-Raphson)  to fit an autonomous dynamical system via basis representation }

\value{
    A list with ten components
    \item{x0.NR}{numeric vector of length N: final estimates of the initial conditions for each curve}
    \item{theta.NR}{numeric vector of length N: final estimates of the scale parameters for each curve}
    \item{beta.NR}{numeric vector of length M: final estimates of the basis coefficients. Note: when method="natural", length(beta.NR)=length(knots.u)+3, since the first three terms of beta.NR corresponds to $(x,x^2,x^3)$}
    \item{vareps.NR}{numeric scalar: final estimate of the error variance}
    \item{gradient.NR}{numeric scalar: the maximum gradient evaluated at the final estimates}
    \item{count.LM}{integer: actual number of LM steps performed}
    \item{count.NR}{integer: actual number of NR steps performed}
    \item{lambda1.NR}{numeric scalar: final value of the penaty parameter for the initial conditions}
    \item{lambda2.NR}{numeric scalar: final value of the penaty parameter for the scale parameters}
    \item{error.c}{logic (TRUE or FALSE): whether an error occurs during the fitting process}
       }

\references{
D. Paul, J. Peng, P. Burman, W. Sacks(2008). Semiparametric modelling of autonomous nonlinear dynamical systems with applications.
}

\author{J. Peng, D. Paul}

\keyword{methods}
\examples{
### (I) Simulation: in this simulated example, the true "base gradient" g is represented by Bspline basis with knots at (0.35 0.60 0.85 1.10) 
##       and basis coefficients being (0.1 1.2 1.6 0.4).
### Data: There are ten clusters (n=10), each having 20 curves (thus total number of curves N=200); 
####The scale parameter for each cluster is randomly sampled from normal with mean zero, sd=0.1; 
###the initial condition for each curve is randomly sampled from chi-square  with mean=0.25, sd=0.05; 
###errors are randomly sampled from normal with mean zero, sd=0.01; 
### For each curve, the number of measurements are uniformly sampled from {5,6,7,8,9,...,20} 
### and measurement times are uniformly sampled from [0,1]; 
### the observed data is then sampled according to the data molde described in Section 2 of the paper. 

data(simu_example)  
####including: data.path -- the list for measurements; obstime -- the list for measurement times; beta.tr: true basis coefficients; knots.tr: true knots for Bsplines 
####subject.id: cluster ID for each curve
summary(data.path)
summary(obstime)
table(subject.id)

###(II) Estimation settings: initial value etc. 
knots.u <- knots.tr   ##beta.ini: use the default (ONE vector); theta.ini: use the default (ZERO vector); x0.ini: use the default (first measurements)

itcountLM <- 200       ### maximum number of iterations for L-M step
itcountNR <- 200        ### maximum number of iterations for N-R step

lambda1.ini<-sigeps.tr^2/sig.x0^2   
lambda2.ini <- sigeps.tr^2/sig.theta^2       ## ratio of variances
lambda3.v <- 1/(1:itcountLM)    #### works better in most simulations

### (III) Estimation: note by default: the LM step is non-adaptive and the NR step is adaptive; also the initial conditions are fitted (ini.fix=FALSE)
result<-dynamics.fit(subid=subject.id,data.path.u=data.path,obstime.u=obstime,method="bspline",
knots.u=knots.u,lambda1.ini=lambda1.ini,lambda2.ini=lambda2.ini,lambda3.v=lambda3.v,niter.LM=itcountLM,niter.NR=itcountNR)


### (IV) Results
summary(result)
x0.est<-result[[1]]
theta.est<-result[[2]]
beta.est<-result[[3]]
sig.theta.est <- sd(result[[2]]) 
sigeps.est <- sqrt(result[[4]])
lambda1.final<-result[[8]]
lambda2.final<-result[[9]]


beta.tr    ##true basis coefficients
beta.est   ##estimated basis coefficients
sig.theta  ##true s.d. of scale parameter
sig.theta.est ## estimated s.d. of scale parameter
sigeps.tr  ## true noise standard deivation
sigeps.est ## estimated noise standard deviation

### CV score
cv<-dynamics.cv(x0.final=x0.est,theta.final=theta.est,beta.final=beta.est,subid=subject.id,data.path.u=data.path,obstime.u=obstime,
method="bspline",knots.u=knots.u,lambda1.final=lambda1.final,lambda2.final=lambda2.final,ini.fix=FALSE)

summary(cv)
print(cv[[1]])    ###CV score based on prediction error
print(cv[[2]])    ###CV score based on loss function
}
