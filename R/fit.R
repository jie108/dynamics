####the final wrapped up function for the package
####12-02-08
dynamics.fit <- function(x0.ini=NULL,theta.ini=NULL,beta.ini=NULL,subid=NULL,data.path.u,obstime.u,N.u=1000,method="bspline",knots.u,lambda1.ini=0,lambda2.ini=1,lambda3.v=NULL,Phi.pen=NULL,
adapt1.LM=FALSE, adapt2.LM=FALSE,adapt1.NR=TRUE, adapt2.NR=TRUE, ini.fix=FALSE, niter.LM=10,thre.LM=1e-3,conv.cri.LM="beta",niter.NR=10,thre.NR=1e-3,conv.cri.NR="grad",
tau.v.LM=NULL,tau.v.NR=NULL,NRadapt=FALSE,kstep=1){
###method: "bspline" or "natural"

####
h.u<-1/N.u
Hbeta.old<-NULL

####default
n.u <- length(obstime.u)                 ##sample size: total number of curves 

if(is.null(x0.ini)){
   for(count in 1:n.u){
    x0.ini[count] <- data.path.u[[count]][1]   ###initial values for x0: the first observations
   }
}

if(is.null(theta.ini)){
theta.ini<-numeric(n.u)                 ##default initial values for theta are all zeros
}

if(is.null(beta.ini)){
M.u<-length(knots.u)
 if (method=="bspline"){
  beta.ini<-rep(1,M.u)                   ##default initial values for theta are all ones
 }
 
 if (method=="natural"){
   M.u<-M.u+3
   beta.ini<-rep(1,M.u)                   ##default initial values for theta are all ones
 }
 
}

if(is.null(subid)){
subid<-1:n.u                           ##default for subid: one curve per cluster, that is no clusters
}




if(is.null(lambda3.v)){
lambda3.v<-0.1/(1:niter.LM)
}

if(is.null(tau.v.LM)){
tau.v.LM<-matrix(1,3,niter.LM)
}

if(is.null(tau.v.NR)){
tau.v.NR<-matrix(1,3,niter.NR)
}



####format: format the measurement times into [0,1]
b <- range(obstime.u)[2]
a <- range(obstime.u)[1]
 for (count in 1:n.u){
  time.c<-obstime.u[[count]]
   for (j in 1:length(time.c)){
     time.c[j]<-(time.c[j]-a)/(b-a) 
     if(time.c[j]==0){
      time.c[j]<-0.000001
     }
   }   
   obstime.u[[count]]<-time.c
 }


####fitting 
if(method=="bspline"){

result<-Fit(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,lambda2.ini,lambda3.v,Phi.pen,
adapt1.LM, adapt2.LM,adapt1.NR, adapt2.NR, ini.fix, niter.LM,thre.LM,conv.cri.LM,niter.NR,thre.NR,conv.cri.NR,
tau.v.LM,tau.v.NR,NRadapt,kstep,Hbeta.old)
}

if(method=="natural"){
result<-Fit.N(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,lambda2.ini,lambda3.v,Phi.pen,
adapt1.LM, adapt2.LM,adapt1.NR, adapt2.NR, ini.fix, niter.LM,thre.LM,conv.cri.LM,niter.NR,thre.NR,conv.cri.NR,
tau.v.LM,tau.v.NR,NRadapt,kstep,Hbeta.old)
}

return(result)

}


