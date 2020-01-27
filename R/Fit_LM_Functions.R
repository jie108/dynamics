###11-16-08: Implementation of the iterative fitting procedure based LM
### this is the final clean up version

####(a)
###penalized least square for updating x0 (initial values) when theta and beta are fixed

PLS.ini <- function(x0.old,theta.old,beta.old,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1,tau=1){

###para: 
###x0.old: vector for current initial conditions, vector of length n (number of curves) 
###theta.old: vector for current random scale parameters, vector of length num.id (number of cluters)
###beta.old: current basis coefficients for g, vector of lenght M; 
###data.path.u : measurements
###obstime.u: measurement times; 
###N.u: number of grid points; h.u: grid step size; knots.u: knots for basis functions
###lambda1: penalty parameters for x0
### tau = dampening factor <= 1
###return: updated estimates of x0 and its gradient

#####
num.id<-length(table(subid))    ##number of clusters
id.list<-unique(subid)   ### cluster indeces
n.u <- length(obstime.u)                 ##sample size: total number of curves 
m.u <- numeric(n.u)                      ## record number of measurements for each curve 
h2 <- 2*h.u

##### evaluation at the current estimates 
y.tilde <- xpath(N.u,h.u,x0.old,theta.old,beta.old,knots.u)        ### sample path evaluated at current estimates
z.tilde <- xderiv.a(N.u,h.u,theta.old,beta.old,knots.u,y.tilde)  ### gradient w.r.t. initial conditions a_il          ######### NEED UPDATE!!!############

                                         
##### COMPUTATION of J_i (Jacobian) and updating for a_il 
JA <- NULL
grad.A <- rep(0,n.u)             
eps.tilde <- NULL
alpha.old <- mean(x0.old)           ##the current estimate for mean of x0 (alpha)
x0.new<-numeric(n.u)

  for(i in 1:n.u){
     m.u[i] <- length(obstime.u[[i]])   ##number of measurements for ith curve 
     grid<-ceiling(obstime.u[[i]]/h2)   ##observation grid
     JA[[i]] <- z.tilde[i,grid]         ## derivative w.r.t. a_i evaluated at the observed times 
     eps.tilde[[i]] <- data.path.u[[i]] - y.tilde[i,ceiling(obstime.u[[i]]/h.u)]   ## current residuals 
     temp1<-sum(JA[[i]]*eps.tilde[[i]])
     J.AA <- sum(JA[[i]]^2)
     udiff.a <- temp1 + lambda1 * (alpha.old - x0.old[i])   ##RHS updates for x0
     x0.new[i] <- x0.old[i] + tau*(J.AA + lambda1)^{-1} * udiff.a          ##new estimates
     grad.A[i] <- -2*temp1 + 2*lambda1*(x0.old[i]-alpha.old) ## gradient of loss function w.r.t. a_{il}
 
  }

return(list(x0.new,grad.A))
}




####(b)
###penalized least square for updating theta (random scale) when x0 and beta are fixed

PLS.theta <- function(x0.new,theta.old,beta.old,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda2,tau=1){
###para: 
###x0.new: vector for current initial conditions, vector of length n (number of curves) 
###theta.old: vector for current random scale parameters, vector of length num.id (number of cluters)
###beta.old: current basis coefficients for g, vector of lenght M; 
###data.path.u : measurements
###obstime.u: measurement times; 
###N.u: number of grid points; h.u: grid step size; knots.u: knots for basis functions
###lambda2: penalty parameters for theta
### tau = dampening factor <= 1
###return: updated estimates of theta and its gradient


#####
num.id<-length(table(subid))    ##number of clusters
id.list<-unique(subid)   ### cluster indeces
n.u <- length(obstime.u)                 ##sample size: total number of curves 
m.u <- numeric(n.u)                      ## record number of measurements for each curve 
h2 <- 2*h.u

#####
y.tilde <- xpath(N.u,h.u,x0.new,theta.old,beta.old,knots.u)        ##  sample path evaluated at current estimates
v.tilde <- xderiv.theta(N.u,h.u,theta.old,beta.old,knots.u,y.tilde)                                ######### NEED UPDATE!!!############

##### COMPUTATION of J_i (Jacobian)  and updating for theta's
JT <- NULL
grad.T <- rep(0,num.id)
eps.tilde <- NULL
##vareps.est <- 0
theta.temp <-numeric(num.id)

  for(i in 1:n.u){
     grid<-ceiling(obstime.u[[i]]/h2)
     JT[[i]] <- v.tilde[i,grid]                              ## derivative w.r.t. theta_i evaluated at the observed times 
     eps.tilde[[i]] <- data.path.u[[i]] - y.tilde[i,ceiling(obstime.u[[i]]/h.u)]   ## current residuals
    ## vareps.est <- vareps.est + sum(eps.tilde[[i]]^2) 
  }
  
  ##vareps.est <- vareps.est/(sum(m.u-1)-length(id.list)-length(knots.u))
  
  for(k in 1:num.id){    #### k = clusterindex,  l = curve index

    locindex <- which(subid ==id.list[k])    ### index of curves which correspond to cluster k
    J.TT <- 0
    udiff.theta <- 0

    for(l in 1:length(locindex)){
      J.TT <- J.TT + sum(JT[[locindex[l]]]^2)
      udiff.theta <- udiff.theta + sum(JT[[locindex[l]]]*eps.tilde[[locindex[l]]]) 
    }
            
    theta.temp[k] <- theta.old[k] + tau * (J.TT + lambda2)^{-1} * (udiff.theta - lambda2 * theta.old[k]) ## dampening factor on the updates
     grad.T[k] <- (-2*udiff.theta + 2*lambda2 * theta.old[k])/length(locindex)  ## (averaged) gradient of loss function w.r.t. theta_{i} ####???#####    
  }
      
  theta.temp <- theta.temp - mean(theta.temp)
  theta.new <- numeric(n.u)
  for(k in 1:num.id){ 
     locindex <- which(subid == id.list[k])
     theta.new[locindex] <- theta.temp[k]
  }

return(list(theta.new,grad.T))
}


#### (c)
####Levenberg-Marquardt for updating beta, when x0 and theta are fixed: with the penalty matrix: Phi.pen
LM.beta <- function(x0.new,theta.new,beta.old,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda3,Phi.pen=NULL,tau=1){

###para:
###x0.new: vector for current initial conditions, vector of length n (number of curves) 
###theta.new: vector for current random scale parameters, vector of length num.id (number of cluters)
###beta.old: current basis coefficients for g, vector of lenght M; 
###data.path.u : measurements
###obstime.u: measurement times; 
###N.u: number of grid points; h.u: grid step size; knots.u: knots for basis functions
###lambda3: penalty parameter for LM when updating beta
###Phi.pen: the matrix in the quadratic penalty term for beta
###tau = dampening factor <= 1 for updating beta
###return:updated beta and its gradient and estimated residual variance (NOT standaridized) at current parameters


#####
id.list<-unique(subid)   ### cluster indecies
M.u <-length(beta.old)                  ###number of knots/coefficients
n.u <-length(obstime.u)                 ##sample size: total number of curves 
h2<-2*h.u


##### evaluation at the current estimates 
y.tilde <- xpath(N.u,h.u,x0.new,theta.new,beta.old,knots.u)         ### sample path 
w.tilde <- xderiv.beta(N.u,h.u,theta.new,beta.old,knots.u,y.tilde)  ### gradient w.r.t. beta     ######### NEED UPDATE!!!############


##### COMPUTATION of J_i (Jacobian) and updating estimates of beta
if(is.null(Phi.pen)){
Phi.pen=matrix(0,M.u,M.u)
}

Jbeta <- NULL
Jbeta.mat <- matrix(0,M.u,M.u)
Jbeta.eps <- rep(0,M.u)
eps.tilde <- NULL
vareps.est<-0 
m.u<-numeric(n.u)              ##record for number of measurements of each curve


for(i in 1:n.u) {
  grid<-ceiling(obstime.u[[i]]/h2)
  m.u[i]<-length(obstime.u[[i]])
  eps.tilde[[i]] <- data.path.u[[i]] - y.tilde[i,ceiling(obstime.u[[i]]/h.u)]   ## current residuals
  
  for(l in 1:M.u){
      Jbeta[[M.u-l+1]][[i]] <- w.tilde[i,grid,M.u-l+1]      ## derivative w.r.t. beta_{M-l+1} 
  }
  
  for(l in 1:M.u){
      for(k in 1:l){
      Jbeta.mat[l,k] <- Jbeta.mat[l,k] + sum(Jbeta[[l]][[i]]*Jbeta[[k]][[i]])
      Jbeta.mat[k,l]<-Jbeta.mat[l,k]  
      }      
     Jbeta.eps[l] <- Jbeta.eps[l] + sum(Jbeta[[l]][[i]]*eps.tilde[[i]])
  }

  vareps.est <- vareps.est + sum(eps.tilde[[i]]^2)
}  

##### update beta  with damping factor tau for beta
  beta.new <- beta.old +  tau * solve(Jbeta.mat + lambda3 * diag(diag(Jbeta.mat)) + Phi.pen) %*% (Jbeta.eps - Phi.pen %*% beta.old)   

#####  
  ##vareps.est<-vareps.est/(sum(m.u-1)-length(id.list)-M.u)
  grad.beta <- -2*(Jbeta.eps - Phi.pen %*% beta.old)/n.u   ## (averaged) gradient of loss function w.r.t. beta with the matrix penalty term  ###????####
  df<-sum(m.u-1)-length(id.list)-M.u
  df.fix<-sum(m.u)-length(id.list)-M.u
 
return(list(beta.new,grad.beta,vareps.est,df,df.fix))

}


###### (d)
####Fit the model by iteration: WITH THE damping factors; and the matrix penalty term for beta

FitLM<-function(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,lambda2.ini,lambda3.v,Phi.pen=NULL,adapt1=FALSE, adapt2=FALSE, ini.fix=FALSE, niter=10,thre=1e-3, conv.cri="beta",tau.v=matrix(1,3,niter)){
#####additional paramters:
##adapt1, adapt2: True (adatively updated lambda1 and/or lambda 2); FALSE (non-adaptive) 
##ini.fix: TRUE (x0 viewed as known); FALSE (need to update x0)
##niter: maximum number of iterations 
##thre: threshold for convergence; 
##conv.cri: "beta" (maximum change of beta and number of iterations);
##or "grad" (maximum gradient and number of iterations)
##tau = dampening factors at each step: three by niter matrix 


#####
x0.u <- NULL
theta.u <-NULL
beta.u <-NULL
vareps <-0


#####update x0
count<-1
print(paste("LM step:",count))
if(ini.fix==FALSE){
temp<-PLS.ini(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,tau=tau.v[1,1])
x0.new<-temp[[1]]
}else{
x0.new<-x0.ini
}


#####update theta
temp <- PLS.theta(x0.new,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda2.ini,tau=tau.v[2,1])
theta.new<-temp[[1]]

#####update beta
temp <- LM.beta(x0.new,theta.new,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda3.v[1],Phi.pen,tau=tau.v[3,1])
beta.new<-temp[[1]]
if(ini.fix==FALSE){
vareps.new<-temp[[3]]/temp[[4]]
}else{
vareps.new<-temp[[3]]/temp[[5]]
}


################updating  steps

beta.old <- beta.new+100
gradsup.norm <- 100

  if(conv.cri=="beta"){
  conv.stat<-(sqrt(sum((beta.old - beta.new)^2)) > thre && (count<=niter))
  }
  
  if(conv.cri=="grad"){
  conv.stat<-(max(gradsup.norm) > thre && (count<=niter))
  }


lambda1 <- lambda1.ini
lambda2 <- lambda2.ini 


while(conv.stat){
 
 x0.u[[count]]<-x0.new
 theta.u[[count]]<-theta.new
 beta.u[[count]]<-beta.new
 beta.old <- beta.new

 vareps[count]<-vareps.new
 vareps.old<-vareps.new

if(adapt1){
 lambda1 <- vareps.new/var(x0.new)
}

if(adapt2){
 lambda2<- vareps.new/var(unique(theta.new))
}



#####update x0
if(ini.fix==FALSE){
temp<-PLS.ini(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1,tau=tau.v[1,count])
x0.new<-temp[[1]]
}else{
x0.new<-x0.ini
}


#####update theta
temp <- PLS.theta(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda2,tau=tau.v[2,count])
theta.new<-temp[[1]]
grad.T <- temp[[2]]


#####update beta
temp <- LM.beta(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda3.v[count],Phi.pen,tau=tau.v[3,count])
beta.new<-temp[[1]]
grad.beta <- temp[[2]] 

if(ini.fix==FALSE){
vareps.new<-temp[[3]]/temp[[4]]
}else{
vareps.new<-temp[[3]]/temp[[5]]
}

##print(paste("beta",beta.new))

######
gradsup.norm <- c(max(abs(grad.T)),max(abs(grad.beta)))
##print(paste("sup norm of gradient of loss function ="))
##print(gradsup.norm)

####
 count<-count+1
 print(paste("LM step:",count))
 if(conv.cri=="beta"){
   conv.stat<-(sqrt(sum((beta.old - beta.new)^2)) > thre && (count<=niter))
   }
   
   if(conv.cri=="grad"){
   conv.stat<-(max(gradsup.norm) > thre && (count<=niter))
   }
 


 }  ## end of while

 
 ###
 x0.u[[count]]<-x0.new
 theta.u[[count]]<-theta.new
 beta.u[[count]]<-beta.new
 vareps[count]<-vareps.new
 
if(conv.cri=="beta"){
   conv.stat<-(sqrt(sum((beta.old - beta.new)^2)) > thre )
   }
   
   if(conv.cri=="grad"){
   conv.stat<-(max(gradsup.norm) > thre)
   }
   
return(list(x0.u,theta.u,beta.u,vareps,count,lambda1,lambda2,1-conv.stat))

}


 
