###11-16-08: Implementation of the iterative fitting procedure (only for theta and beta) based NR
### this is the final clean up version

##### (a)
##### Newton step for theta : SAME AS BEFORE
NR.theta.N <- function(x0.old,theta.old,beta.old,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda2.u,tau=1){
###para: 
###x0.old: vector for current initial conditions, vector of length n (number of curves) 
###theta.old: vector for current random scale parameters, vector of length num.id (number of clusters)
###beta.old: current basis coefficients for g; 
###subid : cluster index
###data.path.u : measurements
###obstime.u: measurement times;
###N.u: number of grid points; h.u: grid step size;  knots.u: vector of knots
###lambda2.u: penalty parameter for theta
###tau: damping factor (<=1)
###return: updated estimates of theta and its gradient

#####
  num.id<-length(table(subid))    ##number of clusters
  id.list<-unique(subid)          ## cluster indices 
  n.u <- length(obstime.u)        ##sample size: total number of curves 
  M.u <- length(beta.old)
  h2 <- 2*h.u
  h4 <- 4*h.u

######
  y.tilde <- xpath.N(N.u,h.u,x0.old,theta.old,beta.old,knots.u)        ### sample path evaluated at current estimates
  v.tilde <- xderiv.theta.N(N.u,h.u,theta.old,beta.old,knots.u,y.tilde)  ### gradient w.r.t. theta                   #####NEED UPDATE!!!!#####
  vv.tilde <- xHess.theta.N(N.u,h.u,theta.old,beta.old,knots.u,y.tilde,v.tilde) ### Hessian w.r.t. theta              #####NEED UPDATE!!!!#####  

  
#####
m.u <- numeric(n.u)
theta.temp <- numeric(num.id)
J.T <- NULL
H.T <- NULL
Jeps.T <- numeric(num.id)
H.TT <- numeric(num.id)
grad.T <- numeric(num.id)             
eps.tilde <- NULL

  for(i in 1:n.u){
   grid<-ceiling(obstime.u[[i]]/h2)
   grid2<-ceiling(obstime.u[[i]]/h4)
   m.u[i] <- length(obstime.u[[i]])     ##number of measurements for ith curve

   eps.tilde[[i]] <- data.path.u[[i]] - y.tilde[i,ceiling(obstime.u[[i]]/h.u)]   ### current residual 
   J.T[[i]] <- v.tilde[i,grid]     ### derivative w.r.t. theta_i evaluated at the observed times 
   H.T[[i]] <- vv.tilde[i,grid2]   ### Hessian w.r.t. theta_i evaluated at the observed times
  }


  for(k in 1:num.id){    #### k = cluster index,  l = curve index
    locindex <- which(subid == id.list[k])    ### index of locations which correspond to kth cluster
    H.TT[k] <- lambda2.u   
    Jeps.T[k] <- -lambda2.u * theta.old[locindex][1] 
    for(l in 1:length(locindex)){
      i <- locindex[l]
      H.TT[k] <- H.TT[k] + sum(J.T[[i]]^2) - sum(eps.tilde[[i]]*H.T[[i]]) 
      Jeps.T[k] <- Jeps.T[k] + sum(eps.tilde[[i]]*J.T[[i]])
    }

    theta.temp[k] <- theta.old[locindex][1] + tau*Jeps.T[k]/H.TT[k]  ## Newton-Raphson update for theta    ###damping 

    grad.T[k] <- -2*Jeps.T[k]/length(locindex)    ### (averaged) gradient of loss function w.r.t. theta         #####???#####
  }


##### generating full n.u x 1 vector theta.new
##theta.temp <- theta.temp - mean(theta.temp)
  theta.new <- numeric(n.u)
  for(k in 1:num.id){ 
     locindex <- which(subid == id.list[k])
     theta.new[locindex] <- theta.temp[k]
  }
 
return(list(theta.new,grad.T))

}



#####(b)
####Newton step for beta: with the matrix penalty term
NR.beta.N<- function(x0.old,theta.old,beta.old,subid,data.path.u,obstime.u,N.u,h.u,knots.u,Phi.pen=NULL,tau=1,NRadapt=FALSE,ncount=1,kstep=1,Hbeta.old=NULL){

###para: 
###x0.old: vector for current initial conditions, vector of length n (number of curves) 
###theta.old: vector for current random scale parameters, vector of length num.id (number of clusters)
###beta.old: current basis coefficients for g; 
###subid : cluster index
###data.path.u : measurements
###obstime.u: measurement times;
###N.u: number of grid points; h.u: grid step size;  knots.u: vector of knots
### Phi.pen: quadratic matrix penalty for beta
###tau: damping factor (<=1)
### NRadapt,ncount,Hbeta.old,kstep: only update  beta by NR if(NRadapt==FALSE  || ncount %% kstep == 1 || is.null(Hbeta.old))
##return: updated estimates of beta and its gradient

#####
  id.list<-unique(subid)                   ## cluster indices
  n.u <- length(obstime.u)                 ##sample size: total number of curves 
  M.u <- length(beta.old)                  ##number of basis  
  h2 <- 2*h.u
  h4 <- 4*h.u

#####
  y.tilde <- xpath.N(N.u,h.u,x0.old,theta.old,beta.old,knots.u)        ## sample path evaluated at current estimates
  w.tilde <- xderiv.beta.N(N.u,h.u,theta.old,beta.old,knots.u,y.tilde)       ## gradient w.r.t. beta                  #####UPDATE NEEDED!!!##########
  
  
#####
if(is.null(Phi.pen)){
Phi.pen=matrix(0,M.u,M.u)
}

m.u <- numeric(n.u)
vareps.est <- 0
Jbeta.eps <- numeric(M.u)
Hbeta.mat <- matrix(0,M.u,M.u)             
eps.tilde <- NULL

  for(i in 1:n.u){
     grid<-ceiling(obstime.u[[i]]/h2)
     grid2<-ceiling(obstime.u[[i]]/h4)
     m.u[i] <- length(obstime.u[[i]])   ##number of measurements for ith curve
     eps.tilde[[i]] <- data.path.u[[i]] - y.tilde[i,ceiling(obstime.u[[i]]/h.u)]   ## current residuals 
     vareps.est <- vareps.est + sum(eps.tilde[[i]]^2)
     
     for(l in 1:M.u){
         Jbeta.eps[l] <- Jbeta.eps[l] + sum(eps.tilde[[i]]*w.tilde[i,grid,l])
     }

  }    #### end of loop for i

    
##### computation of Hessian w.r.t. beta: NRadapt=TRUE: only update Hessian for every ncount steps

if(NRadapt==FALSE  || ncount %% kstep == 1 || is.null(Hbeta.old)){
   
  ww.tilde <- xHess.beta.N(N.u,h.u,theta.old,beta.old,knots.u,y.tilde,w.tilde)     ### Hessian w.r.t. beta;           #####UPDATE NEEDED!!!##########

  for(i in 1:n.u){
     grid<-ceiling(obstime.u[[i]]/h2)
     grid2<-ceiling(obstime.u[[i]]/h4)
     for(l in 1:M.u){
        for(k in 1:l){
     Hbeta.mat[l,k] <- Hbeta.mat[l,k] + sum(w.tilde[i,grid,l]*w.tilde[i,grid,k]) - sum(ww.tilde[i,grid2,l,k]*eps.tilde[[i]])
           Hbeta.mat[k,l]<- Hbeta.mat[l,k]
        }      
     }     
  }   #### end of loof for i
  Hbeta.matinv <- solve(Hbeta.mat + Phi.pen)            ###inverse of Hessian 
}
else{
    Hbeta.matinv <- Hbeta.old
}


  beta.new <- beta.old + tau * Hbeta.matinv %*% (Jbeta.eps - Phi.pen %*% beta.old)  ## Newton-Raphson update for beta with matrix penalty and damping

  grad.beta <- -2*(Jbeta.eps - Phi.pen %*% beta.old)/n.u #### (averaged) gradient of the loss function w.r.t. beta           #####????######
 
  df<-sum(m.u-1)-length(id.list)-M.u
  df.fix<-sum(m.u)-length(id.list)-M.u

  return(list(beta.new,grad.beta,vareps.est,Hbeta.matinv,df,df.fix))

}


######(c)
###### NR algorithm for model fitting 
FitNR.N<-function(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,lambda2.ini,Phi.pen=NULL,adapt1=TRUE, adapt2=TRUE,ini.fix=FALSE,niter=10,thre=1e-3,conv.cri="grad",tau.v=matrix(1,3,niter),NRadapt=FALSE,kstep=1,Hbeta.old=NULL){
#####additional paramters:

##adapt1, adapt2: True (adatively updated lambda1 and/or lambda 2); FALSE (non-adaptive) 
##ini.fix: TRUE (x0 viewed as known); FALSE (need to update x0)
##niter: maximum number of iterations 
##thre: threshold for convergence; 
##conv.cri: "beta" (maximum change of beta and number of iterations);
##or "grad" (maximum gradient and number of iterations)
##tau.v = dampening factors at each step: three by niter matrix 

#####
x0.u <- NULL
theta.u <-NULL
beta.u <-NULL
vareps <-0
gradient<-NULL
num.id <- length(table(subid))


#####update x0 
count<-1
print(paste("NR step:", count))
if(ini.fix==FALSE){
temp<-PLS.ini.N(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,tau=tau.v[1,1])
x0.new<-temp[[1]]
}else{
x0.new<-x0.ini
}


##### update theta
temp <-NR.theta.N(x0.new,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda2.ini,tau=tau.v[2,1])
theta.new <- temp[[1]]

###update beta
temp <-NR.beta.N(x0.new,theta.new,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,Phi.pen,tau=tau.v[3,1],NRadapt,ncount=1,kstep,Hbeta.old)

beta.new <-temp[[1]]
if(ini.fix==FALSE){
vareps.new<-temp[[3]]/temp[[5]]
}else{
vareps.new<-temp[[3]]/temp[[6]]
}
Hbeta.old<-temp[[4]]

################update step
beta.old <- beta.new+100
gradsup.norm <- 100
##Hbeta.est <- NULL           ### Hessian matrix w.r.t. beta

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
temp<-PLS.ini.N(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1,tau=tau.v[1,count])
x0.new<-temp[[1]]
}else{
x0.new<-x0.ini
}


##### update theta
temp <-NR.theta.N(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda2,tau=tau.v[2,count])
theta.new <- temp[[1]]
grad.T <- temp[[2]]


###update beta
temp <-NR.beta.N(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,Phi.pen,tau=tau.v[3,count],NRadapt,ncount=count,kstep,Hbeta.old)

beta.new <-temp[[1]]
if(ini.fix==FALSE){
vareps.new<-temp[[3]]/temp[[5]]
}else{
vareps.new<-temp[[3]]/temp[[6]]
}
Hbeta.old<-temp[[4]]
grad.beta <- temp[[2]]

##### 
 ##print(paste("beta",beta.new))
 gradsup.norm <- c(max(abs(grad.T)),max(abs(grad.beta)))
 gradient[[count]]<- gradsup.norm
 ##print(paste("sup norm of gradient of loss function ="))
 ##print(gradsup.norm)

#####
 count<-count+1
 print(paste("NR step:", count))

 if(conv.cri=="beta"){
    conv.stat<-(sqrt(sum((beta.old - beta.new)^2)) > thre && (count<=niter))
    }
    
    if(conv.cri=="grad"){
    conv.stat<-(max(gradsup.norm) > thre && (count<=niter))
    }
 
 
 }  ## end of while


####
x0.u[[count]]<-x0.new
 theta.u[[count]]<-theta.new
 beta.u[[count]]<-beta.new
 vareps[count]<-vareps.new
 gradient[[count]]<- gradsup.norm
 
if(conv.cri=="beta"){
   conv.stat<-(sqrt(sum((beta.old - beta.new)^2)) > thre )
   }
   
   if(conv.cri=="grad"){
   conv.stat<-(max(gradsup.norm) > thre)
   }

return(list(x0.u,theta.u,beta.u,vareps,count,gradient,lambda1,lambda2,1-conv.stat))

}
