### 07-10-08 
### 10-08-08 (latest)
### for R package: 1-19-09

########################################################################
## approximation to CV score based on leave-one-curve-out cross validation 
#########################################################################

######################################################

CVCUR.N<- function(x0.final,theta.final,beta.final,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.final,lambda2.final,ini.fix=FALSE, maxloop=10,thre=0.0001,Phi.pen=NULL){

##############################################
## subid : index of plant ID;   
## x0.final : final estimates of x0 (or, a);  theta.final : final estimate of theta; beta.final : final estimate of beta
## (lambda1.final,lambda2.final) : final estimates of (lambda1,lambda2)
##  maxloop = maximum number of loops allowed for estimating a_i^{(-il)} : default value = 1 (one-step apprxoimation) 
##  Phi.pen = penalty matrix for beta
## result: two CV scores: one based on prediction error; one based on the loss function
#############################################

  id.list<-unique(subid)
  num.id<-length(id.list)
  alpha.final <- mean(x0.final)           ##final estimate for mean of x0
  
  y.tilde <- xpath.N(N.u,h.u,x0.final,theta.final,beta.final,knots.u)        ### sample path evaluated at current estimates
  temp <- xderiv.N(N.u,h.u,theta.final,beta.final,knots.u,y.tilde)
  z.tilde <- temp[[1]]                                          ### gradient w.r.t. a 
  v.tilde <- temp[[2]]                                          ### gradient w.r.t. theta  
  w.tilde <- xderiv.beta.N(N.u,h.u,theta.final,beta.final,knots.u,y.tilde)       ### gradient w.r.t. beta

  temp2 <- xHess.N(N.u,h.u,theta.final,beta.final,knots.u,y.tilde,z.tilde,v.tilde)
  zz.tilde <- temp2[[1]]                                        ### Hessian w.r.t. a 
  vv.tilde <- temp2[[2]]                                        ### Hessian w.r.t. theta
  ww.tilde <- xHess.beta.N(N.u,h.u,theta.final,beta.final,knots.u,y.tilde,w.tilde)     ### Hessian w.r.t. beta

  n.u <- length(obstime.u)                 ##sample size: number of curves 
  M.u <- length(beta.final)                ##number of basis functions
  h2 <- 2*h.u
  h4 <- 4*h.u


if(is.null(Phi.pen)){
Phi.pen=matrix(0,M.u,M.u)
}

#####
m.u <- numeric(0)
mloc <- numeric(0)
J.T <- NULL
H.T <- NULL
J2.T <- NULL
H.TT <- NULL
Jbeta.eps <- matrix(0,n.u,M.u)
Jbeta2.mat <- matrix(0,M.u,M.u)
Hbeta.mat <- matrix(0,M.u,M.u)
eps.tilde <- NULL

JA.drop <- NULL
J.AA.drop <- 0
eps.drop <- NULL
theta.drop <- numeric(0)
beta.drop <- matrix(0,M.u,n.u)
x0.drop.old <- 1000*x0.final
x0.drop.new <- x0.final
x0.drop <- numeric(0)

  cvcur.score <- 0   ## approximate CV score (based on the loss function, that is including the penalty terms) 
  cvcur.alt <- 0     ## approximate CV score (based on prediction error)

  for(i in 1:n.u){
     m.u[i] <- length(obstime.u[[i]])   ##number of measurements for the ith curve

     grid<-ceiling(obstime.u[[i]]/h2)
     grid2<-ceiling(obstime.u[[i]]/h4)

     eps.tilde[[i]] <- data.path.u[[i]] - y.tilde[i,ceiling(obstime.u[[i]]/h.u)]   ### current residuals 

     J.T[[i]] <- v.tilde[i,grid]   ### derivative w.r.t. theta_i evaluated at the observed times 
     H.T[[i]] <- vv.tilde[i,grid2]   ### Hessian w.r.t. theta_i evaluated at the observed times


     for(l in 1:M.u){
         Jbeta.eps[i,l] <- sum(eps.tilde[[i]]*w.tilde[i,grid,l])
     }

    #### computation of Hessian w.r.t. beta

     for(l in 1:M.u){
        for(k in 1:M.u){
           Jbeta2.mat[l,k] <- Jbeta2.mat[l,k] + sum(eps.tilde[[i]]^2 * w.tilde[i,grid,l] * w.tilde[i,grid,k])
           Hbeta.mat[l,k] <- Hbeta.mat[l,k] + sum(w.tilde[i,grid,l]*w.tilde[i,grid,k]) - sum(ww.tilde[i,grid2,l,k]*eps.tilde[[i]])
        }      
     }     
  }    #### end of for i


##print(Jbeta2.mat)
##print(Hbeta.mat)

  for(k in 1:num.id){    #### k = cluster index,  l = location index

    locindex <- which(subid == id.list[k])    ### index of locations which correspond to cluster k

    mloc[k] <- 0
    for(l in 1:length(locindex)){
        mloc[k] <- mloc[k] + m.u[locindex[l]]        ##### sum_{l=1}^{L_i} m_{il}  ##total number of measurements of the kth cluster
    }
    
    H.TT[[k]] <- lambda2.final   
    J2.T[[k]] <- 0
    for(l in 1:length(locindex)){
      i <- locindex[l]
      H.TT[[k]] <- H.TT[[k]] + sum(J.T[[i]]^2) - sum(eps.tilde[[i]]*H.T[[i]]) 
      J2.T[[k]] <- J2.T[[k]] + sum((eps.tilde[[i]]*J.T[[i]] - lambda2.final * theta.final[k]/mloc[k])^2)
    }
  }


##### estimating a_i for approximate theta_i and beta after dropping  the i-th curve    

  for(i in 1:n.u){
   print(i)
   theta.drop[i] <- theta.final[i] - sum(eps.tilde[[i]]*J.T[[i]])/H.TT[[which(id.list==subid[i])]]
   beta.drop[,i] <- beta.final - solve(Hbeta.mat+Phi.pen) %*% (Jbeta.eps[i,])   
   print(theta.drop[i])
   print(beta.drop[,i])
   
#### esimate x0.drop 
   if(ini.fix==FALSE){
     grid<-ceiling(obstime.u[[i]]/h2)
     loopcount <- 0
     while((abs(x0.drop.old[i] - x0.drop.new[i])/sd(x0.final) > thre)  && (loopcount <= maxloop)){  
       
       loopcount <- loopcount + 1
       x0.drop.old[i] <- x0.drop.new[i]
     
       print(x0.drop.old[i])
       y.drop <- xpath.N(N.u,h.u,x0.drop.old[i],theta.drop[i],beta.drop[,i],knots.u)  
       ###print(summary(y.drop))
       z.drop <- xderiv.a.N(N.u,h.u,theta.drop[i],beta.drop[,i],knots.u,y.drop)  ### gradient w.r.t. a                   
    
       JA.drop[[i]] <- z.drop[grid] 
       eps.drop[[i]] <- data.path.u[[i]] - y.drop[ceiling(obstime.u[[i]]/h.u)]   ### current residuals 
 
       J.AA.drop <- sum(JA.drop[[i]]^2)
       udiff.drop <- sum(JA.drop[[i]]*eps.drop[[i]]) + lambda1.final * (alpha.final - x0.drop.old[i])   ##RHS updates for x0
       x0.drop.new[i] <- x0.drop.old[i] + (J.AA.drop + lambda1.final)^{-1} * udiff.drop          ##new estimates

     }  ## end of while
     
     x0.drop[i] <- x0.drop.new[i] 
    }else{ ###x0 known, so no need to update
    
     x0.drop[i] <- x0.final[i]
    }
    
     y.drop <- xpath.N(N.u,h.u,x0.drop[i],theta.drop[i],beta.drop[,i],knots.u)       
     eps.drop[[i]] <- data.path.u[[i]] - y.drop[ceiling(obstime.u[[i]]/h.u)]  #### residuals with most updated value of a^{-(il)}, theta^{-(il)} and beta^{-(il)}
     print(cvcur.alt) 
     cvcur.alt <- cvcur.alt + sum(eps.drop[[i]]^2)     ###CV score based on prediction error 
  }  ## end of for
  
  cvcur.score <- cvcur.alt + lambda1.final * sum((x0.drop - alpha.final)^2) + lambda2.final * sum(theta.drop^2)   ###CV score based on loss function (without matrix penalty term)

#####

return(list(cvcur.alt,cvcur.score,x0.drop,theta.drop,beta.drop))

}

