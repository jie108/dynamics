### 04-21-08
### Implementation of fourth order Runge-Kutta method for solving an ODE

###R CMD COMPILE RK4_natural.c
###R CMD SHLIB -o RK4_natural.so RK4_natural.o

##########################R wrap functions 
xpath.N <- function(N.u,h.u,x0.u,theta.u,beta.u,knots.u){
  Lx=length(x0.u)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u+1)
  result.v<-as.vector(result)
   ###dyn.load("RK4_natural.so")  
 
   junk=.C("xPathN",
          as.integer(N.u),
          as.double(h.u),
          as.double(x0.u),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

   result=matrix(junk$result.out, nrow=Lx, ncol=N.u+1, byrow=TRUE)
   return(result)


}


###
xderiv.a.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y){

  Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/2+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))

  ###dyn.load("RK4_natural.so")  
 
   junk=.C("xDerivAN",
          as.integer(N.u),
          as.double(h.u),
          as.double(y.v),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

   result=matrix(junk$result.out, nrow=Lx, ncol=N.u/2+1,byrow=TRUE)
   return(result)
}



###
xderiv.theta.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y){
  Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/2+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))

  ###dyn.load("RK4_natural.so")  
 
   junk=.C("xDerivThetaN",
          as.integer(N.u),
          as.double(h.u),
          as.double(y.v),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

   result=matrix(junk$result.out, nrow=Lx, ncol=N.u/2+1,byrow=TRUE)
   return(result)

}


##
xderiv.beta.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y){


  Lx=nrow(y)
  Mk=length(knots.u)
  result=array(0,dim=c(Lx,N.u/2+1,Mk+3))
  result.v<-as.vector(result)
  y.v=as.vector(t(y))

  ###dyn.load("RK4_natural.so")  
 
   junk=.C("xDerivBetaN",
          as.integer(N.u),
          as.double(h.u),
          as.double(y.v),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

  result=array(junk$result.out, dim=c(Lx,N.u/2+1,Mk+3))
  ##result=junk$result.out
  return(result)
}


xderiv.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y){

###para: N.u: number of grid points; h.u: grid step size; 
###theta.u: vector for random scale parameters, vector of length n.u
###beta.u: Bspline coefficients for g; knots.u: Bspline knots for g 
###y: (estimated) sample path X_i(.), i=1...n.u: n.u by N.u+1 matrix from "xpath"

 z<-xderiv.a.N(N.u,h.u,theta.u,beta.u,knots.u,y)
 v<-xderiv.theta.N(N.u,h.u,theta.u,beta.u,knots.u,y)

  return(list(z,v))


}


xHess.theta.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y,v){
 Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/4+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))
  v.v=as.vector(t(v))

  
   ###dyn.load("RK4_natural.so")  
 
   junk=.C("xHessThetaN",
          as.integer(N.u),
          as.double(h.u),
          as.double(y.v),
          as.double(v.v),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

   result=matrix(junk$result.out, nrow=Lx, ncol=N.u/4+1,byrow=TRUE)
   return(result)

}



xHess.A.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y,z){
 Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/4+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))
  z.v=as.vector(t(z))

  
   ###dyn.load("RK4_natural.so")  
 
   junk=.C("xHessAN",
          as.integer(N.u),
          as.double(h.u),
          as.double(y.v),
          as.double(z.v),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

   result=matrix(junk$result.out, nrow=Lx, ncol=N.u/4+1,byrow=TRUE)
   return(result)

}

xHess.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y,z,v){

zz<-xHess.A.N(N.u,h.u,theta.u,beta.u,knots.u,y,z)
vv<-xHess.theta.N (N.u,h.u,theta.u,beta.u,knots.u,y,v)

return(list(zz,vv))
}




xHess.beta.N <- function(N.u,h.u,theta.u,beta.u,knots.u,y,w){
  
  Lx=nrow(y)
  Mk=length(knots.u)
  result=array(0,dim=c(Lx,N.u/4+1,Mk+3,Mk+3))
  result.v<-as.vector(result)
  y.v=as.vector(t(y))
  w.v=as.vector(w)

  
   ###dyn.load("RK4_natural.so")  
 
   junk=.C("xHessBetaN",
          as.integer(N.u),
          as.double(h.u),
          as.double(y.v),
          as.double(w.v),
          as.double(theta.u),
          as.double(beta.u),
          as.double(knots.u),
          as.integer(Lx),
          as.integer(Mk),
          result.out=as.double(result.v)
          )   

   result=array(junk$result.out, dim=c(Lx,N.u/4+1,Mk+3,Mk+3))
   return(result)

}



#####################################################
## Link FUNCTION : g(x) = \sum_{k=1}^M \beta_k \phi_k(x), where \phi_k are 
###                    natural-spline basis functions given knot sequence
######################################################

g1.N <- function(x,beta.u,knots.u){  ### evaluation of g(x)
    eval <- NC.basis(x,knots.u)    
    r <- eval %*% matrix(beta.u)
    return(r)   
}

g2.N <- function(x,beta.u,knots.u){  ### evaluation of g'(x)
    eval <- NPC.basis(x,knots.u)    
    r <- eval %*% matrix(beta.u)
    return(r)   
}

g3.N <- function(x,beta.u,knots.u){  ### evaluation of g''(x)
    eval <- NHC.basis(x,knots.u)    
    r <- eval %*% matrix(beta.u)
    return(r)   
}

#########
## NOTE : first 3 coordinates of the vector beta.u are the coefficients of x, x^2 and x^3, and the rest of  
##        the coefficients correspond to the functions (x-knots.u)_+^3
###      M.u = length(beta.u) = 3 + L.u ; L.u = length(knots.u)
########


##
NC.basis<-function(t,knots.u){

##para:t--evaluation points; knots: equally spaced knots 

L.u<-length(knots.u)
M.u <- L.u + 3

T.mat <- matrix(t,length(t),M.u)
T.mat[,2] <- t^2
T.mat[,3] <- t^3
for(l in 1:L.u){
  T.mat[,(l+3)] <- (t-knots.u[l])^3 * (t > knots.u[l])
}

return(T.mat)
}



NPC.basis<-function(t,knots.u){

##para:t--evaluation points; knots: equally spaced knots 

L.u<-length(knots.u)
M.u <- L.u + 3

T.mat <- matrix(1,length(t),M.u)
T.mat[,2] <- 2*t
T.mat[,3] <- 3*t^2
for(l in 1:L.u){
  T.mat[,(l+3)] <- 3*(t-knots.u[l])^2 * (t > knots.u[l])
}

return(T.mat)
}




##
NHC.basis<-function(t,knots.u){

##para:t--evaluation points; knots: equally spaced knots 

L.u<-length(knots.u)
M.u <- L.u + 3

T.mat <- matrix(0,length(t),M.u)
T.mat[,2] <- 2
T.mat[,3] <- 6*t
for(l in 1:L.u){
  T.mat[,(l+3)] <- 6*(t-knots.u[l]) * (t > knots.u[l])
}

return(T.mat)
}
