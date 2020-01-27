### 04-21-08
### Implementation of fourth order Runge-Kutta method for solving an ODE

###R CMD COMPILE RK4.c
###R CMD SHLIB -o RK4.so RK4.o

##########################R wrap functions 
xpath <- function(N.u,h.u,x0.u,theta.u,beta.u,knots.u){
  Lx=length(x0.u)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u+1)
  result.v<-as.vector(result)
   ##dyn.load("RK4.so")  
 
   junk=.C("xPath",
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
xderiv.a <- function(N.u,h.u,theta.u,beta.u,knots.u,y){

  Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/2+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))

   ##dyn.load("RK4.so")  
 
   junk=.C("xDerivA",
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
xderiv.theta <- function(N.u,h.u,theta.u,beta.u,knots.u,y){
  Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/2+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))

   ##dyn.load("RK4.so")  
 
   junk=.C("xDerivTheta",
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
xderiv.beta <- function(N.u,h.u,theta.u,beta.u,knots.u,y){


  Lx=nrow(y)
  Mk=length(knots.u)
  result=array(0,dim=c(Lx,N.u/2+1,Mk))
  result.v<-as.vector(result)
  y.v=as.vector(t(y))

   ###dyn.load("RK4.so")  
 
   junk=.C("xDerivBeta",
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

  result=array(junk$result.out, dim=c(Lx,N.u/2+1,Mk))
  ##result=junk$result.out
  return(result)
}


xderiv <- function(N.u,h.u,theta.u,beta.u,knots.u,y){

###para: N.u: number of grid points; h.u: grid step size; 
###theta.u: vector for random scale parameters, vector of length n.u
###beta.u: Bspline coefficients for g; knots.u: Bspline knots for g 
###y: (estimated) sample path X_i(.), i=1...n.u: n.u by N.u+1 matrix from "xpath"

 z<-xderiv.a(N.u,h.u,theta.u,beta.u,knots.u,y)
 v<-xderiv.theta(N.u,h.u,theta.u,beta.u,knots.u,y)

  return(list(z,v))


}


xHess.theta <- function(N.u,h.u,theta.u,beta.u,knots.u,y,v){
 Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/4+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))
  v.v=as.vector(t(v))

  
   ###dyn.load("RK4.so")  
 
   junk=.C("xHessTheta",
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



xHess.A <- function(N.u,h.u,theta.u,beta.u,knots.u,y,z){
 Lx=nrow(y)
  Mk=length(knots.u)
  result=matrix(0,Lx,N.u/4+1)
  result.v<-as.vector(result)
  y.v=as.vector(t(y))
  z.v=as.vector(t(z))

  
   ###dyn.load("RK4.so")  
 
   junk=.C("xHessA",
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

xHess <- function(N.u,h.u,theta.u,beta.u,knots.u,y,z,v){

zz<-xHess.A(N.u,h.u,theta.u,beta.u,knots.u,y,z)
vv<-xHess.theta (N.u,h.u,theta.u,beta.u,knots.u,y,v)

return(list(zz,vv))
}




##xHessBeta(int *NN, float *hh, float *y, float *w, float *theta, float *beta, float *knot, int *Lxx, int *Mkk, float *result)
xHess.beta <- function(N.u,h.u,theta.u,beta.u,knots.u,y,w){
  
  Lx=nrow(y)
  Mk=length(knots.u)
  result=array(0,dim=c(Lx,N.u/4+1,Mk,Mk))
  result.v<-as.vector(result)
  y.v=as.vector(t(y))
  w.v=as.vector(w)

  
   ###dyn.load("RK4.so")  
 
   junk=.C("xHessBeta",
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

   result=array(junk$result.out, dim=c(Lx,N.u/4+1,Mk,Mk))
   return(result)

}


#####################################################
## Link FUNCTION : g(x) = \sum_{k=1}^M \beta_k \phi_k(x), where \phi_k are 
###                    cubic B-spline basis functions with equally spaced knots
######################################################

g1 <- function(x,beta.u,knots.u){  ### evaluation of g(x)
    M.u <- length(beta.u)
    eval <- BC.basis(x,knots.u)    
    r <- eval %*% matrix(beta.u)
    return(r)   
}

g2 <- function(x,beta.u,knots.u){  ### evaluation of g'(x)
    M.u <- length(beta.u)
    eval <- BPC.basis(x,knots.u)    
    r <- eval %*% matrix(beta.u)
    return(r)   
}




############################
## cubic b-spline Basis
##
##############################
B.one<-function(x){
return(x^3/6)
}

##
B.two<-function(x){
result<-(-3*x^3+3*x^2+3*x+1)/6
return(result)
}

##
B.three<-function(x){
result<-(3*x^3-6*x^2+4)/6
return(result)
}

##
B.four<-function(x){
result<-(1-x)^3/6
return(result)
}

##
B.cubic<-function(x){
 result<-0

 if(x>=-3&&x<=-2){
 result<-B.one(x+3)
 }
 
 if(x>=-2&&x<=-1){
 result<-B.two(x+2)
 }

 if(x>=-1&&x<=0){
 result<-B.three(x+1)
 }

 if(x>=0&&x<=1){
 result<-B.four(x)
 }

 return(result)
} 

##
BC.basis<-function(t,knots.u){

##para:t--evaluation points; knots: equally spaced knots 

M.u<-length(knots.u)
delta<-knots.u[2]-knots.u[1]
T.mat <- matrix(t,length(t),M.u)
K <- t(matrix(knots.u,M.u,length(t)))
T.v<-(T.mat-K)/delta

result<-t(matrix(apply(matrix(t(T.v)),1,B.cubic),M.u,length(t)))
return(result)
}



######################################
## derivatives of cubic b-spline Basis
##
#################################
BP.one<-function(x){
return(x^2/2)
}

##
BP.two<-function(x){
result<- (-9*x^2+6*x+3)/6
return(result)
}

##
BP.three<-function(x){
result<- (9*x^2-12*x)/6
return(result)
}

##
BP.four<-function(x){
result<- -(1-x)^2/2
return(result)
}

##
BP.cubic<-function(x){
 result<-0

 if(x>=-3&&x<=-2){
 result<-BP.one(x+3)
 }
 
 if(x>=-2&&x<=-1){
 result<-BP.two(x+2)
 }

 if(x>=-1&&x<=0){
 result<-BP.three(x+1)
 }

 if(x>=0&&x<=1){
 result<-BP.four(x)
 }

 return(result)
} 

##
BPC.basis<-function(t,knots.u){

##para:t--evaluation points; knots: equally spaced knots 

M.u<-length(knots.u)
delta<-knots.u[2]-knots.u[1]
T.mat <- matrix(t,length(t),M.u)
K <- t(matrix(knots.u,M.u,length(t)))
T.v<-(T.mat-K)/delta

##result<-t(matrix(apply(matrix(t(T.v)),1,BP.cubic),M.u,length(t)))
result<-t(matrix(apply(matrix(t(T.v)),1,BP.cubic),M.u,length(t)))/delta  ##correct on 04-21-2011!!!
return(result)
}

 