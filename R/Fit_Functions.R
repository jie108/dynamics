####11-16-08
#### Fit the nonlinear dynamics model 
#### using a combination of Levenberg-Marquardt and Newton-Raphson schemes


Fit <- function(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,lambda2.ini,lambda3.v,Phi.pen=NULL,
adapt1.LM=FALSE, adapt2.LM=FALSE,adapt1.NR=TRUE, adapt2.NR=TRUE, ini.fix=FALSE, niter.LM=10,thre.LM=1e-3,conv.cri.LM="beta",niter.NR=10,thre.NR=1e-3,conv.cri.NR="grad",
tau.v.LM=matrix(1,3,niter.LM),tau.v.NR=matrix(1,3,niter.NR),NRadapt=FALSE,kstep=1,Hbeta.old=NULL){

#######(a) LM step
   print(paste("Beginning of Levenberg-Marquardt procedure"))
   result.LM <- try(FitLM(x0.ini,theta.ini,beta.ini,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.ini,lambda2.ini,lambda3.v,Phi.pen,adapt1.LM, adapt2.LM, ini.fix, niter.LM,thre.LM, conv.cri.LM,tau.v.LM))
   error.c <- inherits(result.LM,"try-error")
   if(error.c){
         print("Error in Levenberg-Marquardt step")
         result<- NULL
         result[[10]] <- error.c
         return(result) 
   }
   
   print(paste("End of Levenberg-Marquardt procedure"))
   print(paste("converged?",as.logical(result.LM[[8]])))
   
   
   x0.LM <- result.LM[[1]]
   theta.LM <- result.LM[[2]]
   beta.LM <- result.LM[[3]] 
   vareps.LM <- result.LM[[4]]
   count.LM <- result.LM[[5]]
   
   x0.new<-x0.LM[[count.LM]]
   theta.new<-theta.LM[[count.LM]]
   beta.new<-beta.LM[[count.LM]]
   lambda1.new <- result.LM[[6]]
   lambda2.new <- result.LM[[7]]

######(b) NR step
   print(paste("Beginning of Newton-Raphson procedure"))
   result.NR <- try(FitNR(x0.new,theta.new,beta.new,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.new,lambda2.new,Phi.pen,adapt1.NR, adapt2.NR,ini.fix,niter.NR,thre.NR,conv.cri.NR,tau.v.NR,NRadapt,kstep,Hbeta.old))

   error.c <- inherits(result.NR,"try-error")
      if(error.c){
            print("Error in Newton step")
            result<-result.LM
            ##result[[10]] <- error.c
            return(result) 
      }
   
  
   print(paste("End of Newton-Raphson procedure"))
   print(paste("converged?",as.logical(result.NR[[9]])))
   
   x0.NR <- result.NR[[1]]
   theta.NR <- result.NR[[2]]
   beta.NR <- result.NR[[3]]
   vareps.NR <- result.NR[[4]]
   count.NR <- result.NR[[5]]
   gradient.NR<-result.NR[[6]]
   lambda1.NR <- result.NR[[7]]
   lambda2.NR <- result.NR[[8]]

#####
result.f<-list(x0.NR[[count.NR]],theta.NR[[count.NR]],beta.NR[[count.NR]],vareps.NR[[count.NR]], gradient.NR[[count.NR]], count.LM, count.NR,lambda1.NR,lambda2.NR,error.c)
names(result.f)<-c("initial_condition","scale_parameter","basis_coefficient","error_variance","maximum_gradient","LM_steps","NR_steps","lambda1","lambda2","error_occurence?")
return(result.f)

}
