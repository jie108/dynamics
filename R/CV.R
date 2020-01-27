#####
dynamics.cv<-function(x0.final,theta.final,beta.final,subid=NULL,data.path.u,obstime.u,N.u=1000,method="bspline",knots.u,lambda1.final,lambda2.final,Phi.pen=NULL,ini.fix=FALSE, maxloop=10,thre=0.0001){

####
h.u<-1/N.u

####
n.u <- length(obstime.u)                 ##sample size: total number of curves 
if(is.null(subid)){
subid<-1:n.u                           ##default for subid: one curve per cluster, that is no clusters
}


if(0){
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
}

####CV
if(method=="bspline"){
result<-CVCUR(x0.final,theta.final,beta.final,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.final,lambda2.final,ini.fix, maxloop,thre,Phi.pen)
}

if(method=="natural"){
result<-CVCUR.N(x0.final,theta.final,beta.final,subid,data.path.u,obstime.u,N.u,h.u,knots.u,lambda1.final,lambda2.final,ini.fix, maxloop,thre,Phi.pen)
}


####return
names(result)<-c("cv_prediction_error", "cv_loss_function", "initial_condition", "scale_parameter", "basis_coefficient")
return(result)
}
