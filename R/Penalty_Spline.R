#### 10-05-08 (latest)


### function to compute the penalty weight matrix for beta in the plant dynamics problem
### penalty is of the form
### lambda.u \int_{A.u}^{B.u} \phi'(x) (\phi'(x))^T dx  
### where lambda.u > 0, B.u > A.u > 0 (appropriately chosen)
### and \phi(x) = [\phi_1(x),...,\phi_M(x)]^T  where \phi_j(x) is the j-th basis function


Bspline.penalty <- function(knots.u,A.u,B.u,lambda.u,step=0.0001){

grid.u <- seq(A.u,B.u,step)

eval <- BPC.basis(grid.u,knots.u)

penalty.mat <- lambda.u * t(eval) %*% eval /length(grid.u)*(B.u-A.u)

return(penalty.mat)

}


Nspline.penalty <- function(knots.u,A.u,B.u,lambda.u,step=0.0001){

grid.u <- seq(A.u,B.u, step)

eval <- NPC.basis(grid.u,knots.u)

penalty.mat <- lambda.u * t(eval) %*% eval /length(grid.u)*(B.u-A.u)

return(penalty.mat)

}

