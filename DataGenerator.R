#rm(list=ls())
#n=200
#sigma=1
source("floyd.R")
DataGenerator<-function(n,alpha=0.5){
  u=runif(n,0,1)
  v=runif(n,0,1)
  loc=cbind(u,v)
  s=u^2+v^2
  beta1=sin(pi*s)
  beta2=cos(pi*s)
  beta3=exp(s)
  beta=cbind(beta1,beta2,beta3)
  x1=rnorm(n)
  x2=rnorm(n)
  x3=rnorm(n)
  X=cbind(x1,x2,x3)
  eps=rnorm(n,0,sigma)
  Ystar=x1*beta1+x2*beta2+x3*beta3+eps
  # Calculate the pairwise Euclidean distances and Weights
  
  temp1=t(kronecker(loc,rep(1,n)))
  temp2=t(loc[rep(1:n,n),])
  DE=matrix(sqrt(colSums((temp1-temp2)^2)),n,n)
  temp1=rowSums(exp(-DE))-1
  temp2=matrix(rep(temp1,n),n,n)
  W=exp(-DE)/temp2
  A=diag(n)-alpha*W
  Y=solve(A)%*%Ystar
  list(Y=Y,beta=beta,X=X,Z=loc)
}

# # Boundary
# bnd=rbind(c(0,0),c(1,0),c(1,1),c(0,1))
# bnd=as.data.frame(bnd)
# colnames(bnd)=c("x","y")
# 
# loc=as.data.frame(loc)
# colnames(loc)=c("x","y")
# knots=loc
# D <- msg:::create_distance_matrix(c(loc$x,knots$x),c(loc$y,knots$y),bnd,faster=1)
# # distances from data to knots
# Dg1=D[1:N,(N+1):dim(D)[2]]

