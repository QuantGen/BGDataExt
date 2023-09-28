# A wrappper for ridge regression, possibly shrinking towards b0
 RR <- function(XX, Xy,lambda,p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda0=0, active=1:p, RSS=10000, maxIter=1000, tol=1e-5) {
     
    # adding the shrinkage parameter to the left-hand-side
    diag(XX)=diag(XX)+lambda

     # adding prior mean to the right-hand-side
    if(lambda2>0){
      Xy=Xy+lambda0*b0
    }
    
    ans<-fitSYS(XX, Xy, b, active, RSS, maxIter, tol)[[1]] 
   
    return(ans)
 }


if(FALSE){
library(BGDataExt)

## Ridge-regression with 'prior mean', i am using this to check the internal function in BGDataExt
 RR2=function(XX,Xy,lambda,b0=rep(0,ncol(XX)),lambda0=0){
	diag(XX)=diag(XX)+lambda
	rhs=Xy+lambda0*b0
	sol=solve(XX,rhs)
	return(sol)
 }
 
 

## Data 
 library(BGLR)
 data(wheat)
 X=scale(wheat.X,center=TRUE,scale=FALSE)
 y=wheat.Y[,1]

XX=crossprod(X)
Xy=crossprod(X,y)

lambda=sum(diag(XX))/nrow(X)


# Checking function when shrinking towards 0
 system.time(bHat<-RR(XX,Xy,lambda))

 system.time(bHat2<-RR2(XX,Xy,lambda))
 plot(bHat2,bHat,cex=.5,col=4);abline(a=0,b=1,col=2)

## Clustering 
 SVD=svd(X,nu=5,nv=0)
 group=kmeans(SVD$u,centers=2,nstart=100)$cluster
 plot(SVD$u[,1:2],col=group*2)

## Data split 
 group1=which(group==1)
 group2=which(group==2)

 X1=X[group1,]
 y1=y[group1]
 X2=X[group2,]
 y2=y[group2]
 
 bHat1=RR2(crossprod(X1),crossprod(X1,y1),lambda) # srhinkage towards 0
 

# Now shrinking towards some value different than 0
 XX=crossprod(X2)
 Xy=crossprod(X2,y2)
 
 b0<-RR2(XX,Xy,lambda=lambda)
 system.time(bHat1<-RR( XX,Xy,lambda=lambda,lambda2=lambda,b0=bHat1))
 system.time(bHat2<-RR2(XX,Xy,lambda=lambda,lambda0=lambda,b0=bHat1))

plot(bHat2,bHat1,cex=.5,col=4);abline(a=0,b=1,col=2)




 }
