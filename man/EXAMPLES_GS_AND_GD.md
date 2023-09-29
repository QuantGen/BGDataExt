```r

library(BGDataExt)
```

#### A function for RR based on direct inversion

I use this to check the solution of the functions that work using iterative procedures.

```r
## Ridge-regression with 'prior mean', i am using this to check the internal function in BGDataExt
 RR2=function(XX,Xy,lambda,b0=rep(0,ncol(XX)),lambda0=0){
	diag(XX)=diag(XX)+lambda
	rhs=Xy+lambda0*b0
	sol=solve(XX,rhs)
	return(sol)
 }

```

#### Data 

```r
 library(BGLR)
 data(wheat)
 X=scale(wheat.X,center=TRUE,scale=FALSE)
 y=wheat.Y[,1]

 XX=crossprod(X)
 Xy=crossprod(X,y)

 lambda=sum(diag(XX))/nrow(X)

```

#### Checking the RR (ridge regression) function

```r
 system.time(bHat<-RR(XX,Xy,lambda,tol=1e-8)) # using 1e-5 (the default) renders an algorithm orders of magnitude faster, with good precision.
 system.time(bHat2<-RR2(XX,Xy,lambda))
 plot(bHat2,bHat,cex=.5,col=4);abline(a=0,b=1,col=2)
```


####  Checking the Gradient Descent Function (the defailt rearning rate is 0.2

```r
 PATH=matrix(nrow=ncol(X),ncol=9)
 PATH[,1]=GD(XX,Xy,lambda=lambda,nIter=1) # starting values are zeros
 for(i in 2:9){
 	PATH[,i]<-GD(XX,Xy,lambda=lambda,nIter=2,b=PATH[,i-1])
 } 
 par(mfrow=c(3,3))
 for(i in 1:9){ 
     plot(bHat2,PATH[,i],cex=.5,col=4,ylim=range(bHat2),
             xlim=range(bHat2),main=paste0(i*2, ' iterations.'));
    abline(a=0,b=1,col=2)
 }
```
 
#### Increasing the learning rate

```r
 PATH=matrix(nrow=ncol(X),ncol=9)
 PATH[,1]=GD(XX,Xy,lambda=lambda,nIter=1,learning_rate=.5) # starting values are zeros
 for(i in 2:9){
 	PATH[,i]<-GD(XX,Xy,lambda=lambda,nIter=2,b=PATH[,i-1],learning_rate=.5)
 } 

 par(mfrow=c(3,3))
 for(i in 1:9){ 
     plot(bHat2,PATH[,i],cex=.5,col=4,ylim=range(bHat2),
             xlim=range(bHat2),main=paste0(i*2, ' iterations.'));
    abline(a=0,b=1,col=2)
 }
 
```
