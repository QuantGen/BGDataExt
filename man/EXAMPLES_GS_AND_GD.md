## Srhinkage estimation through penalization and early stopping

### 1) Standard shrinkage twoards zero

```r
 library(BGDataExt)

 # A function for RR using direct inversion (I use this to check the results of the C code impleneted in BGDataExt

 RR2=function(XX,Xy,lambda,b0=rep(0,ncol(XX)),lambda0=0){
	diag(XX)=diag(XX)+lambda
	rhs=Xy+lambda*lambda0*b0
	sol=solve(XX,rhs)
	return(sol)
 }

```

**Data** 

```r
 library(BGLR)
 data(wheat)
 X=scale(wheat.X,center=TRUE,scale=FALSE)
 y=wheat.Y[,1]

 XX=crossprod(X)
 Xy=crossprod(X,y)

 lambda=sum(diag(XX))/nrow(X)

```

**Checking the RR (ridge regression) function**

```r
 system.time(bHat<-RR(XX,Xy,lambda,tol=1e-8)) # using 1e-5 (the default) renders an algorithm orders of magnitude faster, with good precision.
 system.time(bHat2<-RR2(XX,Xy,lambda))
 plot(bHat2,bHat,cex=.5,col=4);abline(a=0,b=1,col=2)
```


**Checking the Gradient Descent Function (the default rearning rate is 0.2)**

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
 
**Increasing the learning rate**

The recommendation is to use learning rates between 0.05 to 0.8 (not much higher)

```r

learning_rate=.5
 PATH=matrix(nrow=ncol(X),ncol=9)
 PATH[,1]=GD(XX,Xy,lambda=lambda,nIter=1,learning_rate=learning_rate) # starting values are zeros
 for(i in 2:9){
 	PATH[,i]<-GD(XX,Xy,lambda=lambda,nIter=2,b=PATH[,i-1],learning_rate=learning_rate)
 } 

 par(mfrow=c(3,3))
 for(i in 1:9){ 
     plot(bHat2,PATH[,i],cex=.5,col=4,ylim=range(bHat2),
             xlim=range(bHat2),main=paste0(i*2, ' iterations.'));
    abline(a=0,b=1,col=2)
 }
 
```


### 2) Shrinkage towards a prior estimate


**Data**

Here we split the wheat data set in two sets (based on a simple clustering of the genotyepes), set 1 is the one we intend to transfer learning from, set 2 is the target data set. Within this set, we create a training and a testing data set.

```r
 set.seed(195021)
 # Clustering 
 SVD=svd(X,nu=5,nv=0)
 group=kmeans(SVD$u,centers=2,nstart=100)$cluster
 plot(SVD$u[,1:2],col=group*2)

 # Data split 
 group1=which(group==1)
 group2=which(group==2)

 X1=X[group1,]
 y1=y[group1]
 X2=X[group2,]
 y2=y[group2]

 tst=sample(1:length(y2),size=100)
 X2.TRN=X2[-tst,]
 y2.TRN=y2[-tst]

 X2.TST=X2[tst,]
 y2.TST=y2[tst]
```

**Cross-group prediction**

```r
 XX1=crossprod(X1)
 Xy1=crossprod(X1,y1)
 bHat1=RR(XX1,Xy1,lambda)
 yHat2_cross=X2.TST%*%bHat1

 COR=c('Across'=cor(yHat2_cross,y2.TST))

```

**Within-group prediction**

```r
 XX2=crossprod(X2.TRN)
 Xy2=crossprod(X2.TRN,y2.TRN)
 bHat2=RR(XX2,Xy2,lambda)
 yHat2_within=X2.TST%*%bHat2

 COR['Within']= cor(yHat2_within,y2.TST)
```


**Transfer learning trhough Gradient Descent with Early Stopping**

```r
  

```

**Transfer learning with srhinkage towards a external estiamte**


