## Transfer Learnin using Early Stopping, Ridge Regression, and a Bayesian method

**Libraries**

These libraries should be installed from GitHub, not CRAN (some updates are not available in CRAN).


```r
 # To install from GitHub
 # library(remotes)
 # install_github('https://github.com/QuantGen/BGDataExt')
 # install_github('https://github.com/gdlc/BGLR-R')
 library(BGDataExt)
 library(BGLR)
```
 
**Data** 

```r
 data(wheat)
 X=scale(wheat.X,center=TRUE,scale=FALSE)
 y=wheat.Y[,1]
```

**Data split**

Here we split the wheat data set in two sets:
  - Set 1 is the one we intend to transfer learning from,
  - Set 2 is the target data set; within this set, we create a training and a testing data set.

We consider **two options**:
  - Cluster the genotypes and use that to split the data; this option lead to a situation with (potential) effect-heterogeneity,
  - Assign genotypes to groups at random, in this case, on average, there is no effect heterogenetiy.

Run only one of the two options:

**Option 1**: Clustering

```r
 set.seed(195021)
 N=nrow(X)
 SVD=svd(X,nu=5,nv=0)
 tmp=kmeans(SVD$u,centers=2,nstart=50)
 group=tmp$cluster
 plot(SVD$u[,1:2],col=group*2)
```

**Option 2:** at random

```r
 set.seed(195021)
 N=nrow(X)
 group=sample(1:2,size=N,replace=TRUE)
 group1=which(group==1)
 group2=which(group==2)
```

 #D1: external data set
 X1=scale(X[group1,],center=TRUE,scale=FALSE)
 y1=scale(y[group1],center=TRUE,scale=FALSE)
 X2=scale(X[group2,],center=TRUE,scale=FALSE)
 y2=y[group2]

 tst=sample(1:length(y2),size=50)
 # D2: training
 X2.TRN=X2[-tst,]
 y2.TRN=scale(y2[-tst],center=TRUE,scale=FALSE)

 # D2: testing
 X2.TST=X2[tst,]
 y2.TST=y2[tst]

```

### Benchmarks

As benchmark we consider:

 - **Cross-group prediction**: fit the model in D1, use the fitted model to predict the testing set of D2.
 - **Within-group prediction**: fit the model to the training set of D2, use themodel to predict the testing set of D2.
 - **Joint Analysis**: fit the model using D1 and the training set of D2. This would be optimal in absence of effect heterogneiety between data set. However, it is not alwasy feasible.

**Cross-group prediction**

```r
 XX1=crossprod(X1)
 Xy1=crossprod(X1,y1)
 lambda=sum(diag(XX1))/nrow(X1)
 bHat1=RR(XX1,Xy1,lambda,tol=1e-5) # ridge regression
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

**Joint Analysis**

```r
 XX=XX1+XX2
 Xy=Xy1+Xy2
 lambda=sum(diag(XX))/(nrow(X1)+nrow(X2.TRN))
 bHatCombined=RR(XX,Xy,lambda,tol=1e-5) # ridge regression
 yHat2_combined=X2.TST%*%bHatCombined
 COR['Combined']=cor(yHat2_combined,y2.TST)

```

**Transfer learning trhough Gradient Descent with Early Stopping**

Here we do Gradient Descent for an OLS ojbective function, using the estiamtes from D1 (`bHat1`) as starting value.

Then, we evaluate prediction accuracy in the testing set of D2 over the solutions obtained through GD.

```r
  # B contains the solutions, B[,1]=bHat1, obtained in the GD path
  B=GD(XX2,Xy2,b=bHat1,lambda=0,nIter=100,returnPath=TRUE,learning_rate=.03)

  COR.GD=rep(NA,ncol(B))
  for(i in 1:ncol(B)){
    COR.GD[i]=cor(y2.TST,X2.TST%*%B[,i])
  }
  plot(COR.GD,type='o');abline(h=COR,col=c(2,4))
  COR['Max GD']=max(COR.GD)
  COR

```

**Transfer learning with srhinkage towards a external estiamte  using Ridge Regression**


The following code fits Ridge Regression over a grid of values of the regularization parameter, shrinking the estimates towards the external estimate (`bHat1`).
```r
 lambda0=seq(from=0,to=1,by=.1)
 B=matrix(nrow=nrow(XX2),ncol=length(lambda0),NA)
 COR.RR=rep(NA,ncol(B))
 for(i in 1:ncol(B)){
    B[,i]=RR(XX2,Xy2, lambda=lambda,lambda0=lambda0[i],b0=bHat1)
    COR.RR[i]=cor(y2.TST,X2.TST%*%B[,i])
 }
 plot(COR.RR,type='o');abline(h=COR,col=c(2,4))
 COR['Max RR']=max(COR.RR)
 COR

```

**Bayesian method**

The following approach uses a Bayesian mixture model with two components, one centered in 0 and the other one centered in bHat1

```r
source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/misc/mixturesWithNonZeroPriorMeans.R')
my=mean(y2.TRN)
vy=var(y2.TRN)
n=nrow(X2.TRN)

B0=as.matrix(bHat1)
bayes=BMM(C=XX2,rhs=Xy2,my=my,vy=vy,n=n,nIter=1500,burnIn=500,B0=B0,verbose=FALSE)
COR['Bayes_1_comp.']= cor(y2.TST,X2.TST%*%bayes$b)

B0=cbind(0,bHat1)
bayes=BMM(C=XX2,rhs=Xy2,my=my,vy=vy,n=n,nIter=1500,burnIn=500,B0=B0,verbose=FALSE)
COR['Bayes_2_comp']= cor(y2.TST,X2.TST%*%bayes$b)

```

**Plot with results**

```r
p=ggplot(data=data.frame(model=names(COR),RSq=COR^2),aes(y=RSq,x=model,fill=model))+geom_bar(stat='identity')
plot(p)
```
# End
## Old code I used to test some of the functions

**A function for RR using direct inversion**

I use this to check the results of the C code impleneted in BGDataExt.

```r
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




