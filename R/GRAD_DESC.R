# A wrapper to the C-function that performs Gradient Descent in a system of linear equations

GD<- function(XX,Xy,b=rep(0,ncol(XX)),active=1:ncol(C), RSS=1,nIter=10,learning_rate=mean(diag(XX))/10,lambda=0,b0=rep(0,ncol(XX)),lambda0=1){
    #C, rhs, b, active=1:ncol(C), RSS=1, nIter, learning_rate) {
    diag(XX)=diag(XX)+lambda
    Xy=Xy+lambda*lambda0*b0
    active <- active - 1L # for the 0-based index
    ans <- .Call("GRAD_DESC", C, rhs, b, active, nIter, learning_rate)
    return(ans[[1]])
}
