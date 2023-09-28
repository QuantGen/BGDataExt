# A wrappper for ridge regression, possibly shrinking towards b0
RR <- function(XX, Xy,lambda,p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda2=0, active=1:p, RSS=1, maxIter=1000, tol=1e-5) {
    active <- active - 1L # for the 0-based index
    diag(XX)=diag(XX)+lambda
    Xy=Xy+lambda*b0
    ans <- .Call("fitLSYS", XX, Xy, b, active, RSS, maxIter, tol)
    return(list(b = ans[[1]], RSS = ans[[2]]))
}
