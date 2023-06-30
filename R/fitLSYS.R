
fitLSYS=function(C,rhs,p=nrow(C),b=rep(0,p),active=1:p,RSS=1e5,maxIter=10,tol=1e-4){          
  
    fm <- fitSYS(
                C = C,
                rhs = rhs,
                b = b,
                active = c(inactiveSet[i], activeSet),
                RSS = RSS,
                maxIter = maxIter,
                tol = tol
            )
      return(fm)
}
