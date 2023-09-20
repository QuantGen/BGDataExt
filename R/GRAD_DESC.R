# A wrapper to the C-function that performs Gradient Descent in a system of linear equations

GRAD_DESC<- function(C, rhs, b, active, RSS, nIter, learning_rate) {
    active <- active - 1L # for the 0-based index
    ans <- .Call("fitLSYS", C, rhs, b, active, RSS, maxIter, tol)
    return(list(b = ans[[1]], RSS = ans[[2]]))
}
