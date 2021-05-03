FWD <- function(y, X, df = 20, tol = 1e-7, maxIter = 1000, centerImpute = TRUE, verbose = TRUE) {
    if (is.null(colnames(X))) {
        colnames(X) <- paste0("X", 1:ncol(X))
    }
    colNames <- colnames(X)
    y <- y - mean(y)
    if (centerImpute) {
        X <- BGData::preprocess(X, center = TRUE, impute = TRUE)
    }
    X <- cbind(1, X)
    df <- df + 1
    colnames(X) <- c("Int", colNames)
    colNames <- colnames(X)
    C <- crossprod(X)
    rhs <- crossprod(X, y)
    n <- length(y)
    p <- ncol(X)
    B <- matrix(data = 0, nrow = p, ncol = df)
    rownames(B) <- colNames
    B[1, 1] <- mean(y)
    RSS <- rep(NA_real_, df)
    LogLik <- rep(NA_real_, df)
    VARE <- rep(NA_real_, df)
    AIC <- rep(NA_real_, df)
    DF <- rep(NA_real_, df)
    BIC <- rep(NA_real_, df)
    path <- rep(NA_character_, df)
    RSS[1] <- sum((y - B[1, 1])^2)
    tol <- tol * RSS[1]
    DF[1] <- 1
    VARE[1] <- RSS[1] / (n - DF[1])
    LogLik[1] <- -(n / 2) * log(2 * pi * VARE[1]) - RSS[1] / (2 * VARE[1])
    AIC[1] <- -2 * LogLik[1] + 2 * DF[1]
    BIC[1] <- -2 * LogLik[1] + log(n) * (DF[1] + 1)
    path[1] <- colNames[1]
    for (i in 2:df) {
        tmp <- addOne(C = C, rhs = rhs, b = B[, i - 1], RSS = RSS[i - 1], maxIter = maxIter, tol = tol)
        B[, i] <- tmp[["b"]]
        if (length(tmp[["newPred"]]) > 0) {
            path[i] <- colNames[tmp[["newPred"]]]
        } else {
            path[i] <- NA
        }
        RSS[i] <- tmp[["RSS"]]
        DF[i] <- sum(tmp[["b"]] != 0)
        VARE[i] <- RSS[i] / (n - DF[i])
        LogLik[i] <- -(n / 2) * log(2 * pi * VARE[i]) - RSS[i] / VARE[i] / 2
        AIC[i] <- -2 * LogLik[i] + 2 * (DF[i] + 1)
        BIC[i] <- -2 * LogLik[i] + log(n) * (DF[i] + 1)
        if (verbose) {
            message("  ", DF[i] - 1, " predictors, AIC=", round(AIC[i], 2))
        }
    }
    OUT <- list(
        B = B,
        path = data.frame(
            variable = path,
            RSS = RSS,
            LogLik = LogLik,
            VARE = VARE,
            DF = DF,
            AIC = AIC,
            BIC = BIC
        )
    )
    return(OUT)
}

addOne <- function(C, rhs, b, RSS, maxIter = 100, tol = 1e-5) {
    inactive <- which(b == 0)
    active <- which(b != 0)
    nInactive <- length(inactive)
    nActive <- length(active)
    # if model is null
    if (nActive == 0) {
        bOLS <- rhs / diag(C)
        dRSS <- diag(C) * bOLS^2
        k <- which.max(dRSS)
        b[k] <- bOLS[k]
        RSS <- RSS - bOLS^2 * C[k, k]
        ans <- list(b = b, newPred = inactive[k], RSS = RSS)
    # when model is not null
    } else {
        RSSNew <- rep(NA_real_, nInactive)
        for (i in 1:nInactive) {
            fm <- fitSYS(C = C, rhs = rhs, b = b, active = c(inactive[i], active), RSS = RSS, maxIter = maxIter, tol = tol)
            RSSNew[i] <- fm[["RSS"]]
        }
        k <- which.min(RSSNew)
        fm <- fitSYS(C = C, rhs = rhs, b = b, active = c(inactive[k], active), RSS = RSS, maxIter = maxIter, tol = tol)
        ans <- list(b = fm[["b"]], newPred = inactive[k], RSS = fm[["RSS"]])
    }
    return(ans)
}

fitSYS <- function(C, rhs, b, active, RSS, maxIter, tol) {
    active <- active - 1L # for the 0-based index
    ans <- .Call("fitLSYS", C, rhs, b, active, RSS, maxIter, tol)
    return(list(b = ans[[1]], RSS = ans[[2]]))
}
