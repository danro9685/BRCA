# perform the assignment of K somatic mutational signatures provided as input to samples given a set of observed counts x
"signaturesAssignment" <- function( x, beta ) {
    
    # initialize alpha with an empty matrix
    alpha <- array(NA,c(nrow(x),nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    
    # perform signatures assignments
    for(j in 1:nrow(alpha)) {
        res <- cv.glmnet(t(beta),as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=100,family="gaussian",lower.limits=0.00)
        alpha[j,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
    }

    # return the assigned signatures
    results <- list(alpha=alpha,beta=beta)
    return(results)
    
}
