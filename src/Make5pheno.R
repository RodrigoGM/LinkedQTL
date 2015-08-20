#Generation of 10 correlated phenotypes using a Cholesky Decomposition Matrix

make.5.pheno <- 
    function(n.ind=10) {
        mm<- matrix(data=c(  1, 0.9, 0.6, 0.3,
                        0.9,   1, 0.2, 0.2,
                        0.6, 0.2,   1, 0.2,
                        0.3, 0.2, 0.2,   1),
                    nrow=4, ncol=4)
        q <- chol(mm)
        
        pheno <- matrix(NA, nrow=n.ind, ncol=4)
        sapply(1:4, function(x) pheno[,x] <<- rnorm(n.ind, 0, sqrt(100)))
        pheno <- cbind(pheno %*% q, rnorm(n.ind, 0, sqrt(100)))
        pheno
    }
