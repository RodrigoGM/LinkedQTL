sim.cross.f2.npheno=
function (map, model, n.ind, error.prob, missing.prob, partial.missing.prob,
    keep.errorind, m, p, map.function, five.phenos=TRUE, ...)
{
    if (map.function == "kosambi")
        mf <- mf.k
    else if (map.function == "c-f")
        mf <- mf.cf
    else if (map.function == "morgan")
        mf <- mf.m
    else mf <- mf.h
    if (any(sapply(map, is.matrix)))
        stop("Map must not be sex-specific.")
    chr.type <- sapply(map, function(a) if (is.null(class(a)))
        return("A")
    else return(class(a)))
    n.chr <- length(map)
    if (is.null(model))
        n.qtl <- 0
    else {
        if (!((!is.matrix(model) && length(model) == 4) || (is.matrix(model) &&
            ncol(model) == 4))) {
            stop("Model must be a matrix with 4 columns (chr, pos and effects).")
        }
        if (!is.matrix(model))
            model <- rbind(model)
        n.qtl <- nrow(model)
        if (any(model[, 1] < 0 | model[, 1] > n.chr))
            stop("Chromosome indicators in model matrix out of range.")
        model[, 2] <- model[, 2] + 1e-14
    }
    if (n.qtl > 0) {
        for (i in 1:n.qtl) {
            temp <- map[[model[i, 1]]]
            if (model[i, 2] < min(temp)) {
                temp <- c(model[i, 2], temp)
                names(temp)[1] <- paste("QTL", i, sep = "")
            }
            else if (model[i, 2] > max(temp)) {
                temp <- c(temp, model[i, 2])
                names(temp)[length(temp)] <- paste("QTL", i,
                  sep = "")
            }
            else {
                j <- max((seq(along = temp))[temp < model[i,
                  2]])
                temp <- c(temp[1:j], model[i, 2], temp[(j + 1):length(temp)])
                names(temp)[j + 1] <- paste("QTL", i, sep = "")
            }
            map[[model[i, 1]]] <- temp
        }
    }
    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    n.mar <- sapply(map, length)
    mar.names <- lapply(map, names)
    for (i in 1:n.chr) {
        thedata <- sim.bcg(n.ind, map[[i]], m, p, map.function)
        dimnames(thedata) <- list(NULL, mar.names[[i]])
        if (chr.type[i] != "X")
            thedata <- thedata + sim.bcg(n.ind, map[[i]], m,
                p, map.function) - 1
        geno[[i]] <- list(data = thedata, map = map[[i]])
        class(geno[[i]]) <- chr.type[i]
        class(geno[[i]]$map) <- NULL
    }

    if(five.phenos == FALSE) pheno <- rnorm(n.ind, mean = 0, sd = sqrt(100))
    if(five.phenos == TRUE) pheno <- make.5.pheno(n.ind)
    
    if (n.qtl > 0) {
        QTL.chr <- QTL.loc <- NULL
        for (i in 1:n.chr) {
            o <- grep("^QTL[0-9]+", mar.names[[i]])
            if (length(o) > 0) {
                QTL.chr <- c(QTL.chr, rep(i, length(o)))
                QTL.loc <- c(QTL.loc, o)
            }
        }
        for (i in 1:n.qtl) {
            QTL.geno <- geno[[QTL.chr[i]]]$data[, QTL.loc[i]]
            pheno[QTL.geno == 1] <- pheno[QTL.geno == 1] - model[i,
                3]
            pheno[QTL.geno == 2] <- pheno[QTL.geno == 2] + model[i,
                4]
            pheno[QTL.geno == 3] <- pheno[QTL.geno == 3] + model[i,
                3]
        }
    }
    n.mar <- sapply(geno, function(a) length(a$map))
    if (error.prob > 0) {
        for (i in 1:n.chr) {
            if (chr.type[i] == "X") {
                a <- sample(0:1, n.mar[i] * n.ind, repl = TRUE,
                  prob = c(1 - error.prob, error.prob))
                geno[[i]]$data[a == 1] <- 3 - geno[[i]]$data[a ==
                  1]
            }
            else {
                a <- sample(0:2, n.mar[i] * n.ind, repl = TRUE,
                  prob = c(1 - error.prob, error.prob/2, error.prob/2))
                if (any(a > 0 & geno[[i]]$data == 1))
                  geno[[i]]$data[a > 0 & geno[[i]]$data == 1] <- (geno[[i]]$data +
                    a)[a > 0 & geno[[i]]$data == 1]
                if (any(a > 0 & geno[[i]]$data == 2)) {
                  geno[[i]]$data[a > 0 & geno[[i]]$data == 2] <- (geno[[i]]$data +
                    a)[a > 0 & geno[[i]]$data == 2]
                  geno[[i]]$data[geno[[i]]$data > 3] <- 1
                }
                if (any(a > 0 & geno[[i]]$data == 3))
                  geno[[i]]$data[a > 0 & geno[[i]]$data == 3] <- (geno[[i]]$data -
                    a)[a > 0 & geno[[i]]$data == 3]
            }
            if (keep.errorind) {
                errors <- matrix(0, n.ind, n.mar[i])
                errors[a > 0] <- 1
                colnames(errors) <- colnames(geno[[i]]$data)
                geno[[i]]$errors <- errors
            }
        }
    }
    if (partial.missing.prob > 0) {
        for (i in 1:n.chr) {
            if (chr.type[i] != "X") {
                o <- sample(c(TRUE, FALSE), n.mar[i], repl = TRUE,
                  prob = c(partial.missing.prob, 1 - partial.missing.prob))
                if (any(o)) {
                  o2 <- grep("^QTL[0-9]+", mar.names[[i]])
                  if (length(o2) > 0)
                    x <- geno[[i]]$data[, o2]
                  m <- (1:n.mar[i])[o]
                  for (j in m) {
                    if (runif(1) < 0.5)
                      geno[[i]]$data[geno[[i]]$data[, j] == 1 |
                        geno[[i]]$data[, j] == 2, j] <- 4
                    else geno[[i]]$data[geno[[i]]$data[, j] ==
                      3 | geno[[i]]$data[, j] == 2, j] <- 5
                  }
                  if (length(o2) > 0)
                    geno[[i]]$data[, o2] <- x
                }
            }
        }
    }
    if (missing.prob > 0) {
        for (i in 1:n.chr) {
            o <- grep("^QTL[0-9]+", mar.names[[i]])
            if (length(o) > 0)
                x <- geno[[i]]$data[, o]
            geno[[i]]$data[sample(c(TRUE, FALSE), n.mar[i] *
                n.ind, repl = TRUE, prob = c(missing.prob, 1 -
                missing.prob))] <- NA
            if (length(o) > 0)
                geno[[i]]$data[, o] <- x
        }
    }

    pheno <- data.frame(phenotype = pheno)
    cross <- list(geno = geno, pheno = pheno)
    class(cross) <- c("f2", "cross")
    cross
}
