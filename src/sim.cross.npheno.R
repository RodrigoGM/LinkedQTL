## R/qtl 1.36-6
sim.cross.npheno <- 
    function (map, model = NULL, n.ind = 100, type = c("f2", "f2p5", "bc", 
                 "4way", "risib", "riself", "ri4sib", "ri4self", "ri8sib", 
                 "ri8self", "bcsft"), error.prob = 0, missing.prob = 0, partial.missing.prob = 0, 
              keep.qtlgeno = TRUE, keep.errorind = TRUE, m = 0, p = 0, 
              map.function = c("haldane", "kosambi", "c-f", "morgan"), 
              founderGeno, random.cross = TRUE, ...) 
{
    type <- match.arg(type)
    map.function <- match.arg(map.function)
    if (error.prob < 1e-50) 
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    if (type == "risib" || type == "riself") {
        if (type == "risib") 
            type <- "sibmating"
        else type <- "selfing"
        cross <- sim.ril(map, n.ind, type, "2", m = m, p = p, 
                         error.prob = error.prob, missing.prob = missing.prob)
        cross$cross <- NULL
        return(cross)
    }
    if (type == "ri4sib" || type == "ri4self" || type == "ri8sib" || 
        type == "ri8self") {
        if (substr(type, 4, nchar(type)) == "self") 
            crosstype <- "selfing"
        else crosstype <- "sibmating"
        n.str <- substr(type, 3, 3)
        cross <- sim.ril(map, n.ind, crosstype, n.str, m = m, 
                         p = p, random.cross = random.cross, error.prob = 0, 
                         missing.prob = missing.prob)
        rcross <- convertMWril(cross, founderGeno, error.prob = error.prob)
        for (i in names(cross$geno)) if (!("truegeno" %in% names(rcross$geno[[i]]))) 
                                         rcross$geno[[i]]$truegeno <- cross$geno[[i]]$data
        class(rcross)[1] <- substr(class(cross)[1], 1, nchar(class(cross)[1]) - 
                                       2)
        fg <- t(founderGeno[[1]])
        if (length(founderGeno) > 1) 
            for (i in 2:length(founderGeno)) fg <- cbind(fg, 
                                                         t(founderGeno[[i]]))
        colnames(fg) <- markernames(rcross)
        rcross$founderGeno <- fg
        return(rcross)
    }
    if (!is.null(model) && is.matrix(model)) 
        model <- model[order(model[, 1], model[, 2]), ]
    if (type == "bc") 
        cross <- sim.cross.bc(map, model, n.ind, error.prob, 
                              missing.prob, keep.errorind, m, p, map.function)
    else if (type == "f2p5") 
        cross <- sim.cross.f2.npheno(map, model, n.ind, error.prob, 
                                     missing.prob, partial.missing.prob, keep.errorind, 
                                     m, p, map.function)
    else if (type == "f2") 
        cross <- sim.cross.f2(map, model, n.ind, error.prob, 
                              missing.prob, partial.missing.prob, keep.errorind, 
                              m, p, map.function)
    else if (type == "bcsft") 
        cross <- sim.cross.bcsft(map, model, n.ind, error.prob, 
                                 missing.prob, partial.missing.prob, keep.errorind, 
                                 m, p, map.function, ...)
    else cross <- sim.cross.4way(map, model, n.ind, error.prob, 
                                 missing.prob, partial.missing.prob, keep.errorind, m, 
                                 p, map.function)
    qtlgeno <- NULL
    for (i in 1:nchr(cross)) {
        o <- grep("^QTL[0-9]+", colnames(cross$geno[[i]]$data))
        if (length(o) != 0) {
            qtlgeno <- cbind(qtlgeno, cross$geno[[i]]$data[, 
                                                           o, drop = FALSE])
            cross$geno[[i]]$data <- cross$geno[[i]]$data[, -o, 
                                                         drop = FALSE]
            if (is.matrix(cross$geno[[i]]$map)) 
                cross$geno[[i]]$map <- cross$geno[[i]]$map[, 
                                                           -o, drop = FALSE]
            else cross$geno[[i]]$map <- cross$geno[[i]]$map[-o]
        }
    }
    if (keep.qtlgeno) 
        cross$qtlgeno <- qtlgeno
    for (i in 1:nchr(cross)) storage.mode(cross$geno[[i]]$data) <- "integer"
    if (is.null(names(cross$geno))) 
        names(cross$geno) <- 1:length(cross$geno)
    cross
}

