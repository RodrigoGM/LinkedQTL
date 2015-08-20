sim.cross.npheno=
function (map, model = NULL, n.ind = 100, type = c("f2", "bc",
    "4way"), error.prob = 0, missing.prob = 0, partial.missing.prob = 0,
    keep.qtlgeno = TRUE, keep.errorind = TRUE, m = 0, p = 0,
    map.function = c("haldane", "kosambi", "c-f", "morgan"), ...)
{
    type <- match.arg(type)
    map.function <- match.arg(map.function)
    if (error.prob < 1e-50)
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    if (!is.null(model) && is.matrix(model))
        model <- model[order(model[, 1], model[, 2]), ]
    if (type == "bc")
        cross <- sim.cross.bc(map, model, n.ind, error.prob,
            missing.prob, keep.errorind, m, p, map.function)
    else if (type == "f2")
        cross <- sim.cross.f2.npheno(map, model, n.ind, error.prob,
            missing.prob, partial.missing.prob, keep.errorind,
            m, p, map.function)
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
    cross
}
