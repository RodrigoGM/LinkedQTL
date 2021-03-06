plotmedians <-
function (formula, data = NULL, subset, na.action, bars = TRUE,
    p = c(.25, 0.75), minsd = 0, minbar = NULL, maxbar = NULL, xlab = names(mf)[2],
    ylab = names(mf)[1], mean.labels = FALSE, ci.label = FALSE,
    n.label = TRUE, digits = getOption("digits"), col = "black",
    barwidth = 1, barcol = "blue", connect = TRUE, ccol = col,
    legends = names(medians), xaxt, use.t = TRUE, ...)
{
    is.R <- get("is.R")
    if (is.null(is.R))
        is.R <- function(x) FALSE
    if (!is.R()) {
        if (col == "black")
            col <- 1
        if (barcol == "blue")
            barcol <- 2
    }
    if (invalid(formula) || (length(formula) != 3))
        stop("formula missing or incorrect")
    if (invalid(na.action))
        na.action <- options("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.R()) {
        if (is.matrix(eval(m$data, parent.frame())))
            m$data <- as.data.frame(data)
    }
    else {
        if (is.matrix(eval(m$data, FALSE)))
            m$data <- as.data.frame(data)
    }
    m$... <- m$bars <- m$barcol <- m$p <- NULL
    m$minsd <- m$minbar <- m$maxbar <- NULL
    m$xlab <- m$ylab <- NULL
    m$col <- m$barwidth <- NULL
    m$digits <- m$mean.labels <- m$ci.label <- m$n.label <- NULL
    m$connect <- m$ccol <- m$legends <- m$labels <- NULL
    m$xaxt <- m$use.t <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    wFact <- which(attr(attr(mf, "terms"), "dataClasses") ==
        "factor")
    for (i in wFact) mf[, i] <- factor(mf[, i])
    medians <- sapply(split(mf[[response]], mf[[-response]]), median,
        na.rm = TRUE)
    ns <- sapply(sapply(split(mf[[response]], mf[[-response]]),
        na.omit, simplify = FALSE), length)
    xlim <- c(0.5, length(medians) + 0.5)
    if (!bars) {
        plot(medians, ..., col = col, xlim = xlim)
    }
    else {
        myq25 <- function(x) quantile(x[!is.na(x)], probs = p[1])
        myq75 <- function(x) quantile(x[!is.na(x)], probs = p[2])
        q25s <- sapply(split(mf[[response]], mf[[-response]]),
            myq25)
        q75s <- sapply(split(mf[[response]], mf[[-response]]),
            myq75)
#                vars <- ifelse(vars < (minsd^2), (minsd^2), vars)
#        if (use.t)
#            ci.width <- qt((1 + p)/2, ns - 1) * sqrt(vars/ns)
#        else ci.width <- qnorm((1 + p)/2) * sqrt(vars/ns)
        if (length(mean.labels) == 1) {
            if (mean.labels == TRUE)
                mean.labels <- format(round(medians, digits = digits))
            else if (mean.labels == FALSE)
                mean.labels <- NULL
        }

        plotCI(x = 1:length(medians), y = medians, uiw = q75s, liw = q25s,
            xaxt = "n", xlab = xlab, ylab = ylab, labels = mean.labels,
            col = col, xlim = xlim, lwd = barwidth, barcol = barcol,
            minbar = minbar, maxbar = maxbar, ...)

        if (invalid(xaxt) || xaxt != "n")
            axis(1, at = 1:length(medians), labels = legends)

        if (ci.label) {
            ci.lower <- medians - q25s
            ci.upper <- medians + q75s
            if (!invalid(minbar))
                ci.lower <- ifelse(ci.lower < minbar, minbar,
                  ci.lower)
            if (!invalid(maxbar))
                ci.upper <- ifelse(ci.upper > maxbar, maxbar,
                  ci.upper)
            labels.lower <- paste(" \n", format(round(ci.lower,
                digits = digits)), sep = "")
            labels.upper <- paste(format(round(ci.upper, digits = digits)),
                "\n ", sep = "")
            text(x = 1:length(medians), y = ci.lower, labels = labels.lower,
                col = col)
            text(x = 1:length(medians), y = ci.upper, labels = labels.upper,
                col = col)
        }
    }
    if (n.label)
        if (is.R())
            text(x = 1:length(medians), y = par("usr")[3], labels = paste("n=",
                ns, "\n", sep = ""))
        else {
            axisadj <- (par("usr")[4] - (par("usr")[3]))/75
            text(x = 1:length(medians), y = par("usr")[3] + axisadj,
                labels = paste("n=", ns, "\n", sep = ""))
        }
    if (!invalid(connect) & !identical(connect, FALSE)) {
        if (is.list(connect)) {
            if (length(ccol) == 1)
                ccol <- rep(ccol, length(connect))
            for (which in 1:length(connect)) lines(x = connect[[which]],
                y = medians[connect[[which]]], col = ccol[which])
        }
        else lines(medians, ..., col = ccol)
    }
}

