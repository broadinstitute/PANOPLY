### Based on heatmap2 3.0.1 from the gplots package [https://cran.r-project.org/package=gplots]
### Fix: Avoid calculating dendrogram when Rowv or Colv are set to NA, FALSE, or NULL.
### License: GPL-2

### Package: gplots
### Title: Various R Programming Tools for Plotting Data
### Description: Various R programming tools for plotting data, including:
###   - calculating and plotting locally smoothed summary function as
###     ('bandplot', 'wapply'),
###   - enhanced versions of standard plots ('barplot2', 'boxplot2',
###     'heatmap.2', 'smartlegend'),
###   - manipulating colors ('col2hex', 'colorpanel', 'redgreen',
###     'greenred', 'bluered', 'redblue', 'rich.colors'),
###   - calculating and plotting two-dimensional data summaries ('ci2d',
###     'hist2d'),
###   - enhanced regression diagnostic plots ('lmplot2', 'residplot'),
###   - formula-enabled interface to 'stats::lowess' function ('lowess'),
###   - displaying textual data in plots ('textplot', 'sinkplot'),
###   - plotting a matrix where each cell contains a dot whose size
###     reflects the relative magnitude of the elements ('balloonplot'),
###   - plotting "Venn" diagrams ('venn'),
###   - displaying Open-Office style plots ('ooplot'),
###   - plotting multiple data on same region, with separate axes
###     ('overplot'),
###   - plotting means and confidence intervals ('plotCI', 'plotmeans'),
###   - spacing points in an x-y plot so they don't overlap ('space').
### Depends: R (>= 3.0)
### Imports: gtools, gdata, stats, caTools, KernSmooth
### Suggests: grid, MASS
### Version: 3.0.1
### Date: 2016-03-30
### Author: Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert
### 	Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin
### 	Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz, Bill
### 	Venables
### Maintainer: Gregory R. Warnes <greg@warnes.net>
### License: GPL-2
### NeedsCompilation: no
### Packaged: 2016-03-30 17:32:42 UTC; gwarnes1
### Repository: CRAN
### Date/Publication: 2016-03-30 21:04:18


invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0) 
        return(TRUE)
    if (is.list(x)) 
        return(all(sapply(x, invalid)))
    else if (is.vector(x)) 
        return(all(is.na(x)))
    else return(FALSE)
}

heatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
        w), symm = FALSE, scale = c("none", "row", "column"), 
    na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
    symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
    col = "heat.colors", colsep, rowsep, sepcolor = "white", 
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
    na.color = par("bg"), trace = c("column", "row", "both", 
        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
    linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
        NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
    colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5, 
    density.info = c("histogram", "density", "none"), denscol = tracecol, 
    symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
    key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, 
    key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, 
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && any(duplicated(breaks))) 
        stop("breaks may not contian duplicate values")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (all(Colv == "Rowv")) 
        Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && 
            (dendrogram %in% c("both", "row"))) {
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
            if (dendrogram == "both") 
                dendrogram <- "column"
            else dendrogram <- "none"
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && 
            (dendrogram %in% c("both", "column"))) {
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
            if (dendrogram == "both") 
                dendrogram <- "row"
            else dendrogram <- "none"
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
            nr)) 
            stop("Rowv dendrogram doesn't match size of x")
        if (length(rowInd) < nr) 
            nr <- length(rowInd)
    }
    else if (is.integer(Rowv)) {
        distr <- distfun(x)
        hcr <- hclustfun(distr)
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        distr <- distfun(x)
        hcr <- hclustfun(distr)
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    #else if (!isTRUE(Rowv)) {
    #    rowInd <- nr:1
    #    ddr <- as.dendrogram(hclust(dist(diag(nr))))
    #}
    else {
        rowInd <- nr:1
        #ddr <- as.dendrogram(Rowv)
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd > 
            nc)) 
            stop("Colv dendrogram doesn't match size of x")
        if (length(colInd) < nc) 
            nc <- length(colInd)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        distc <- distfun(if (symm) 
            x
        else t(x))
        hcc <- hclustfun(distc)
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        distc <- distfun(if (symm) 
            x
        else t(x))
        hcc <- hclustfun(distc)
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    #else if (!isTRUE(Colv)) {
    #    colInd <- 1:nc
    #    ddc <- as.dendrogram(hclust(dist(diag(nc))))
    #}
    else {
        colInd <- 1:nc
        #ddc <- as.dendrogram(Colv)
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (!is.null(colRow)) 
        colRow <- colRow[rowInd]
    if (!is.null(colCol)) 
        colCol <- colCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    plot.index <- 1
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        plot.index <- plot.index + 1
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        plot.index <- plot.index + 1
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol) && is.null(colCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.null(srtCol) || is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            if (is.null(srtCol)) 
                srtCol <- 90
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol, col = colCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow) && is.null(colRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.null(srtRow) || is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                srt = srtRow, col = colRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in 1:length(colInd)) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in 1:length(rowInd)) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol, 
                  lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    plot.index <- plot.index + 1
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, 
            yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i", 
            leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab)) 
            mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab)) 
            mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title)) 
            mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0) 
            do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min.breaks
            max.raw <- max.breaks
        }
        z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row") 
                key.xlab <- "Row Z-Score"
            else if (scale == "column") 
                key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
                cex = par("cex") * par("cex.lab"))
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE, 
                from = min.scale, to = max.scale)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[!omit]
            dens$y <- dens$y[!omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                  0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Density"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5, cex = par("cex") * par("cex.lab"))
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                  labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Count"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5, cex = par("cex") * par("cex.lab"))
        }
        else if (is.null(key.title)) 
            title("Color Key")
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else {
        par(mar = c(0, 0, 0, 0))
        plot.new()
    }
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    retval$layout <- list(lmat = lmat, lhei = lhei, lwid = lwid)
    if (!is.null(extrafun)) 
        extrafun()
    invisible(retval)
}