### NOTE: Sourced from https://github.com/Displayr/flipPlots/. Needed to run SankeyDiagram(). DO NOT EDIT.

#' Sankey Diagram
#'
#' Creates a sankey diagram of the relationship between different variables.
#' @param data A \code{\link{data.frame}} of variables.
#' @param output.data.only Only outputs the link and node structure without printing the sankey diagram.
#' @param links.and.nodes Input a list structure from a previous call of SankeyDiagram
#'  with \code{output.data.only = TRUE}.
#' @param max.categories When the number of unique values
#' of numeric data exceeds this value, the variable is quantized.
#' @param subset An optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights An optional vector of weights. If not provided, the width of
#' the links is proportional to the number of rows in the data frame that correspond
#' to the path (flow).
#' @param font.size Font size of node labels.
#' @param font.family Font family of node labels.
#' @param font.unit One of "px" of "pt". By default all font sizes are specified in terms of
#'  pixels ("px"). But changing this to "pt" will mean that the font sizes will be in terms
#'  points ("pt"), which will be consistent with font sizes in text boxes.
#' @param colors Colors of the nodes, supplied as a vector of hex colors.
#' @param node.width Width of the width.
#' @param node.padding Vertical space between nodes.
#' @param link.color One of \code{"None", "Source", "Target", "First variable",
#'    "Last variable"}. This specifies whether the links are shown in grey (None);
#'    the same color as the source node; the same color as the target node; or
#'    the same color as node in the first or last variable.
#' @param node.position.automatic Position nodes automatically to reduce
#'    the number of overlapping links. This is turned on by default, and when
#'    turned off the nodes will be positioned in order they are given in the data.
#' @param variables.share.values If \code{TRUE}, and \code{link.color = "Source"}
#'    or \code{"Target"}, then the same set colors will be used for each variable.
#' @param label.max.length Maximum number of characters in each node label.
#' @param label.show.varname Show variable name in node label.
#' @param label.show.counts Append count data to node labels.
#' @param label.show.percentages Append percentages to node labels.
#' @param hovertext.show.percentages Show percentages instead of counts in the hovertext.
#' @param sinks.right Logical indicating whether to move the last nodes to the right border of the plot.
#' @importFrom networkD3 sankeyNetwork JS
#' @importFrom grDevices col2rgb rgb
#' @return A sankey diagram (HTMLwidget).
#' @details Text variables are grouped as having text or not having text.
#' To see patterns with text variables, they should first be turned into
#' factors.
#' @export
SankeyDiagram <- function(data = NULL, links.and.nodes = NULL, output.data.only = FALSE,
                          max.categories = 8, subset = NULL, weights = NULL,
                          font.size = 12, font.family = "Times New Roman",
                          font.unit = "px", colors = NULL,
                          link.color = c("None", "Source", "Target", "First variable",
                          "Last variable")[1], variables.share.values = FALSE,
                          label.show.varname = TRUE, label.max.length = 100,
                          label.show.counts = FALSE, label.show.percentages = FALSE,
                          node.position.automatic = TRUE,
                          node.width = 30, node.padding = 10, sinks.right = TRUE,
                          hovertext.show.percentages = FALSE)
{
    if (!is.null(links.and.nodes))
    {
        link.color <- links.and.nodes$link.color
        links <- links.and.nodes$links
        nodes <- links.and.nodes$nodes
        hovertext.show.percentages <- links.and.nodes$hovertext.show.percentages

    } else
    {
        if (!is.data.frame(data))
            data <- as.data.frame(data)
        else if (!is.null(attr(data[[1]], "questiontype")))
        {
            # From R 4.0, stringsAsFactors defaults to FALSE.
            # Because of this, a specific call to as.data.frame() in the
            # Visualization - Sankey Diagram wiki code no longer behaves as
            # before. This only affects users who supplied variables as data,
            # which is why we check for the questiontype attribute.
            # The code below converts any strings to factors, which would have
            # previously been done by the call to as.data.frame().
            data <- as.data.frame(as.list(data),
                                  stringsAsFactors = TRUE,
                                  check.names = FALSE)
        }

        if (nrow(data) < 2)
            stop(paste0(nrow(data), "observations: more data is required to create a Sankey diagram."))
        if (max.categories < 2)
            stop("Maximum number of categories must be greater than 1.")

        if (!is.null(weights) & length(weights) != nrow(data))
            stop("'weights' and 'data' are required to have the same number of observations. They do not.")
        if (!is.null(subset) & length(subset) > 1 & length(subset) != nrow(data))
            stop("'subset' and 'data' are required to have the same number of observations. They do not.")

        # Take subset and resample to generate a weighted sample, if necessary.
        if (is.null(subset))
            subset.data <- data
        else {
            subset.data <- data[subset, , drop = FALSE]
            weights <- weights[subset]
        }

        tmp.is.numeric <- sapply(subset.data, is.numeric)
        if (link.color %in% c("Source", "Target") && variables.share.values &&
            any(tmp.is.numeric) && !all(tmp.is.numeric))
            stop("'Variables share common values' has been set to true so variables must be of the same type.")
        variables <- categorizeData(subset.data, weights, max.categories,
            variables.share.values, label.show.varname, label.max.length)
        links <- computeLinks(variables, weights, hovertext.show.percentages)
        nodes <- nodeDictionary(variables, weights, label.show.counts, label.show.percentages)

        # Determine color of nodes
        grps <- 0:(nrow(nodes)-1)
        if (variables.share.values && link.color %in% c("Source", "Target"))
        {
            all.levels <- attr(variables, "all.levels")
            tmp.var.names <- paste(names(variables), ": ", sep = "", collapse = "|")
            tmp.node.names <- gsub(tmp.var.names, "", nodes$name)
            if (label.show.counts || label.show.percentages)
                tmp.node.names <- gsub(" \\([a-zA-Z0-9 =%,]*\\)$", "", tmp.node.names, perl = TRUE)
            nodes$group <- factor(match(tmp.node.names, all.levels))
        }
        else if (link.color == "None")      # nodes at each level to be the same
             nodes$group <- factor(rep(1:ncol(variables), sapply(variables, function(x){length(unique(x))})))
        else if (link.color == "Last variable" || link.color == "First variable")
            nodes$group <- factor(getNodeGroups(link.color, links))
        else
            nodes$group <- factor(grps) # all nodes different colors

        # Determine colors of links
        if (link.color == "Source")
            links$group <- as.factor(nodes$group[links$source+1])
        else if (link.color == "Target")
            links$group <- as.factor(nodes$group[links$target+1])
        else if (link.color == "Last variable")
            links$group <- as.factor(nodes$group[links$target+1])
        else if (link.color == "First variable")
            links$group <- as.factor(nodes$group[links$source+1])
    }

    if (output.data.only)
        return(list(link.color = link.color, links = links, nodes = nodes,
                    hovertext.show.percentages = hovertext.show.percentages))

    if (is.null(colors))
        color.str <- "d3.scaleOrdinal(d3.schemeCategory20);"
    else
        color.str <- paste0('d3.scaleOrdinal() .domain([\'',
                     paste(levels(nodes$group), collapse = '\',\''), '\']) .range([\'',
                     paste(colorsToHex(colors), collapse = '\',\''), '\']);')

    # For the other chart types, the font size conversion
    # happens inside flipChart::CChart but SankeyDiagram is called separately.
    if (tolower(font.unit) %in% c("pt", "point", "points"))
    {
        fsc <- 1.3333
        font.size = round(fsc * font.size, 0)
    }

    res <- sankeyNetwork(Links = links, LinkGroup = if (link.color == "None") NULL else 'group',
                Nodes = nodes, NodeID = 'name', NodeGroup = 'group', nodeWidth = node.width,
                Source = "source", Target = "target", Value = "value", nodePadding = node.padding,
                fontSize = font.size, fontFamily = font.family, colourScale = JS(color.str),
                iterations = if (node.position.automatic) 32 else 0,
                units = if (hovertext.show.percentages) "%" else "", sinksRight = sinks.right)
    class(res) <- c(class(res), "visualization-selector")
    return(res)
}

getNodeGroups <- function(type, links)
{
    num.nodes <- length(unique(unlist(links$source, links$target)))
    grps <- rep(NA, num.nodes)
    gval <- rep(NA, num.nodes)

    indexes <- 1:nrow(links)
    src.ind <- 1
    tgt.ind <- 2
    if (type == "Last variable")
    {
        indexes <- rev(indexes)
        src.ind <- 2
        tgt.ind <- 1
    }

    for (i in indexes)
    {
        tmp.node <- links[i,src.ind]
        if (is.na(grps[tmp.node+1]))
            grps[tmp.node+1] <- tmp.node
        else
            tmp.node <- grps[tmp.node+1]

        if (is.na(gval[links[i,tgt.ind]+1]) ||
            gval[links[i,tgt.ind]+1] < links$value[i])
        {
            grps[links[i,tgt.ind]+1] <- tmp.node
            gval[links[i,tgt.ind]+1] <- links$value[i]
        }
    }
    return(grps)
}


# Ensures it is always a 6-digit hex
colorsToHex <- function(xx)
{
    res <- sapply(xx, function(x){tryCatch({ tmp <- col2rgb(x, alpha = TRUE);
        return(rgb(t(tmp[1:3]), alpha = if (tmp[4] == 255) NULL else tmp[4], maxColorValue = 255)) },
        error=function(cond){NA})})
    ind <- which(is.na(res))
    if (length(ind) > 0)
    {
        res[ind] <- "#000000"
        for (i in ind)
            warning("Invalid color '", names(res)[i], "' replaced with '#000000'")
    }
    return(res)
}

#' computeLinks
#'
#' Computes the links between the nodes, so that can be expressed as a network.
#' @param data A \code{\link{data.frame}} or \code{\link{list}} of variables.
#' @param weights A numeric vector with length equal to the number of rows in
#'  \code{data}. This is used to adjust the width of the links.
#' @param show.percentages Whether to show percentages or counts
#' @importFrom stats xtabs
computeLinks <- function(data, weights, show.percentages = FALSE)
{
    links <- NULL
    counter <- 0
    n <- length(data)
    tmp.total <- nrow(data)

    for (i in 1:(n - 1)){
        x <- data[[i]]
        y <- data[[i + 1]]
        x.names <- levels(x)
        n.x <- length(x.names)
        n.y <- nlevels(y)
        x.y <- if (is.null(weights)) xtabs(~ x + y)
                  else               xtabs(weights ~ x + y)
        for (i.x in 1:n.x) {
            row.node <- counter + i.x - 1
            for (i.y in 1:n.y) {
                value <- x.y[i.x, i.y]
                if (value > 0)
                {
                    if (show.percentages)
                        value <- value / tmp.total * 100

                    column.node <- counter + n.x + i.y - 1
                    links <- rbind(links,c(row.node,column.node, value))
                }
            }
        }
        counter <- counter + n.x
    }
    links <- as.data.frame(links)
    names(links) <- c("source", "target", "value")
    links
}

#' @importFrom verbs Sum
nodeDictionary <- function(list.of.factors, weights, show.counts, show.percentages)
{
    nodes <- NULL
    for (vr in list.of.factors)
    {
        if (show.counts || show.percentages)
        {
            tmp.info <- if (is.null(weights)) xtabs(~vr)
                        else                  xtabs(weights~vr)

            suffix <- ""
            if (show.counts)
                suffix <- paste("n =", tmp.info)
            if (show.percentages)
            {
                denom <- if (!is.null(weights)) Sum(weights) else length(vr)
                suffix <- paste(suffix, sprintf("%.0f%%", tmp.info/denom*100),
                    sep = if (show.counts) ", " else "")
            }
            nodes <- c(nodes, paste0(levels(vr), " (", suffix, ")"))
        }
        else
            nodes <- c(nodes, levels(vr))
    }
    data.frame(name = nodes)
}

#' categorizeData
#'
#' Creates a dictionary of the nodes.
#' Quantizes numeric variables.
#' @param data A \code{\link{data.frame}} or \code{\link{list}} of variables.
#' @param weights A numeric vector giving the weight of each row in \code{data}.
#' @param max.categories When the number of unique values exceeds this value,
#' the variable is quantized.
#' @param share.values Whether the values in each variable are expected to be the same.
#' @param label.show.varname Whether to include the variable name in the node label.
#' @param label.max.length Maximum number of characters in the node label.
categorizeData <- function(data, weights, max.categories, share.values,
                           label.show.varname, label.max.length)
{
    var.names <- names(data)
    n <- length(var.names)
    breaks <- NULL
    .truncate <- function(labels, n) ifelse(nchar(labels) > label.max.length, 
        paste0(substr(labels, 1, label.max.length), " ..."), labels)

    if (share.values)
    {
        if (is.numeric(data[[1]]) && length(unique(unlist(data))) > max.categories)
        {
            n.cuts <- max.categories - if (any(is.na(unlist(data)))) 1 else 0
            tmp.range <- range(unlist(data), na.rm = TRUE)
            breaks <- seq(from = tmp.range[1], to = tmp.range[2], length = n.cuts+1)
            breaks[1] <- (1 - 0.001) * breaks[1]
            breaks[n.cuts+1] <- (1 + 0.001) * breaks[n.cuts+1]
        }
        all.dat <- quantizeVariable(unlist(data), max.categories, breaks = breaks)
    }
    num.values <- sapply(data, function(vv) length(unique(vv)))
    var.ordered <- order(num.values, decreasing = TRUE)

    # Numeric breaks can be identified without any interaction between variables
    # Factor is used if variable is non-numeric or if there are only a few values
    for (i in var.ordered)
        data[[i]] <- quantizeVariable(data[[i]], max.categories, breaks = breaks)

    for (i in var.ordered)
    {
        while (length(unique(data[[i]])) > max.categories)
        {
            node.change <- findNodesToMerge(data, i, weights)
            data[[i]] <- mergeNodes(data[[i]], node.change)
            if (share.values)
            {
               for (j in 1:n)
                    data[[i]] <- mergeNodes(data[[i]], node.change)
                all.dat <- mergeNodes(all.dat, node.change)
            }
        }
        data[[i]] <- addNA(data[[i]], ifany = TRUE)
        if (label.show.varname)
            levels(data[[i]]) <- paste(var.names[i], .truncate(levels(data[[i]])), sep = ": ")
        else
            levels(data[[i]]) <- paste0("", .truncate(levels(data[[i]])))
    }
    if (share.values)
        attr(data, "all.levels") <- levels(all.dat)
    data
}

#' findNodesToMerge
#'
#' Identifies the pair of nodes to merge which will minimize
#' distinct links. Note this uses the weights, so adding multiple
#' small links may be preferred to adding a single heavy link.
#' Note also that the link to the merged node i.e. (A,B+C) is the
#' has a weight that is the sum of (A,B) and (A,C). The relative
#' proportions are not taken into account.
#" If there are multiple node-pairs with the same number of distinct
#  links, then the pair to give the smallest merged node is used.
#' @param df dataframe containing factor variables only.
#' @param column index of the column which we are looking to merge nodes for
#' @param weights numeric vector containing weights for each row of \code{df}.
#' @importFrom verbs Sum SumEmptyHandling
findNodesToMerge <- function(df, column, weights = NULL)
{
    lvls <- levels(df[[column]])
    n <- length(lvls)

    # If there are many categories ignore linkage patterns
    # and merge the smallest categories in the old way
    if (n > 50)
    {
        tb <- sort(table(df[[column]]))
        return(names(tb)[1:2])
    }

    profile <- do.call("paste", c(df[-column], sep = "\\r"))
    if (is.null(weights))
        weights <- rep(1, nrow(df))

    # Minimize number fo links that are not shared and
    # If tied, minimize size of merged nodes
    m.diff <- matrix(NA, n, n) # number of differences
    m.size <- matrix(NA, n, n) # weight of merged node
    for (i in 1:(n-1))
    {
        for (j in (i+1):n)
        {
            ind.i  <- which(df[[column]] == lvls[i])
            ind.j  <- which(df[[column]] == lvls[j])
            ind.ij <- which(!profile[ind.i] %in% profile[ind.j])
            ind.ji <- which(!profile[ind.j] %in% profile[ind.i])
            m.diff[i,j] <- SumEmptyHandling(weights[c(ind.i[ind.ij], ind.j[ind.ji])], remove.missing = FALSE)
            m.size[i,j] <- SumEmptyHandling(weights[c(ind.i, ind.j)])
        }
    }
    min.diff <- min(m.diff, na.rm = TRUE)
    ind.min <- which(m.diff == min.diff)
    if (length(ind.min) > 1)
    {
        ind.min.sz <- which.min(m.size[ind.min])
        ind.min <- ind.min[ind.min.sz]
    }
    ind.min <- as.numeric(arrayInd(ind.min, .dim = dim(m.diff)))
    return(lvls[ind.min])
}

#' Merge levels in a factor
#'
#' @param x factor variable
#' @param old.nodes a vector containing factor levels to merge.
#'   Note level names are used rather than the numeric representation
#'   is used so it can be applied to multiple factors
#' @param nchar Maximum number of characters in each label.
mergeNodes <- function(x, old.nodes)
{
    new.node <- paste(old.nodes, collapse = ", ")
    ind <- match(old.nodes, levels(x))
    levels(x)[ind] <- new.node
    return(x)
}

#' @importFrom flipTransformations Factor
quantizeVariable <- function(x, max.categories, breaks = NULL)
{
    if (is.null(breaks) && (length(unique(x)) <= max.categories || is.factor(x)))
    {
        x.fac <- Factor(x)

    } else if (is.numeric(x)) # if breaks is supplied, this clause is always used
    {
        n.cuts <- max.categories - if(any(is.na(x))) 1 else 0
        x.fac <- if (is.null(breaks)) cut(x, max(2, n.cuts)) else cut(x, breaks)

    } else
    {
        x.fac <- rep("Text", length(x))
        x.fac[x == ""] <- "BLANK"
    }
    x.fac <- addNA(x.fac, ifany = TRUE)
    return(x.fac)
}

