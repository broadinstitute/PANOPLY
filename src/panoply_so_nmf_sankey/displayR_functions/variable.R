### NOTE: Sourced from https://github.com/Displayr/flipTransformations/. Needed to run SankeyDiagram(). DO NOT EDIT.

#' \code{Factor}
#' @description Converts to a factor, but retains the label attribute. Automatically combines duplicate
#' levels, with a warning (\code{\link{factor}} instead throws an error).
#' @param x A vector of data, usually taking a small number of distinct values.
#' @param ... Further arguments passed to \code{\link{factor}}.
#' @details If the variable is already a factor, removes any empty levels
#' @importFrom flipU CopyAttributes
#' @importFrom flipFormat Names
#' @export
Factor <- function(x, ...)
{

    if ("QDate" %in% class(x))
        # When the input is a Date from Q (class QDate)
        # extract the date categories which tell us
        # about the aggregation that has been used in Q.
        result <- attr(x, "QDate")
    else if ("factor" %in% class(x))
    {
        # Warning and addressing duplicate factor levels (this throws an error in 'factor' since R4.3)
        level.count <- table(levels(x))
        if (max(level.count) > 1)
        {
            duplicates <- paste0(names(level.count[level.count > 1]), collapse = ", ")
            warning(paste0(Names(x), " contains duplicate factor levels: ", duplicates, "."))
            result <- droplevels(factor(x, levels <- unique(levels(x)), ...))  # bug in 3.4.0
        }
        else
            result <- factor(x, ...)
    }
    else if (is.character(x)) # retain ordering in which it appears
        result <- factor(x, levels = unique(x))
    else
        result <- factor(x, ...)
    return(CopyAttributes(result, x))
}

#' \code{Ordered}
#' @description Converts to an ordered, but retains the label attribute.
#' @param x A vector of data, usually taking a small number of distinct values.
#' @param ... Further arguments passed to \code{factor}.
#' @details If the variable is already a factor, nothing is done to it.
#' @export
Ordered <- function(x, ...)
{
    if ("QDate" %in% class(x))
        # When the input is a Date from Q (class QDate)
        # extract the date categories which tell us
        # about the aggregation that has been used in Q.
        result <- attr(x, "QDate")
    else
        result <- ordered(x, ...)
    return(CopyAttributes(result, x))
}

#' Unclasses and removes the levels of a factor.
#'
#' @param x A factor.
#' @param warn Show warning.
#' @param use.numeric.levels If \code{TRUE}, numeric factor levels are used as the values of the output,
#'    else levels are output as 1, 2, 3, ...
#' @importFrom stats model.matrix
#' @importFrom flipFormat Labels
#' @export
Unclass <- function(x, warn = TRUE, use.numeric.levels = TRUE)
{
    in.loop <- !is.null(attr(x, "InLoop"))
    labels <- Labels(x, show.name = TRUE)
    ints <- suppressWarnings(as.numeric(as.character(x)))

    if (use.numeric.levels && !any(is.na(ints)))
    {
        to.factor.levels <- TRUE # Factor level labels have been used as values
        result <- ints
    }
    else
    {
        to.factor.levels <- FALSE # Levels converted to 1, 2, 3 ... ignoring the level labels
        result <- unclass(x)
        attr(result, "levels") <- NULL # As some functions get confused by levels and use as a proxy for factor.
    }
    if (in.loop && !to.factor.levels)
        attr(result, "Unclassed") <- labels
    else if (warn)
        warning(asNumericWarning(labels, to.factor.levels = to.factor.levels))
    result
}

#' \code{UnclassIfNecessary}
#' @description Unclasses a variable if it is a factor. Otherwise, returns x.
#' @param x A vector.
#' @param warn Show warning.
#' @return A vector
#' @export
UnclassIfNecessary <- function(x, warn = TRUE)
{
    if(is.factor(x))
        return(Unclass(x, warn));
    return(x);
}

asNumericWarning <- function(variables, to.factor.levels = FALSE)
{
    value.message <- if (to.factor.levels)
        "Values are assigned according to the labels of the categories."
    else
        "Values are assigned in the order of the categories: 1, 2, 3, ...;"
    if ((n <- length(variables)) > 1)
        full.variable.string <- paste0(paste0(variables[1:(n - 1)], collapse = ", "),
                                       " and ", variables[n])
    else
        full.variable.string <- variables
    if (!is.null(variables))
        variable.message <- ngettext(length(variables),
                                     paste0("The variable ", full.variable.string,
                                            " has been converted."),
                                     paste0("The variables ", full.variable.string,
                                            " have been converted."))
    else
        variable.message <- NULL
    paste("Data has been automatically converted to numeric.", value.message,
          "To use alternative numeric values, transform the data prior including it in this analysis",
          "(e.g. by changing its structure).",
          variable.message)
}

#' Convert an ordered factor to a numeric vector.
#'
#' @param x An ordered factor.
#' @details It is first checked if \code{x} is a categorical variable
#'     set (question) from Displayr (Q), in which case their value
#'     attributes are used.  To check this, (\code{x} is searched for
#'     attributes "questiontype", and "sourcevalues", "values",
#'     "codeframe" (for "PickOne" questiontype/vector \code{x}) and
#'     "sourcevariablevalues", "variablevalues", and "codeframe"
#'     attributes, in the case of "PickOneMulti"
#'     questiontype/data.frame \code{x}). See the Examples.
#'
#' If \code{x} is missing these attributes, \code{\link{Unclass}} is
#' used. If all labels of \code{x} (or columns of \code{x}) can be
#' coerced to numeric, these values will be used to create the numeric
#' variable; otherwise, the integers 1, 2, ... as per the ordering of
#' the factors.
#' @seealso \code{\link{numbersFromCategoricalVariableSets}}, \code{\link{Unclass}}
#' @importFrom stats model.matrix
#' @examples
#' file <- system.file("extdata", "variable.sets.rda", package = "flipTransformations")
#' vs.env <- new.env()
#' load(file, vs.env)
#' ord.num <- AsNumeric(vs.env$ordinal, TRUE)
#' ## Compare
#' table(ord.num)
#' levels(vs.env$ordinal)
#' @export
OrderedToNumeric <- function(x)
{
    is.variable.set <-isVariableSet(x)
    if (is.variable.set)
        out <- numbersFromCategoricalVariableSets(x)
    else
        out <- Unclass(x)
    return(CopyAttributes(out, x))
}

#' Convert a factor variable to a numeric vector
#'
#' Convert a factor variable to a numeric vector (when the factor is ordered),
#' or a matrix of indicator variables (when the factor is not ordered).
#' @param x A factor or ordered factor.
#' @param binary Returns the factor as binary variables. If \code{FALSE}, \code{\link{OrderedToNumeric}}.
#' @param name The name of the variable.
#' @param remove.first Remove the first binary variable, if a binary variable is being created.
#' @importFrom flipFormat RemoveParentName Names
#' @importFrom flipU CopyAttributes
#' @seealso \code{\link[stats]{model.matrix}}, \code{\link{AsNumeric}}, \code{\link{OrderedToNumeric}}
#' @return a vector of length \code{x}, if \code{x} is an ordered factor; otherwise, a 0-1 matrix
#' with number of rows equal to the length of \code{x} an number of columns equal to the number of
#' levels in \code{x}
#' @export
FactorToNumeric <- function(x, binary = TRUE, name = NULL, remove.first = TRUE)
{
    if (!binary)#(is.ordered(x))
        return(OrderedToNumeric(x))
    if (is.null(name))
        name = RemoveParentName(Names(x))
    indicators <- FactorToIndicators(x, name)
    if (nrow(indicators) < length(x))
    {
        new.indicators <- matrix(NA, length(x), ncol(indicators))
        row.names <- as.numeric(dimnames(indicators)[[1]])
        colnames(new.indicators) <- colnames(indicators)
        new.indicators[row.names, ] <- as.matrix(indicators)
        new.indicators <- as.data.frame(new.indicators)
        indicators <- CopyAttributes(new.indicators, indicators)
    }
    if (remove.first)
    {
        new.indicators <- indicators[, -1, drop = FALSE]
        indicators <- CopyAttributes(new.indicators, indicators)
    }
    return(indicators)
}

#' \code{FactorToIndicators}
#' @description Convert a factor variable to a matrix whose columns are binary variables
#' representing one of the levels from the factor variable.
#' @param variable The factor variable to convert.
#' @param name The name of the input variable.
#' @importFrom stats model.matrix
#' @export
FactorToIndicators <- function(variable, name = NULL)
{
    if (is.null(name))
        name <- Names(variable)
    if (nlevels(variable) == 1) { 
        # When only a single level, model.matrix will fail.
        # Instead, create a single-column matrix where all
        # non-na values are set to 1.
        missings <- is.na(variable)  
        result <- matrix(as.numeric(!missings), ncol = 1)
        result[missings] <- NA

    } else {
        result <- stats::model.matrix( ~ variable - 1)    
    }
    levs <- levels(variable)
    colnames(result) <- paste0(name, ".", 1:nlevels(variable))
    result <- as.data.frame(result)
    result <- CopyAttributes(result, variable)
    label <- Labels(variable)
    if (!is.null(label))
    {
        labels <- paste0(label, ": ", levs)
        for (i in 1:nlevels(variable))
            if (!is.null(label))
                attr(result[, i], "label") <- labels[i]
    }
    result
}

#' Converts a list of variable or data frames into a data.frame.
#'
#' @param variable A variable in a DataSet or data.frame.
#' @param cutoff The cutoff point to split the variable into.
#' @param warning If TRUE, raise a warning showing the new levels.
#' @param name An alternate name to show instead of the deparsed variable name.
#' @importFrom flipFormat RemoveParentName
#' @export
DichotomizeFactor <- function(variable, cutoff = 0.5, warning = FALSE,
                              name = RemoveParentName(deparse(substitute(variable)))) {
    label <- attr(variable, "label")
    if (is.null(label))
        label <- name
    if (!is.factor(variable))
        variable <- factor(variable)
    if (nlevels(variable) == 1)
        stop(paste(deparse(substitute(variable)), "cannot be dichotimized as it only contains one level."))
    else if (nlevels(variable) == 2)
        return(variable)
    cumulative.probs <- cumsum(prop.table(table(variable)))
    cut.point <- match(TRUE, cumulative.probs > cutoff)
    if (cut.point == 1)
        stop(paste(name, "cannot be dichotimized (e.g., perhaps only has 1 value)."))
    new.factor <- factor(suppressWarnings(Unclass(variable, use.numeric.levels = FALSE)) >= cut.point)
    levels(new.factor) <- paste(c("<=", ">="), levels(variable)[c(cut.point - 1, cut.point )])
    if (warning)
        warning(paste(name, "has been dichotimized into", paste(levels(new.factor), collapse = " & ")))
    attr(new.factor, "label") <- paste(label, levels(new.factor)[2])
    new.factor
}



#' Turn the outcome variable from a formula into a factor
#'
#'   if not already a factor.
#' @param formula A formula.
#' @param data A data.frame containing the variable.
#' @importFrom flipU OutcomeName
#' @export
CreatingFactorDependentVariableIfNecessary <- function(formula, data)
{
    outcome.name <- OutcomeName(formula, data)
    data[, outcome.name] <- Factor(data[[outcome.name]])
    data
}


#' \code{CreatingBinaryDependentVariableIfNecessary}
#' @description Dichotomizes the dependent variable in a data.frame if not already dichotomized.
#' @param formula A formula.
#' @param data A data.frame
#' @importFrom flipU OutcomeName
#' @export
CreatingBinaryDependentVariableIfNecessary <- function(formula, data)
{
    outcome.name <- OutcomeName(formula, data)
    data[, outcome.name] <- CreatingBinaryVariableIfNecessary(data, outcome.name)
    data
}

#' \code{CreatingBinaryVariableIfNecessary}
#' @description Dichotomizes a variable.
#' @param data A data.frame
#' @param name The name of the variable in the data.frame.
#' @export
CreatingBinaryVariableIfNecessary <- function(data, name)
{
    variable <- data[[name]]
    n.unique <- length(unique(variable))
    if (n.unique < 2)
        warning("The Outcome variable needs to contain two or more categories. It does not.")
    else
    {
        if(!is.factor(variable))
        {
            variable <- Factor(variable)
        }
        if (nlevels(variable) > 2)
            variable <- DichotomizeFactor(variable, warning = TRUE, name = name)
    }
    variable
}

#' @title ProcessQVariables
#' @description Processes Q variables, e.g.: converting date variables to be categorical
#' based on their period. This function should be called to process Q variables before they are used.
#' @param x A Q variable or a data frame containing Q variables.
#' @importFrom flipU CopyAttributes
#' @export
ProcessQVariables <- function(x)
{
    .processQVariable <- function(v)
    {
        if ("QDate" %in% class(v))
            CopyAttributes(attr(v, "QDate"), v, c("dimnames", "names",
                                                "dim", "class", "QDate"))
        else
            v
    }

    if (is.null(x))
        NULL
    else if (is.data.frame(x))
        data.frame(lapply(x, function(v) .processQVariable(v)), check.names = FALSE,
                   stringsAsFactors = FALSE, row.names = attr(x, "row.names"))
    else
        .processQVariable(x)
}
