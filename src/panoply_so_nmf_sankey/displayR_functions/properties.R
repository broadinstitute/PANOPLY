### NOTE: Sourced from https://github.com/Displayr/flipU/. Needed to run SankeyDiagram(). DO NOT EDIT.

#' Test if all items in a variable are integers
#'
#' @param x A vector.
#' @return logical.
#' @export
AllIntegers <- function(x)
{
    all(x == floor(x))
}

#' Find the names of the variables in a formula
#'
#' Handles \code{.} on right hand side of formula and \code{$} within backticks in
#' variable names.
#' @param formula A \code{\link{formula}}.
#' @param data A \code{\link{data.frame}} from which to extract variable
#' names if \code{.} is used in the formula.
#' @return A character vector of variable names appearing in \code{formula}.
#' @export
#' @importFrom stats terms
#' @examples
#' dat <- data.frame("dat$Var$y" = 1, x = 2, "`dat$Var$z`" = 3,
#'                                check.names = FALSE)
#' AllVariablesNames(`dat$Var$y` ~ ., data = dat)
#' AllVariablesNames(`dat$Var$y` ~ `dat$Var$z`, data = dat)
#' AllVariablesNames(`dat$Var$y` ~ `dat$Var$z`*x)
#' AllVariablesNames(log(y)~I(x^2))
AllVariablesNames <- function(formula, data = NULL)
{
    terms <- stats::terms(formula, data = data, keep.order = TRUE)
    vars <- attr(terms, "variables")

    out <- parseVars(vars)

    ##     ## add backticks for any var with non-syntactic names
    ##     idx1 <- make.names(out, unique = FALSE) != out
    ##     ## AND doesn't already have backticks
    ##     idx2 <- !grepl("^`.*`$", out)
    ##     ## AND doesnt have $ outside backticks
    ##     idx3 <- !grepl("[$](?=([^`]*`[^`]*`)*[^`]*$)", out, perl = TRUE)
    ##     idx <- idx1 & idx2 & idx3
    ## browser()
    ##     if (any(idx))
    ##         out[idx] <- paste0("`", out[idx], "`")

    ## need unique() below because of strange behaviour where
    ## response sometimes appears twice in output (once with
    ## backticks, once without) when dot on RHS
    unique(out)
}

#' @noRd
#' @param vars A call; the "variables" attribute of a
#' \code{\link{terms.object}}
#' @return A character vector of variable names
parseVars <- function(vars)
{
    vars <- vars[-1]  # remove "list" call
    out <- character(length(vars))
    ## Need for loop and not lapply to allow check for backticks in parseVar
    ## The check needs the call (`[`) not the call component (`[[`)
    for (i in seq_along(vars))
        out[i] <- parseVar(vars[i])
    out
}

#' Parses a component of the variables attribute
#' of a \code{\link{terms.object}}, which is a call to 'list'
#' of the variables in the model
#' @param var A call, a component of the "variables", attribute
#' of a call to \code{\link[stats]{terms}}
## @param label A character string, variable label
#' @return A character string containing the variable name
#' @importFrom stats as.formula
#' @noRd
parseVar <- function(var)
{
    ## var[[1L]] can have one of three classes
    ## 1) in simplest case, var has class "AsIs"; e.g. "x"
    ## 2)a) if var is backtick'd in formula, it has class "name"; e.g. "`x$y`"
    ## 2)b) if var is included b/c of "." on RHS of formula, it has class "name"
    ## 3) if var involves fun. call (including $), it has class "call";
    ## e.g, "log(x)" or "I(x^2)" or "dat$x" or df[["x"]] or a[1, "var", ]
    v1 <- var[[1L]]
    out <- as.character(v1)
    ## check for backticks
    if (inherits(v1, "name"))
    {
        ## check 1) don't backtick if variable enters through ".",
        ## but doesn't need backtick
        ## check 2) don't backtick if already backtick'd
        if (grepl("^`.*`[(][)]$", deparse(var)) && !grepl("^`.*`$", out))
            return(paste0("`", out, "`"))
    }
    else if (inherits(v1, "call"))
    {
        if (length(out) == 3L){
            out <- if (out[1L] == "[")
                paste0(out[2L], "[", addQuotesOrComma(var, out[3L]), "]")
            else if (out[1L] == "[[")
                paste0(out[2L], "[[", addQuotesOrComma(var, out[3L]), "]]")
            else  # $; hopefully, nothing else
                paste0(out[2L], out[1L], out[3L])
        }
        else if (length(out) > 1L)
        {
            ## deal with e.g. I(log(x)), I(x^2), I(dat$x); strip I and re-parse
            out <- if (out[1L] == "I")
                parseVar(as.formula(paste0("~", out[2L]))[2L])
            else if (out[1L] == "[")
                paste0(out[2L], "[", paste(vapply(out[3:length(out)],
                                                  addQuotesOrComma, "", v = var),
                                           collapse = ","), "]")
            else if (out[1L] == "[[")  ## unlikely, list extract with multiple indices
                paste0(out[2L], "[[", paste(vapply(out[3:length(out)],
                                                   addQuotesOrComma, "", v = var),
                                            collapse = ","), "]]")
            else  # e.g. log(x) with no I
                out[2L]
        }
    }
    out
}

#' For parsing terms involving extraction using [ or [[
#'
#' If the term is e.g. df[1, "y", ], then as.character(var[[1L]])
#' in parseVar is c("df", "[", "1", "y", "")
#' @noRd
addQuotesOrComma <- function(v, term)
{
    if (grepl("^$", term))  # if empty string, assume term is
        return("")
    out <- term
    if (grepl(paste0("(?:\\\"|\\')", term), deparse(v), perl = TRUE))
        out <- paste0("'", out, "'")
    out
}
#' Copy attributes from one object to another
#'
#' Copies attributes such as "label", "name", "question" and "questiontype"
#' from one object to another. If both objects are lists (data frames),
#' elements (columns) in the recipient list (data frame) will also receive
#' attributes from elements (columns) with the same name in the donor list
#' (data frame). The function is recursive, and will copy attributes when the
#' inputs are nested lists.
#' @param data.without.attributes an object to receive attributes from, such as
#' a data.frame, list, or matrix
#' @param data.with.attributes an object to copy attributes from
#' @param attr.to.not.copy character vector of attribute names appearing in
#' \code{data.with.attributes} that should not be copied
#' @return A copy of \code{data.without.attributes} with all the attributes
#' of \code{data.with.attributes}.
#' @details In the case when both arguments are data frames, any attributes
#' in the columns of \code{data.with.attributes} will also be copied to
#' \code{data.without.attributes} excluding \code{class} and \code{levels}
#' Names are used when copying attributes in each component. Nothing will be
#' copied for the case of lists with \code{NULL} names attribute.
#'
#' In the case when \code{data.without.attributes} is not a data frame
#' and \code{data.with.attributes} is a data frame, the
#' attributes of \code{data.with.attributes} are copied over to
#' \code{data.without.attributes} and the columns are not changed.
#'
#' Similarly, if \code{data.without.attributes} is a data frame and
#' and \code{data.with.attributes} is not a data frame the
#' attributes of \code{data.with.attributes} are copied over to
#' \code{data.without.attributes} and the columns are also not changed.
#' @examples
#' v1 <- 1:10
#' v2 <- 11:20
#' attr(v1, "label") <- "label for v1"
#' CopyAttributes(v2, v1) # returns v2 with the label attribute from v1
#'
#' df1 <- data.frame(a = 1:10, b = 11:20)
#' df2 <- data.frame(a = 21:30, b = 31:40)
#' attr(df1, "label") <- "label for df1"
#' attr(df1$a, "label") <- "label for column a of df1"
#' attr(df1$b, "label") <- "label for column b of df1"
#' CopyAttributes(df2, df1) # returns df2 with label attributes (including column attributes) from df1
#'
#' CopyAttributes(df2, v1) # returns df2 with label attribute from v1 (columns not changed)
#' CopyAttributes(v2, df1) # returns v2 with label attribute from df1 (columns not changed)
#' @export
CopyAttributes <- function(data.without.attributes, data.with.attributes,
                           attr.to.not.copy = c("dimnames", "names", "row.names",
                                                "dim", "class", "levels"))
{
    ## for data.frame recursion when arg1 has columns arg2 does not,
    ## atts.to.copy will be NULL and data.without.attributes is returned
    atts.to.copy <- names(attributes(data.with.attributes))
    atts.to.copy <- atts.to.copy[!atts.to.copy %in% attr.to.not.copy]
    for (a in atts.to.copy)
        attr(data.without.attributes, a) <- attr(data.with.attributes, a)

    if (is.list(data.without.attributes) && is.list(data.with.attributes))
        for (n in names(data.without.attributes))
            data.without.attributes[[n]] <- CopyAttributes(data.without.attributes[[n]],
                                                           data.with.attributes[[n]])

    data.without.attributes
}


#' \code{CopyAttributes}
#' @description Copies the "label", "name", "question" and "questiontype" attributes
#' for each for variable in a \code{\link{data.frame}}.
#' @param data.without.attributes A \code{\link{data.frame}}.
#' @param data.with.attributes A \code{\link{data.frame}}.
#' @return A \code{\link{data.frame}}.
#' @noRd
copyAttributesOld <- function(data.without.attributes, data.with.attributes)
{
    if (is.list(data.without.attributes))
    {
        for (i in seq_along(data.without.attributes))
            data.without.attributes[[i]] <- CopyAttributes(data.without.attributes[[i]],
                                                           data.with.attributes[[i]])
        return(data.without.attributes)
    }
    # Attention: the vector of attribute names below should be kept up to date
    # to match the attributes being assigned to Q variables.
    for (a in c("name", "label", "question", "questiontype"))
        attr(data.without.attributes, a) <- attr(data.with.attributes, a)
    data.without.attributes
}

#' Find the name of the outcome variable from a formula
#'
#' @param formula A \code{\link{formula}} or a \code{\link[stats]{terms}} object.
#' @param data A \code{\link{data.frame}} containing the variables in the formula.
#' Currently, ignored.
#' @return Character string giving the response variable name, or \code{NULL} if
#' no response is present in \code{formula}.
#' @export
OutcomeName <- function(formula, data = NULL)
{
    if (HasOutcome(formula))
    {
        if(inherits(formula, "terms"))
            class(formula) <- formula
        return(parseVar(formula[2]))
    }
    return(NULL)
}


#' Check if a formula contains an outcome variable
#'
#' @param formula A \code{\link{formula}}.
#' @return logical
#' @export
HasOutcome <- function(formula)
{
    length(formula) == 3
}

#' A print function for error checking
#'
#' Prints its name and a \code{\link{summary}}.
#' @param x Something to be printed.
#' @export
PrintDetails <- function(x)
{
    cat(paste0(deparse(substitute(x)), " n:", length(x), " valid:", sum(!is.na(x)), " missing:",sum(is.na(x)), "\n"))
    print(summary(x))
    cat("\n")
}


#' Test whether a vector contains any negative values
#'
#' @param x A vector.
#' @return logical.
#' @export
AnyNegative <- function(x)
{
    min(c(x, NA), na.rm = TRUE) < 0
}

#' Check if data, or, a model description, counts or represents counts
#'
#' @param x A variable or text string describing a family (e.g., "Poisson").
#' @return logical.
#' @export
IsCount <- function(x) {
    if(is.factor(x))
        return(FALSE)
    else if (is.logical(x))
        return (FALSE)
    if(!is.numeric(x)) {
        if (!is.character(x))
            x <- x$type
        return(x == "Poisson" | x == "Quasi-Poisson" | x == "NBD")
    }
    x <- x[!is.na(x)]
    if (length(x) == 0)
        stop("No data.")
    u = unique(x)
    if (min(u, na.rm = TRUE) < 0)
        return(FALSE)
    sum(as.integer(u) != u, na.rm = TRUE) == 0}


#' Return the outcome variable from a model
#'
#' @param formula A \code{\link{formula}}.
#' @param data A \code{\link{data.frame}} from which to extract the variable.
#' @return A vector of data.
#' @export
OutcomeVariable <- function(formula, data)
{
    data[[OutcomeName(formula, data)]]
}

#' Check that a subset contains data
#'
#' @param subset The filter used to filter data in a model.
#' @return true if the subset contains information
#' @export
HasSubset <- function(subset)
{
    !is.null(subset) & length(subset) != 1
}


#' Check if there are any NAs in a data frame
#'
#' @param data A \code{\link{data.frame}}.
#' @param formula A \code{\link{formula}}. Where supplied, only variables in the formula are checked.
#' @return logical.
#' @export
AnyNA <- function(data, formula = NULL)
{
    if (!is.null(formula))
    {
        data <- data[, AllVariablesNames(formula, data)]
    }
    any(is.na(data))
}
