dist <-
function(x, y = NULL, method = NULL, ...,
         diag = FALSE, upper = FALSE,
         by_rows = TRUE, auto_convert = TRUE)
{

### PARAMETER HANDLING
    ## convenience hack to allow dist(x, "method")
    if ((is.function(y) || is.character(y)) && is.null(method)) {
        method <- y
        y <- NULL
    }

    ## method lookup
    reg_entry <- NULL
    if (is.null(method))
        method <- if (is.data.frame(x))
            "Gower"
        else if (is.logical(x))
            "Jaccard"
        else
            "Euclidean"
    if (!is.function(method))
        reg_entry <- pr_DB$get_entry(method)

    ## vector handling
    if (is.vector(x) && is.atomic(x))
        x <- as.matrix(x)
    if (!is.null(y) && is.vector(y) && is.atomic(y))
        y <- as.matrix(y)

    ## some checks
    if (!is.data.frame(x) && !is.matrix(x) && !is.list(x))
        stop("Can only handle data frames, vectors, matrices, and lists!")
    if (!is.null(y)) {
        if (is.data.frame(x) && !is.data.frame(y)
            || is.matrix(x) && !is.matrix(y)
            || is.list(x) && !is.list(y))
            stop("x and y must be of same type.")
        if (is.matrix(x) && is.matrix(y) || is.data.frame(x) && is.data.frame(y))
            if (by_rows && (ncol(x) != ncol(y)))
                stop("x and y must be conform in columns.")
            else if (!by_rows && (nrow(x) != nrow(y)))
                stop("x and y must be conform in rows.")
    }

### PREPROCESS
    params <- list(...)
    if (!is.null(reg_entry)) {
        if(!is.na(reg_entry$PREFUN)) {
            tmp <- do.call(reg_entry$PREFUN, c(list(x, y, params, reg_entry)))
            if (!is.null(tmp)) {
                x <- tmp$x
                y <- tmp$y
                params <- tmp$p
                reg_entry <- tmp$reg_entry
            }
        }
        method <- reg_entry$FUN
    }

    ## helper function for calling the C-level loops
    .proxy_external <- function(CFUN, x, y)
        do.call(".External",
                c(list(CFUN, x, y,
                       if (!is.function(method)) get(method) else method),
                  params))

    result <-
### PASS-THROUGH-cases
        if (!is.null(reg_entry) && !reg_entry$loop) {
            if (reg_entry$C_FUN)
                do.call(".Call", c(list(method), list(x), list(y), params))
            else
                do.call(method, c(list(x), list(y), params))
        } else if (is.null(y)) {
### LOOP WORKHORSE for auto-proximities
            ## transpose data for column-wise loop
            if (!by_rows && !is.list(x))
                x <- t(x)
            if (is.matrix(x) && !is.null(reg_entry) && reg_entry$abcd)
                ## binary matrix
                .proxy_external("R_apply_dist_binary_matrix", x != 0, NULL)
            else if (is.matrix(x))
                ## real, integer matrix
                .proxy_external("R_apply_dist_matrix", x, NULL)
            else if (is.list(x) && !(is.data.frame(x) && by_rows))
                ## list
                .proxy_external("R_apply_dist_list", x, NULL)
            else ## data frame (by rows)
                .proxy_external("R_apply_dist_data_frame", x, NULL)

        } else {
### LOOP WORKHORSE for cross-proximities
            ## transpose data for column-wise loop
            if (!by_rows && !is.list(x)) {
                x <- t(x)
                y <- t(y)
            }
            if (is.matrix(x) && !is.null(reg_entry) && reg_entry$abcd)
                ## binary matrices
                .proxy_external("R_apply_dist_binary_matrix", x != 0, y != 0)
            else if (is.matrix(x))
                ## real, integer matrices
                .proxy_external("R_apply_dist_matrix", x, y)
            else if (is.list(x) && !(is.data.frame(x) && by_rows))
                ## lists
                .proxy_external("R_apply_dist_list", x, y)
            else ## data frames (by rows)
                .proxy_external("R_apply_dist_data_frame", x, y)

        }

### set col/rownames for cross-proximity-objects (if needed)
    if (is.matrix(result) && is.null(dimnames(result)))
        if (is.list(x) && !is.data.frame(x)) {
            rownames(result) <- names(x)
            colnames(result) <- names(y)
        } else if (by_rows) {
            rownames(result) <- rownames(x)
            colnames(result) <- rownames(y)
        } else {
            rownames(result) <- colnames(x)
            colnames(result) <- colnames(y)
        }

### POSTPROCESS
    if (!is.null(reg_entry)) {
        if (!is.na(reg_entry$POSTFUN))
            result <- do.call(reg_entry$POSTFUN, c(list(result, params)))
        if (!reg_entry$distance && !(is.logical(auto_convert) && !auto_convert)) {
            result <- if (is.function(auto_convert) || is.character(auto_convert))
                do.call(auto_convert, list(result))
            else if (is.null(reg_entry$convert))
                pr_simil2dist(result)
            else
                do.call(reg_entry$convert, list(result))
        }
        method <- reg_entry$names[1]
    }

### RETURN DIST-OBJECT
    result <-
        if (is.matrix(result))
            structure(result, class = "crossdist")
        else {
            if (!inherits(result, "dist"))
                stop("debug missing class")
            structure(result, Diag = diag,
                      Upper = upper)
        }
    structure(result,
              method = if (is.character(method)) method else deparse(substitute(method)),
              call = match.call())
}

simil <-
function(x, y = NULL, method = NULL, ...,
         diag = FALSE, upper = FALSE,
         by_rows = TRUE, auto_convert = TRUE)
{
    ## convenience to allow dists(x, "method")
    if ((is.function(y) || is.character(y)) && is.null(method)) {
        method <- y
        y <- NULL
    }
    if (is.null(method))
        method <- if (is.data.frame(x))
            "Gower"
        else if (is.logical(x))
            "Jaccard"
        else
            "correlation"

    ret <- dist(x, y, method, ..., diag = diag, upper = upper,
                by_rows = by_rows, auto_convert = FALSE)

    ## possibly convert to similarity
    reg_entry <- pr_DB$get_entry(attr(ret, "method"), stop_if_missing = FALSE)
    if (!is.null(reg_entry)) {
        if (reg_entry$distance && !(is.logical(auto_convert) && !auto_convert)) {
            ret <- if (is.function(auto_convert) || is.character(auto_convert))
                do.call(auto_convert, list(ret))
            else if (is.null(reg_entry$convert))
                pr_simil2dist(ret)
            else
                do.call(reg_entry$convert, list(ret))
        }
    }

    class(ret) <- c(if (inherits(ret, "crossdist")) "crosssimil" else "simil",
                    class(ret))
    ret
}

as.matrix.simil <-
function (x, ...)
{
    ret <- NextMethod(x, ...)
    diag(ret) <- 1
    ret
}

as.simil <-
function(x, FUN = NULL)
{
    if (inherits(x, c("dist", "crossdist"))) {
        class(x) <- if (inherits(x, "dist"))
            c("simil", class(x))
        else
            c("crosssimil", class(x))
        if (!is.null(FUN))
            FUN(x)
        else {
            reg_entry <- pr_DB$get_entry(attr(x, "method"), stop_if_missing = FALSE)
            if (!is.null(reg_entry) && !is.null(reg_entry$convert))
                do.call(reg_entry$convert, list(x))
            else
                pr_dist2simil(x)
        }
    } else if (inherits(x, "simil"))
        x
    else
        structure(stats::as.dist(x), class = c("simil", "dist"))
}

as.dist <-
function(x, FUN = NULL)
{
    if (inherits(x, c("simil","crosssimil"))) {
        class(x) <- if (inherits(x, "simil"))
            c("dist", class(x))
        else
            c("crossdist", class(x))
        if (!is.null(FUN))
            FUN(x)
        else {
            reg_entry <- pr_DB$get_entry(attr(x, "method"), stop_if_missing = FALSE)
            if (!is.null(reg_entry) && !is.null(reg_entry$convert))
                do.call(reg_entry$convert, list(x))
            else
                pr_simil2dist(x)
        }
    } else if (inherits(x, "dist"))
        x
    else
        stats::as.dist(x)
}

print.crossdist <-
print.crosssimil <-
function (x, digits = getOption("digits"),
          justify = "none", right = TRUE, ...)
{
    if (length(x) > 0) {
        m <- as.matrix(x)
        cf <- format(m, digits = digits, justify = justify)
        print(cf, quote = FALSE, right = right, ...)
    } else {
        cat(data.class(x), "(0)\n", sep = "")
    }
    invisible(x)
}

pr_simil2dist <-
function(x)
    1 - x

pr_dist2simil <-
function(x)
    1 / (1 - x)

###
