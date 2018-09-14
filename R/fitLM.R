#' @export
#' @importFrom BiocParallel bplapply bpnworkers bpparam
fitLM <- function(x, design, rows=NULL, BPPARAM=bpparam())
# Fits a linear model in 'design' to each row of 'x'.
# 
# written by Aaron Lun
# created 14 September 2018
{
    if (!is.null(dim(design))) {
        if (nrow(design)==0L) {
            qr.out <- list(qr=matrix(0, nrow(design), ncol(design)), qraux=numeric(ncol(design)), pivot=seq_len(ncol(design)))
        } else {
            qr.out <- qr(design, LAPACK=TRUE)
        }

        diags <- abs(diag(qr.out$qr))
        mat.rank <- sum(diags >= .Machine$double.eps * max(c(diags, -Inf))) # -Inf to avoid warnings.
        if (mat.rank < ncol(design)) { 
            stop("design matrix is not of full rank")
        }

        groups <- NULL
        coef.names <- colnames(design)
    } else {
        qr.out <- list()
        groups <- droplevels(factor(design))
        coef.names <- levels(groups)
        groups <- as.integer(groups) - 1L
    }

    # Splitting the specified rows into jobs (note the zero indexing for C++).
    if (!is.null(rows)) {
        if (is.logical(rows)) { rows <- which(rows) }
        else if (is.character(rows)) { rows <- match(rows, rownames(x)) }
        else { rows <- as.integer(rows) }
    } else {
        rows <- seq_len(nrow(x))
    }

    n.cores <- bpnworkers(BPPARAM)
    if (n.cores > 1L && length(rows)) {
        by.core <- split(rows - 1L, cut(seq_along(rows), n.cores, labels=FALSE))
    } else {
        by.core <- list(rows - 1L)
    }

    # Collating the results across all cores.
    results <- bplapply(by.core, FUN=.fit_lm_internal, x=x, qr.out=qr.out, groups=groups)
    all.coef <- lapply(results, "[[", i=1)
    all.coef <- t(do.call(cbind, all.coef))
    all.var <- unlist(lapply(results, "[[", i=2))

    results <- list(coefficients=all.coef, variance=all.var)
    dimnames(results$coefficients) <- list(rownames(x)[rows], coef.names)
    names(results$variance) <- rownames(x)[rows]
    results
}

.fit_lm_internal <- function(x, qr.out, groups, subset) {
    .Call(cxx_fit_lm, x, qr.out, groups, subset)
}
