#' @export
fitLM <- function(x, design) 
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

    results <- .Call(cxx_fit_lm, x, qr.out, groups)
    names(results) <- c("coefficients", "variance")
    results$coefficients <- t(results$coefficients)

    dimnames(results$coefficients) <- list(rownames(x), coef.names)
    names(results$variance) <- rownames(x)
    results
}
