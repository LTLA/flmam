#' @export
fitLM <- function(x, design) 
# Fits a linear model in 'design' to each row of 'x'.
# 
# written by Aaron Lun
# created 14 September 2018
{
    if (!is.null(dim(design))) {
        qr.out <- qr(design, LAPACK=TRUE)
        d <- diag(qr.out$qr)

        if (!all(abs(d) > 1e-8)) { 
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

    results <- .Call(cxx_fit_lm, x, qr.out$qr, qr.out$qraux, groups)
    names(results) <- c("coefficients", "variance")
    results$coefficients <- t(results$coefficients)

    dimnames(results$coefficients) <- list(rownames(x), coef.names)
    names(results$variance) <- rownames(x)
    results
}
