\name{fitLM}
\alias{fitLM}

\title{Fit a linear model}
\description{Fit a linear model to any matrix representation.}

\usage{
fitLM(x, design) 
}

\arguments{
    \code{x}{A matrix-like object where each row is a feature and columns are samples.}
    \code{design}{A numeric design matrix with number of rows equal to \code{ncol(x)}.
        Alternatively, a factor of length \code{ncol(x)} for one-way layouts.}
}

\details{
This fits a linear model to each row of observations in \code{x}.
With the default settings, it effectively does the same calculations as:
\preformatted{    lm.fit(x=design, y=as.matrix(t(x)))
}
but without requiring the expensive realization of sparse or file-backed matrix representations.

If \code{design} is a matrix, the linear model is fitted by a QR decomposition and backsolving with LAPACK.
If \code{design} is a factor for a one-way layout, the coefficients are estimated directly as sample means.
}

\value{
A list containing:
\describe{
\item{\code{coefficients}:}{A numeric matrix with one row per row of \code{x} and one column per coefficient.
Contains the estimated values for each coefficient.}
\item{\code{variance}:}{A numeric vector of length equal to the number of rows of \code{x}.
Contains the residual variance for each row.}
}
}

\author{
Aaron Lun
}

\seealso{
\code{\link{lm.fit}}
}

\examples{
library(Matrix)
y <- rsparsematrix(10000, 100, 0.1)

# Using a design matrix.
cov <- runif(ncol(y))
X <- model.matrix(~cov)
fit <- fitLM(y, X)

# Using a one-way layout.
g <- gl(20, 5)
fit.o <- fitLM(y, g)
}
