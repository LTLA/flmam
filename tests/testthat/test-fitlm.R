# Tests the fitLM function.
# library(flmam); library(testthat); source("test-fitlm.R")

set.seed(1000)
test_that("fitLM works for design matrices", {
    library(Matrix)
    y <- rsparsematrix(1000, 100, 0.1)

    cov <- runif(ncol(y))
    X <- model.matrix(~cov)
    fit <- fitLM(y, X)

    ref <- lm.fit(x=X, y=as.matrix(t(y)))
    expect_equal(fit$coefficients, t(ref$coefficients))
    var.est <- colMeans(ref$effects[-seq_len(ref$rank),]^2)
    expect_equal(fit$variance, var.est)

    # Same results for a dense array.
    fit.dense <- fitLM(as.matrix(y), X)
    expect_equal(fit, fit.dense)

    # Behaves sensibly with just one group.
    X.solo <- cbind(rep(1, ncol(y)))
    fit.solo <- fitLM(y, X.solo)
    expect_equal(fit.solo$coefficients[,1], Matrix::rowMeans(y))
    expect_equal(fit.solo$variance, DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(y)))
   
    # Handles more complex design matrices (with non-trivial pivoting).
    f1 <- factor(rep(1:10, each=10))
    f2 <- factor(rep(1:10, 10))
    X.add <- model.matrix(~f1 + f2)
    expect_false(identical(qr(X.add, LAPACK=TRUE)$pivot, 1:ncol(X.add)))

    fit <- fitLM(y, X.add)
    ref <- lm.fit(x=X.add, y=as.matrix(t(y)))
    expect_equal(fit$coefficients, t(ref$coefficients))
    var.est <- colMeans(ref$effects[-seq_len(ref$rank),]^2)
    expect_equal(fit$variance, var.est)
})

set.seed(1001)
test_that("fitLM works for one-way layouts", {
    y <- rsparsematrix(1000, 100, 0.1)
    
    grouping <- factor(sample(10, ncol(y), replace=TRUE))
    X <- model.matrix(~0 + grouping)
    colnames(X) <- levels(grouping)
    fit <- fitLM(y, grouping)

    ref <- lm.fit(x=X, y=as.matrix(t(y)))
    expect_equal(fit$coefficients, t(ref$coefficients))
    var.est <- colMeans(ref$effects[-seq_len(ref$rank),]^2)
    expect_equal(fit$variance, var.est)

    # Behaves sensibly with just one group.
    solo <- integer(ncol(y))
    fit.solo <- fitLM(y, solo)
    expect_equal(fit.solo$coefficients[,1], Matrix::rowMeans(y))
    expect_equal(fit.solo$variance, DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(y)))

    # Works properly for single-sample groups.
    outlier <- rep(1:2, c(1, 99))
    fit.out <- fitLM(y, outlier)
    expect_equal(fit.out$coefficients[,1], y[,1])
    expect_equal(fit.out$coefficients[,2], Matrix::rowMeans(y[,-1]))
    expect_equal(fit.out$variance, DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(y[,-1])))

    # Same results for a dense array.
    fit.dense <- fitLM(as.matrix(y), X)
    expect_equal(fit, fit.dense)
})

set.seed(1002)
test_that("fitLM fails gracefully for silly inputs", {
    y <- matrix(rnorm(1000), ncol=10)
    design <- cbind(1, rep(1, ncol(y)))
    expect_error(fitLM(y, design), "full rank")

    # No features.
    zero.des <- fitLM(y[0,], design[,1,drop=FALSE])
    expect_identical(dim(zero.des$coefficients), c(0L, 1L))
    expect_identical(length(zero.des$variance), 0L)

    zero.group <- fitLM(y[0,], gl(5, 2))
    expect_identical(dim(zero.group$coefficients), c(0L, 5L))
    expect_identical(length(zero.group$variance), 0L)
    
    # No libraries.
    expect_error(empty.des <- fitLM(y[,0], design[0,]), "not of full rank")
    empty.des <- fitLM(y[,0], design[0,0])
    expect_identical(dim(empty.des$coefficients), c(nrow(y), 0L))
    expect_identical(empty.des$variance, rep(NaN, nrow(y)))

    empty.group <- fitLM(y[,0], integer(0))
    expect_identical(dim(empty.group$coefficients), c(nrow(y), 0L))
    expect_identical(empty.group$variance, rep(NaN, nrow(y)))

    # No coefficients.
    free.des <- fitLM(y, design[,0])
    expect_identical(dim(free.des$coefficients), c(nrow(y), 0L))
    expect_equal(rowMeans(y^2), free.des$variance)

    # No residual d.f.
    full.des <- fitLM(y, diag(ncol(y)))
    expect_equivalent(full.des$coefficients, y)
    expect_true(all(is.na(full.des$variance)))

    full.group <- fitLM(y, seq_len(ncol(y)))
    expect_equivalent(full.group$coefficients, y)
    expect_true(all(is.na(full.group$variance)))
})

