#include "flmam.h"

template<class M>
SEXP fit_lm_internal (M emat, SEXP qr, SEXP groups, SEXP subset) {
    auto ptr=flmam::dispatcher(qr, groups);

    const int ncoefs=ptr->get_ncoefs();
    const int ncells=ptr->get_nobs();
    if (ncells!=int(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    }

    // Checking the subset input.
    Rcpp::IntegerVector rows(subset);
    const int all_genes=emat->get_nrow();
    const size_t ngenes=rows.size();
    for (auto r : rows) {
        if (r == NA_INTEGER || r< 0 || r >= all_genes) {
            throw std::runtime_error("row index is NA or out of range");
        }
    }

    // Setting up output objects (note 'coefs' is transposed).
    Rcpp::NumericMatrix coefs(ncoefs, ngenes);
    Rcpp::NumericVector vars(ngenes);

    // Running through each gene and fitting the linear model.
    Rcpp::NumericVector tmp(ncells);
    auto vIt=vars.begin();
    auto cIt=coefs.begin();
    
    for (auto rIt=rows.begin(); rIt<rows.end(); ++rIt, ++vIt, cIt+=ncoefs) {
        emat->get_row(*rIt, tmp.begin());
        (*vIt)=ptr->fit(tmp.begin(), cIt);
    }
    
    return Rcpp::List::create(coefs, vars);
}

SEXP fit_lm (SEXP exprs, SEXP qr, SEXP groups, SEXP subset) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto emat=beachmat::create_integer_matrix(exprs);
        return fit_lm_internal(emat.get(), qr, groups, subset);
    } else {
        auto emat=beachmat::create_numeric_matrix(exprs);
        return fit_lm_internal(emat.get(), qr, groups, subset);
    }
    END_RCPP
}


