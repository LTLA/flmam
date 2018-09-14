#include "flmam.h"

template<class M>
SEXP fit_lm_internal (M emat, SEXP qr, SEXP groups) {
    auto ptr=flmam::dispatcher(qr, groups);

    const int ncoefs=ptr->get_ncoefs();
    const int ncells=ptr->get_nobs();
    if (ncells!=int(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    }
    
    // Setting up output objects (note coefs is transposed).
    const size_t ngenes=emat->get_nrow();
    Rcpp::NumericMatrix coefs(ncoefs, ngenes);
    Rcpp::NumericVector vars(ngenes);

    // Running through each gene and fitting the linear model.
    Rcpp::NumericVector tmp(ncells);
    auto vIt=vars.begin();
    auto cIt=coefs.begin();
    
    for (size_t s=0; s<ngenes; ++s, ++vIt, cIt+=ncoefs) {
        emat->get_row(s, tmp.begin());
        (*vIt)=ptr->fit(tmp.begin(), cIt);
    }
    
    return Rcpp::List::create(coefs, vars);
}

SEXP fit_lm (SEXP exprs, SEXP qr, SEXP groups) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto emat=beachmat::create_integer_matrix(exprs);
        return fit_lm_internal(emat.get(), qr, groups);
    } else {
        auto emat=beachmat::create_numeric_matrix(exprs);
        return fit_lm_internal(emat.get(), qr, groups);
    }
    END_RCPP
}


