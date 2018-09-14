#ifndef QR_MANAGER_H
#define QR_MANAGER_H

#include "Rcpp.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include <sstream>

namespace flmam {

/* Class to work with the QR decomposition from base::qr().
 * Options are provided for multiplication (QY or QtY),
 * or solution of a linear system based on back-solving.
 */

class QR_manager {
public:
    QR_manager(SEXP qr, SEXP qraux, SEXP pivot) : QR(qr), AUX(qraux), PIVOT(pivot), nobs(QR.nrow()), ncoef(QR.ncol()), info(0), lwork(-1) {
        if (AUX.size()!=ncoef) { 
            throw std::runtime_error("QR auxiliary vector should be of length 'ncol(Q)'"); 
        }
        if (PIVOT.size()!=ncoef) {
            throw std::runtime_error("pivoting vector should be of length 'ncol(Q)'");
        }
        if (ncoef && (*std::max_element(PIVOT.begin(), PIVOT.end()) < 1 || *std::min_element(PIVOT.begin(), PIVOT.end()) > ncoef)) {
            throw std::runtime_error("pivoting vector should contain values in [1, ncoef]");                        
        }
    
        // Workspace query (should only depend on 'side', 'nobs', and 'ncoef').
        double tmpwork=0;
        work.resize(nobs); // using 'work' as a dummy for 'rhs' in this call, just in case.
        constexpr char side='L', trans='T';
        constexpr int ncol=1;

        if (nobs) { 
            F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef,
                    QR.begin(), &nobs, AUX.begin(), work.data(), &nobs,
                    &tmpwork, &lwork, &info); 
            if (info) { 
                throw std::runtime_error("workspace query failed for 'dormqr'");
            }
        }
    
        lwork=static_cast<int>(tmpwork+0.5);
        work.resize(std::max(ncoef, lwork));
        return;
    } 
       
    int get_nobs() const {
        return nobs;
    }

    int get_ncoefs() const { 
        return ncoef;
    }

    void multiply(double* rhs, const char trans='T') {
        if (!nobs) { return; }

        constexpr char side='L';
        constexpr int ncol=1;
        F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef, QR.begin(), &nobs, AUX.begin(), rhs, &nobs, work.data(), &lwork, &info); 
        if (info) { 
            std::stringstream err;
            err << "Q matrix multiplication failed (error code " << info << ")";
            throw std::runtime_error(err.str().c_str());
        }
        return;
    }

    void backsolve(double* rhs) {
        if (!nobs) { return; }

        constexpr char uplo='U', diag='N', trans='N';
        constexpr int ncol=1;
        F77_CALL(dtrtrs)(&uplo, &trans, &diag, &ncoef, &ncol, QR.begin(), &nobs, rhs, &nobs, &info);
        if (info) { 
            std::stringstream err;
            err << "coefficient calculations failed (error code " << info << ")";
            throw std::runtime_error(err.str().c_str());
        }

        // Reordering the coefficients.
        for (size_t p=0; p<ncoef; ++p) {
            work[PIVOT[p]-1]=rhs[p];
        }
        std::copy(work.begin(), work.begin()+ncoef, rhs);
        return;
    }

private:
    Rcpp::NumericMatrix QR;
    Rcpp::NumericVector AUX;
    Rcpp::IntegerVector PIVOT;

    const int nobs, ncoef;
    int info, lwork;
    std::vector<double> work;
};

}

#endif
