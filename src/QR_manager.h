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
    QR_manager(SEXP qr, SEXP qraux) : QR(qr), AUX(qraux), nobs(QR.nrow()), ncoef(QR.ncol()), info(0), lwork(-1) {
        if (AUX.size()!=ncoef) { 
            throw std::runtime_error("QR auxiliary vector should be of length 'ncol(Q)'"); 
        }
    
        // Workspace query (should only depend on 'side', 'nobs', and 'ncoef').
        double tmpwork=0;
        const char dummytrans='T';
        work.resize(nobs); // using 'work' as a dummy for 'rhs' in this call, just in case.
    
        F77_CALL(dormqr)(&side, &dummytrans, &nobs, &ncol, &ncoef,
                QR.begin(), &nobs, AUX.begin(), work.data(), &nobs,
                &tmpwork, &lwork, &info); 
        if (info) { 
            throw std::runtime_error("workspace query failed for 'dormqr'");
        }
    
        lwork=static_cast<int>(tmpwork+0.5);
        work.resize(lwork);
        return;
    } 
       
    int get_nobs() const {
        return nobs;
    }
    
    int get_ncoefs() const { 
        return ncoef;
    }

    void multiply(double* rhs, const char trans='T') {
        F77_CALL(dormqr)(&side, &trans, &nobs, &ncol, &ncoef, QR.begin(), &nobs, AUX.begin(), rhs, &nobs, work.data(), &lwork, &info); 
        if (info) { 
            std::stringstream err;
            err << "Q matrix multiplication failed (error code " << info << ")";
            throw std::runtime_error(err.str().c_str());
        }
        return;
    }

    void backsolve(double* rhs) {
        F77_CALL(dtrtrs)(&uplo, &xtrans, &diag, &ncoef, &ncol, QR.begin(), &nobs, rhs, &nobs, &info);
        if (info) { 
            std::stringstream err;
            err << "coefficient calculations failed (error code " << info << ")";
            throw std::runtime_error(err.str().c_str());
        }
        return;
    }

private:
    Rcpp::NumericMatrix QR;
    Rcpp::NumericVector AUX;

    const int nobs, ncoef;
    int info, lwork;
    std::vector<double> work;

    static const int ncol=1;
    static const char side='L', uplo='U', xtrans='N', diag='N';
};

}

#endif
