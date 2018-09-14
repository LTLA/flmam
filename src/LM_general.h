#ifndef LM_GENERAL_H
#define LM_GENERAL_H

#include "LM_base.h"
#include "QR_manager.hpp"

namespace flmam {

class LM_general : public LM_base {
public:
    LM_general(SEXP, SEXP);
    ~LM_general() = default;
    LM_general(const LM_general&) = default;
    LM_general& operator=(const LM_general&) = default;
    LM_general(LM_general&&) = default;
    LM_general& operator=(LM_general&&) = default;

    // Virtual methods
    double fit (const double*, double*);
    int get_nobs() const;
    int get_ncoefs() const;

    double fit0 (double*, bool=false);
private:
    QR_manager qr;
};

LM_general::LM_general (SEXP Qr, SEXP Aux) : qr(Qr, Aux) {
    workspace.resize(get_nobs());
}

double LM_general::fit(const double* values, double* out_coefs) {
    std::copy(values, values+get_nobs(), workspace.begin());
    double resid_var=fit0(workspace.data());
    std::copy(workspace.begin(), workspace.begin()+get_ncoefs(), out_coefs);
    return resid_var; 
}

double LM_general::fit0 (double* in, bool var_only) {
    qr.multiply(in);
    const int nobs=qr.get_nobs(), ncoefs=qr.get_ncoefs();
    
    double curvar=0;
    for (int i=ncoefs; i<nobs; ++i) {
        const auto& curval=in[i];
        curvar += curval * curval;
    }
    curvar /= nobs - ncoefs;
    
    if (!var_only) {
        qr.solve(in);
    }
    return curvar;
}

int LM_general::get_nobs() const{
    return qr.get_nobs();
}

int LM_general::get_ncoefs() const{
    return qr.get_ncoefs();
}

}

#endif
