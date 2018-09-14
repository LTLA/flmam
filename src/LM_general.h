#ifndef LM_GENERAL_H
#define LM_GENERAL_H

#include "LM_base.h"
#include "QR_manager.h"

namespace flmam {

class LM_general : public LM_base {
public:
    LM_general(SEXP Qr, SEXP Aux) : qr(Qr, Aux) {
        workspace.resize(get_nobs());
    }

    // Rule of 5.
    ~LM_general() = default;
    LM_general(const LM_general&) = default;
    LM_general& operator=(const LM_general&) = default;
    LM_general(LM_general&&) = default;
    LM_general& operator=(LM_general&&) = default;

    // Subclassed virtual methods
    double fit (const double* values, double* out_coefs) {
        std::copy(values, values+get_nobs(), workspace.begin());
        double resid_var=fit0(workspace.data());
        std::copy(workspace.begin(), workspace.begin()+get_ncoefs(), out_coefs);
        return resid_var; 
    }
    
    int get_nobs() const{
        return qr.get_nobs();
    }
    
    int get_ncoefs() const{
        return qr.get_ncoefs();
    }

    double fit0(double* in, bool var_only=false) {
        qr.multiply(in);
        const int nobs=qr.get_nobs(), ncoefs=qr.get_ncoefs();
        
        double curvar=0;
        for (int i=ncoefs; i<nobs; ++i) {
            const auto& curval=in[i];
            curvar += curval * curval;
        }
        curvar /= nobs - ncoefs;
        
        if (!var_only) {
            qr.backsolve(in);
        }
        return curvar;
    }
private:
    QR_manager qr;
};

}

#endif
