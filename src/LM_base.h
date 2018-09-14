#include "Rcpp.h"

namespace flmam {

class LM_base {
public:
    LM_base() = default;
    virtual ~LM_base() = default;
    LM_base(const LM_base&) = default;
    LM_base& operator=(const LM_base&) = default;
    LM_base(LM_base&&) = default;
    LM_base& operator=(LM_base&&) = default;

    virtual double fit(const double*, double*)=0;
    virtual int get_nobs() const=0;
    virtual int get_ncoefs() const=0;
protected:
    std::vector<double> workspace;
};

};
