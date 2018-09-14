#ifndef LM_ONEWAY_H
#define LM_ONEWAY_H

#include "LM_base.h"

namespace flmam {

class LM_oneway : public LM_base {
public:
    LM_oneway(SEXP ids) : assignments(ids), ngroups(*std::max_element(assignments.begin(), assignments.end()) + 1), resid_df(0) {
        npergroup.resize(ngroups);
        for (auto a : assignments) {
            if (a<0 || a==NA_INTEGER) {
                throw std::runtime_error("group IDs must be non-negative integers");
            }
            ++npergroup[a];
        }
        
        for (auto n : npergroup) {
            if (n==0) {
                throw std::runtime_error("group IDs must be consecutive integers");
            }
        }
        
        resid_df=std::accumulate(npergroup.begin(), npergroup.end(), -static_cast<int>(ngroups));
        workspace.resize(ngroups);
    }

    // Rule of 5.
    ~LM_oneway() = default;
    LM_oneway(const LM_oneway&) = default;
    LM_oneway& operator=(const LM_oneway&) = default;
    LM_oneway(LM_oneway&&) = default;
    LM_oneway& operator=(LM_oneway&&) = default;

    // Subclassed virtual methods.
    double fit(const double* values, double* out_means) {
        fit0(values, out_means, workspace.data());
        
        double output=0;
        for (size_t i=0; i<ngroups; ++i) {
            if (npergroup[i]>1) {
                output+=workspace[i]/(npergroup[i]-1);
            }
        }
        
        return output;
    }

    int get_nobs() const {
        return assignments.size();
    }
    
    int get_ncoefs() const {
        return ngroups;    
    }

    // Special methods.
    void fit0(const double* values, double* out_means, double* out_vars, bool mean_only=false) {
        // Computing the mean for each group.
        auto copy=values;
        std::fill(out_means, out_means+npergroup.size(), 0);
    
        for (auto a : assignments) {
            const double& val=*(copy++);
            out_means[a]+=val;
        }
        for (size_t i=0; i<npergroup.size(); ++i) {
            out_means[i]/=npergroup[i];
        }
    
        if (!mean_only) {
            // Computing the RSS for each group (numerically stable).
            std::fill(out_vars, out_vars+npergroup.size(), 0);
            for (auto a : assignments) {
                const double& val = *(values++);
                const double tmp = val - out_means[a];
                out_vars[a] += tmp * tmp;
            }
    
            for (size_t i=0; i<npergroup.size(); ++i) {
                if (npergroup[i]>1) { 
                    out_vars[i]/=npergroup[i] - 1;
                } else {
                    out_vars[i]=R_NaReal;
                }
            }
        }
    
        return;
    }

    const std::vector<int>& get_group_sizes() const {
        return npergroup;
    }
private:
    Rcpp::IntegerVector assignments;
    const size_t ngroups;
    std::vector<int> npergroup;
    int resid_df;
};

}

#endif
