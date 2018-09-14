#ifndef LM_ONEWAY_H
#define LM_ONEWAY_H

#include "LM_base.h"

namespace flmam {

class LM_oneway : public LM_base {
public:
    LM_oneway(Rcpp::IntegerVector ids) : assignments(ids), ngroups(0), resid_df(0) {
        if (!assignments.size()) {
            return;
        }

        // Checking the validity of all indices.
        for (auto a : assignments) {
            if (a<0 || a==NA_INTEGER) {
                throw std::runtime_error("group IDs must be non-negative integers");
            }
        }

        // Counting the number in each group. 
        ngroups=*std::max_element(assignments.begin(), assignments.end()) + 1;
        npergroup.resize(ngroups);
        for (auto a : assignments) {
            ++npergroup[a];
        }
        
        for (auto n : npergroup) {
            if (n==0) {
                throw std::runtime_error("group IDs must be consecutive integers");
            }
        }
           
        resid_df=std::accumulate(npergroup.begin(), npergroup.end(), -static_cast<int>(ngroups));
        workspace.resize(ngroups);
        return;
    }

    // Rule of 5.
    ~LM_oneway() = default;
    LM_oneway(const LM_oneway&) = default;
    LM_oneway& operator=(const LM_oneway&) = default;
    LM_oneway(LM_oneway&&) = default;
    LM_oneway& operator=(LM_oneway&&) = default;

    // Subclassed virtual methods.
    double fit(const double* values, double* out_means) {
        fit0(values, out_means, workspace.data(), false, true);

        if (get_nobs()==ngroups) {
            return R_NaN;
        } else {
            return std::accumulate(workspace.begin(), workspace.end(), 0.0) / (get_nobs() - ngroups);
        }
    }

    int get_nobs() const {
        return assignments.size();
    }
    
    int get_ncoefs() const {
        return ngroups;    
    }

    // Special methods.
    void fit0(const double* values, double* out_means, double* out_vars, bool mean_only=false, bool return_rss=false) {
        std::fill(out_means, out_means+npergroup.size(), 0);

        // Computing the mean for each group.
        auto copy=values;
        for (auto a : assignments) {
            const double& val=*(copy++);
            out_means[a]+=val;
        }
        for (size_t i=0; i<ngroups; ++i) {
            out_means[i]/=npergroup[i];
        }
    
        if (!mean_only) {
            std::fill(out_vars, out_vars+npergroup.size(), 0);

            // Computing the RSS for each group (numerically stable).
            for (auto a : assignments) {
                const double& val = *(values++);
                const double tmp = val - out_means[a];
                out_vars[a] += tmp * tmp;
            }
   
            if (!return_rss) {
                for (size_t i=0; i<npergroup.size(); ++i) {
                    if (npergroup[i]>1) { 
                        out_vars[i]/=npergroup[i] - 1;
                    } else {
                        out_vars[i]=R_NaN;
                    }
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
    size_t ngroups;
    std::vector<int> npergroup;
    int resid_df;
};

}

#endif
