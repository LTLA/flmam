#ifndef DISPATCHER_H
#define DISPATCHER_H

#include "LM_base.h"
#include "LM_oneway.h"
#include "LM_general.h"

namespace flmam {

inline std::unique_ptr<LM_base> dispatcher (SEXP qr, SEXP groups) {
    if (groups==R_NilValue) {
        return std::unique_ptr<LM_base>(new LM_general(qr));
    } else {
        return std::unique_ptr<LM_base>(new LM_oneway(groups));
    }
}

}

#endif
