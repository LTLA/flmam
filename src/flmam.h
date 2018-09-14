#ifndef FLMAM_H
#define FLMAM_H

#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "dispatcher.h"

extern "C" {

SEXP fit_lm (SEXP, SEXP, SEXP);

}

#endif
