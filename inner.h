/* ============================================================================
 File   : inner.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file defines operations specific to Kaucher arithmetic for computing inner ranges
 ============================================================================ */

#ifndef INNER_H
#define INNER_H

#include "filib_interval.h"
#include "fadbad_aa.h"
#include "matrix.h"

// Kaucher multiplication Jacobian * (z-z0) where we consider the improper dual(eps=(z-z0))
interval Kaucher_multeps(interval Jac, interval eps);
// Kaucher addition where we consider pro as proper and impro as improper and we want an improper result
interval Kaucher_add_pro_impro(interval pro, interval impro);
// Kaucher addition where we consider pro as proper and impro as improper and we want an proper result
interval Kaucher_add_pro_impro_resultpro(interval pro, interval impro);
// computer inner and outer-approx by mean-value theorem
// using special addition corresponding to addition of proper and improper intervals: the result is an inner range for the solution of ODE system
void InnerOuter(vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xinner_minimal, vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xouter_minimal, vector<AAF> &x0p1, vector<vector<AAF>> &Jtau, vector<interval> &eps, double tnp1);

#endif
