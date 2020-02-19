/* ============================================================================
 File   : discrete system.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file computes ranges for functions and discrete dynamical systems
 ============================================================================*/


#ifndef DISCRETE_SYSTEM_H
#define DISCRETE_SYSTEM_H

#include "filib_interval.h"
#include "tadiff.h"
#include "fadiff.h"
#include "fadbad_aa.h"

class DiscreteFunc {
public:
    template <class C>
    vector<C> operator()(vector<C> x) {
        vector<C> z(x.size());
        z[0] = x[0]*x[0] - x[0];
        return z;
    }
};

void range_discrete_system(void);


#endif
