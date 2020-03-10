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
#include "ode_def.h"

class DiscreteFunc {
public:
    template <class C>
    vector<C> operator()(vector<C> x) {
        vector<C> z(sysdim);
        
        if (syschoice == 1)
            z[0] = x[0]*x[0] - x[0];
        else if (syschoice == 2)
            z[0] = x[1]*x[1] - 2.0*x[0];
        else if (syschoice == 3) // example 3.5 Goldztejn
        {
            z[0] = 81.0*x[0]*x[0] + x[1]*x[1] + 18.0*x[0]*x[1] - 100.0;
            z[1] = x[0]*x[0] + 81.0*x[1]*x[1] + 18.0*x[0]*x[1] - 100.0;
        }
        else if (syschoice == 4) // example 3 CDC
        {
            z[0] = 5.0*x[0]*x[0] + x[1]*x[1]  - 2.0*x[0]*x[1] - 4.0;
            z[1] = x[0]*x[0] + 5.0*x[1]*x[1] - 2.0*x[0]*x[1] - 4.0;
        }
        else if (syschoice == 5) { // example 5.1 Goldstzjen - surprising that square joint inner-approx is empty !?
            z[0] = x[0]*x[0]*x[0]*x[0]*x[0]*x[0] + x[1]*x[1]*x[1]*x[1]*x[1]*x[1] + x[0]*x[1] - 3.0;
            z[1] = x[0]*x[0]*x[0]*x[0]*x[0]*x[0] - x[1]*x[1]*x[1]*x[1]*x[1]*x[1] - x[0]*x[1] + 1.0;
        }
        else if (syschoice == 6) { //
            z[0] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1] + 2.0*x[0]*x[1] - 4.0;
            z[1] = x[0]*x[0]*x[0] - x[1]*x[1]*x[1] - 2.0*x[0]*x[1] + 2.0;
        }
        else if (syschoice == 7) { // skewed inner-approx is not very good - try to see with zonotope ?
            z[0] = 2.0*x[0]*x[0] + 2.0*x[1]*x[1] - 2.0*x[0]*x[1] - 2.0;
            z[1] = x[0]*x[0] - x[1]*x[1] + 4.0*x[0]*x[1] - 4.0;
        }
        else if (syschoice == 8) { // skewed inner-approx is empty - try to see with zonotope ?
            z[0] = 2.0*x[0]*x[0] + 1.0*x[1]*x[1] - 2.0*x[0]*x[1] - 1.0;
            z[1] = x[0]*x[0] - x[1]*x[1] + 3.0*x[0]*x[1] - 3.0;
        }
        else if (syschoice == 9) { //
            z[0] = 2.0*x[0]*x[0] - x[0]*x[1] - 1.0;
            z[1] = x[0]*x[0] + x[1]*x[1]  - 2.0;
        }
        else if (syschoice == 10) { // (z0,z1) empty ???
            z[0] = 2.0*x[0]*x[0] - x[0]*x[1] + x[0]*x[2] + x[2]*x[2] - 3.0;
            z[1] = x[0]*x[0] + x[1]*x[1] - x[2]*x[2] - 1.0;
            z[2] = x[0] + x[1] - 2.0;
        }
        else if (syschoice == 11) { // ok
            z[0] = 2.0*x[0]*x[0] - x[0]*x[1] - 1.0;
            z[1] = x[0]*x[0] + x[1]*x[1]  - 2.0;
            z[2] = x[2] - 1.0;
        }
        else if (syschoice == 12) { // ok
            z[0] = x[0];
            z[1] = 0.707*x[1] + 0.707*x[1]*x[1]  - 2.0;
            z[2] = x[2] - 1.0;
        }
        return z;
    }
};


extern vector<vector<vector<interval>>> constr_eps;  // constraints on noise symbols
extern vector<vector<vector<interval>>> eps_loc;     // consequence on eps=[x]-x0

void range_discrete_system(void);

void evaluate_projections(vector<interval> &z0,  vector<vector<interval>> &Jacf);
void evaluate_projections_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int index1, int n1, int index2, int n2);
void evaluate_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, vector<bool> &is_existential);

// for only one component: used for joint range
interval evaluate_innerrange_x_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, vector<int> &exist_quantified, int i, int index1, int index2);
interval evaluate_innerrange_x(vector<interval> &z0,  vector<vector<interval>> &Jacf, vector<bool> &is_existential, int i);

void joint_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, int varx, int vary);
void joint_ranges_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int varx, int vary);

void preconditioned_joint_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, int varx, int vary);
void preconditioned_joint_ranges_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAfff, int varx, int vary);
#endif
