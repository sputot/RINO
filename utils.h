/* ============================================================================
 File   : utils.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 ============================================================================ */

#ifndef UTILS_H
#define UTILS_H

#include "ode_def.h"

extern vector<ofstream> outFile_outer_minimal;   //  minimal outer-approximated range for each variable of the system
extern vector<ofstream> outFile_outer;   // output outer-approximated range for each variable of the system
extern vector<ofstream> outFile_outer_robust;
extern vector<ofstream> outFile_inner_minimal;   //  minimal inner-approximated range for each variable of the system
extern vector<ofstream> outFile_inner;   // output inner-approximated range for each variable of the system
extern vector<ofstream> outFile_inner_robust;   // output inner-approximated range for each variable of the system
extern vector<ofstream> outFile_center;

extern ofstream outFile_width_ratio;     // min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
extern ofstream outFile_meanerror_outer; // mean on xi of error between outer-approx and analytical solution if any
extern ofstream outFile_meanerror_inner; // mean on xi of error between inner-approx and analytical solution if any
extern ofstream outFile_meanerror_diff;  // mean on xi of error between outer-approx and inner-approx
extern ofstream outFile_relmeanerror_outer; // mean on xi of error between outer-approx and analytical solution if any, over width of exact tube
extern ofstream outFile_relmeanerror_inner; // mean on xi of error between inner-approx and analytical solution if any, over width of exact tube
extern ofstream outFile_relmeanerror_diff;  // mean on xi of error between outer-approx and inner-approx, over width of over-approx tube

#define tol_noise 0.00001


#endif
