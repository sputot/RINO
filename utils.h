/* ============================================================================
 File   : utils.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 ============================================================================ */

#ifndef UTILS_H
#define UTILS_H

//#include "ode_def.h"

void open_outputfiles();

extern vector<ofstream> outFile_outer_minimal;   //  minimal outer-approximated range for each variable of the system
extern vector<ofstream> outFile_outer;   // output outer-approximated range for each variable of the system
extern vector<ofstream> outFile_exact;
extern vector<ofstream> outFile_outer_robust;
extern vector<ofstream> outFile_inner_minimal;   //  minimal inner-approximated range for each variable of the system
extern vector<ofstream> outFile_inner;   // output inner-approximated range for each variable of the system
extern vector<vector<ofstream>> outFile_joint_inner;   // output inner-approximated range for each couple of variables of the system
extern vector<vector<vector<ofstream>>> outFile_joint_inner3d;   // output inner-approximated range for each triple of variables of the system
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

extern int interactive_visualization; // 0 or 1
extern vector<bool> variables_to_display;

void print_finalstats(clock_t begin);

void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs);

void run_pythonscript_visualization();

void print_ErrorMeasures(int current_iteration, double d0);
void print_finalsolution(int max_it, double d0);

std::ostream& operator<<(std::ostream& os, const std::vector<double> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<AAF> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<vector<double>> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<vector<vector<double>>> &input);

#endif
