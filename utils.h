/* ============================================================================
 File   : utils.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 ============================================================================ */

#ifndef UTILS_H
#define UTILS_H

#include "filib_interval.h"
#include "tadiff.h"
#include "fadiff.h"
#include "fadbad_aa.h"

#include "yaml-cpp/yaml.h"

//#include "ode_def.h"

void readfromfile_syschoice(const char * params_filename, char* sfx_filename, char* onnx_filename, int &nb_sample_per_dim);

void open_outputfiles();

extern ofstream approxreachsetfile;
extern YAML::Emitter out_approx;


#define tol_noise 0.00001

extern int interactive_visualization; // 0 or 1
extern vector<bool> variables_to_display;

// close XML and print final time: common to ODEs / discrete systems
void print_finalstats(clock_t begin);

// print initial conditions and init XML in the discrete systems case
void print_init_discrete(vector<interval> &x, bool skew);
// print initial conditions and init XML in the ODE case
void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs);

void run_pythonscript_visualization();

void print_ErrorMeasures(int current_iteration, double d0);
void print_finalsolution(int max_it, double d0);

std::ostream& operator<<(std::ostream& os, const std::vector<double> &input);
//std::ostream& operator<<(std::ostream& os, const std::vector<interval> &input);
//std::ostream& operator<<(std::ostream& os, const interval &input);
std::ostream& operator<<(std::ostream& os, const std::vector<AAF> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<T<AAF>> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<T<F<AAF>>> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<vector<double>> &input);
std::ostream& operator<<(std::ostream& os, const std::vector<vector<vector<double>>> &input);

#endif
