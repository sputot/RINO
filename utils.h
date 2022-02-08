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

extern int systype;    // systype = 0 (ODE) or 1 (DDE) -- initialized in main.cpp / command line
extern int syschoice;  // to choose among the predefined systems of ODE or DDE -- initialized in main.cpp / command line

extern bool nn_analysis; // whether there is a network in the loop or not

extern int sysdim; // dimension of the system of ODE/DDE to analyze
extern int inputsdim; // dimension of the uncertain inputs and parameters of the system
extern int fullinputsdim; // full dimension of the uncertain inputs and parameters of the system: taking into account variable inputs
extern int jacdim; // Jacobian will be dimension jacdim = sysdim + fullinputsdim
extern int paramsdim;  // dimension of the vector of parameters params that do not appear in Jacobian
extern int nncontroldim;  // dimension of the neural network control - does not appear in Jacobian

extern bool create_png; // whether or not the python script is called to create the .png result files

class ReachSet
{
public:
    vector<interval> Xsampled;
    vector<interval> Xouter;
    vector<interval> Xinner;
    
    ReachSet() : Xsampled(sysdim), Xouter(sysdim), Xinner(sysdim) {}
    ReachSet(vector<interval> &Xs, vector<interval> &Xo, vector<interval> &Xi) : Xsampled(Xs), Xouter(Xo), Xinner(Xi) {}  // ou recopie ?
};


void readfromfile_syschoice(const char * params_filename, char* sfx_filename, char* onnx_filename, int &nb_sample_per_dim);

void open_outputfiles();

extern ofstream approxreachsetfile;
extern YAML::Emitter out_approx;

extern int points_per_graph; // number of steps in time printed on each graphs (and in yaml files)
extern int printing_period; // deduced from nb of steps and points per graph

#define tol_noise 0.00001

extern int interactive_visualization; // 0 or 1
extern vector<bool> variables_to_display;


// print initial conditions and init XML in the discrete systems case
void print_init_discrete(vector<interval> &x, bool skew);
// print initial conditions and init XML in the ODE case
void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs);

void run_pythonscript_visualization();

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
