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

extern int nb_sample_per_dim; // for range estimation by sampling: # of samples per dimension

extern int points_per_graph; // number of steps in time printed on each graphs (and in yaml files)
extern int printing_period; // deduced from nb of steps and points per graph

#define tol_noise 0.00001

extern int interactive_visualization; // 0 or 1
extern vector<bool> variables_to_display;

/*----- begin parameters for discrete systems -----*/
extern int nb_steps;  // nb of discrete time steps
extern int AEextension_order;
extern int iter_method;
extern bool skewing; // compute skewboxes or regular boxes approx
/*----- end parameters for discrete systems -----*/

extern const int LINESZ;
extern char sfx_filename[];
extern char onnx_filename[];




// for subdivisions of the initial domain to refine precision
extern int nb_subdiv_init; // number of subdivisiions
extern int component_to_subdiv, component_to_subdiv2;

extern double recovering; // percentage of recovering between subdivisions
extern vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xinner_print, Xinner_joint_print, Xinner_robust_print, Xexact_print; // store results of subdivision
extern vector<double> t_print; // times where results are stored
extern int current_subdiv;
extern int current_iteration;

// for robust inner-approximations
extern int uncontrolled;  // number of uncontrolled parameters (forall params)
extern int controlled;  // number of controlled parameters (forall params)
extern vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
//extern vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust inner-approx)
//extern int variable;  // number of non constant parameters
//extern vector<bool> is_variable; // for each parameter, constant or variable

extern vector<interval> target_set;
extern vector<interval> unsafe_set;

extern bool refined_mean_value;

extern bool print_debug;

extern bool recompute_control;




char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);

//void read_system(int argc, char* argv[]);

void readfromfile_syschoice(const char * params_filename, char* sfx_filename, char* onnx_filename);

void open_outputfiles();

extern ofstream approxreachsetfile;
extern YAML::Emitter out_approx;




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



class ReachSet
{
public:
    vector<interval> Xsampled;
    vector<interval> Xouter;        // maximal outer-approx
    vector<interval> Xouter_rob;    // robust outer-approx
    vector<interval> Xinner;        // maximal inner-approx
    vector<interval> Xinner_rob;    // robust inner-approx
    
    
    ReachSet() : Xsampled(sysdim), Xouter(sysdim), Xouter_rob(sysdim), Xinner(sysdim), Xinner_rob(sysdim) {}
    ReachSet(vector<interval> &Xs, vector<interval> &Xo, vector<interval> &Xor, vector<interval> &Xi, vector<interval> &Xir) : Xsampled(Xs), Xouter(Xo), Xouter_rob(Xor), Xinner(Xi), Xinner_rob(Xir) {}
};


class NQuant
{
public:
    vector<int> var_id;   // variable id (from 0 to jacdim-1)
    vector<bool> exists;  // true for existential quantification, false for universal
    
    int uncontrolled; // number of uncontrolled parameters (forall params)
    
    NQuant() {}
     NQuant(int n) : var_id(n), exists(n) {}
};

extern NQuant nquant;

#endif
