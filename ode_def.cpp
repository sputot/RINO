
/* ============================================================================
 File   : ode_def.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
The place to define initial conditions and parameters the systems of ODEs or DDEs on which to perform reachability
============================================================================ */

#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"
#include "ode_def.h"
#include "matrix.h"

using namespace std;

int sysdim; // dimension of system of ODE/DDE
int inputsdim; // dimension of the uncertain inputs and parameters of the system
int fullinputsdim; // full dimension of the uncertain inputs and parameters of the system: taking into account variable inputs
int jacdim;  //  Jacobian will be dimension sysdim * jacdim, for ODEs jacdim = sysdim + fullinputsdim
int sysdim_params;

double t_end; // ending time of integration
double t_begin; // starting time of initialization
double control_period = 0.0;

// parameters of the system of the ODEs
vector<AAF> params;  // params of the ODE that don't appear in the Jcaobian such as control given by NN output (or nondeterministic disturbances ?)
vector<vector<AAF>> Jac_params;   // (\partial u) / (partial x)
vector<vector<AAF>> Jac_params_order2;   // (\partial u) / (partial x)

vector<AAF> initial_values;
vector<AAF> center_initial_values;

vector<AAF> inputs; // uncertain inputs and parameters
vector<AAF> fullinputs; // uncertain inputs and parameters
vector<int> nb_inputs; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
vector<AAF> center_fullinputs;
vector<int> index_param;
vector<int> index_param_inv;

vector<interval> eps;

vector<vector<interval>> Jac_param_inputs; // for inputs defined as g(x1,...xn): we give the jacobian

//vector<F<AAF>> nn_outputs; // result of NN evaluation

// for subdivisions of the initial domain to refine precision
int nb_subdiv_init = 1; // number of subdivisiions
int component_to_subdiv = -1;
int component_to_subdiv2 = -1;

double recovering = 0.0; // percentage of recovering between subdivisions
vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print, Xinner_joint_print, Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
vector<double> t_print; // times where results are stored
int current_subdiv;
int current_iteration;

// for robust inner-approximations
int uncontrolled; // number of uncontrolled parameters (forall params)
int controlled; // number of controlled parameters (forall params)
vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
//vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust outer-approx)
//int variable; // number of non constant parameters
//vector<bool> is_variable;  // for each parameter, constant or variable

vector<interval> target_set;
vector<interval> unsafe_set;

bool refined_mean_value;

bool print_debug = true;

bool recompute_control = true;

// define the dimensions of your system (ODE or DDE) and if we want initial subdivisions
void define_system_dim(const char * params_filename)
{
    /*************************************************************************** ODE ************************************************************/
    
    sysdim_params = 0;
    inputsdim = 0;
    nb_subdiv_init = 1; // nb of initial subdivisions of the input range

    if (systype == 0) // ODE
    {
      
        if (syschoice == 1)  // running example
        {
            sysdim = 1;
            sysdim_params = 1;
        }
        else if (syschoice == 2)  // Brusselator
        {
            sysdim = 2;
            sysdim_params = 2;
        }
        else if (syschoice == 3)  // ballistic
        {
            sysdim = 4;
        }
        else if (syschoice == 4)  // ballistic linearise + masse incertaine
        {
            sysdim = 4;
            inputsdim = 1;
        }
        else if (syschoice == 5)  // self-driving car
        {
            sysdim = 2;
            sysdim_params = 2;
        }
        else if (syschoice == 6)  //  self-driving car
        {
            sysdim = 2;
            inputsdim = 2;
            //     sysdim_params = 2;
        }
        else if (syschoice == 7)  //  self-driving car
        {
            sysdim = 4;
            //     sysdim_params = 2;
        }
        else if (syschoice == 8)
        {
            sysdim = 1;
        }
        else if (syschoice == 9)  // Acrobatic quadcopter
        {
            sysdim = 6;
            inputsdim = 2;
        }
        else if (syschoice == 99)  // Acrobatic quadcopter  with m et Iyy as disturbances
        {
            sysdim = 6;
            inputsdim = 4;
        }
        else if (syschoice == 10)  // 10-D near-hover quadrotor
        {
            sysdim = 10;
            inputsdim = 6;
        }
        else if (syschoice == 11)  // Dubbins vehicle
        {
            sysdim = 3;
            inputsdim = 4;
        }
        else if (syschoice == 12)  // academic example to investigate time-varying parameters
        {
            sysdim = 2;
        }
        else if (syschoice == 13)  // Laub-Loomis Benchmark [Arch 2019]
        {
            sysdim = 7;
            sysdim_params = 1;
        }
        else if (syschoice == 14) // Van der Pol oscillator [Arch 2019]
        {
            sysdim = 2;
        }
        else if (syschoice == 15) // Van der Pol oscillator [Arch 2018 and Sparse Polynomial zonotopes]
        {
            sysdim = 2;
        }
        else if(syschoice == 17) // quadrotor model [Arch 2019]
        {
            sysdim = 12;
        }
        else if (syschoice == 18) // HSCC 2019 paper crazyflie example
        {
            sysdim = 14;
            // sysdim_params = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 181) // HSCC 2019 paper crazyflie example with neural network controller
        {
            sysdim = 14;
            inputsdim = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 182) // HSCC 2019 paper crazyflie example with neural network controller and agressive altitude PID
        {
            sysdim = 14;
            inputsdim = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 183) // HSCC 2019 paper crazyflie example with neural network controller and agressive altitude PID
        {
            sysdim = 14;
            inputsdim = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 19) {  // academic example, time-varying (piecewise constant) parameters
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 21) {  // academic example, time-varying (piecewise constant) parameters
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 22) {  // academic example, time-varying (piecewise constant) parameters
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 23) {   // pursuer-evader example Mitchell
            sysdim = 3;
            inputsdim = 2;
        }
        else if (syschoice == 24) { // [Franzle et al.]
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 25) { // [Franzle et al.] reversed time van der pol oscillator with uncertainty
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 26) { // [Franzle et al.] 7-d biological system
            sysdim = 7;
            inputsdim = 1;
        }
        else if (syschoice == 27) { // [Franzle et al.] 7-d biological system but with sharing
            sysdim = 7;
            inputsdim = 1;
        }
        else if (syschoice == 28) { // [Franzle et al.] 7-d biological system but with sharing
            sysdim = 7;
            inputsdim = 0;
        }
        else if (syschoice == 29) { // EX_10 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 4;
            inputsdim = 2;
        }
        else if (syschoice == 30) { // EX_1 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 301) { // EX_1 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 2;
            inputsdim = 0;
        }
        else if (syschoice == 31) { // Quadcopter Mikhail Bessa
            sysdim = 14;
            inputsdim = 3;
        }
        else if (syschoice == 311) { // Quadcopter Mikhail Bessa avec 3 composantes en plus
            sysdim = 11;
            inputsdim = 3;
        }
        else if (syschoice == 32) { // EX_2 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 33) { // EX_3 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 34) { // EX_4 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            inputsdim = 1;
        }
        else if (syschoice == 35) { // EX_5 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            inputsdim = 1;
        }
        else if (syschoice == 36) { // EX_6 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            inputsdim = 1;
        }
        else if (syschoice == 37) { // EX_7 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            inputsdim = 1;
        }
        else if (syschoice == 38) { // EX_8 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 4;
            inputsdim = 1;
        }
        else if (syschoice == 381) { // EX_8 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 4;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 382) {
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 383) { // EX_2 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 2;                // EX13  sherlock/systems_with_networks
            sysdim_params = 1;
        }
        else if (syschoice == 384) { // EX_3 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 2;              // EX14  sherlock/systems_with_networks
            sysdim_params = 1;      // EX15  sherlock/systems_with_networks
        }
        else if (syschoice == 385) { // EX_4 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            sysdim_params = 1;
        }
        else if (syschoice == 386) { // EX_5 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            sysdim_params = 1;
        }
        else if (syschoice == 387) { // EX_6 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            sysdim_params = 1;
        }
        else if (syschoice == 388) { // EX_7 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 3;
            sysdim_params = 1;
        }
        else if (syschoice == 39) { // EX_9 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 4;
            inputsdim = 1;
        }
        else if (syschoice == 40) { // EX_10 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            sysdim = 4;
            inputsdim = 2;
        }
        else if (syschoice == 41) { // essai
            sysdim = 2;
            inputsdim = 1;
        }
        else if (syschoice == 42) // HSCC 2019 paper crazyflie example+aerodynamical effects
        {
            sysdim = 14;
            // sysdim_params = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 43) // HSCC 2019 paper crazyflie example - new PID, no aerodynamical effects
        {
            sysdim = 14;
            // sysdim_params = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 44) // HSCC 2019 paper crazyflie example+ new PID, with aerodynamical effects
        {
            sysdim = 14;
            // sysdim_params = 3;
            // 0 for sysdim params
        }
        else if (syschoice == 45) { // Mountain car Verisig
            sysdim = 2;
            inputsdim = 0;
        }
        else if (syschoice == 451) { // Mountain car Verisig
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 46) {  // Tora Heterogeneous ARCH-COMP 2020 - NNV (avec RNN format sfx obtenu a partir du .mat),
            sysdim = 4;
            inputsdim = 0;
        }
        else if (syschoice == 461) {  // Tora Heterogeneous ARCH-COMP 2020 - NNV (avec RNN format sfx obtenu a partir du .mat),
            sysdim = 4;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 471) {  // // Ex 1 ReachNNstar (avec RNN format sfx obtenu a partir du .txt),
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 1111) {  // toy example
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 1113) {  // toy example
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 0;
        }
        else if (syschoice == 481) {  // // Ex 2 ReachNNstar (avec RNN format sfx obtenu a partir du .txt),
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 482) {  // // Ex 3 ReachNNstar (avec RNN format sfx obtenu a partir du .txt),
            sysdim = 2;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 483) {  // // Ex 4 ReachNNstar (avec RNN format sfx obtenu a partir du .txt),
            sysdim = 3;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 484|| syschoice == 4841) {  // // Ex 5 ReachNNstar (avec RNN format sfx obtenu a partir du .txt),
            sysdim = 3;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 491) // Ex ACC de Verisig (avec nn obtenu a partir du yaml)
        {
            sysdim = 6;
            inputsdim = 0;
            sysdim_params = 1;
        }
        else if (syschoice == 493) // Ex QMPC de Verisig (avec nn obtenu a partir du yaml)
        {
            sysdim = 6;
            inputsdim = 0;
            sysdim_params = 3;
        }
    }
    /*************************************************************************** DDE ************************************************************/
    else if (systype == 1) // DDE
    {
        
        if (syschoice == 1)  // running example
        {
            sysdim = 1;
      //      nb_subdiv_init = 2;
      //      nb_subdiv_init = 10;
      //      recovering = 0.1; // recovering between 2 subdivisions when subdividing initial parameters
            
        }
        if (syschoice == 2)  //
        {
            sysdim = 2;
        }
        else if (syschoice == 3) // Xue 2017 (Ex 3)
        {
            sysdim = 7;
        }
        else if (syschoice == 4) // Szczelina 1 and 2 2014
        {
            sysdim = 1;
        }
        else if (syschoice == 5) // Szczelina 2 2014
        {
            sysdim = 1;
        }
        else if (syschoice == 6) // self-driving car
        {
            sysdim = 2;
            sysdim_params = 2;
        }
        else if (syschoice == 8) // self-driving car but with coeff in interv
        {
            sysdim = 2;
            inputsdim = 2;
        }
        else if (syschoice == 9) 
        {
            sysdim = 1;
        }
        else if (syschoice == 10) // platoon of 5 vehicles
        {
            sysdim = 9;
        }
        else if (syschoice == 11) // platoon of 7 vehicles
        {
            sysdim = 19;
        }
    }
    
    if (params_filename) // called with configuration file: we overwrite the initialization of init_system
        readfromfile_nbsubdiv(params_filename, nb_subdiv_init);
    
}


// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &param_inputs, vector<AAF> &param_inputs_center, vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J)
{
    assert(systype == 0); // ODE
   
    // par défaut
    for (int i=0 ; i<sysdim ; i++) {
        x[i] = initial_values[i];
        xcenter[i] = center_initial_values[i];
    }
    
    // param_inputs between 0 and jacdim-sysdim, inputs between 0 and inputsdim
        for (int j = 0; j < fullinputsdim ; j++) {
            
            param_inputs[j] = fullinputs[j]; // inputs[sysdim+i].convert_int(); // create a new independent AAF ?
            if (innerapprox == 1)
                param_inputs_center[j] = center_fullinputs[j]; // center_inputs[sysdim+i].convert_int();
            else
                param_inputs_center[j] = fullinputs[j]; // inputs[sysdim+i].convert_int();
        }
    
    setId(J);
    
}

void readfromfile_nbsubdiv(const char * params_filename, int &nb_subdiv_init)
{
    const int LINESZ = 2048;
    char buff[LINESZ];
    
    cout << "****** Reading system dimensions from file " <<  params_filename << " ******" << endl;
    FILE *params_file = fopen(params_filename,"r");
    if (params_file == NULL)
        cout << "Error reading " << params_filename << ": file not found" << endl;
    while (fgets(buff,LINESZ,params_file)) {
//        sscanf(buff, "system-dimension = %d\n", &sysdim);
//        sscanf(buff, "inputs-dimension = %d\n", &inputsdim);
//        sscanf(buff, "sys-parameters-dimension = %d\n", &sysdim_params);
        sscanf(buff, "nb-initial-subdivisions = %d\n", &nb_subdiv_init);
    }
    fclose(params_file);
}

// d0 and t_begin are for DDEs only, rest are common to ODE and DDE
void read_parameters(const char * params_filename, double &tau, double &t_end, double &d0, double &t_begin, int &order, int &nb_subdiv)
{
    const int LINESZ = 2048;
    char buff[LINESZ];
    char initialcondition[LINESZ];
    const char space[2] = " ";
    double a, b;
    int c;
    
    cout << "****** Reading system parameter from file " <<  params_filename << " ******" << endl;
    FILE *params_file = fopen(params_filename,"r");
    if (params_file == NULL)
        cout << "Error reading " << params_filename << ": file not found" << endl;
    
    
    
    while (fgets(buff,LINESZ,params_file)) {
        //     sscanf(buff, "system = %s\n", sys_name);
        //      sscanf(buff, "initially = %[^\n]\n", initial_condition);   // tell separator is newline, otherwise by default it is white space
        sscanf(buff, "time-horizon = %lf\n", &t_end);
        sscanf(buff, "integration-step = %lf\n", &tau);
        sscanf(buff, "control-step = %lf\n", &control_period);
        //      sscanf(buff, "output-variables = %[^\n]\n", output_variables);
        sscanf(buff, "delay = %lf\n", &d0);               // for DDEs
        sscanf(buff, "starting-time = %lf\n", &t_begin);  // for DDEs
        sscanf(buff, "nb-time-subdivisions = %d\n", &nb_subdiv); // for DDEs : subdiv of time interval d0 : tau is deduced
        sscanf(buff, "order = %d\n", &order);
        sscanf(buff, "interactive-visualization = %d\n", &interactive_visualization);
        sscanf(buff, "refined-mean-value = %d\n", &refined_mean_value);
        if (sscanf(buff, "variables-to-display = %s\n", initialcondition) == 1)
        {
            for (int i=0; i< sysdim; i++)
                variables_to_display[i] = false;
            
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i;
            while( token != NULL ) {
                sscanf(token,"%d",&i);
                variables_to_display[i-1] = true;
             //   cout <<"input="<<inputs[i].convert_int()<<endl;
                token = strtok(NULL,space);
            }
        }
        if (sscanf(buff, "initial-values = %s\n", initialcondition) == 1)
        {
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i = 0;
            nb_subdiv_init = 1;
            component_to_subdiv = -1;
            component_to_subdiv2 = -1;
            // only one component can be subdivided for now (we keep the last if several...)
            while( token != NULL ) {
                if (sscanf(token,"([%lf,%lf],%d)",&a,&b,&c) == 3) {
                    if (component_to_subdiv > -1) // we already have one component to subdivide
                    {
                        component_to_subdiv2 = i;
                        cout << "component "<< i << " should also be subdivided " << endl;
                    }
                    else
                    {
                    nb_subdiv_init = c;
                    component_to_subdiv = i;
                    cout << "component "<< i << " should be dubdivided " << c << "times" << endl;
                    }
                }
                else
                    sscanf(token,"[%lf,%lf]",&a,&b);
                initial_values[i] = interval(a,b);
                cout <<"intial_value="<<interval(a,b)<<endl;
                i++;
                token = strtok(NULL,space);
            }
        }
        if (sscanf(buff, "inputs = %s\n", initialcondition) == 1)
        {
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i = 0;
            while( token != NULL ) {
                if (sscanf(token,"([%lf,%lf],%d)",&a,&b,&c) == 3)
                    nb_inputs[i] = c;
                else {
                    nb_inputs[i] = 1;
                    sscanf(token,"[%lf,%lf]",&a,&b);
                }
                inputs[i] = interval(a,b);
           //     cout <<"input="<<interval(a,b)<<endl;
           //     cout <<"nb_inputs["<<i<<"]="<<nb_inputs[i]<<endl;
                i++;
                token = strtok(NULL,space);
            }
        }
        if (sscanf(buff, "params = %s\n", initialcondition) == 1)
        {
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i = 0;
            while( token != NULL ) {
                sscanf(token,"[%lf,%lf]",&a,&b);
                params[i] = interval(a,b);
                  //       cout <<"params="<<params[i].convert_int()<<endl;
                i++;
                token = strtok(NULL,space);
            }
        }
        if (sscanf(buff, "uncontrolled = %s\n", initialcondition) == 1)
        {
            for (int i=0 ; i<inputsdim; i++)
                is_uncontrolled[i] = false;
                
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i;
            while( token != NULL ) {
                sscanf(token,"%d",&i);
                is_uncontrolled[i-1] = true;
          //      cout <<"is_uncontrolled="<<i<<endl;
                token = strtok(NULL,space);
            }
        }
     
    }
    fclose(params_file);
        
    // cout << "system name = " << sys_name << endl;
      cout << "****** End parameter reading ******" << endl << endl;
}



// maps the ODE sys coordinates to the neural network inputs when they do not correspond
//vector<F<AAF>> syst_to_nn(vector<F<AAF>> &sysval){
    template <class C> vector<C> syst_to_nn(vector<C> &sysval) {
    vector<C> res; //= vector<AAF>(NH.n_inputs);
                                  
    if ((syschoice == 491) || (syschoice == 492)) // Ex ACC de Verisig (avec nn obtenu a partir du yaml)
     {
         res = vector<C>(NH.n_inputs);
         res[0] = 30.0;
         res[1] = 1.4;
         res[2] = sysval[4]; // v_ego
         res[3] = sysval[0]-sysval[3]; // x_lead-x_ego
         res[4] = sysval[1]-sysval[4]; // v_lead-v_ego
    }
    else
        res = sysval;
    return res;
}

// maps the neural network output to the control commands when they do not correspond
//vector<AAF> nn_to_control(vector<AAF> &nnoutput)
template <class C> vector<C> nn_to_control(vector<C> nnoutput) {
    
    if (syschoice == 493) // QMPC
    {
        //vector<AAF> res = vector<AAF>(sysdim_params);
        vector<double> actions(NH.n_outputs);
        for (int i=0 ; i<NH.n_outputs ; i++)
            actions[i] = sup(nnoutput[i].convert_int());
        // unsound implementation for now: should do the union over all argmax for the interval outputs of the NN...
        int index = argmax(actions.begin(),actions.end());
        cout << "index=" << index;
        if (index == 0)
            return {interval(-0.1,-0.1), interval(-0.1,-0.1), interval(7.81,7.81)};
        else if (index == 1)
            return {interval(-0.1,-0.1), interval(-0.1,-0.1), interval(11.81,11.81)};
        else if (index == 2)
            return {interval(-0.1,-0.1), interval(0.1,0.1), interval(7.81,7.81)};
        else if (index == 3)
            return {interval(-0.1,-0.1), interval(0.1,0.1), interval(11.81,11.81)};
        else if (index == 4)
            return {interval(0.1,0.1), interval(-0.1,-0.1), interval(7.81,7.81)};
        else if (index == 5)
            return {interval(0.1,0.1), interval(-0.1,-0.1), interval(11.81,11.81)};
        else if (index == 6)
            return {interval(0.1,0.1), interval(0.1,0.1), interval(7.81,7.81)};
        else
            return {interval(0.1,0.1), interval(0.1,0.1), interval(11.81,11.81)};
    }
    else
        return nnoutput;
}






// the main function to define the system
// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(const char * params_filename, double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order /*, vector<interval> &ix*/)
{
    define_system_dim(params_filename); // defines value of sysdim: depends on syschoice -- reads from file if input at command-line
    
    // inputs
    cout << "inputsdim=" << inputsdim << "sysdim=" << sysdim << endl;
    
    inputs = vector<AAF>(inputsdim);      // bounds
    nb_inputs = vector<int>(inputsdim);   // number of instances for each input
    for (int i=0; i<inputsdim; i++)
        nb_inputs[i] = 1;
    uncontrolled = 0;
    controlled = 0;
    is_uncontrolled = vector<bool>(inputsdim);
    for (int i=0 ; i<inputsdim; i++)
        is_uncontrolled[i] = false;  // controlled input or parameter
    
    target_set = vector<interval>(sysdim);
    unsafe_set = vector<interval>(sysdim);
    
    // initial values
    initial_values = vector<AAF>(sysdim);
    
    // parameters not part of the jacobian
    if (sysdim_params > 0) {
        params = vector<AAF>(sysdim_params);
        Jac_params = vector<vector<AAF>>(sysdim, vector<AAF>(sysdim+inputsdim));  // should probably be sysdim \times jacdim but jacdim not yet defined ?
        Jac_params_order2 = vector<vector<AAF>>(sysdim, vector<AAF>(sysdim+inputsdim));  // should probably be sysdim \times jacdim but jacdim not yet defined ?
    }

    refined_mean_value = false;
    
    if (systype == 0) // ODE
    {
        t_begin = 0;
        if (syschoice == 1)  // running example
        {
            tau = 0.01;
            t_end = 2.;
            order = 3;
            
            initial_values[0] = interval(0.9,1);
            params[0] = 1.0;
         //   nb_subdiv_init = 2;
         //   component_to_subdiv = 0;
        }
        else if (syschoice == 2) // Brusselator
        {
            tau = 0.05;
            t_end = 10.;
            order = 4;
            
            initial_values[0] = interval(0.9,1);
            initial_values[1] = interval(0,0.1);
            
            params[0] = 1;
            params[1] = 1.5;
        }
        else if (syschoice == 3) // ballistic
        {
            tau = 0.1;
            t_end = 4.;
            order = 3;
            
            initial_values[0] = interval(181.,185.); // velocity // interval(175.0,190.0); pour Eric
            // ix[1] = 3.14159/180*interval(2.5,3.5);  // angle   // interval(0,5) pour Eric
            //  ix[1] = mid(3.14159/180*interval(2.5,3.5)) + interval(-0.00872664, -0.00497644); // almost fault trajectories
            initial_values[1] = interval(0.0436,0.0611); // 3.14159/180*interval(2.5,3.5); // mid(3.14159/180*interval(2.5,3.5)) + interval( -0.00497644,0.00872664); // complement = safe trajectories
            initial_values[2] = interval(0.0,0.01);
            initial_values[3] = interval(0.0,0.01);
        }
        else if (syschoice == 4) // ballistic linearise
        {
            tau = 0.1;
            t_end = 1.4;
            order = 3;
            
            initial_values[0] = interval(181.,185.); // velocity // interval(175.0,190.0); pour Eric
            // ix[1] = 3.14159/180*interval(2.5,3.5);  // angle   // interval(0,5) pour Eric
            //  ix[1] = mid(3.14159/180*interval(2.5,3.5)) + interval(-0.00872664, -0.00497644); // almost fault trajectories
            initial_values[1] =  interval(0.0436,0.0611); // 3.14159/180*interval(2.5,3.5); // mid(3.14159/180*interval(2.5,3.5)); // + interval( -0.00497644,0.00872664); // complement = safe trajectories
            initial_values[2] = interval(0.0,0.01);
            initial_values[3] = interval(0.0,0.25); // interval(0.0,0.01);
            inputs[0]= interval(11.,15.); // 14.... la masse (incontrollable)
            is_uncontrolled[0] = true;
  //          is_variable[4] = true;
            nb_subdiv_init = 1;
            component_to_subdiv = 3;
            component_to_subdiv2 = 4;
        }
        else if (syschoice == 5) // self-driving car; sysdim = 2, jacdim = 2
        {
            tau = 0.05;
            t_end = 5.;
            order = 3;  // order of Taylor Models
            // initial condition
            initial_values[0] = interval(-0.1,0.1);
            initial_values[1] = interval(0,0.1);
            //  uncertain parameter 
            params[0] =  interval(1.9,2.1);  // Kp
            params[1] =  interval(2.9,3.1);    // Kd
        }
        else if (syschoice == 6) // self-driving car with piecewise constant parameters; sysdim = 4, jacdim = 4
        {
            tau = 0.02;
            t_end = 5.;
            order = 3;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(-0.1,0.1);
            initial_values[1] = interval(0,0.1);
            inputs[0] =  interval(1.9,2.1);  // Kp
            inputs[1] =  interval(2.9,3.1);    // Kd
            is_uncontrolled[1] = true; // Kd uncontrolled
         //   is_variable[2] = true;  // piecewise constant
          //  is_uncontrolled[2] = true;  // Kp uncontrolled
          //   is_variable[3] = true; // piecewise constant
        }
  // REFLECHIR COMMENT GERER CA DIFFEREMMENT
        else if (syschoice == 7) // self-driving car with time varying parameters; sysdim = 4, jacdim = 4
        {
            tau = 0.02;
            t_end = 5.;
            order = 3;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(-0.1,0.1);
            initial_values[1] = interval(0,0.1);
            initial_values[2] =  interval(1.9,2.1);  // Kp
            initial_values[3] =  interval(2.9,3.1);    // Kd
         //   is_uncontrolled[3] = true; // Kd uncontrolled
            // REFLECHIR COMMENt GERER is_variable[2] et is_variable[3]
         //   is_variable[2] = true;  // attention, when changing from const to time-varying the differential system must also be modified in ode_def.h
          //  is_uncontrolled[2] = true;  // Kp uncontrolled
          //   is_variable[3] = true; // attention the differential system must also be modified in ode_def.h
        }
        else if (syschoice == 8)
        {
            tau = 0.01;
            t_end = 5.;
            order = 3;
            initial_values[0] = interval(0.4,0.5);
        }
        else if (syschoice == 9) // acrobatic quadrotor
        {
            tau = 0.01;
            t_end = 0.5;
            order = 4;
            initial_values[0] = interval(-1.,1.);       // px
            initial_values[1] = interval(-0.1,0.1);        // vx
            initial_values[2] = interval(-1.,1.);       // py
            initial_values[3] = interval(-0.1,0.1);        // vy
            initial_values[4] = interval(-0.1,0.1);     // phi
            initial_values[4] = interval(-0.1,0.1);     // omega
            inputs[0] = interval(9,9.5125); // interval(0,18.39375);           // T1
            inputs[1] = interval(9,9.5125); // interval(0,18.39375);           // T2
        }
        else if (syschoice == 99)  // acrobatic quadrotor  with m et Iyy as disturbances
        {
            tau = 0.01;
            t_end = 0.5;
            order = 4;
            initial_values[0] = interval(-1.,1.);       // px
            initial_values[1] = interval(-0.1,0.1);        // vx
            initial_values[2] = interval(-1.,1.);       // py
            initial_values[3] = interval(-0.1,0.1);        // vy
            initial_values[4] = interval(-0.1,0.1);     // phi
            initial_values[4] = interval(-0.1,0.1);     // omega
            inputs[0] = interval(9,9.5125); // interval(0,18.39375);           // T1
            inputs[1] = interval(9,9.5125); // interval(0,18.39375);           // T2
            inputs[2] = interval(1.25,1.25);     // m
            is_uncontrolled[2] = true;
            inputs[3] = interval(0.03,0.03);     // Iyy
            is_uncontrolled[3] = true; 
        }
        else if (syschoice == 10) // 10-D near-hover quadrotor
        {
            tau = 0.01;
            t_end = 0.3;
            order = 3;
            initial_values[0] = interval(-1.,1.);       // px
            initial_values[1] = interval(-0.1,0.1);     // vx
            initial_values[2] = interval(-0.1,0.1);     // thetax
            initial_values[3] = interval(-0.1,0.1);     // omegax
            initial_values[4] = interval(-1.,1.);       // py
            initial_values[5] = interval(-0.1,0.1);     // vy
            initial_values[6] = interval(-0.1,0.1);     // thetay
            initial_values[7] = interval(-0.1,0.1);     // omegay
            initial_values[8] = interval(-2.5,2.5);     // pz
            initial_values[9] = interval(-0.1,0.1);     // vz
            inputs[0] = interval(-0.5,0.5); is_uncontrolled[0] = true;  // disturbance dx
            inputs[1] = interval(-0.5,0.5); is_uncontrolled[1] = true;  // disturbance dy
            inputs[1] = interval(-0.5,0.5); is_uncontrolled[2] = true;  // disturbance dz
            inputs[3] = interval(-0.17453,0.17453);     // control input Sx - desired pitch angle (+/-Pi/18)
            inputs[4] = interval(-0.17453,0.17453);     // control input Sy - desired roll angle
            inputs[5] = interval(0,19.62);              // control input Sz - vertical thrust <= 2g
        }
        else if (syschoice == 11)   // Dubbins vehicle
        {
            tau = 0.01;
            t_end = 1.;
            order = 3;
            initial_values[0] = interval(-0.5,0.5);       // px
            initial_values[1] = interval(-0.5,0.5);     // py
            initial_values[2] = interval(-0.1,0.1);     // theta
            inputs[0] = interval(-1,1); is_uncontrolled[0] = true;  // disturbance b1
            inputs[1] = interval(-1,1); is_uncontrolled[1] = true;  // disturbance b2
            inputs[1] = interval(-5,5); is_uncontrolled[2] = true;  // disturbance b3
            inputs[3] = interval(-1,1);                  // control a
        }
        else if (syschoice == 12)
        {
            tau = 1.;
            t_end = 2.;
            order = 4;
            initial_values[0] = 1;
            initial_values[1] = 0;
            inputs[0] = interval(0,0.1);
            inputs[1] = interval(0,0.1);
        //    is_variable[2] = true;
        //    is_variable[3] = true;
        }
        else if (syschoice == 13)  // Laub-Loomis Benchmark [Arch 2019]
        {
            tau = 0.1;
            t_end = 20.;
            order = 3;
            interval W = interval(-0.05,0.05);
            // to express that it is the same interval: use params[0] = interval(-0.1,0.1) and inputs[0] = 1.2 + params[0]; etc
            initial_values[0] = 1.2 + W;
            initial_values[1] = 1.05 + W;
            initial_values[2] = 1.5 + W;
            initial_values[3] = 2.4 + W;
            initial_values[4] = 1. + W;
            initial_values[5] = 0.1 + W;
            initial_values[6] = 0.45 + W;
        }
        else if (syschoice == 14) // Van der Pol oscillator [Arch 2019]
        {
            tau = 0.01;
            t_end = 6.;
            order = 10;
            // for mu = 1
            initial_values[0] = interval(1.25,1.55);
            initial_values[1] = interval(2.35,2.45);
           // for mu = 2
            // inputs[0] = interval(1.55,1.85);
           // inputs[1] = interval(2.35,2.45);
        }
        else if (syschoice == 15) // Van der Pol oscillator [Arch 2018 and Sparse Polynomial zonotopes]
        {
            tau = 0.005;
            t_end = 3.15;
            order = 10;
            // for mu = 1
            initial_values[0] = interval(1.23,1.57);
            initial_values[1] = interval(2.34,2.46);
            // for mu = 2
            // inputs[0] = interval(1.55,1.85);
            // inputs[1] = interval(2.35,2.45);
        }
        else if(syschoice == 17) // quadrotor model [Arch 2019]
        {
            tau = 0.1;
            t_end = 5;
            order = 3;
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            // positions
            initial_values[0] = interval(-0.4,0.4);
            initial_values[1] = interval(-0.4,0.4);
            initial_values[2] = interval(-0.4,0.4);
            // velocities
            initial_values[3] = interval(-0.4,0.4);
            initial_values[4] = interval(-0.4,0.4);
            initial_values[5] = interval(-0.4,0.4);
        }
        else if (syschoice == 18) // crazyflie HSCC 2019 paper
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.01;
            t_end = 5.;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
           initial_values[3] = interval(-0.01,0.01); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[5] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] = interval(-0.05,0.05); // * M_PI/180.0;  // z ?
            
	    //            initial_values[3] = //interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
	    //	      initial_values[4] = //interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
	    //	      initial_values[12] = //interval(-0.2,0.2); // * M_PI/180.0;  // z ?
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
           // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
           // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
           // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
          //  inputs[3] = interval(-0.05,0.05);
          //  inputs[4] = interval(-0.05,0.05);
          //  inputs[5] = interval(-0.01,0.01);;
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
          //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
        }
        else if (syschoice == 181) // crazyflie HSCC 2019 paper with neural network controoller
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.02;
            t_end = 2.;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
         //   initial_values[3] = 0; // interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
         //   initial_values[4] = 0; //interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
         //   initial_values[12] = interval(-0.001,0.001); // interval(-0.2,0.2); // * M_PI/180.0;  // z ?
       
            initial_values[3] =  interval(-0.001,0.001); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.001,0.001); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] =  interval(-0.01,0.01); // * M_PI/180.0;  // z ?
            
          //  initial_values[3] =  interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
          //  initial_values[4] = interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
          //  initial_values[12] =  interval(-0.2,0.2); // * M_PI/180.0;  // z ?
            
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
            // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
            //  inputs[3] = interval(-0.05,0.05);
            //  inputs[4] = interval(-0.05,0.05);
            //  inputs[5] = interval(-0.01,0.01);;
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
            //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
            
            inputs[0] = interval(0.0,0.0); // cmd_phi (initial value)
            inputs[1] = interval(0.0,0.0); // cmd_theta
            inputs[2] = interval(0.0,0.0); // cmd_psi
            is_uncontrolled[0] = true;
            is_uncontrolled[1] = true;
            is_uncontrolled[2] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[1] = 50; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[2] = 50; // control is constant for each step of the control loop: will take 67 different values overall
        }
        else if (syschoice == 182) // crazyflie HSCC 2019 paper with neural network controller and agressive PID
        {   // do not forget to initialize the setpoints in the ode_def.h file...
	  tau = 0.03; // PB : dt_commands=0.03s it this OK still for integration?
            t_end = 1.;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
         //   initial_values[3] = 0; // interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
         //   initial_values[4] = 0; //interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
         //   initial_values[12] = interval(-0.001,0.001); // interval(-0.2,0.2); // * M_PI/180.0;  // z ?
       
          //  initial_values[3] =  interval(-0.01,0.01); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
          //  initial_values[4] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
          //  initial_values[12] = interval(-0.01,0.01); // interval(-0.01,0.01); // * M_PI/180.0;  // z ?
            
          //  initial_values[3] =  interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
          //  initial_values[4] = interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
          //  initial_values[12] =  interval(-0.2,0.2); // * M_PI/180.0;  // z ?
            
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
            // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
            //  inputs[3] = interval(-0.05,0.05);
            //  inputs[4] = interval(-0.05,0.05);
            //  inputs[5] = interval(-0.01,0.01);;
            
           // initial_values[3] = interval(-0.001,0.001);
            initial_values[3] = interval(-0.101,-0.099); // -0.1; // interval(-0.101,-0.099); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.001,0.001);  // 0; // interval(-0.001,0.001); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] = interval(-0.01,0.01); // 0; // interval(-0.01,0.01); // * M_PI/180.0;  // z ?
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
            //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
            
            inputs[0] = interval(0.0,0.0); // cmd_phi (initial value)
            inputs[1] = interval(0.0,0.0); // cmd_theta
            inputs[2] = interval(0.0,0.0); // cmd_psi
            is_uncontrolled[0] = true;
            is_uncontrolled[1] = true;
            is_uncontrolled[2] = true;
            nb_inputs[0] = 33; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[1] = 33; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[2] = 33; // control is constant for each step of the control loop: will take 67 different values overall
	    // PB HERE for 2 seconds, this should be 100?
        }
        else if (syschoice == 183) // crazyflie HSCC 2019 paper with neural network controller and agressive PID
        {   // do not forget to initialize the setpoints in the ode_def.h file...
	  tau = 0.02; // PB : dt_commands=0.03s it this OK still for integration?
            t_end = 2.;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
         //   initial_values[3] = 0; // interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
         //   initial_values[4] = 0; //interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
         //   initial_values[12] = interval(-0.001,0.001); // interval(-0.2,0.2); // * M_PI/180.0;  // z ?
       
            initial_values[3] =  interval(-0.001,0.001); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.001,0.001); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] =  interval(-0.01,0.01); // * M_PI/180.0;  // z ?
            
          //  initial_values[3] =  interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
          //  initial_values[4] = interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
          //  initial_values[12] =  interval(-0.2,0.2); // * M_PI/180.0;  // z ?
            
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
            // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
            //  inputs[3] = interval(-0.05,0.05);
            //  inputs[4] = interval(-0.05,0.05);
            //  inputs[5] = interval(-0.01,0.01);;
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
            //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
            
            inputs[0] = interval(0.0,0.0); // cmd_phi (initial value)
            inputs[1] = interval(0.0,0.0); // cmd_theta
            inputs[2] = interval(0.0,0.0); // cmd_psi
            is_uncontrolled[0] = true;
            is_uncontrolled[1] = true;
            is_uncontrolled[2] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[1] = 50; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[2] = 50; // control is constant for each step of the control loop: will take 67 different values overall
	    // PB HERE for 2 seconds, this should be 100?
        }
        else if (syschoice == 19) {  // academic example, time-varying (piecewise constant) parameters
            tau = 1.0;
            t_end = 2;
            order = 3;
            initial_values[0] = 0;
            initial_values[1] = 0;
            inputs[0] = interval(0,1.);
            // nb_inputs[2] = 1; // default value: constant coeff
            // solution at t=2 is 6 + u1 - u2 (u being the piecewise constant value of param_inputs[1] on each time interval)
        }
        else if (syschoice == 21) {  // academic example, time-varying (piecewise constant) parameters
            tau = 1.0;
            t_end = 2;
            order = 3;
            initial_values[0] = 0;
            initial_values[1] = 0;
            inputs[0] = interval(0,1.);
            nb_inputs[0] = 2; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            // solution at t=2 is 6 + u1 - u2 (u being the piecewise constant value of param_inputs[1] on each time interval)
        }
        else if (syschoice == 22) {  // academic example, time-varying (piecewise constant) parameters
            tau = 1.;
            t_end = 2;
            order = 3;
            initial_values[0] = 0;
            initial_values[1] = 0;
            inputs[0] = interval(0,1.);
            nb_inputs[0] = 2; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            //  nb_inputs[2] = 1; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            // solution at t=2 is 6 + u1 - u2 (u being the piecewise constant value of param_inputs[1] on each time interval)
        }
        else if (syschoice == 23) {   // pursuer-evader example Mitchell
            tau = 0.1;
            t_end = 1;
            order = 3;
            initial_values[0] = interval(-5,5);
            initial_values[1] = interval(-5,5);
            initial_values[2] = interval(-1,1);
            inputs[0] = interval(-1,1.);    // angular velocity a of the evader (control)
            inputs[1] = interval(-1,1.);    // angular velocity b of the pursuer (disturbance)
            is_uncontrolled[1] = true;  // disturbance
        }
        else if (syschoice == 24) { // [Franzle et al.]
            tau = 0.1;
            t_end = 1;
            order = 3;
            initial_values[0] = interval(-0.8,0.8);
            initial_values[1] = interval(-0.8,0.8);
            inputs[0] = interval(-0.01,0.01);
            nb_inputs[0] = 10; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            is_uncontrolled[0] = true;  // disturbance
        }
        else if (syschoice == 25) { // [Franzle et al.] reversed time van der pol oscillator with uncertainty
            tau = 0.01;
            t_end = 1;
            order = 3;
            initial_values[0] = interval(-0.5,0.5);
            initial_values[1] = interval(-0.5,0.5);
            inputs[0] = interval(-0.01,0.01);
            nb_inputs[0] = 10; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            is_uncontrolled[0] = true;  // disturbance
        }
        else if (syschoice == 26) { // [Franzle et al.] 7-d biological system
            tau = 0.1;
            t_end = 1;
            order = 3;
            initial_values[0] = interval(-0.25,0.25);
            initial_values[1] = interval(-0.45,0.05);
            initial_values[2] = interval(-0.25,0.25);
            initial_values[3] = interval(-0.25,0.25);
            initial_values[4] = interval(-0.25,0.25);
            initial_values[5] = interval(-0.25,0.25);
            initial_values[6] = interval(-0.25,0.25);
            inputs[0] = interval(-0.1,0.1);
            nb_inputs[0] = 10; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            is_uncontrolled[0] = true;  // disturbance
        }
        else if (syschoice == 27) { // [Franzle et al.] 7-d biological system but with sharing
            tau = 0.01;
            t_end = 1;
            order = 3;
            initial_values[0] = interval(-0.25,0.25);
            initial_values[1] = interval(-0.45,0.05);
            initial_values[2] = interval(-0.25,0.25);
            initial_values[3] = interval(-0.25,0.25);
            initial_values[4] = interval(-0.25,0.25);
            initial_values[5] = interval(-0.25,0.25);
            initial_values[6] = interval(-0.25,0.25);
            inputs[0] = interval(-0.1,0.1);
            nb_inputs[0] = 10; // piecewise constant input changes value every t_end/nb_inputs[i] seconds
            is_uncontrolled[0] = true;  // disturbance
        }
        else if (syschoice == 28) { // [Franzle et al.] 7-d biological system without disturbance
            tau = 0.01;
            t_end = 1;
            order = 3;
            initial_values[0] = interval(-0.25,0.25);
            initial_values[1] = interval(-0.45,0.05);
            initial_values[2] = interval(-0.25,0.25);
            initial_values[3] = interval(-0.25,0.25);
            initial_values[4] = interval(-0.25,0.25);
            initial_values[5] = interval(-0.25,0.25);
            initial_values[6] = interval(-0.25,0.25);
        }
        else if (syschoice == 29) { // EX_10 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.01;
            t_end = 0.2;
            order = 3;
            initial_values[0] = interval(9.5,9.55);
            initial_values[1] = interval(-4.5,-4.45);
            initial_values[2] = interval(2.1,2.11);
            initial_values[3] = interval(1.5,1.51);
            inputs[0] = interval(0.0,0.0);
            inputs[1] = interval(0.0,0.0);
        }
        else if (syschoice == 30) { // EX_1 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.05;
            t_end = 0.05*30;
            order = 3;
            initial_values[0] = interval(0.5,0.9);
            initial_values[1] = interval(0.5,0.9);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 30; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 301) { // EX_1 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.05;
            t_end = 0.05*30;
            order = 3;
            initial_values[0] = interval(0.5,0.9);
            initial_values[1] = interval(0.5,0.9);
            //inputs[0] = interval(0.0,0.0);
            //is_uncontrolled[0] = true;
            //nb_inputs[0] = 30; // control is constant for each step of the control loop: will take 30 different values overall
            
            NH = network_handler("networks/sherlock_nn_ex1.sfx");
       //     NH = network_handler("dnnuoa/networks/network_1.sfx");
            //        L.reserve(NH.n_hidden_layers+1);
           // for (int i=0 ; i<NH.n_hidden_layers+1 ; i++ ) {
           //     cout << "before layer" << endl;
           //     L[i] = Layer(NH,i);
          //      cout << "after layer" << endl;
          //  }
            
        }
        else if (syschoice == 31) // Quadcopter MB vec 3 composantes en plus
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.0005;
            t_end = 20*tau*10;
            order = 3;
            
           for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
            initial_values[7] = interval(-0.00872,0.00872); // p
            initial_values[8] = interval(-0.00872,0.00872); // q
            initial_values[0] = interval(-0.2,0.2); // z
            inputs[0] = interval(0.0,0.0); // cmd_phi
            inputs[1] = interval(0.0,0.0); // cmd_theta
            inputs[2] = interval(0.0,0.0); // cmd_psi
            is_uncontrolled[0] = true;
            is_uncontrolled[1] = true;
            is_uncontrolled[2] = true;
            nb_inputs[0] = 10; // control is constant for each step of the control loop
            nb_inputs[1] = 10; // control is constant for each step of the control loop
            nb_inputs[2] = 10; // control is constant for each step of the control loop
        }
        else if (syschoice == 311) // Quadcopter MB a
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.0005;
            t_end = 20*tau*10;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
            initial_values[7] = interval(-0.008,0.008); // p
            initial_values[8] = interval(0.0,0.0); // q
            initial_values[0] = interval(-0.01,0.01); // z
            inputs[0] = interval(0.0,0.0); // cmd_phi
            inputs[1] = interval(0.0,0.0); // cmd_theta
            inputs[2] = interval(0.0,0.0); // cmd_psi
            is_uncontrolled[0] = true;
            is_uncontrolled[1] = true;
            is_uncontrolled[2] = true;
            nb_inputs[0] = 10; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[1] = 10; // control is constant for each step of the control loop: will take 67 different values overall
            nb_inputs[2] = 10; // control is constant for each step of the control loop: will take 67 different values overall
           // cout << "syschoice 311" << endl;
        }
        else if (syschoice == 32) { // EX_2 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.7,0.9);
            initial_values[1] = interval(0.42,0.58);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 33) { // EX_3 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;
            t_end = 0.1*100;
            order = 3;
            initial_values[0] = interval(0.8,0.9);
            initial_values[1] = interval(0.4,0.5);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 100; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 34) { // EX_4 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.35,0.45);
            initial_values[1] = interval(0.25,0.35);
            initial_values[2] = interval(0.35,0.45);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 35) { // EX_5 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.3,0.4);
            initial_values[1] = interval(0.3,0.4);
            initial_values[2] = interval(-0.4,-0.3);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 36) { // EX_6 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.35,0.4);
            initial_values[1] = interval(-0.35,-0.3);
            initial_values[2] = interval(0.35,0.4);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 37) { // EX_7 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.05;
            t_end = 0.5*20;
            order = 3;
            initial_values[0] = interval(0.35,0.45);
            initial_values[1] = interval(0.45,0.55);
            initial_values[2] = interval(0.25,0.35);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 20; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 38) { // EX_8 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;
            t_end = 0.2*25;
            order = 3;
            initial_values[0] = interval(0.5,0.6);
            initial_values[1] = interval(0.5,0.6);
            initial_values[2] = interval(0.5,0.6);
            initial_values[3] = interval(0.5,0.6);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 25; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 381) { // EX_8 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference - EX 19 dnas le repository
            tau = 0.02;
            t_end = 0.2*25;
            order = 3;
            initial_values[0] = interval(0.5,0.501); //interval(0.5,0.6);
            initial_values[1] = interval(0.5,0.501); //interval(0.5,0.6);
            initial_values[2] = interval(0.5,0.501); //interval(0.5,0.6);
            initial_values[3] = interval(0.5,0.501); //interval(0.5,0.6);
            params= NH.eval_network(initial_values);
            control_period = 0.2;
//            inputs[0] = interval(0.0,0.0);
//            is_uncontrolled[0] = true;
//            nb_inputs[0] = 25; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 382) { // Ex 12 in sherlock/systems_with_networks
            tau = 0.01;
            t_end = 6;
            order = 3;
            initial_values[0] = interval(0.5,0.55);  // interval(0.5,0.9);
            initial_values[1] = interval(0.5,0.55);   // interval(0.5,0.9);
            params= NH.eval_network(initial_values);
            control_period = 0.2; // a verifier
        }
        else if (syschoice == 383) { // EX_2 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;                 // Ex 13 in sherlock/systems_with_networks
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.7,0.9); // interval(0.7,0.9);
            initial_values[1] = interval(0.42,0.5); // interval(0.42,0.58);
            params= NH.eval_network(initial_values);
            control_period = 0.2;
        }
        else if (syschoice == 384) { // EX_3 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;              // Ex 14 in sherlock/systems_with_networks
            t_end = 0.1*100;
            order = 3;
            initial_values[0] = interval(0.8,0.9); // interval(0.8,0.9);
            initial_values[1] = interval(0.4,0.5); // interval(0.4,0.5);
            params= NH.eval_network(initial_values);
            control_period = 0.1;
        }
        else if (syschoice == 385) { // EX_4 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;                 // Ex 15 in sherlock/systems_with_networks
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.35,0.45);
            initial_values[1] = interval(0.25,0.35);
            initial_values[2] = interval(0.35,0.45);
            params= NH.eval_network(initial_values);
            control_period = 0.2;
        }
        else if (syschoice == 386) { // EX_5 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;             // Ex 16 in sherlock/systems_with_networks
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.3,0.4);
            initial_values[1] = interval(0.3,0.4);
            initial_values[2] = interval(-0.4,-0.3);
            params= NH.eval_network(initial_values);
            control_period = 0.2;
        }
        else if (syschoice == 387) { // EX_6 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.02;             // Ex 17 in sherlock/systems_with_networks
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(0.35,0.4);
            initial_values[1] = interval(-0.35,-0.3);
            initial_values[2] = interval(0.35,0.4);
            params= NH.eval_network(initial_values);
            control_period = 0.2;
        }
        else if (syschoice == 388) { // EX_7 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.05;
            t_end = 0.5*20;
            order = 3;
            initial_values[0] = interval(0.39,0.41); // interval(0.35,0.45);
            initial_values[1] = interval(0.49,0.51); // interval(0.45,0.55);
            initial_values[2] = interval(0.29,0.31); // interval(0.25,0.35);
            params= NH.eval_network(initial_values);
            control_period = 0.5;
        }
        else if (syschoice == 39) { // EX_9 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.1;
            t_end = 20;
            order = 3;
            initial_values[0] = interval(0.6,0.7);
            initial_values[1] = interval(-0.7,-0.6);
            initial_values[2] = interval(-0.4,-0.3);
            initial_values[3] = interval(0.5,0.6);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 20; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 40) { // EX_10 Reachability for Neural Feedback Systems using Regressive Polynomial Rule Inference
            tau = 0.01;
            t_end = 0.2*50;
            order = 3;
            initial_values[0] = interval(9.5,9.55);
            initial_values[1] = interval(-4.5,-4.45);
            initial_values[2] = interval(2.1,2.11);
            initial_values[3] = interval(1.5,1.51);
            inputs[0] = interval(0.0,0.0);
            inputs[1] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            is_uncontrolled[1] = true;
            nb_inputs[0] = 50; // control is constant for each step of the control loop: will take 30 different values overall
            nb_inputs[1] = 50; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 41) { // essai sys couple
            tau = 0.1;
            t_end = 0.1*5;
            order = 3;
            initial_values[0] = interval(0.,1.);
            initial_values[1] = interval(1.,2.);
            inputs[0] = interval(0.0,0.0);
            is_uncontrolled[0] = true;
            nb_inputs[0] = 5; // control is constant for each step of the control loop: will take 30 different values overall
        }
        else if (syschoice == 42) // crazyflie HSCC 2019 paper
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.01;
            t_end = 2.5;
            order = 5;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
            initial_values[3] = interval(-0.00872,0.00872); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.00872,0.00872); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] = interval(-0.2,0.2); // * M_PI/180.0;  // z ?
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
            // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
            //  inputs[3] = interval(-0.05,0.05);
            //  inputs[4] = interval(-0.05,0.05);
            //  inputs[5] = interval(-0.01,0.01);;
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
            //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
        }
        else if (syschoice == 43) // crazyflie HSCC 2019 paper
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.01;
            t_end = 5;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
            initial_values[3] = interval(-0.01,0.01); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[5] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] = interval(-0.05,0.05); // * M_PI/180.0;  // z ?
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
            // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
            //  inputs[3] = interval(-0.05,0.05);
            //  inputs[4] = interval(-0.05,0.05);
            //  inputs[5] = interval(-0.01,0.01);;
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
            //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
        }
       else if (syschoice == 44) // crazyflie HSCC 2019 paper
        {   // do not forget to initialize the setpoints in the ode_def.h file...
            tau = 0.01;
            t_end = 3;
            order = 3;
            
            for (int j=0 ; j<sysdim; j++)
                initial_values[j] = 0;
            
            initial_values[3] = interval(-0.01,0.01); // = interval(-0.5,0.5) * M_PI/180.0;  // p ?
            initial_values[4] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[5] = interval(-0.01,0.01); //interval(-0.5,0.5) * M_PI/180.0;  // q ?
            initial_values[12] = interval(-0.05,0.05); // * M_PI/180.0;  // z ?
            
            // roll yaw pitch (degree) inputs value (here we consider input as initial)
            // inputs[0] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[1] = interval(3.0 , 5.0) * M_PI/180.0;
            // inputs[2] = 0.0 * M_PI/180.0;
            
            // p , q , r in rad/s -> the value here is an upper bound of the gyro noise of crazyflie
            //  inputs[3] = interval(-0.05,0.05);
            //  inputs[4] = interval(-0.05,0.05);
            //  inputs[5] = interval(-0.01,0.01);;
            
            // err_p , err_q , err_r
            initial_values[6] = 0.0;
            initial_values[7] = 0.0;
            initial_values[8] = 0.0;
            
            // body speed u , v and w -> for embedded verif we instead use world speed
            initial_values[9] = 0.0;
            initial_values[10] = 0.0;
            initial_values[11] = 0.0;
            
            // Z and err_Vz
            //  inputs[12] = interval(-0.1 , 0.1);
            initial_values[13] = 0.0;
        }
       else if (syschoice == 45) { // Mountain car Verisig
           tau = 0.1;
           t_end = 15.;
           order = 3;
           initial_values[0] = interval(-0.5,-0.48);
           initial_values[1] = interval(0.,0.001);
           //inputs[0] = interval(0.0,0.0);
           //is_uncontrolled[0] = true;
           //nb_inputs[0] = 30; // control is constant for each step of the control loop: will take 30 different values overall
           
           NH = network_handler("dnnuoa/networks/verisig_mc_16x16.sfx");
          // for (int i=0 ; i<NH.n_hidden_layers+1 ; i++ )
          //     L[i] = Layer(NH,i);
       }
       else if (syschoice == 451) { // Mountain car Verisig
           tau = 0.5;
           t_end = 71.; // t_end = 115.;
           order = 3;
      //     initial_values[0] = interval(-0.41,-0.4);  // interval(-0.53,-0.5); // interval(-0.53,-0.5);  // interval(-0.5,-0.48);
           initial_values[0] = interval(-0.53,-0.5);
           initial_values[1] = interval(0.,0.);   // interval(0.,0.001);
           
           params= NH.eval_network(initial_values);
           //inputs[0] = interval(0.0,0.0);
           //is_uncontrolled[0] = true;
           //nb_inputs[0] = 30; // control is constant for each step of the control loop: will take 30 different values overall
           
         //  NH = network_handler("dnnuoa/networks/verisig_mc_16x16.sfx");
           // for (int i=0 ; i<NH.n_hidden_layers+1 ; i++ )
           //     L[i] = Layer(NH,i);
           control_period = 1.;
       }
       else if (syschoice == 46) {  // Tora Heterogeneous ARCH-COMP 2020 - NNV (avec RNN format sfx obtenu a partir du .mat),
           tau = 0.1;
           t_end = 5.;
           order = 3;
           initial_values[0] = interval(-0.77,-0.75);
           initial_values[1] = interval(-0.45,-0.43);
           initial_values[2] = interval(0.51,0.54);
           initial_values[3] = interval(-0.3,-0.28);
           // lb = [-0.77; -0.45; 0.51; -0.3];
           // ub = [-0.75; -0.43; 0.54; -0.28];
           // reachStep = 0.01; controlPeriod = 0.5;
          // goal = Box([-0.1;-0.9],[0.2;-0.6]);
       }
       else if (syschoice == 461) {  // Tora Heterogeneous ARCH-COMP 2020 - NNV (avec RNN format sfx obtenu a partir du .mat), // ReachNN avec version tanh
           // values coming from https://github.com/verivital/ARCH-COMP2020/blob/master/benchmarks/Tora_Heterogeneous/reachTora_sigmoid.m
           tau = 0.05;
           t_end = 5.;
           order = 3;
           initial_values[0] = interval(-0.77,-0.75);
           initial_values[1] = interval(-0.45,-0.43);
           initial_values[2] = interval(0.51,0.54);
           initial_values[3] = interval(-0.3,-0.28);
           // lb = [-0.77; -0.45; 0.51; -0.3];
           // ub = [-0.75; -0.43; 0.54; -0.28];
           // reachStep = 0.01; controlPeriod = 0.5;
           // goal = Box([-0.1;0.2 ],[-0.9;-0.6]);
            params= NH.eval_network(initial_values);
           
            control_period = 0.1;
       }
       else if (syschoice == 471) {  // Ex 1 ReachNNstar  (avec RNN format sfx obtenu a partir du .txt)
           tau = 0.05; // a bit strange, is more precise with 0.01 or 0.02 this way (and satisfies the spec) than when tau = 0.005
           t_end = 7.;  // 35 control steps
          //  t_end = 6.;
           order = 3;
           initial_values[0] = interval(0.8,0.9);
           initial_values[1] = interval(0.5,0.6);
           params= NH.eval_network(initial_values);
           // goal [0,0.2],[0.05,0.3]
           target_set[0] = interval(0,0.2);
           target_set[1] = interval(0.05,0.3);
           unsafe_set[0] = interval(0.3,0.8);
           unsafe_set[1] = interval(-0.1,0.4);
           control_period = 0.2;
       }
       else if (syschoice == 1111) {  //toy example
           tau = 0.01; //
           t_end = 1.0;
           order = 3;
           initial_values[0] = interval(0.8,0.9);
           initial_values[1] = interval(0.5,0.6);
           params= NH.eval_network(initial_values);
           control_period = 0.2; // 0.01; // control_period = 0.2;
       }
       else if (syschoice == 1113) {  //toy example
           tau = 0.01; //
           t_end = 1.0;
           order = 3;
           initial_values[0] = interval(0.8,0.9);
           initial_values[1] = interval(0.5,0.6);
           control_period = 0.2; // 0.01; // control_period = 0.2;
       }
       else if (syschoice == 481) {  // Ex 2 ReachNNstar  (avec RNN format sfx obtenu a partir du .txt)
           tau = 0.05; //
           t_end = 1.8;  // 9 control steps
           //  t_end = 6.;
           order = 3;
           initial_values[0] = interval(0.7,0.9);
           initial_values[1] = interval(0.7,0.9);
           params= NH.eval_network(initial_values);
           // goal [-0.3,0.1],[-0.35,0.5]
           control_period = 0.2;
       }
       else if (syschoice == 482) {  // Ex 3 ReachNNstar  (avec RNN format sfx obtenu a partir du .txt)
           tau = 0.05; //
           t_end = 6.0;  // 60 control steps
           //  t_end = 6.;
           order = 3;
           initial_values[0] = interval(0.8,0.9);
           initial_values[1] = interval(0.4,0.5);
           params= NH.eval_network(initial_values);
           // goal [0.2,0.3],[-0.3,-0.05]
           control_period = 0.1;
       }
       else if (syschoice == 483) {  // Ex 4 ReachNNstar  (avec RNN format sfx obtenu a partir du .txt)
           tau = 0.05; //
           t_end = 1.0;  // 10 control steps
           //  t_end = 6.;
           order = 3;
           initial_values[0] = interval(0.25,0.27);
           initial_values[1] = interval(0.08,0.1);
           initial_values[2] = interval(0.25,0.27);
           params= NH.eval_network(initial_values);
           // goal [-0.05,0.05],[-0.05,-0.]
           control_period = 0.1;
       }
       else if (syschoice == 484) {  // Ex 5 ReachNNstar  (avec RNN format sfx obtenu a partir du .txt)
           tau = 0.05; //
           t_end = 2.0;  // 10 control steps
           //  t_end = 6.;
           order = 3;
           initial_values[0] = interval(0.38,0.4);
           initial_values[1] = interval(0.45,0.47);
           initial_values[2] = interval(0.25,0.27);
           params= NH.eval_network(initial_values);
           // goal [-0.4,-0.28],[0.05,0.22]
           control_period = 0.2;
       }
       else if (syschoice == 491) // Ex ACC de Verisig (avec nn obtenu a partir du yaml)
       {
           tau = 0.05; //
           t_end = 5.0;
           order = 3;
           initial_values[0] = interval(90.0,91.0);
           initial_values[1] = interval(32,32.05);
           initial_values[2] = interval(0,0);
           initial_values[3] = interval(10,11);
           initial_values[4] = interval(30,30.05);
           initial_values[5] = interval(0,0);
           params= NH.eval_network(syst_to_nn(initial_values));
         //  cout << params;
         //  params= NH.eval_network(initial_values);
           control_period = 0.1;
       }
       else if (syschoice == 493) // Ex QMPC de Verisig (avec nn obtenu a partir du yaml)
       {
           tau = 0.01; //
           t_end = 5.0; // ? a voir
           order = 3;
           initial_values[0] = interval(0.025,0.05);
           initial_values[1] = interval(0,0.025);
           initial_values[2] = interval(0,0);
           initial_values[3] = interval(0,0);
           initial_values[4] = interval(0,0);
           initial_values[5] = interval(0,0);
           // A REPRENDRE UN JOUR
       //    params= nn_to_control(NH.eval_network(syst_to_nn(initial_values)));
         //  cout << params;
         //  params= NH.eval_network(initial_values);
           control_period = 0.1; // ?  a verifier
       }
//        vector<vector<AAF>> J(sysdim, vector<AAF>(sysdim+inputsdim));  // should be jacdim but not yet defined ?
 //       for (int i=0; i<sysdim; i++)
 //           J[i][i] = 1.0;
    //    eval_valandJacobian_nn(initial_values,inputs,0,tau,J);  // remplace l'evaluation params= NH.eval_network(syst_to_nn(initial_values)); ... (a supprimer ensuite de chaque systeme)
    }
    if (systype == 1) // DDE
    {
        if (syschoice == 1)  // running example
        {
            d0 = 1; // delay in DDE
            // nb_subdiv = 50;
             nb_subdiv = 20;  // number of Taylor models on [0,d0] - defines the timestep here
            t_begin = -d0;  // starting time is -d0 (delay)
            // t_end = 15;
            t_end = 2.;  // final time
            // order = 3;
            order = 2;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(0.33,1.0);
            // nb_subdiv_init = 5;
        }
        else if (syschoice == 2)
        {
            d0 = 1; // delay in DDE
            nb_subdiv = 33;  // number of Taylor models on [0,d0]
            t_begin = -d0;  // starting time is -d0 (delay)
            t_end = 10.;
            order = 3;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(0.9,1.1);
            initial_values[1] = interval(0.9,1.1);
        }
        else if (syschoice == 3)  // Xue 2017 (Ex 3)
        {
            d0 = 0.01; // delay in DDE
            nb_subdiv = 1;  // number of Taylor models on [0,d0]
            t_begin = 0.0;  // starting time is 0 (delay)
            t_end = 0.1;
            order = 2;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(1.1,1.2); // et non (1.1,1.3) comme indique dans le papier
            initial_values[1] = interval(0.95,1.15);
            initial_values[2] = interval(1.4,1.6);
            initial_values[3] = interval(2.3,2.5);
            initial_values[4] = interval(0.9,1.1);
            initial_values[5] = interval(0.0,0.2);
            initial_values[6] = interval(0.35,0.55);
        }
        else if (syschoice == 4)  // Szczelina 1  2014
        {
            d0 = 1.0;
            nb_subdiv = 10;  // number of Taylor models on [0,d0]
            t_begin = -d0;  // starting time is -d0 (delay)
            t_end = 2.;
            order = 2;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(0.9,1.1);
        }
        else if (syschoice == 5) // Szczelina 2  2014
        {
            d0 = 1.0;
            nb_subdiv = 10;  // number of Taylor models on [0,d0]
            t_begin = -d0;  // starting time is -d0 (delay)
            t_end = 2.;
            order = 3;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(0.9,1.1);
        }
        else if (syschoice == 6) // self-driving car; sysdim = 2, jacdim = 2
        {
            d0 = 0.2;
            nb_subdiv = 5;  // number of Taylor models on [0,d0]
            t_begin = -d0;  // starting time is -d0 (delay)
            t_end = 5.;
            order = 3;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(-0.1,0.1);
            initial_values[1] = interval(0,0.1);
            params[0] =  interval(1.9,2.1); // 2;  // Kp
            params[1] =  interval(2.9,3.1);  // 3;   // Kd
        }
        else if (syschoice == 8) // self-driving car bt with coeff in interv; sysdim = 2; jacdim = 4
        {
            d0 = 0.2;
            nb_subdiv = 5;  // number of Taylor models on [0,d0]
            t_begin = -d0;  // starting time is -d0 (delay)
            t_end = 5.;
            order = 3;  // order of Taylor Models
            // uncertain parameter occurring in initial condition
            initial_values[0] = interval(-0.1,0.1);
            initial_values[1] = interval(0,0.1);
            inputs[0] = interval(1.9,2.1);  // Kp
            inputs[1] = interval(2.9,3.1);   // Kd
        //    is_uncontrolled[2] = true;
            is_uncontrolled[1] = true;
        }
        else if (syschoice == 9) // Zou CAV 2015
        {
            d0 = 1.;
            nb_subdiv = 20;
            t_begin = 0;
            t_end = 2;
            order = 5;
            initial_values[0] = interval(3.,6.);
        }
        else if (syschoice == 10)  // platoon of 5 vehicles
        {
            d0 = 0.3;
            t_begin = -d0;
            t_end = 10;
            nb_subdiv = 3;
            order = 3;
            initial_values[0] = interval(-0.01,0.01);
            initial_values[1] = interval(-1.2,-0.8);
            initial_values[2] = interval(1.99,2.01);
            initial_values[3] = interval(-2.2,-1.8);
            initial_values[4] = interval(1.99,2.01);
            initial_values[5] = interval(-3.2,-2.8);
            initial_values[6] = interval(1.99,2.01);
            initial_values[7] = interval(-4.2,-3.8);
            initial_values[8] = interval(1.99,2.01);
        }
        else if (syschoice == 11)  // platoon of 10 vehicles
        {
            d0 = 0.3;
            t_begin = -d0;
            t_end = 10;
            nb_subdiv = 3;
            order = 3;
            initial_values[0] = interval(-0.01,0.01);
            initial_values[1] = interval(-1.2,-0.8);
            initial_values[2] = interval(1.99,2.01);
            initial_values[3] = interval(-2.2,-1.8);
            initial_values[4] = interval(1.99,2.01);
            initial_values[5] = interval(-3.2,-2.8);
            initial_values[6] = interval(1.99,2.01);
            initial_values[7] = interval(-4.2,-3.8);
            initial_values[8] = interval(1.99,2.01);
            initial_values[9] = interval(-5.2,-4.8);
            initial_values[10] = interval(1.99,2.01);
            initial_values[11] = interval(-6.2,-5.8);
            initial_values[12] = interval(1.99,2.01);
            initial_values[13] = interval(-7.2,-6.8);
            initial_values[14] = interval(1.99,2.01);
            initial_values[15] = interval(-8.2,-7.8);
            initial_values[16] = interval(1.99,2.01);
            initial_values[17] = interval(-9.2,-8.8);
            initial_values[18] = interval(1.99,2.01);
        }
    }
    
    variables_to_display = vector<bool>(sysdim);
    for (int i=0; i< sysdim; i++)
        variables_to_display[i] = true;
    
    
        
}




void init_utils_inputs(double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv)
{
    interval temp;
    int nb_points;
    
    // ******************* for piecewise constant inputs
    fullinputsdim = 0;
    for (int i=0; i<inputsdim; i++) {
        cout << "nb_inputs=" <<nb_inputs[i]<<endl;
        fullinputsdim += nb_inputs[i];
    }
    cout << "fullinputsdim=" << fullinputsdim << endl;
    
    // correspondance between variable inputs wwhich sucessive bounds are stored in inputs of size [jacdim]
    index_param = vector<int>(fullinputsdim);
    index_param_inv = vector<int>(inputsdim);
    int j = 0;
    for (int i = 0; i < inputsdim ; i++) {
        for (int k=0; k<nb_inputs[i] ; k++) {
            index_param[j] = i;
            if (k == 0)
                index_param_inv[i] = j;
            j++;
        }
    }
    
    fullinputs = vector<AAF>(fullinputsdim);
    for (int i=0; i<inputsdim; i++) {
        for (int j=0; j<nb_inputs[i]; j++)
            fullinputs[index_param_inv[i]+j] = inputs[i].convert_int();
    }
    // ****************** end  for piecewise constant inputs
    
    if (systype == 0)
    {
        nb_points = ceil((t_end-t_begin)/tau)+2;
    }
    else // systype == 1
    {
        tau = d0/nb_subdiv;
        nb_points = (ceil((t_end-t_begin)/d0+2))*(nb_subdiv+1);
    }
    
    
    // common to EDO and DDE
    center_initial_values = vector<AAF>(sysdim);
    center_fullinputs = vector<AAF>(fullinputsdim);
    
    Jac_param_inputs = vector<vector<interval>>(inputsdim,vector<interval>(sysdim));
    
    jacdim = sysdim+fullinputsdim;
    
    eps = vector<interval>(jacdim);
    for (int i=0 ; i<sysdim ; i++)
    {
        temp = initial_values[i].convert_int();
        center_initial_values[i] = mid(temp);
        eps[i] = temp-mid(temp);
        //   cout << "initial_values[i] = " << initial_values[i] << endl;
    }
    for (int i=0 ; i<inputsdim ; i++)
    {
        if (is_uncontrolled[i])
            uncontrolled ++;
        if (!is_uncontrolled[i])
            controlled++;
    }
    for (int i=0 ; i<fullinputsdim ; i++)
    {
        temp = fullinputs[i].convert_int();
        center_fullinputs[i] = mid(temp);
        eps[i+sysdim] = temp-mid(temp);
    }
    
    //  cout << "controlled=" << controlled  << " uncontrolled=" << uncontrolled << endl;
    
    
    
    // for saving results
    //  cout << "(t_end-t_begin)*nb_subdiv/d0+1=" << ((t_end-t_begin)/d0+1)*(nb_subdiv+1) << endl;
    Xouter_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xouter_robust_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xouter_minimal_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_joint_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_robust_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_minimal_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xexact_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    t_print = vector<double>(nb_points);
}



void init_subdiv(int current_subdiv, vector<AAF> initial_values_save, vector<AAF> inputs_save, int param_to_subdivide)
{
 //   center_inputs = vector<AAF>(jacdim);
 //   eps = vector<interval>(jacdim);
    
    if (param_to_subdivide < sysdim)
    {
        interval save = initial_values_save[param_to_subdivide].convert_int();
        double delta = (save.sup()-save.inf())/nb_subdiv_init;
        if ((current_subdiv > 1) && (current_subdiv < nb_subdiv_init))
            initial_values[param_to_subdivide] = interval(save.inf()+delta*(current_subdiv-1-recovering),save.inf()+delta*(current_subdiv+recovering));
        else if (current_subdiv == 1)
            initial_values[param_to_subdivide] = interval(save.inf()+delta*(current_subdiv-1),save.inf()+delta*(current_subdiv+recovering));
        else if (current_subdiv == nb_subdiv_init)
            initial_values[param_to_subdivide] = interval(save.inf()+delta*(current_subdiv-1-recovering),save.inf()+delta*(current_subdiv));
 //       cout << "initial_values_save[param_to_subdivide] " << inputs[param_to_subdivide] << endl;
        
        interval temp = initial_values[param_to_subdivide].convert_int();
        center_initial_values[param_to_subdivide] = mid(temp);
        eps[param_to_subdivide] = temp-mid(temp);
    }
    else
    {
        interval save = inputs_save[param_to_subdivide-sysdim].convert_int();
        double delta = (save.sup()-save.inf())/nb_subdiv_init;
        if ((current_subdiv > 1) && (current_subdiv < nb_subdiv_init))
            fullinputs[param_to_subdivide-sysdim] = interval(save.inf()+delta*(current_subdiv-1-recovering),save.inf()+delta*(current_subdiv+recovering));
        else if (current_subdiv == 1)
            fullinputs[param_to_subdivide-sysdim] = interval(save.inf()+delta*(current_subdiv-1),save.inf()+delta*(current_subdiv+recovering));
        else if (current_subdiv == nb_subdiv_init)
            fullinputs[param_to_subdivide-sysdim] = interval(save.inf()+delta*(current_subdiv-1-recovering),save.inf()+delta*(current_subdiv));
      //  cout << "inputs[param_to_subdivide] " << inputs[param_to_subdivide-sysdim] << endl;
    
         interval temp = fullinputs[param_to_subdivide-sysdim].convert_int();
            center_fullinputs[param_to_subdivide-sysdim] = mid(temp);
            eps[param_to_subdivide-sysdim] = temp-mid(temp);
    }
}


// initial condition on [t0,t0+d0] given as x = g(t)
vector <T<AAF>> Initfunc(const  T<AAF> &t, vector<AAF> &beta_initial, vector<AAF> &beta_inputs)
{
    vector<T<AAF>> res(sysdim);
    
    // by default
    for (int i=0 ; i<sysdim; i++)
    res[i] = beta_initial[i];
    
    if (syschoice == 1)  // running example
        res[0] = (1+beta_initial[0]*t)*(1+beta_initial[0]*t);  // ix[0] = beta
    else if (syschoice == 2) //  example 5.15
    {
        res[0] = beta_initial[0]*exp(t);
        res[1] = beta_initial[1]*(1-exp(-1));
    }
    else if (syschoice == 4) // Szczelina_1 2014
    {
        res[0] = beta_initial[0]*sin(M_PI*t/2.0);
    }
    else if (syschoice == 5) // Szczelina_2 2014
    {
        res[0] = beta_initial[0]*sin(M_PI*t/2.0);
    }
    
    return res;
}



// initial condition on [-d0,0] given as x = g(t)
vector <T<F<AAF>>> Initfunc(const  T<F<AAF>> &t, vector<T<F<AAF>>> &beta_initial, vector<T<F<AAF>>> &beta_inputs)
{
    vector<T<F<AAF>>> res(sysdim);
    
    // by default
    for (int i=0 ; i<sysdim; i++)
    res[i] = beta_initial[i];
    
    if (syschoice == 1)  // running example
        res[0] = (1+beta_initial[0]*t)*(1+beta_initial[0]*t);  // ix[0] = beta
    else if (syschoice == 2) //  example 5.15
    {
        res[0] = beta_initial[0]*exp(t);
        res[1] = beta_initial[1]*(1-exp(-1));
    }
    else if (syschoice == 4) // Szczelina_1 2014
    {
        res[0] = beta_initial[0]*sin(M_PI*t/2.0);
    }
    else if (syschoice == 5) // Szczelina_2 2014
    {
        res[0] = beta_initial[0]*sin(M_PI*t/2.0);
    }
    return res;
}



// analytical solution if any (for comparison purpose)
void AnalyticalSol(int current_iteration, double d0)
{
  //  vector<interval> res(sysdim);
    vector<interval> Xouter_min(sysdim), Xouter_max(sysdim);
    
    double t = t_print[current_iteration];
    double beta_inf = initial_values[0].convert_int().inf();
    double beta_sup = initial_values[0].convert_int().sup();
    
    // running example
    //res[0] = (1+beta[0]*t)*(1+beta[0]*t);  // ix[0] = beta
    
    for (int i=0 ; i<sysdim ; i++)
        Xexact_print[0][current_iteration][i] = interval(1,-1); // no analytical solution : bot
    
    if (systype == 1)
    {
    if ((syschoice == 1) && beta_sup <= 1)  // running example
    {
        if (t <= 0)   // on [-d0,0], solution is defined by init function
        {
            Xouter_min[0] = ((1.+beta_inf*t)*(1.+beta_inf*t));
            Xouter_max[0] = ((1.+beta_sup*t)*(1.+beta_sup*t));
            Xexact_print[0][current_iteration][0] = hull(Xouter_min[0],Xouter_max[0]);
        }
        else if (t >= 0 && t <= d0)
        {
            Xouter_min[0] = exp((-1./(3.*beta_inf)*(pow(1.+(t-1.)*beta_inf,3)-pow(1.-beta_inf,3))));
            Xouter_max[0] = exp((-1./(3.*beta_sup)*(pow(1.+(t-1.)*beta_sup,3)-pow(1.-beta_sup,3))));
            Xexact_print[0][current_iteration][0] = hull(Xouter_min[0],Xouter_max[0]);
            
        }
        else if (t >= d0 && t <= 2*d0)
        {
            double aux, temp1, temp2;
            aux = exp((-1./(3.*beta_inf)*(pow(1.+(d0-1.)*beta_inf,3)-pow(1.-beta_inf,3))));
            temp1 = pow(1+(t-2)*beta_inf,3)/(3*beta_inf);
            temp2 = pow(1-beta_inf,3)/(3*beta_inf);
            Xouter_min[0] =  aux * exp(-exp(temp2)*pow(3*beta_inf,-2/3.0)*2.6789385347077476337*(gsl_sf_gamma_inc_P(1/3.0,temp1)-gsl_sf_gamma_inc_P(1/3.0,temp2)));
            aux = exp((-1./(3.*beta_sup)*(pow(1.+(d0-1.)*beta_inf,3)-pow(1.-beta_sup,3))));
            temp1 = pow(1+(t-2)*beta_sup,3)/(3*beta_sup);
            temp2 = pow(1-beta_sup,3)/(3*beta_sup);
            Xouter_max[0] =  aux * exp(-exp(temp2)*pow(3*beta_sup,-2/3.0)*2.6789385347077476337*(gsl_sf_gamma_inc_P(1/3.0,temp1)-gsl_sf_gamma_inc_P(1/3.0,temp2)));
            Xexact_print[0][current_iteration][0] = hull(Xouter_min[0],Xouter_max[0]);
            // res[0] = interval(min(inf(Xouter_min[0]),inf(Xouter_max[0])),max(sup(Xouter_min[0]),sup(Xouter_max[0])));
          //  cout << "beta_inf = " << beta_inf << "beta_sup=" << beta_sup << " res[0] = " << res[0] << endl;
        }
        else
        {
            Xexact_print[0][current_iteration][0] = interval(1,-1); // no analytical solution : bot
        }
        
    }
    else if (syschoice == 2) //  example 5.15
    {
        // example 5.15
        if (t <0)
        {
            Xexact_print[0][current_iteration][0] = initial_values[0].convert_int()*exp(t);
            Xexact_print[0][current_iteration][1] = initial_values[1].convert_int()*(1-exp(-1.));
        }
        else
        {
            Xexact_print[0][current_iteration][0] = initial_values[0].convert_int()*exp(t);
            Xexact_print[0][current_iteration][1] = initial_values[1].convert_int()*(exp(t)-exp(t-1.));
        }
    }
    }
    
    // iterative definition
    // mmmh idealementil faudrait calculer en AAF ici ???
    else if ((systype == 0) && ((syschoice == 19) || (syschoice == 20)|| (syschoice == 21)))
    {
        if (current_iteration == 0)
            Xexact_print[0][current_iteration][0] = initial_values[0].convert_int();
        else {
            if (jacdim>sysdim+inputsdim) {
                double delta_t = t_print[current_iteration]-t_print[current_iteration-1];
                double delta_t_sq = t_print[current_iteration]*t_print[current_iteration]-t_print[current_iteration-1]*t_print[current_iteration-1];
                Xexact_print[0][current_iteration][0] = Xexact_print[0][current_iteration-1][0] + inputs[0].convert_int()*(2*delta_t-delta_t_sq) + 2*delta_t + delta_t_sq/2;
            }
            else {
                double delta_t = t_print[current_iteration]-t_print[0];
                double delta_t_sq = t_print[current_iteration]*t_print[current_iteration]-t_print[0]*t_print[0];
                Xexact_print[0][current_iteration][0] = initial_values[0].convert_int() + inputs[0].convert_int()*(2*delta_t-delta_t_sq) + 2*delta_t + delta_t_sq/2;
            }
        }
    }
    else if ((systype == 0) && ((syschoice == 22)))
    {
        if (current_iteration == 0)
            Xexact_print[0][current_iteration][0] = initial_values[0].convert_int();
        else {
            if (jacdim>sysdim+inputsdim) {
                double delta_t = t_print[current_iteration]-t_print[current_iteration-1];
                double delta_t_sq = t_print[current_iteration]*t_print[current_iteration]-t_print[current_iteration-1]*t_print[current_iteration-1];
                Xexact_print[0][current_iteration][0] = Xexact_print[0][current_iteration-1][0] + (inputs[0].convert_int()+inputs[0].convert_int()*inputs[0].convert_int())*(2*delta_t-delta_t_sq) + 2*delta_t + delta_t_sq/2;
            }
            else {
                double delta_t = t_print[current_iteration]-t_print[0];
                double delta_t_sq = t_print[current_iteration]*t_print[current_iteration]-t_print[0]*t_print[0];
                Xexact_print[0][current_iteration][0] = initial_values[0].convert_int() + (inputs[0].convert_int()+inputs[0].convert_int()*inputs[0].convert_int())*(2*delta_t-delta_t_sq) + 2*delta_t + delta_t_sq/2;
            }
        }
    }
    
}

