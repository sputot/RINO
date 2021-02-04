/* ============================================================================
 File   : main.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This is the main file:
    - the system to analyze is chosen (among those defined in file ode_def.h/cpp)
    - subdivision of the input domain may be asked
    - the actual reachability analyis on the system of ODE or DDE is called
    - results are output (output files + gnuplot script)
 ============================================================================ */

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"
#include "utils.h"
#include "ode_def.h"
#include "matrix.h"
#include "network_handler.h"
//#include "taylor.h"
#include "inner.h"
#include "ode_integr.h"
#include "dde_integr.h"
#include "discrete_system.h"
//#include "tinyxml2.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <ctime>
#include <array>

#include <stdlib.h>
#include <algorithm>

using namespace std;
//using namespace tinyxml2;




//using namespace fadbad;



// void read_system(const char * system_filename, char *sys_name);

//int systype; // 0 is ODE, 1 is DDE
//int syschoice; // choice of system to analyze

// void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs);


// for command line options
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

// for command line options
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


int main(int argc, char* argv[])
{
    // all these parameters in initialization functions defined in ode_def.cpp
    double tn;    // current time
    double tau;   // integration time step (fixed step for now)
    int order;    // order of Taylor expansion
    
    double d0; // = 1;   // delay in DDE
    int nb_subdiv; // = 10;   // number of Taylor models on [0,d0]
    const char* network_filename;
    
    vector<interval> Xouter(sysdim);
   
    /********* DEFINING SYSTEM *******************/
    // default is running example of CAV 2018 paper
    systype = 1; // EDO = 0, DDE = 1, discrete systems = 2, DNN = 3
    syschoice = 1;
    
    char * str_systype = getCmdOption(argv, argv + argc, "-systype");
    if (str_systype)
    {
        if (strcmp(str_systype,"ode")==0)
            systype = 0;
        else if (strcmp(str_systype,"dde")==0)
            systype = 1;
        else if (strcmp(str_systype,"discrete")==0)
            systype = 2;
        else if (strcmp(str_systype,"nn")==0)
            systype = 3;
   //     cout << "systype= " << systype << endl;
    }
    
    char * str_syschoice = getCmdOption(argv, argv + argc, "-syschoice");
    if (str_syschoice)
        syschoice = atoi(str_syschoice);
  //  cout << "syschoice= " << syschoice << endl;
    
    network_filename = getCmdOption(argv, argv + argc, "-nnfile");
    if (network_filename)
    {
        NH = network_handler(network_filename);
       // cout << NH.n_hidden_layers+1 << endl;
        for (int i=0 ; i<NH.n_hidden_layers+1 ; i++ )
            L[i] = Layer(NH,i);
    }
    
    char * config_filename = getCmdOption(argv, argv + argc, "-configfile");
    
    if ((!str_systype) || cmdOptionExists(argv, argv + argc, "-help")  || cmdOptionExists(argv, argv + argc, "-h"))
    {
        cout << "Usage is: " << endl;
        cout << "-systype xxx: ode for ODE, dde for DDE, discrete for discrete-time systems, nn for neural networks" << endl;
        cout << "-syschoice x: integer setting the predefined system ID" << endl;
        cout << "-nbsteps x: for discrete systems only, (optional) number of steps - default value is 1" << endl;
        cout << "-AEextension_order x: for discrete systems only, (optional) order of AE-extensions (1 or 2) - default value is 1" << endl;
        cout << "-configfile xxx: name of the (optional) config file - for ODE and DDE" << endl;
       // cout << "For ODEs: 1 for 1D running example, 2 for Brusselator, 3 for ballistic," << endl;
       // cout << "4 for linearized ballitsic, 5 for self driving car with jacdim = 2, 6 for self driving car with jacdim = 4,"<<endl;
      //  cout << "For DDEs: 1 for 1D running example, 6 for self driving car with jacdim = 2, 7 for self driving car with sysdim = 4,"<<endl;
     //   cout << "8 for selfdriving with sysdim=2, jacdim=4"<<endl;
        exit(0);
    }
    
    /*******************************************************************************************/
   /************* Discrete Systems ************/
    if (systype == 2) {
        
        int nb_steps = 1;
        int AEextension_order = 1;
    
        char * str_nbsteps = getCmdOption(argv, argv + argc, "-nbsteps");
        if (str_nbsteps)
            nb_steps = atoi(str_nbsteps);
        char * str_AEextension_order = getCmdOption(argv, argv + argc, "-AEextension_order");
        if (str_AEextension_order)
            nb_steps = atoi(str_AEextension_order);
        
        
        if (syschoice == 23)
            discrete_dynamical_preconditioned(nb_steps,AEextension_order);
           // discrete_dynamical_method2_preconditioned(nb_steps);
        
        // if (syschoice == 15)
       //     discrete_dynamical(nb_steps,order);
         if (syschoice == 15)
            discrete_dynamical_preconditioned(nb_steps,AEextension_order);
         else if (syschoice == 16)
             discrete_dynamical_preconditioned_3d(nb_steps,AEextension_order);
         // discrete_dynamical_preconditioned(nb_steps,order);
        else if (syschoice == 16)
            discrete_dynamical_method2(nb_steps);
        else if (syschoice == 18)
//            function_range();
            // discrete_dynamical_preconditioned_3d(nb_steps,order);
            discrete_dynamical_method2_preconditioned(nb_steps);
        else if (syschoice == 17 ||  syschoice == 21 ) {
           discrete_dynamical_method2(nb_steps);
         // discrete_dynamical(nb_steps,order);
        }
        else if (syschoice == 20 ||  syschoice == 18)
         discrete_dynamical_method2_preconditioned(nb_steps);
        else if (syschoice == 15 || syschoice == 18 || syschoice == 19  ||  syschoice == 16 || syschoice == 20 )
            discrete_dynamical_preconditioned(nb_steps,AEextension_order);
        else
            function_range();
        
 //       discrete_dynamical_preconditioned(nb_steps);
      //  discrete_dynamical_method2(nb_steps);
        
     /*   if (syschoice == 15 || syschoice == 18 || syschoice == 19  ||  syschoice == 16)
            discrete_dynamical_preconditioned(nb_steps);
        else if (syschoice == 17)
            discrete_dynamical_method2(nb_steps);
    //    else if (syschoice == 16)
     //       discrete_dynamical_preconditioned_3d(nb_steps);
            // discrete_dynamical();
        else
            function_range(); */
        
        return 0;
    }
    else if (systype == 3) {
        // neural network range evaluation
        
        vector<vector<AAF>> net_outputs(NH.n_hidden_layers+2); // outputs for each layer
        
        vector<AAF> net_inputs(NH.n_inputs);
        net_inputs[0] = interval(0,0.1);
        net_inputs[1] = interval(0,0.1);
        
        net_outputs[0] = net_inputs;
        for (int i=0 ; i<NH.n_hidden_layers+1 ; i++ ) {
            net_outputs[i+1] = eval_layer(L[i],net_outputs[i]);
        }
        for (int i=-1 ; i<NH.n_hidden_layers+1 ; i++ )
            cout << "output layer " << i << " is " << net_outputs[i+1] << endl;
        return 0;
    }
    
    
   
    
    clock_t begin = clock();

    init_system(config_filename, t_begin,t_end,tau,d0,nb_subdiv,order); // reads from file if input at command-line
    
    if (config_filename) // called with configuration file: we overwrite the initialization of init_system
        read_parameters(config_filename, tau, t_end, d0, t_begin, order, nb_subdiv);
    
    init_utils_inputs(t_begin,t_end,tau,d0,nb_subdiv);
    
    vector<AAF> initial_values_save(sysdim);
    vector<AAF> fullinputs_save(fullinputsdim);
    for (int i=0 ; i<sysdim; i++)
        initial_values_save[i] = initial_values[i];
    for (int i=0 ; i<fullinputsdim; i++)
        fullinputs_save[i] = fullinputs[i];
    
    DdeFunc bf;    // contains the differential system - defined in ode_def.h
    DdeJacFunc bbf;
    OdeFunc obf;    // contains the differential system - defined in ode_def.h
    
    /*************************************************************************** DDE ************************************************************/
    if (systype == 1) // DDE
    {
        // printing exact solution if any for comparison
      //  print_exactsolutiondde(t_begin, d0, tau, t_end, nb_subdiv/*, ip*/);
        
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, initial_values_save, fullinputs_save, component_to_subdiv);
            
            // a changer de 0 en 1 comme pour ODE (l'iteration 0 servant a stocker l'instant initial non stocke sinon - c'est du moins le cas en ODE, a verifier pour DDE)
            current_iteration = 0;
            
            tn = t_begin;
            
            HybridStep_dde prev_step = HybridStep_dde(bf,bbf,order,tn,tau,d0,nb_subdiv);
            prev_step.init_dde();
            
            HybridStep_dde cur_step = prev_step.init_nextbigstep(tau);
            
            //******** Integration loop ************************
            
            while (cur_step.tn+d0 <= t_end)
            {
                // build Taylor model for Value and Jacobian and deduce guards for each active mode
                for (int j=0 ; j<nb_subdiv ; j++)
                {
                    cur_step.TM_build(j);  // j-th Taylor model valid on [tn+j*d0,tn+(j+1)*d0]
                    // eval at tn+tau
                    cur_step.TM_evalandprint_solutionstep(j,eps); // cur_step.tn+(j+1)*tau/nb_subdiv);
                    cout << endl;
                    
                    if (j < nb_subdiv-1)
                        cur_step.init_nextsmallstep(j);
                }
                cur_step = cur_step.init_nextbigstep(tau);
                
                
            }
            if (current_subdiv<nb_subdiv_init) {
                for (int i=0 ; i<sysdim ; i++)
                    outFile_inner[i] << "\n";
                if (uncontrolled > 0) {
                    for (int i=0 ; i<sysdim ; i++) {
                        outFile_inner_robust[i] << "\n";
                        outFile_outer_robust[i] << "\n";
                    }
                }
                if (controlled > 0 || uncontrolled > 0) {
                    for (int i=0 ; i<sysdim ; i++){
                        outFile_inner_minimal[i] << "\n";
                        outFile_outer_minimal[i] << "\n";
                    }
                    
                }
            }
        }
        
        print_finalsolution(ceil((t_end-t_begin)*nb_subdiv/d0), d0);  // a verifier cf def nb_points dans ode_def.cpp
    }
     /*************************************************************************** EDO ************************************************************/
    else // systype == 0: EDO
    {
        vector<vector<AAF>> J(jacdim, vector<AAF>(jacdim));
        vector<AAF> x(sysdim);
        vector<AAF> param_inputs(jacdim-sysdim);
        vector<AAF> param_inputs_center(jacdim-sysdim);
        vector<AAF> xcenter(sysdim);
        
        
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, initial_values_save, fullinputs_save, component_to_subdiv);
            
            set_initialconditions(param_inputs,param_inputs_center,x,xcenter,J);  //            setId(J0);
            
         //   for (int i=0 ; i<sysdim ; i++)
         //       cout << "x[i]=" << x[i] << endl;
            
            
            tn = t_begin;
            print_initstats(initial_values,param_inputs);
            
            // pb sur les fonctiosn trigos en formes affines a corriger a l'occasion ci-dessous:
        //    vector<AAF> yp(sysdim);
        //    yp[0] = -(5.0*cos(x[2].convert_int()) + param_inputs[0].convert_int());        // px' = v.cos(theta) + b1
        //    yp[1] = -(5.0*sin(x[2].convert_int()) + param_inputs[1].convert_int());        // py' = v.sin(theta) + b2
        //    yp[2] = -(param_inputs[3] + param_inputs[2]);    // theta' = a + b3
        //    for (int i=0 ; i<sysdim ; i++)
        //        cout << "yp[i]=" <<yp[i].convert_int() << " ";
        //    cout << endl;
            
            current_iteration = 1;
            HybridStep_ode cur_step = init_ode(obf,xcenter,x,J,tn,tau,order);
            
            while (cur_step.tn < t_end)
            {
                // build Taylor model for Value and Jacobian and deduce guards for each active mode
                cur_step.TM_build(param_inputs,param_inputs_center);
                Xouter = cur_step.TM_evalandprint_solutionstep(eps,cur_step.tn+tau);
                cur_step.init_nextstep(tau);
                
                if ((Xouter[0] == interval::EMPTY()) || (Xouter[0] == interval::ENTIRE())) {
                    printf("Terminated due to too large overestimation.\n");
                    break;
                }
            }
            // adding a white line separator between subdivisions in the output result (except for maximal outer approx which is computed by union of all subdivisions)
            if (current_subdiv<nb_subdiv_init) {
                for (int i=0 ; i<sysdim ; i++)
                    outFile_inner[i] << "\n";
                if (uncontrolled > 0) {
                    for (int i=0 ; i<sysdim ; i++) {
                        outFile_inner_robust[i] << "\n";
                        outFile_outer_robust[i] << "\n";
                    }
                }
                if (controlled > 0 || uncontrolled > 0) {
                    for (int i=0 ; i<sysdim ; i++){
                        outFile_inner_minimal[i] << "\n";
                        outFile_outer_minimal[i] << "\n";
                    }
                    
                }
            }
            
            
        }
        print_finalsolution(ceil((t_end-t_begin)/tau), d0);
    }
    
    print_finalstats(begin);
    if (print_debug)
        run_pythonscript_visualization();
}




















