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

using namespace std;
//using namespace tinyxml2;




//using namespace fadbad;



// void read_system(const char * system_filename, char *sys_name);

//int systype; // 0 is ODE, 1 is DDE
//int syschoice; // choice of system to analyze

// void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs);





int main(int argc, char* argv[])
{
    // all these parameters in initialization functions defined in ode_def.cpp
    double tn;    // current time
    double tau;   // integration time step (fixed step for now)
    int order;    // order of Taylor expansion
    
    double d0; // = 1;   // delay in DDE
    int nb_subdiv; // = 10;   // number of Taylor models on [0,d0]
   
    /********* DEFINING SYSTEM *******************/
    // default is running example of CAV 2018 paper
    systype = 1; // EDO = 0, DDE = 1, discrete systems = 2
    syschoice = 1;
    
    if (argc >= 3)
    {
        systype = atoi(argv[1]);
        syschoice = atoi(argv[2]);
    }
    else
    {
        cout << "You can choose your system: " << endl;
        cout << "First argument is equation type : 0 for ODE, 1 for DDE" << endl;
        cout << "Second argument is system number :" << endl;
        cout << "For discrete systems, 3rd (optional) argument is number of steps, 4th (optional) argument is order of approx (1 or 2)" << endl;
       // cout << "For ODEs: 1 for 1D running example, 2 for Brusselator, 3 for ballistic," << endl;
       // cout << "4 for linearized ballitsic, 5 for self driving car with jacdim = 2, 6 for self driving car with jacdim = 4,"<<endl;
      //  cout << "For DDEs: 1 for 1D running example, 6 for self driving car with jacdim = 2, 7 for self driving car with sysdim = 4,"<<endl;
     //   cout << "8 for selfdriving with sysdim=2, jacdim=4"<<endl;
    }
    
    /*******************************************************************************************/
   
    if (systype == 2) {
        int nb_steps = 1;
        int order = 1;
        if (argc >= 4)
          nb_steps = atoi(argv[3]);
        if (argc >= 5)
            order = atoi(argv[4]);
        
       // if (syschoice == 15)
       //     discrete_dynamical(nb_steps,order);
         if (syschoice == 15)
            discrete_dynamical_preconditioned(nb_steps,order);
         else if (syschoice == 16)
             discrete_dynamical_preconditioned_3d(nb_steps,order);
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
            discrete_dynamical_preconditioned(nb_steps,order);
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
    
    
   
    
    clock_t begin = clock();

    init_system(argc, argv, t_begin,t_end,tau,d0,nb_subdiv,order); // reads from file if input at command-line
    
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
                cur_step.TM_evalandprint_solutionstep(eps,cur_step.tn+tau);
                cur_step.init_nextstep(tau);
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




















