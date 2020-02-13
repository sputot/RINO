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

void print_finalstats(clock_t begin);
void run_pythonscript_visualization();

void read_system(const char * system_filename, char *sys_name);

int systype; // 0 is ODE, 1 is DDE
int syschoice; // choice of system to analyze

void print_initstats(vector<AAF> &x);

void print_ErrorMeasures(int current_iteration, double d0);
void print_finalsolution(int max_it, double d0);



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
    systype = 1; // EDO = 0, DDE = 1
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
        cout << "First argument is system numbe :" << endl;
       // cout << "For ODEs: 1 for 1D running example, 2 for Brusselator, 3 for ballistic," << endl;
       // cout << "4 for linearized ballitsic, 5 for self driving car with jacdim = 2, 6 for self driving car with jacdim = 4,"<<endl;
      //  cout << "For DDEs: 1 for 1D running example, 6 for self driving car with jacdim = 2, 7 for self driving car with sysdim = 4,"<<endl;
     //   cout << "8 for selfdriving with sysdim=2, jacdim=4"<<endl;
    }
    
    /*******************************************************************************************/
   
   
    
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
            
            tn = t_begin;
            print_initstats(initial_values);
            
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
    run_pythonscript_visualization();
}

// printig solution from all stored results of subdivisiions
void print_finalsolution(int max_it, double d0)
{
    // print final solution + output some stats / error information
    bool no_hole = true;
    
    // starts at delta_t (initial condition is not stored)
    for (current_iteration = 0; current_iteration <= max_it; current_iteration++)
    {
        for (int i=0 ; i<sysdim ; i++)
        {
            Xouter_print[0][current_iteration][i] = Xouter_print[1][current_iteration][i];
            Xouter_robust_print[0][current_iteration][i] = Xouter_robust_print[1][current_iteration][i];
            Xinner_print[0][current_iteration][i] = Xinner_print[1][current_iteration][i];
            Xinner_joint_print[0][current_iteration][i] = Xinner_joint_print[1][current_iteration][i];
            Xinner_robust_print[0][current_iteration][i] = Xinner_robust_print[1][current_iteration][i];
            for (int j=2; j <=nb_subdiv_init; j++)
            {
                Xouter_print[0][current_iteration][i] = hull(Xouter_print[0][current_iteration][i],Xouter_print[j][current_iteration][i]);
                Xouter_robust_print[0][current_iteration][i] = hull(Xouter_robust_print[0][current_iteration][i],Xouter_robust_print[j][current_iteration][i]);
                // verifying there is no hole in the inner-approx
                if (!((Xinner_print[j][current_iteration][i].sup() >= Xinner_print[0][current_iteration][i].inf()) &&
                      (Xinner_print[0][current_iteration][i].sup() >= Xinner_print[j][current_iteration][i].inf())))
                {
                    //  cout << "coucou " << Xinner_print[j][current_iteration][i] << " " << Xinner_print[j][current_iteration][i] << endl;
                    no_hole = false;
                }
                // Warning: joining tubes is correct for inner approx only if no hole
                Xinner_print[0][current_iteration][i] = hull(Xinner_print[0][current_iteration][i],Xinner_print[j][current_iteration][i]);
                Xinner_joint_print[0][current_iteration][i] = hull(Xinner_joint_print[0][current_iteration][i],Xinner_joint_print[j][current_iteration][i]);
                Xinner_robust_print[0][current_iteration][i] = hull(Xinner_robust_print[0][current_iteration][i],Xinner_robust_print[j][current_iteration][i]);
            }
            if (nb_subdiv_init > 1)
                outFile_outer[i] << t_print[current_iteration] << "\t" << inf(Xouter_print[0][current_iteration][i]) << "\t" << sup(Xouter_print[0][current_iteration][i]) << endl;
            
            //  outFile_inner[i] << t_print[current_iteration] << "\t" << inf(Xinner_print[0][current_iteration][i]) << "\t" << sup(Xinner_print[0][current_iteration][i]) << endl;
        }
        
        print_ErrorMeasures(current_iteration,d0);
    }
    
    if (no_hole)
        cout << "NO HOLE when joining the inner-approx tubes";
}

void print_ErrorMeasures(int current_iteration, double d0)
{
    double aux, minwidth_ratio, sum, rel_sum;
  //  vector<interval> Xexact(sysdim);
    
    
    // correct only of no hole
    // min on xi of ratio of width of inner / width of outer
    minwidth_ratio = (sup(Xinner_print[0][current_iteration][0])-inf(Xinner_print[0][current_iteration][0]))/(sup(Xouter_print[0][current_iteration][0])-inf(Xouter_print[0][current_iteration][0]));
    for (int i=1 ; i<sysdim ; i++) {
        aux = (sup(Xinner_print[0][current_iteration][i])-inf(Xinner_print[0][current_iteration][i]))/(sup(Xouter_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
        if (minwidth_ratio > aux)
            minwidth_ratio = aux;
    }
    if (t_print[current_iteration] != 0)
        outFile_width_ratio << t_print[current_iteration] << "\t" << minwidth_ratio << endl;
    
    
    // fills Xexact_print with analytical solution
    AnalyticalSol(current_iteration, d0);
   // cout << "before testing x_exact, current_iteration=" << current_iteration << "t_print[current_iteration] " << t_print[current_iteration]  << endl;
    if (sup(Xexact_print[0][current_iteration][0]) >= inf(Xexact_print[0][current_iteration][0])) // an analytical solution exists
    {
        for (int i=0 ; i<sysdim ; i++)
            outFile_exact[i] << t_print[current_iteration] << "\t" << inf(Xexact_print[0][current_iteration][i]) << "\t" << sup(Xexact_print[0][current_iteration][i]) << endl;
        // mean over the xi of the error between over-approx and exact solution
        sum = 0;
        rel_sum = 0;;
        for (int i=0 ; i<sysdim ; i++)
        {
            aux = max(sup(Xouter_print[0][current_iteration][i])-sup(Xexact_print[0][current_iteration][i]),inf(Xexact_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
            sum += aux;
            rel_sum += aux / (sup(Xexact_print[0][current_iteration][i])-inf(Xexact_print[0][current_iteration][i]));
        }
        sum = sum/sysdim;
        rel_sum = rel_sum/sysdim;
        outFile_meanerror_outer << t_print[current_iteration] << "\t" << sum << endl;
        outFile_relmeanerror_outer << t_print[current_iteration] << "\t" << rel_sum << endl;
        
        // mean over the xi of the error between inner-approx and exact solution
        sum = 0;
        rel_sum = 0;
        for (int i=0 ; i<sysdim ; i++)
        {
            aux = max(sup(Xexact_print[0][current_iteration][i])-sup(Xinner_print[0][current_iteration][i]),inf(Xinner_print[0][current_iteration][i])-inf(Xexact_print[0][current_iteration][i]));
            sum += aux;
            rel_sum += aux / (sup(Xexact_print[0][current_iteration][i])-inf(Xexact_print[0][current_iteration][i]));
        }
        sum = sum/sysdim;
        rel_sum = rel_sum/sysdim;
        outFile_meanerror_inner << t_print[current_iteration] << "\t" << sum << endl;
        outFile_relmeanerror_inner << t_print[current_iteration] << "\t" << rel_sum << endl;
    }
  
    
    // mean over the xi of the error between over-approx and inner-approx
    sum = 0;
    rel_sum = 0;
    for (int i=0 ; i<sysdim ; i++)
    {
        aux = max(sup(Xouter_print[0][current_iteration][i])-sup(Xinner_print[0][current_iteration][i]),inf(Xinner_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
        sum += aux;
        rel_sum += aux / (sup(Xouter_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
    }
    sum = sum/sysdim;
    rel_sum = rel_sum/sysdim;
    outFile_meanerror_diff << t_print[current_iteration] << "\t" << sum << endl;
    outFile_relmeanerror_diff << t_print[current_iteration] << "\t" << rel_sum << endl;
}







void run_pythonscript_visualization()
{
    char command_line[1000];
    cout << "......" << endl;
    cout << "Result files are in the output directory" << endl;
    cout << "Visualize them with by : cd GUI; python3 Visu_output.py ; cd .." << endl;
    cout << "The python script will save files as png files (in same directory output)" << endl;
    char displayed_variables[100];
    bool init=true;
    int j = 0;
    for (int i=0; i<sysdim; i++) {
        if (variables_to_display[i]) {
            j++;
            if (init) {
                sprintf(displayed_variables,"-%d",i+1);
                init = false;
            }
            else
                sprintf(displayed_variables,"%s-%d",displayed_variables,i+1);
        }
    }
    if (!init)
        sprintf(displayed_variables,"%s-",displayed_variables);
    cout << displayed_variables << endl;
    if (j == sysdim)
        sprintf(displayed_variables,"%s","all");
    cout << displayed_variables << endl;
    sprintf(command_line,"cd GUI; python3 Visu_output.py --interactive=%d --printvar=%s; cd ..",interactive_visualization,displayed_variables);
    cout << command_line << endl;
    system(command_line);
}






void print_initstats(vector<AAF> &x)
{
   
    // print initial conditions of the ODE
    for (int i=0 ; i<sysdim ; i++) {
        // print in exit files
        outFile_outer[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        if (uncontrolled > 0) {
            outFile_inner_robust[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
            outFile_outer_robust[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        }
        if (controlled > 0 || uncontrolled > 0) {
            outFile_inner_minimal[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
            outFile_outer_minimal[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        }
        outFile_center[i] << 0<< "\t" << mid(x[i].convert_int()) << "\t" << mid(x[i].convert_int()) << endl;
        outFile_inner[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        outFile_inner_joint[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        
        Xouter_print[current_subdiv][0][i] = x[i].convert_int();
        Xouter_robust_print[current_subdiv][0][i] = x[i].convert_int();
        Xouter_minimal_print[current_subdiv][0][i] = x[i].convert_int();
        Xinner_print[current_subdiv][0][i] = x[i].convert_int();
        Xinner_joint_print[current_subdiv][0][i] = x[i].convert_int();
        Xinner_robust_print[current_subdiv][0][i] = x[i].convert_int();
        Xinner_minimal_print[current_subdiv][0][i] = x[i].convert_int();
        Xexact_print[current_subdiv][0][i] = x[i].convert_int();
    }
    outFile_width_ratio << 0 << "\t" << 1.0 << endl;
    
    // a changer un jour pour t_begin (notamment pour DDE)?
    t_print[0] = 0;
    
    
    cout << "printing at t=0 : x=" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x[" << i <<"]=" << x[i] << "\t";
    cout << endl;
    cout << "x0=" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x0[" << i <<"]=" << mid(x[i].convert_int()) << "\t";
    cout << endl;
    
    
}




// print after the end of the analysis
void print_finalstats(clock_t begin)
{
    clock_t end = clock();
    // double end_time = getTime ( );
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; //getTotalTime (start_time , end_time );
    cout << "elpased time (sec) =" << elapsed_secs << endl;
    
    
    for (int i=0 ; i<sysdim ; i++) {
        
        if (uncontrolled > 0) {
            outFile_inner_robust[i].close();
            outFile_outer_robust[i].close();
        }
        if (uncontrolled > 0 || controlled > 0) {
            outFile_inner_minimal[i].close();
            outFile_outer_minimal[i].close();
        }
        outFile_outer[i].close();
        outFile_exact[i].close();
        outFile_inner[i].close();
        outFile_inner_joint[i].close();
        outFile_center[i].close();
    }
    outFile_width_ratio.close();
    outFile_meanerror_outer.close();
    outFile_meanerror_inner.close();
    outFile_meanerror_diff.close();
    outFile_relmeanerror_outer.close();
    outFile_relmeanerror_inner.close();
    outFile_relmeanerror_diff.close();
}

