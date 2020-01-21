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
void generate_gnuplot_script();


void read_system(const char * system_filename, char *sys_name);

int systype; // 0 is ODE, 1 is DDE
int syschoice; // choice of system to analyze

void print_initstats(vector<AAF> &x);

void print_ErrorMeasures(int current_iteration, vector<AAF> inputs_save, double d0);
void print_finalsolution(vector<AAF> inputs_save, int max_it, double d0);



int main(int argc, char* argv[])
{
    // all these parameters in initialization functions defined in ode_def.cpp
    double tn;    // current time
    double tau;   // integration time step (fixed step for now)
    double t_begin; // starting time of initialization
    double t_end; // ending time of integration
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
    
    vector<AAF> inputs_save(jacdim);
    for (int i=0 ; i<jacdim; i++)
        inputs_save[i] = inputs[i];
    
    DdeFunc bf;    // contains the differential system - defined in ode_def.h
    DdeJacFunc bbf;
    OdeFunc obf;    // contains the differential system - defined in ode_def.h
    
    /*************************************************************************** DDE ************************************************************/
    if (systype == 1) // DDE
    {
        // printing exact solution if any for comparison
        print_exactsolutiondde(t_begin, d0, tau, t_end, nb_subdiv/*, ip*/);
        
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, inputs_save, 0);
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
            for (int i=0 ; i<sysdim ; i++) {
                outFile_inner[i] << "\n \n";
            }
        }
        
        print_finalsolution(inputs_save, (t_end-t_begin)*nb_subdiv/d0, d0);
    }
     /*************************************************************************** EDO ************************************************************/
    else // systype == 0: EDO
    {
        vector<vector<AAF>> J(sysdim, vector<AAF>(jacdim));
        vector<AAF> x(sysdim);
        vector<AAF> xcenter(sysdim);
        
        
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, inputs_save, 0);
            current_iteration = 0;
            
            
    //        cout << "center_inputs[0]" << center_inputs[0] << endl;
    //        cout << inputs[0] << endl;
            
            set_initialconditions(x,xcenter,J);  //            setId(J0);
            
            tn = t_begin;
            print_initstats(inputs);
            
            HybridStep_ode cur_step = init_ode(obf,xcenter,x,J,tn,tau,order);
            
            while (cur_step.tn < t_end)
            {
                // build Taylor model for Value and Jacobian and deduce guards for each active mode
                cur_step.TM_build();
                cur_step.TM_evalandprint_solutionstep(eps,cur_step.tn+tau);
                cur_step.init_nextstep(tau);
            }
        }
        print_finalsolution(inputs_save, (t_end-t_begin)/tau, d0);
    }
    
    print_finalstats(begin);

    generate_gnuplot_script();
}

// printig solution from all stored results of subdivisiions
void print_finalsolution(vector<AAF> inputs_save, int max_it, double d0)
{
    // print final solution + output some stats / error information
    bool no_hole = true;
    
    for (current_iteration = 0; current_iteration <= max_it; current_iteration++)
    {
        for (int i=0 ; i<sysdim ; i++)
        {
            Xouter_print[0][current_iteration][i] = Xouter_print[1][current_iteration][i];
            Xouter_robust_print[0][current_iteration][i] = Xouter_robust_print[1][current_iteration][i];
            Xinner_print[0][current_iteration][i] = Xinner_print[1][current_iteration][i];
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
                Xinner_robust_print[0][current_iteration][i] = hull(Xinner_robust_print[0][current_iteration][i],Xinner_robust_print[j][current_iteration][i]);
            }
            if (nb_subdiv_init > 1)
                outFile_outer[i] << t_print[current_iteration] << "\t" << inf(Xouter_print[0][current_iteration][i]) << "\t" << sup(Xouter_print[0][current_iteration][i]) << endl;
            
            //  outFile_inner[i] << t_print[current_iteration] << "\t" << inf(Xinner_print[0][current_iteration][i]) << "\t" << sup(Xinner_print[0][current_iteration][i]) << endl;
        }
        
        print_ErrorMeasures(current_iteration,inputs_save,d0);
    }
    
    if (no_hole)
        cout << "NO HOLE when joining the inner-approx tubes";
}

void print_ErrorMeasures(int current_iteration, vector<AAF> inputs_save, double d0)
{
    double aux, minwidth_ratio, sum, rel_sum;
    vector<interval> Xexact(sysdim);
    
    
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
    
    
 /*
    Xexact = AnalyticalSol(t_print[current_iteration],inputs_save,d0);
    
    if (sup(Xexact[0]) >= inf(Xexact[0])) // an analytical solution exists
    {
        // mean over the xi of the error between over-approx and exact solution
        sum = 0;
        rel_sum = 0;;
        for (int i=0 ; i<sysdim ; i++)
        {
            aux = max(sup(Xouter_print[0][current_iteration][i])-sup(Xexact[i]),inf(Xexact[i])-inf(Xouter_print[0][current_iteration][i]));
            sum += aux;
            rel_sum += aux / (sup(Xexact[i])-inf(Xexact[i]));
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
            aux = max(sup(Xexact[i])-sup(Xinner_print[0][current_iteration][i]),inf(Xinner_print[0][current_iteration][i])-inf(Xexact[i]));
            sum += aux;
            rel_sum += aux / (sup(Xexact[i])-inf(Xexact[i]));
        }
        sum = sum/sysdim;
        rel_sum = rel_sum/sysdim;
        outFile_meanerror_inner << t_print[current_iteration] << "\t" << sum << endl;
        outFile_relmeanerror_inner << t_print[current_iteration] << "\t" << rel_sum << endl;
    }
  */
    
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







void generate_gnuplot_script()
{
    ofstream gnuplot_script;
    gnuplot_script.open("output/gnuplot_script.gp");
    
    // ---------- plotting variables ----------
    int nb_lignes, nb_colonnes;
    nb_lignes = sqrt(sysdim);
    nb_colonnes = nb_lignes;
    if (nb_lignes * nb_colonnes < sysdim)
        nb_lignes++;
    if (nb_lignes * nb_colonnes < sysdim)
        nb_colonnes++;
    
    // to screen
/*    gnuplot_script << "set terminal aqua 0" << endl;
    gnuplot_script << "set multiplot layout " << nb_lignes << ", " << nb_colonnes << " title \"System solutions \" font \"arial,18\""  << endl;
    for (int i=0 ; i<sysdim ; i++) {
        gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w lp, 'x"<<i+1<<"outer.out' using 1:3 w lp, ";
        gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2 w lp, 'x"<<i+1<<"inner.out' using 1:3 w lp" << endl;
    }
    gnuplot_script << "unset multiplot" << endl; */
    
    // to file
 /*   gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"solutions.png\"" << endl;
    gnuplot_script << "set multiplot layout " << nb_lignes << ", " << nb_colonnes << " title \"System solutions \" font \"arial,18\""  << endl;
    for (int i=0 ; i<sysdim ; i++) {
        gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w lp, 'x"<<i+1<<"outer.out' using 1:3 w lp, ";
 //       gnuplot_script << "'x"<<i+1<<"center.out' using 1:2 w lp, 'x"<<i+1<<"center.out' using 1:3 w lp, ";
        gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2 w lp, 'x"<<i+1<<"inner.out' using 1:3 w lp" << endl;
        
    }
    gnuplot_script << "unset multiplot" << endl;
    gnuplot_script << "unset output" << endl; */
    
    // plotting indiviudually each solution to a file
    gnuplot_script << "set term pngcairo font \"Helvetica,18\"" << endl;
    
    for (int i=0 ; i<sysdim ; i++) {
        if (sysdim == 1)
            gnuplot_script << "set output \"x.png\"" << endl;
        else
            gnuplot_script << "set output \"x"<<i+1<<".png\"" << endl;
        if (nb_subdiv_init > 1)
            gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
        gnuplot_script << "set xlabel 't (seconds)'" << endl;
        if ((systype == 0) && (syschoice == 4) && (i==3)) // ballistic example
        {
            gnuplot_script << "set ylabel 'y(t)'" << endl;
             gnuplot_script << "set arrow 2 from 0,0 to 2,0 nohead  front" << endl;  // black axis
            gnuplot_script << "set object 1 poly from 0,6 to 2,6 to 2,3 to 0,5 to 0,6 fs solid fc rgb \"red\"" << endl;  // red box to avoid
            gnuplot_script << "set key left bottom" << endl;
        //    gnuplot_script << "set object 1 rect from 0,5 to 2,3 fc rgb \"red\"" << endl;  // boite rouge a eviter
        }
        else
            gnuplot_script << "set ylabel 'x(t)'" << endl;
        if ((systype == 0) && (syschoice == 4) && (i==3)) // ballistic example, y coordinate
        {
            gnuplot_script << "set yrange [0:6]" << endl;
            gnuplot_script << "set object 2 rect from 0.75,3 to 1.,4 fc rgb \"green\" front" << endl;  // green box to reach
            gnuplot_script << "set object 1 poly from 0.,2 to 1.4,4.5 to 1.4,5.5 to 0,4 to 0,2 fs solid  fc rgb \"cyan\"" << endl;  // green box to reach
            gnuplot_script << "set object 3 rect from 0.4,0. to 0.6,2 fc rgb \"red\"" << endl;  // red box to avoid
            gnuplot_script << "LABEL1 = \"U\"" << endl;
            gnuplot_script << "set label at 0.5,1.5 LABEL1 front center" << endl;
            gnuplot_script << "LABEL2 = \"T1\"" << endl;
            gnuplot_script << "set label at 0.9,3.5 LABEL2 front center" << endl;
            gnuplot_script << "LABEL3 = \"T2\"" << endl;
            gnuplot_script << "set label at 0.1,3 LABEL3 front center" << endl;
            gnuplot_script << "set key left top" << endl;
        }
        
     /*   if ((uncontrolled > 0)) // && (controlled > 0))
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 1  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 1 notitle, ";
            //       gnuplot_script << "'x"<<i+1<<"exact.out' using 1:2 w l lt 3  dashtype 2 title \"analytical solution\", 'x"<<i+1<<"exact.out' using 1:3 w l lt 3 dashtype 2 title \"\", ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu title \"maximal inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_robust.out' using 1:2:3 w filledcu title \"robust outer flowpipe\", 'x"<<i+1<<"outer_robust.out' using 1:2 w l lt 2 notitle, 'x"<<i+1<<"outer_robust.out' using 1:3 w l lt 2 notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_robust.out' using 1:2:3 w filledcu title \"robust inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_minimal.out' using 1:2:3 w filledcu title \"minimal outer flowpipe\", 'x"<<i+1<<"outer_minimal.out' using 1:2 w l lt 2 notitle, 'x"<<i+1<<"outer_minimal.out' using 1:3 w l lt 2 notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_minimal.out' using 1:2:3 w filledcu title \"minimal inner flowpipe\", " << endl;
            //    gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu title \"maximal inner flowpipe\", 'x"<<i+1<<"inner.out' using 1:2 w l lt 2 notitle, 'x"<<i+1<<"inner.out' using 1:3 w l lt 2 notitle" << endl;
        } */
        if ((uncontrolled > 0)) // && (controlled > 0))
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
            //       gnuplot_script << "'x"<<i+1<<"exact.out' using 1:2 w l lt 3  dashtype 2 title \"analytical solution\", 'x"<<i+1<<"exact.out' using 1:3 w l lt 3 dashtype 2 title \"\", ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_robust.out' using 1:2 w l lt 5  lc rgb \"red\"  lw 2 title \"robust outer flowpipe\", 'x"<<i+1<<"outer_robust.out' using 1:3 w l lt 5  lc rgb \"red\"  lw 2notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_robust.out' using 1:2:3 w filledcu lc rgb \"red\" title \"robust inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe\", 'x"<<i+1<<"outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe\" " << endl;
            //    gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu title \"maximal inner flowpipe\", 'x"<<i+1<<"inner.out' using 1:2 w l lt 2 notitle, 'x"<<i+1<<"inner.out' using 1:3 w l lt 2 notitle" << endl;
        }
        else if (controlled > 0)
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
            //       gnuplot_script << "'x"<<i+1<<"exact.out' using 1:2 w l lt 3  dashtype 2 title \"analytical solution\", 'x"<<i+1<<"exact.out' using 1:3 w l lt 3 dashtype 2 title \"\", ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe\", 'x"<<i+1<<"outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe\" " << endl;
            
        }
        else // no uncertain parameters, only innitial conditions
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
            //       gnuplot_script << "'x"<<i+1<<"exact.out' using 1:2 w l lt 3  dashtype 2 title \"analytical solution\", 'x"<<i+1<<"exact.out' using 1:3 w l lt 3 dashtype 2 title \"\", ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\" " << endl;
            
        }
        
        
        gnuplot_script << "unset output" << endl;
    }
    
    
    gnuplot_script << "set term pngcairo font \"Helvetica,18\"" << endl;
//    gnuplot_script << "set terminal png size 1280, 480 font \"Helvetica,30\"" << endl;
    gnuplot_script << "set output \"xi.png\"" << endl;
    if (nb_subdiv_init > 1)
        gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
    if (syschoice == 6 || syschoice == 7)
    {
        gnuplot_script << "set xrange [0:5]" << endl;
        gnuplot_script << "set yrange [0:1]" << endl;
        gnuplot_script << "set ylabel 'x(t), v(t)'" << endl;
        gnuplot_script << "set key right center" << endl;
    }
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    
    if (uncontrolled > 0)
    {
        gnuplot_script << "set style fill noborder"<<endl;
        
        gnuplot_script << "plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\", ";
        gnuplot_script << "'x1outer_robust.out' using 1:2 w l lt 5  lc rgb \"red\"  lw 2 title \"robust outer flowpipe\", 'x1outer_robust.out' using 1:3 w l lt 5  lc rgb \"red\"  lw 2 notitle,";
        gnuplot_script << "'x1inner_robust.out' using 1:2:3 w filledcu lc rgb \"red\" title \"robust inner flowpipe \", ";
        gnuplot_script << "'x1outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe\", 'x1outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x1inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe\", ";
        
        gnuplot_script << " 'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle , ";
        gnuplot_script << "'x2outer_robust.out' using 1:2 w l lt 5  lc rgb \"red\"  lw 2 notitle, 'x2outer_robust.out' using 1:3 w l lt 5  lc rgb \"red\"  lw 2 notitle,";
        gnuplot_script << "'x2inner_robust.out' using 1:2:3 w filledcu lc rgb \"red\" notitle, ";
        gnuplot_script << "'x2outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle, 'x2outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x2inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' notitle " << endl;
    
    }
    else if (controlled > 0)
    {
        gnuplot_script << "set style fill noborder"<<endl;
        
        gnuplot_script << "plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe \", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe \", ";
        gnuplot_script << "'x1outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe \", 'x1outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x1inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe \", ";
        
        gnuplot_script << " 'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle, 'x2outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x2inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' notitle " << endl;
    }
    else // only initial conditions => only maximal flowpipes
    {
        gnuplot_script << "set style fill noborder"<<endl;
        
        gnuplot_script << "plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe \", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe \", ";
        
        gnuplot_script << " 'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle " << endl;
    }
    gnuplot_script << "unset output" << endl;
    
    
    // plotting the width ration
    gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"width_ratio.png\"" << endl;
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    gnuplot_script << "set ylabel 'min over x_i of width ratios'" << endl;
    gnuplot_script << "plot 'width_ratio.out' w l title \"width(inner-approx) / width (outer-approx)" << endl;
    gnuplot_script << "unset output" << endl;
    
    // plotting  mean on xi of error between outer-approx and analytical solution if any
    gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"meanerror.png\"" << endl;
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    gnuplot_script << "set ylabel 'mean over x_i of max error'" << endl;
    if (nb_subdiv_init > 1)
        gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
    gnuplot_script << "plot 'meanerror_outer.out' w l lt 1 title \"outer-approx error\", 'meanerror_inner.out' w l lt 1 dashtype 2 title \"inner-approx error\", 'meanerror_diff.out' w l lt 2 title \"max distance between inner and outer-approx\"" << endl;
    gnuplot_script << "unset output" << endl;
    
    // plotting  mean on xi of relative error between outer-approx and analytical solution if any (relative because divided by width of exact or over-approx solution)
    gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"relmeanerror.png\"" << endl;
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    gnuplot_script << "set ylabel 'mean over x_i of max relative error'" << endl;
    if (nb_subdiv_init > 1)
        gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
    gnuplot_script << "plot 'relmeanerror_outer.out' w l lt 1 title \"outer-approx error\", 'relmeanerror_inner.out' w l lt 1 dashtype 2 title \"inner-approx error\", 'relmeanerror_diff.out' w l lt 2 title \"max distance between inner and outer-approx\"" << endl;
    gnuplot_script << "unset output" << endl;
    
    
    // gnuplot_script << "set term wxt" << endl;  // vue par defaut
    // gnuplot_script << "set term aqua // pour avoir des chiffres plus gros sur les axes pour les papiers..." << endl;
    
    gnuplot_script << "set terminal aqua" << endl;
    
    gnuplot_script.close();
//    system("cd output; gnuplot -p 'gnuplot_script.gp'"); // marche plus quand je l'appelle de mon programme alors que ca a eu fonctionné ???
    
    cout << "......" << endl;
    cout << "Result files are in the output directory" << endl;
    cout << "Visualize them with gnuplot by : cd output; gnuplot -p 'gnuplot_script.gp' ; cd .." << endl;
    cout << "The gnuplot script will also save files as png files (in same directory)" << endl;
    system("cd output; gnuplot -p 'gnuplot_script.gp' ; cd ..");
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
        
    }
    outFile_width_ratio << 0 << "\t" << 1.0 << endl;
    
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
        outFile_inner[i].close();
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

