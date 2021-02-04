/* ============================================================================
 File   : utils.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
============================================================================ */

#include <assert.h>
#include <math.h>

//#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_gamma.h>


//#include "filib_interval.h"
//#include "tadiff.h" 
//#include "fadiff.h"
#include "fadbad_aa.h"
#include "utils.h"
#include "ode_def.h"
//#include "matrix.h"

vector<ofstream> outFile_outer_minimal;   //  minimal outer-approximated range for each variable of the system
vector<ofstream> outFile_outer_robust;  // robust outer-approximated range for each variable of the system
vector<ofstream> outFile_outer;   //  maximal outer-approximated range for each variable of the system
vector<ofstream> outFile_inner_minimal;   //  minimal inner-approximated range for each variable of the system
vector<ofstream> outFile_inner;   //  maximal inner-approximated range for each variable of the system

vector<vector<ofstream>> outFile_joint_inner;   // output inner-approximated range for each couple of variables of the system
vector<vector<vector<ofstream>>> outFile_joint_inner3d;   // output inner-approximated range for each triple of variables of the system

vector<ofstream> outFile_inner_robust;   // robust inner-approximated range for each variable of the system
vector<ofstream> outFile_center;
vector<ofstream> outFile_exact; // analytical solution if any

ofstream outFile_width_ratio;     //  min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
ofstream outFile_meanerror_outer; // mean on xi of error between outer-approx and analytical solution if any
ofstream outFile_meanerror_inner; // mean on xi of error between inner-approx and analytical solution if any
ofstream outFile_meanerror_diff;  // mean on xi of error between outer-approx and inner-approx
ofstream outFile_relmeanerror_outer; // mean on xi of error between outer-approx and analytical solution if any, over width of exact tube
ofstream outFile_relmeanerror_inner; // mean on xi of error between inner-approx and analytical solution if any, over width of exact tube
ofstream outFile_relmeanerror_diff;  // mean on xi of error between outer-approx and inner-approx, over width of over-approx tube

int systype; // 0 is ODE, 1 is DDE
int syschoice; // choice of system to analyze


int interactive_visualization = 0; // 0 or 1
vector<bool> variables_to_display;

using namespace std;

void open_outputfiles()
{
    system("mv output output_sauv");
    system("rm -r output");
    system("mkdir output");
    
    outFile_outer = vector<ofstream>(sysdim);   // output outer-approximated range for each variable of the system
    outFile_exact = vector<ofstream>(sysdim);
    if (uncontrolled > 0) {
        outFile_outer_robust = vector<ofstream>(sysdim);
        outFile_inner_robust = vector<ofstream>(sysdim);
    }
    outFile_inner = vector<ofstream>(sysdim); // vector<ofstream> outFile_inner(sysdim);   // output inner-approximated range for each variable of the system
    if (uncontrolled > 0 || controlled > 0) {
        outFile_outer_minimal = vector<ofstream>(sysdim);
        outFile_inner_minimal = vector<ofstream>(sysdim);
    }
    outFile_center  = vector<ofstream>(sysdim);
    
    outFile_joint_inner = vector<vector<ofstream>>(sysdim);
    for (int i=0 ; i<sysdim ; i++)
        outFile_joint_inner[i] = vector<ofstream>(sysdim);
    
    outFile_joint_inner3d = vector<vector<vector<ofstream>>>(sysdim);
    for (int i=0 ; i<sysdim ; i++) {
        outFile_joint_inner3d[i] = vector<vector<ofstream>>(sysdim);
        for (int j=0 ; j<sysdim ; j++)
            outFile_joint_inner3d[i][j] = vector<ofstream>(sysdim);
    }
    
 
    char file_name[1028];
    
    for (int i=0 ; i<sysdim ; i++) {
        sprintf(file_name,"output/x%douter.out",i+1);
        outFile_outer[i].open(file_name);
        sprintf(file_name,"output/x%dexact.out",i+1);
        outFile_exact[i].open(file_name);
        if (uncontrolled > 0) {
            sprintf(file_name,"output/x%douter_robust.out",i+1);
            outFile_outer_robust[i].open(file_name);
            sprintf(file_name,"output/x%dinner_robust.out",i+1);
            outFile_inner_robust[i].open(file_name);
        }
        if (uncontrolled > 0 || controlled > 0) {
            sprintf(file_name,"output/x%douter_minimal.out",i+1);
            outFile_outer_minimal[i].open(file_name);
            sprintf(file_name,"output/x%dinner_minimal.out",i+1);
            outFile_inner_minimal[i].open(file_name);
        }
        sprintf(file_name,"output/x%dcenter.out",i+1);
        outFile_center[i].open(file_name);
        sprintf(file_name,"output/x%dinner.out",i+1);
        outFile_inner[i].open(file_name);
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            sprintf(file_name,"output/x%dx%dinner_joint.out",i+1,j+1);
            outFile_joint_inner[i][j].open(file_name);
        }
    }
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            for (int k=j+1 ; k < sysdim ; k++) {
                sprintf(file_name,"output/x%dx%dx*dinner_joint3d.out",i+1,j+1,k+1);
                outFile_joint_inner3d[i][j][k].open(file_name);
            }
        }
    }
  
 
    

  
 /*
    for (int i=0 ; i<sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer.out";
        outFile_outer[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x"<<i+1<<"exact.out";
        outFile_exact[i].open(file_name.str().c_str());
        if (uncontrolled > 0) {
            file_name.str("");
            file_name << "output/x" << i+1 << "outer_robust.out";
            outFile_outer_robust[i].open(file_name.str().c_str());
            file_name.str("");
            file_name << "output/x" << i+1 << "inner_robust.out";
            outFile_inner_robust[i].open(file_name.str().c_str());
        }
        if (uncontrolled > 0 || controlled > 0) {
            file_name.str("");
            file_name << "output/x" << i+1 << "outer_minimal.out";
            outFile_outer_minimal[i].open(file_name.str().c_str());
            file_name.str("");
            file_name << "output/x" << i+1 << "inner_minimal.out";
            outFile_inner_minimal[i].open(file_name.str().c_str());
        }
        file_name.str("");
        file_name << "output/x" << i+1 << "center.out";
        outFile_center[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_inner[i].open(file_name.str().c_str());
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            file_name.str("");
            file_name << "output/x" << i+1 << "x" << j+1 << "inner_joint.out";
            outFile_joint_inner[i][j].open(file_name.str().c_str());
        }
    }
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            for (int k=j+1 ; k < sysdim ; k++) {
                file_name.str("");
                file_name << "output/x" << i+1 << "x" << j+1 <<  "x" << k+1 << "inner_joint3d.out";
                outFile_joint_inner3d[i][j][k].open(file_name.str().c_str());
            }
        }
    }
    */
    
    
    outFile_width_ratio.open("output/width_ratio.out");
    outFile_meanerror_outer.open("output/meanerror_outer.out");
    outFile_meanerror_inner.open("output/meanerror_inner.out");
    outFile_meanerror_diff.open("output/meanerror_diff.out");
    outFile_relmeanerror_outer.open("output/relmeanerror_outer.out");
    outFile_relmeanerror_inner.open("output/relmeanerror_inner.out");
    outFile_relmeanerror_diff.open("output/relmeanerror_diff.out");
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
        outFile_center[i].close();
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            outFile_joint_inner[i][j].close();
        }
    }
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            for (int k=j+1 ; k < sysdim ; k++) {
                outFile_joint_inner3d[i][j][k].close();
            }
        }
    }
    
    
    outFile_width_ratio.close();
    outFile_meanerror_outer.close();
    outFile_meanerror_inner.close();
    outFile_meanerror_diff.close();
    outFile_relmeanerror_outer.close();
    outFile_relmeanerror_inner.close();
    outFile_relmeanerror_diff.close();
}



void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs)
{
    interval range_x;
    
    
    // print initial conditions of the ODE
    for (int i=0 ; i<sysdim ; i++) {
        range_x = x[i].convert_int();
        // print in exit files
        outFile_outer[i] << 0 << "\t" << inf(range_x) << "\t" << sup(range_x) << endl;
        if (uncontrolled > 0) {
            outFile_inner_robust[i] << 0 << "\t" << inf(range_x) << "\t" << sup(range_x) << endl;
            outFile_outer_robust[i] << 0 << "\t" << inf(range_x) << "\t" << sup(range_x) << endl;
        }
        if (controlled > 0 || uncontrolled > 0) {
            outFile_inner_minimal[i] << 0 << "\t" << inf(range_x) << "\t" << sup(range_x) << endl;
            outFile_outer_minimal[i] << 0 << "\t" << inf(range_x) << "\t" << sup(range_x) << endl;
        }
        outFile_center[i] << 0<< "\t" << mid(range_x) << "\t" << mid(range_x) << endl;
        outFile_inner[i] << 0 << "\t" << inf(range_x) << "\t" << sup(range_x) << endl;
        
        Xouter_print[current_subdiv][0][i] = range_x;
        Xouter_robust_print[current_subdiv][0][i] = range_x;
        Xouter_minimal_print[current_subdiv][0][i] = range_x;
        Xinner_print[current_subdiv][0][i] = range_x;
        Xinner_robust_print[current_subdiv][0][i] = range_x;
        Xinner_minimal_print[current_subdiv][0][i] = range_x;
        Xexact_print[current_subdiv][0][i] = range_x;
    }
    outFile_width_ratio << 0 << "\t" << 1.0 << endl;
    
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            outFile_joint_inner[i][j] << "maximal" << "\t" << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << "\t" << inf(x[j].convert_int()) << "\t" << sup(x[j].convert_int()) << endl;
            if (uncontrolled > 0)
                outFile_joint_inner[i][j] << "robust" << "\t" << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << "\t" << inf(x[j].convert_int()) << "\t" << sup(x[j].convert_int()) << endl;
            if (uncontrolled > 0 || controlled > 0)
                outFile_joint_inner[i][j] << "minimal"  << "\t" << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << "\t" << inf(x[j].convert_int()) << "\t" << sup(x[j].convert_int()) << endl;
        }
    }
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            for (int k=j+1 ; k < sysdim ; k++) {
                outFile_joint_inner3d[i][j][k] << "maximal"  << "\t" << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << "\t" << inf(x[j].convert_int()) << "\t" << sup(x[j].convert_int()) << "\t" << inf(x[k].convert_int()) << "\t" << sup(x[k].convert_int()) << endl;
                if (uncontrolled > 0)
                    outFile_joint_inner3d[i][j][k] << "robust" << "\t" << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << "\t" << inf(x[j].convert_int()) << "\t" << sup(x[j].convert_int()) << "\t" << inf(x[k].convert_int()) << "\t" << sup(x[k].convert_int()) << endl;
                if (uncontrolled > 0 || controlled > 0)
                    outFile_joint_inner3d[i][j][k] << "minimal" << "\t" << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << "\t" << inf(x[j].convert_int()) << "\t" << sup(x[j].convert_int()) << "\t" << inf(x[k].convert_int()) << "\t" << sup(x[k].convert_int()) << endl;
            }
        }
    }
    
    
    // a changer un jour pour t_begin (notamment pour DDE)?
    t_print[0] = 0;
    
    
    cout << "At t=0 :" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x[" << i <<"]=" << x[i].convert_int() << "\t";
    cout << endl;
    for (int i=0 ; i<fullinputsdim ; i++)
        cout << "param_inputs[" << i <<"]=" << param_inputs[i].convert_int() << "\t";
    cout << endl;
     
     
    /*
     cout << "x0=" << endl;
     for (int i=0 ; i<sysdim ; i++)
     cout << "x0[" << i <<"]=" << mid(x[i].convert_int()) << "\t";
     cout << endl;
     */
    
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


std::ostream& operator<<(std::ostream& os, const std::vector<double> &input)
{
    for (auto const& i: input) {
        os << i << " ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<AAF> &input)
{
    for (auto const& i: input) {
        os << i.convert_int() << " ";
    }
    return os;
}


std::ostream& operator<<(std::ostream& os, const std::vector<vector<double>> &input)
{
    for (int i=0 ; i<input.size(); i++) {
        for (int j=0; j<input[i].size(); j++) {
            os << input[i][j] << " ";
        }
        os << "\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<vector<vector<double>>> &input)
{
    for (int i=0 ; i<input.size(); i++) {
        for (int j=0; j<input[i].size(); j++) {
            for (int k=0; k<input[i][j].size(); k++) {
                os << input[i][j][k] << " ";
            }
            os << "\n";
        }
        os << "\n";
    }
    return os;
}
