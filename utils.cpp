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
vector<ofstream> outFile_inner_robust;   // robust inner-approximated range for each variable of the system
vector<ofstream> outFile_center;

ofstream outFile_width_ratio;     //  min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
ofstream outFile_meanerror_outer; // mean on xi of error between outer-approx and analytical solution if any
ofstream outFile_meanerror_inner; // mean on xi of error between inner-approx and analytical solution if any
ofstream outFile_meanerror_diff;  // mean on xi of error between outer-approx and inner-approx
ofstream outFile_relmeanerror_outer; // mean on xi of error between outer-approx and analytical solution if any, over width of exact tube
ofstream outFile_relmeanerror_inner; // mean on xi of error between inner-approx and analytical solution if any, over width of exact tube
ofstream outFile_relmeanerror_diff;  // mean on xi of error between outer-approx and inner-approx, over width of over-approx tube




using namespace std;

void open_outputfiles()
{
    system("mv output output_sauv");
    system("rm -r output");
    system("mkdir output");
    
    outFile_outer = vector<ofstream>(sysdim);   // output outer-approximated range for each variable of the system
    if (uncontrolled > 0) {
        outFile_outer_robust = vector<ofstream>(sysdim);
        outFile_inner_robust = vector<ofstream>(sysdim);
    }
    outFile_inner = vector<ofstream>(sysdim); // vector<ofstream> outFile_inner(sysdim);   // output inner-approximated range for each variable of the system
    outFile_outer_minimal = vector<ofstream>(sysdim);
    outFile_inner_minimal = vector<ofstream>(sysdim);
    outFile_center  = vector<ofstream>(sysdim);
    
    stringstream file_name;
    
    for (int i=0 ; i<sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer.out";
        outFile_outer[i].open(file_name.str().c_str());
        if (uncontrolled > 0) {
            file_name.str("");
            file_name << "output/x" << i+1 << "outer_robust.out";
            outFile_outer_robust[i].open(file_name.str().c_str());
            file_name.str("");
            file_name << "output/x" << i+1 << "inner_robust.out";
            outFile_inner_robust[i].open(file_name.str().c_str());
        }
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_minimal.out";
        outFile_outer_minimal[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "center.out";
        outFile_center[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_inner[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "inner_minimal.out";
        outFile_inner_minimal[i].open(file_name.str().c_str());
    }
    
    
    outFile_width_ratio.open("output/width_ratio.out");
    outFile_meanerror_outer.open("output/meanerror_outer.out");
    outFile_meanerror_inner.open("output/meanerror_inner.out");
    outFile_meanerror_diff.open("output/meanerror_diff.out");
    outFile_relmeanerror_outer.open("output/relmeanerror_outer.out");
    outFile_relmeanerror_inner.open("output/relmeanerror_inner.out");
    outFile_relmeanerror_diff.open("output/relmeanerror_diff.out");
}

