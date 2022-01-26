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

ofstream approxreachsetfile;
YAML::Emitter out_approx;



int systype; // 0 is ODE, 1 is DDE
int syschoice; // choice of system to analyze


int interactive_visualization = 0; // 0 or 1
vector<bool> variables_to_display;

using namespace std;


// reading system choice and neural network from file
void readfromfile_syschoice(const char * params_filename, char* sfx_filename, char* onnx_filename)
{
    const int LINESZ = 2048;
    char buff[LINESZ];
    char str_systype[LINESZ];
    
    cout << "****** Reading system choice from file " <<  params_filename << " ******" << endl;
    FILE *params_file = fopen(params_filename,"r");
    if (params_file == NULL)
        cout << "Error reading " << params_filename << ": file not found" << endl;
    while (fgets(buff,LINESZ,params_file)) {
        sscanf(buff, "systype = %s\n", str_systype);
        sscanf(buff, "syschoice = %d\n", &syschoice);
        sscanf(buff, "nnfile-sfx = %s\n", sfx_filename);
        sscanf(buff, "nnfile-onnx = %s\n", onnx_filename);
    }
    
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
        // cout << "systype in readfromfile= " << systype << endl;
    }
    if (sfx_filename)
        cout << "sfx =" << sfx_filename;
    
    fclose(params_file);
}


void open_outputfiles()
{
    system("mv output output_sauv");
    system("rm -r output");
    system("mkdir output");
    
    approxreachsetfile.open("output/approxreachset.yaml");
    
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "approx";
    out_approx << YAML::Value << YAML::BeginSeq;
 
    char file_name[1028];
    
    

}

// print after the end of the analysis
void print_finalstats(clock_t begin)
{
    clock_t end = clock();
    // double end_time = getTime ( );
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; //getTotalTime (start_time , end_time );
    cout << "elapsed time (sec) =" << elapsed_secs << endl;
    
    out_approx << YAML::EndSeq;
    out_approx << YAML::EndMap;
    approxreachsetfile << out_approx.c_str();
    approxreachsetfile.close();
    
    
}


// print initial conditions and init XML in the discrete systems case
void print_init_discrete(vector<interval> &x, bool skew)
{
    
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "tn";
    out_approx << YAML::Value << 0;
    
    vector<double> temp(2*sysdim);
    
    
    for (int i=0 ; i<sysdim ; i++) {
        temp[2*i] = x[i].inf();
        temp[2*i+1] = x[i].sup();
    }
        
    out_approx << YAML::Key << "outer";
    out_approx << YAML::Value << temp; //
    out_approx << YAML::Key << "inner";
    out_approx << YAML::Value << temp; //
    
/*    if (uncontrolled > 0) {
        out_approx << YAML::Key << "innerrobust";
        out_approx << YAML::Value << temp; //
        out_approx << YAML::Key << "outerrobust";
        out_approx << YAML::Value << temp; //
    }
    if (controlled > 0 || uncontrolled > 0) {
        out_approx << YAML::Key << "innerminimal";
        out_approx << YAML::Value << temp;
        out_approx << YAML::Key << "outerminimal";
        out_approx << YAML::Value << temp;
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        range_x = x[i].convert_int();
        temp[2*i] = range_x.mid();
        temp[2*i+1] = temp[2*i];
    }
    out_approx << YAML::Key << "center";
    out_approx << YAML::Value << temp;
   */
    
    // error measures
    vector<double> temp2(sysdim);
    for (int i=0 ; i<sysdim ; i++)
        temp2[i] = 1.0;
    out_approx << YAML::Key << "etaouter";
    out_approx << YAML::Value << temp2;
    out_approx << YAML::Key << "etainner";
    out_approx << YAML::Value << temp2;
    out_approx << YAML::Key << "gamma";
    out_approx << YAML::Value << temp2;
    
    

    out_approx << YAML::Key << "inner2d";
    out_approx << YAML::Value << YAML::BeginSeq;
    
    vector<double> temp2d(4);
    vector<double> temp3d(6);
    vector<double> tempskew(8);
        
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            
            
            out_approx << YAML::BeginMap;
            
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << i;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << j;
            
            temp2d[0] = inf(x[i]); temp2d[1] = sup(x[i]); temp2d[2] = inf(x[j]); temp2d[3] = sup(x[j]);
            if (!skew) {
                out_approx << YAML::Key << "maxbox";
                out_approx << YAML::Value << temp2d;
            }
            else {
        
                tempskew[0] = temp2d[0];
                tempskew[1] = temp2d[2];
                tempskew[2] = temp2d[0];
                tempskew[3] = temp2d[3];
                tempskew[4] = temp2d[1];
                tempskew[5] = temp2d[3];
                tempskew[6] = temp2d[1];
                tempskew[7] = temp2d[2];
            
                out_approx << YAML::Key << "maxskew";
                out_approx << YAML::Value << tempskew;
            }
            
            /*    if (uncontrolled > 0 || controlled > 0) {
                    out_approx << YAML::Key << "minbox";
                    out_approx << YAML::Value << temp2d;
                }
                if (uncontrolled > 0) {
                    out_approx << YAML::Key << "robbox";
                    out_approx << YAML::Value << temp2d;
                }
                */
            
            
      /*      if (uncontrolled > 0 || controlled > 0) {
                out_approx << YAML::Key << "minskew";
                out_approx << YAML::Value << tempskew;
            }
            if (uncontrolled > 0) {
                out_approx << YAML::Key << "robskew";
                out_approx << YAML::Value << tempskew;
            }
            */
            out_approx << YAML::EndMap;
        }
    }
    
    out_approx << YAML::EndSeq;
    
    /* TODO. A AJOUTER PLUS TARD
    
    out_approx << YAML::Key << "inner3d";
    out_approx << YAML::Value << YAML::BeginSeq;
    
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            for (int k=j+1 ; k < sysdim ; k++) {
               
                out_approx << YAML::BeginMap;
                
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << i;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << j;
                out_approx << YAML::Key << "x3";
                out_approx << YAML::Value << k;
                
                temp3d[0] = inf(x[i]); temp3d[1] = sup(x[i]);
                temp3d[2] = inf(x[j]); temp3d[3] = sup(x[j]);
                temp3d[4] = inf(x[k]); temp3d[5] = sup(x[k]);
                
                out_approx << YAML::Key << "maxbox";
                out_approx << YAML::Value << temp3d;
         //       if (uncontrolled > 0 || controlled > 0)
         //       {
         //           out_approx << YAML::Key << "minbox";
         //           out_approx << YAML::Value << temp3d;
         //       }
         //       if (uncontrolled > 0) {
         //           out_approx << YAML::Key << "robbox";
         //           out_approx << YAML::Value << temp3d;
         //       }
         //
                out_approx << YAML::EndMap;

            }
        }
    }
    
    out_approx << YAML::EndSeq;
    
     FIN A AJOUTER PLUS TARD */
     
     
    out_approx << YAML::EndMap;
    
     
    
    
 /*   cout << "At t=0 :" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x[" << i <<"]=" << x[i].convert_int() << "\t";
    cout << endl;
   */
     
}

// print initial conditions and init XML in the ODE case
void print_initstats(vector<AAF> &x, vector<AAF> &param_inputs)
{
    interval range_x;
    
    
    // print initial conditions of the ODE
    
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "tn";
    out_approx << YAML::Value << 0;
    
    vector<double> temp(2*sysdim);
    
    
    for (int i=0 ; i<sysdim ; i++) {
        range_x = x[i].convert_int();
        temp[2*i] = range_x.inf();
        temp[2*i+1] = range_x.sup();
        
        Xouter_print[current_subdiv][0][i] = range_x;
        Xouter_robust_print[current_subdiv][0][i] = range_x;
        Xouter_minimal_print[current_subdiv][0][i] = range_x;
        Xinner_print[current_subdiv][0][i] = range_x;
        Xinner_robust_print[current_subdiv][0][i] = range_x;
        Xinner_minimal_print[current_subdiv][0][i] = range_x;
        Xexact_print[current_subdiv][0][i] = range_x;
    }
        
    out_approx << YAML::Key << "outer";
    out_approx << YAML::Value << temp; //
    out_approx << YAML::Key << "inner";
    out_approx << YAML::Value << temp; //
    
    if (uncontrolled > 0) {
        out_approx << YAML::Key << "innerrobust";
        out_approx << YAML::Value << temp; //
        out_approx << YAML::Key << "outerrobust";
        out_approx << YAML::Value << temp; //
    }
    if (controlled > 0 || uncontrolled > 0) {
        out_approx << YAML::Key << "innerminimal";
        out_approx << YAML::Value << temp;
        out_approx << YAML::Key << "outerminimal";
        out_approx << YAML::Value << temp;
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        range_x = x[i].convert_int();
        temp[2*i] = range_x.mid();
        temp[2*i+1] = temp[2*i];
    }
    out_approx << YAML::Key << "center";
    out_approx << YAML::Value << temp;
    
    // error measures
    vector<double> temp2(sysdim);
    for (int i=0 ; i<sysdim ; i++)
        temp2[i] = 1.0;
    out_approx << YAML::Key << "etaouter";
    out_approx << YAML::Value << temp2;
    out_approx << YAML::Key << "etainner";
    out_approx << YAML::Value << temp2;
    out_approx << YAML::Key << "gamma";
    out_approx << YAML::Value << temp2;
    
    

    out_approx << YAML::Key << "inner2d";
    out_approx << YAML::Value << YAML::BeginSeq;
    
    vector<double> temp2d(4);
    vector<double> temp3d(6);
    vector<double> tempskew(8);
        
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            
            
            out_approx << YAML::BeginMap;
            
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << i;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << j;
            
            out_approx << YAML::Key << "maxbox";
            temp2d[0] = inf(x[i].convert_int()); temp2d[1] = sup(x[i].convert_int()); temp2d[2] = inf(x[j].convert_int()); temp2d[3] = sup(x[j].convert_int());
            out_approx << YAML::Value << temp2d;
            
            if (uncontrolled > 0 || controlled > 0) {
                out_approx << YAML::Key << "minbox";
                out_approx << YAML::Value << temp2d;
            }
            if (uncontrolled > 0) {
                out_approx << YAML::Key << "robbox";
                out_approx << YAML::Value << temp2d;
            }
            
            tempskew[0] = temp2d[0];
            tempskew[1] = temp2d[2];
            tempskew[2] = temp2d[0];
            tempskew[3] = temp2d[3];
            tempskew[4] = temp2d[1];
            tempskew[5] = temp2d[3];
            tempskew[6] = temp2d[1];
            tempskew[7] = temp2d[2];
            
            out_approx << YAML::Key << "maxskew";
            out_approx << YAML::Value << tempskew;
            
            if (uncontrolled > 0 || controlled > 0) {
                out_approx << YAML::Key << "minskew";
                out_approx << YAML::Value << tempskew;
            }
            if (uncontrolled > 0) {
                out_approx << YAML::Key << "robskew";
                out_approx << YAML::Value << tempskew;
            }
            
            out_approx << YAML::EndMap;
        }
    }
    
    out_approx << YAML::EndSeq;
    
    out_approx << YAML::Key << "inner3d";
    out_approx << YAML::Value << YAML::BeginSeq;
    
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            for (int k=j+1 ; k < sysdim ; k++) {
               
                out_approx << YAML::BeginMap;
                
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << i;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << j;
                out_approx << YAML::Key << "x3";
                out_approx << YAML::Value << k;
                
                temp3d[0] = inf(x[i].convert_int()); temp3d[1] = sup(x[i].convert_int());
                temp3d[2] = inf(x[j].convert_int()); temp3d[3] = sup(x[j].convert_int());
                temp3d[4] = inf(x[k].convert_int()); temp3d[5] = sup(x[k].convert_int());
                
                out_approx << YAML::Key << "maxbox";
                out_approx << YAML::Value << temp3d;
                if (uncontrolled > 0 || controlled > 0)
                {
                    out_approx << YAML::Key << "minbox";
                    out_approx << YAML::Value << temp3d;
                }
                if (uncontrolled > 0) {
                    out_approx << YAML::Key << "robbox";
                    out_approx << YAML::Value << temp3d;
                }
                
                out_approx << YAML::EndMap;

            }
        }
    }
    
    out_approx << YAML::EndSeq;
    
    out_approx << YAML::EndMap;
    
    // a changer un jour pour t_begin (notamment pour DDE)?
    t_print[0] = 0;
    
    
    cout << "At t=0 :" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x[" << i <<"]=" << x[i].convert_int() << "\t";
    cout << endl;
    for (int i=0 ; i<fullinputsdim ; i++)
        cout << "param_inputs[" << i <<"]=" << param_inputs[i].convert_int() << "\t";
    cout << endl;
     
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
    if (systype == 2) // discrete systems
        sprintf(displayed_variables,"%s","all");
    else {
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
    }
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
         //   if (nb_subdiv_init > 1)
         //       outFile_outer[i] << t_print[current_iteration] << "\t" << inf(Xouter_print[0][current_iteration][i]) << "\t" << sup(Xouter_print[0][current_iteration][i]) << endl;
            
            //  outFile_inner[i] << t_print[current_iteration] << "\t" << inf(Xinner_print[0][current_iteration][i]) << "\t" << sup(Xinner_print[0][current_iteration][i]) << endl;
        }
        
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "tn";
        out_approx << YAML::Value << t_print[current_iteration];
        out_approx << YAML::Key << "outer";
        vector<double> temp(2*sysdim);
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xouter_print[0][current_iteration][i].inf();
            temp[2*i+1] = Xouter_print[0][current_iteration][i].sup();
        }
        out_approx << YAML::Value << temp;
        out_approx << YAML::Key << "currentsubdiv";
        out_approx << YAML::Value << 0; // to indicate union of all subdivisions
        
        
        print_ErrorMeasures(current_iteration,d0);
        
        out_approx << YAML::EndMap;
    }
    
    if (no_hole)
        cout << "NO HOLE when joining the inner-approx tubes";
}

void print_ErrorMeasures(int current_iteration, double d0)
{
    double aux, sum, rel_sum;
    //  vector<interval> Xexact(sysdim);

    out_approx << YAML::Key << "tn";
    out_approx << YAML::Value << t_print[current_iteration];
    
    // fills Xexact_print with analytical solution
    AnalyticalSol(current_iteration, d0);
    // cout << "before testing x_exact, current_iteration=" << current_iteration << "t_print[current_iteration] " << t_print[current_iteration]  << endl;
    if (sup(Xexact_print[0][current_iteration][0]) >= inf(Xexact_print[0][current_iteration][0])) // an analytical solution exists
    {
        out_approx << YAML::Key << "exact";
        vector<double> temp(2*sysdim);
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xexact_print[0][current_iteration][i].inf();
            temp[2*i+1] = Xexact_print[0][current_iteration][i].sup();
        }
        out_approx << YAML::Value << temp;
        
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
        out_approx << YAML::Key << "meanerrorouter"; // mean on xi of error between outer-approx and analytical solution if any
        out_approx << YAML::Value << sum;
        out_approx << YAML::Key << "relmeanerrorouter"; // mean on xi of error between outer-approx and analytical solution if any, over width of exact tube
        out_approx << YAML::Value << rel_sum;
        
        
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
        out_approx << YAML::Key << "meanerrorinner";  // mean on xi of error between inner-approx and analytical solution if any
        out_approx << YAML::Value << sum;
        out_approx << YAML::Key << "relmeanerrorinner"; // mean on xi of error between inner-approx and analytical solution if any, over width of exact tube
        out_approx << YAML::Value << rel_sum;
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
    
    out_approx << YAML::Key << "meanerrordiff";  // mean on xi of error between outer-approx and inner-approx
    out_approx << YAML::Value << sum;
    out_approx << YAML::Key << "relmeanerrordiff"; // mean on xi of error between outer-approx and inner-approx, over width of over-approx tube
    out_approx << YAML::Value << rel_sum;
    
}


std::ostream& operator<<(std::ostream& os, const std::vector<double> &input)
{
    for (auto const& i: input) {
        os << i << " ";
    }
    return os;
}

/*std::ostream& operator<<(std::ostream& os, const std::vector<interval> &input)
{
    interval temp;
    for (auto const& i: input) {
        //temp = i.convert_int();
        os << "[" << i.inf() << ", " << i.sup() << "] ";
       // os << i.convert_int() << " ";
    }
    return os;
}*/


/*std::ostream& operator<<(std::ostream& os, const interval &input)
{
    os.precision(6);
    os << "[" << input.inf() << ", " << input.sup() << "] ";
}*/


std::ostream& operator<<(std::ostream& os, const std::vector<AAF> &input)
{
    interval temp;
    for (auto const& i: input) {
        temp = i.convert_int();
        os << "[" << temp.inf() << ", " << temp.sup() << "] ";
       // os << i.convert_int() << " ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<T<AAF>> &input)
{
    interval temp;
    for (auto const& i: input) {
        temp = i.val().convert_int();
        os << "[" << temp.inf() << ", " << temp.sup() << "] ";
      //  os << i.val().convert_int() << " ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<T<F<AAF>>> &input)
{
    
    
    
 //   for (auto const& i: input) {
 //       os << i.x().convert_int() << " ";
 //   }
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
