/* ============================================================================
 File   : utils.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
============================================================================ */

#include <assert.h>
#include <math.h>
#include <cstring>

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

bool nn_analysis = false; // whether there is a network in the loop or not

int sysdim; // dimension of system of ODE/DDE
int inputsdim; // dimension of the uncertain inputs and parameters of the system
int fullinputsdim; // full dimension of the uncertain inputs and parameters of the system: taking into account variable inputs
int jacdim;  //  Jacobian will be dimension sysdim * jacdim, for ODEs jacdim = sysdim + fullinputsdim
int paramsdim; // dimension of the vector of parameters params that do not appear in Jacobian
int nncontroldim;  // dimension of the neural network control - does not appear in Jacobian

bool create_png = 0;

int nb_sample_per_dim = 20; // for range estimation by sampling: # of samples per dimension

int points_per_graph = 50;
int printing_period = 1;

int interactive_visualization = 0; // 0 or 1
vector<bool> variables_to_display;



/*----- begin parameters for discrete systems -----*/
int nb_steps = 1; // nb of discrete time steps
int AEextension_order = 1;
int iter_method = 1;
bool skewing = true;
/*----- end parameters for discrete systems -----*/


const int LINESZ = 2048;
char sfx_filename[LINESZ]={0};
char onnx_filename[LINESZ]={0};



// for subdivisions of the initial domain to refine precision
int nb_subdiv_init = 1; // number of subdivisiions
int component_to_subdiv = -1;
int component_to_subdiv2 = -1;
double recovering = 0.0; // percentage of recovering between subdivisions

vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xinner_print, Xinner_joint_print, Xinner_robust_print, Xexact_print; // store results of subdivision
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





using namespace std;


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



void read_system(int argc, char* argv[])
{
  //  char sfx_filename[LINESZ]={0};
  //  char onnx_filename[LINESZ]={0};
    
    systype = 1; // EDO = 0, DDE = 1, discrete systems = 2, DNN = 3
    syschoice = 1;
    nb_sample_per_dim = 20;
    
    char * str_systype = getCmdOption(argv, argv + argc, "-systype");
    if (str_systype)
    {
        if (strcmp(str_systype,"ode")==0)
            systype = 0;
        else if (strcmp(str_systype,"dde")==0)
            systype = 1;
        else if (strcmp(str_systype,"discrete")==0)
            systype = 2;
    }
    
    char * str_syschoice = getCmdOption(argv, argv + argc, "-syschoice");
    if (str_syschoice)
        syschoice = atoi(str_syschoice);
  //  cout << "syschoice= " << syschoice << endl;
    
    
    
    
    char* sfx_filename_temp = getCmdOption(argv, argv + argc, "-nnfile-sfx");
    char* onnx_filename_temp = getCmdOption(argv, argv + argc, "-nnfile-onnx");
    
    char * config_filename = getCmdOption(argv, argv + argc, "-configfile");
    
    if (config_filename)
        readfromfile_syschoice(config_filename,sfx_filename,onnx_filename);
    
    if ((!str_systype && !config_filename) || cmdOptionExists(argv, argv + argc, "-help")  || cmdOptionExists(argv, argv + argc, "-h"))
    {
        cout << "Usage is: " << endl;
        cout << "-systype xxx: ode for ODE, dde for DDE, discrete for discrete-time systems" << endl;
        cout << "-nnfile-sfx xxx or -nnfile-onnx xxx: in case systype is nn, then you should enter the name of file that contains the network, either in old sherlock-like format (sfx) or onnx format" << endl;
        cout << "-syschoice x: integer setting the predefined system ID" << endl;
        cout << "-nbsteps x: for discrete systems only, (optional) number of steps - default value is 1" << endl;
        cout << "-AEextension_order x: for discrete systems only, (optional) order of AE-extensions (1 or 2) - default value is 1" << endl;
        cout << "-iter_method x: for discrete systems only, (optional) choice of iterating method (1 or 2) - default value is 1" << endl;
        cout << "-skew x: for discrete systems only, (optional) skewing/preconditioning for joint range (0 is false or 1 is true) - default value is 1" << endl;
        cout << "-configfile xxx: name of the (optional) config file " << endl;
        exit(0);
    }
    
    
    // parsing neural network
    if (sfx_filename && (sfx_filename[0] != 0)) {// read from config file
        cout << "reading network from" << sfx_filename << endl;
        NH = network_handler(sfx_filename);
        nn_analysis = true;
    }
    else if (sfx_filename_temp) { // read from command-line
        cout << "reading network from" << sfx_filename_temp << endl;
        NH = network_handler(sfx_filename_temp);
        nn_analysis = true;
    }
    else if (onnx_filename && (onnx_filename[0] != 0)) // read from config file
    {
#if ONNX_active
        onnx_parser my_parser(onnx_filename);
        map<string, ParameterValues <uint32_t> > tensor_mapping;
        my_parser.build_graph(CG, tensor_mapping);
        nn_analysis = true;
#endif
    }
    else if (onnx_filename_temp) // read from command-line
    {
#if ONNX_active
        onnx_parser my_parser(onnx_filename_temp);
        map<string, ParameterValues <uint32_t> > tensor_mapping;
        my_parser.build_graph(CG, tensor_mapping);
        nn_analysis = true;
#endif
    }
    
    
    /*----- begin parameters for discrete systems -----*/
    char * str_nbsteps = getCmdOption(argv, argv + argc, "-nbsteps");
    if (str_nbsteps)
        nb_steps = atoi(str_nbsteps);
    
    char * str_AEextension_order = getCmdOption(argv, argv + argc, "-AEextension_order");
    if (str_AEextension_order)
        AEextension_order = atoi(str_AEextension_order);
    
    char * str_method = getCmdOption(argv, argv + argc, "-iter_method");
    if (str_method)
        iter_method = atoi(str_method); //
    char * str_skew = getCmdOption(argv, argv + argc, "-skew");
    if (str_skew)
        skewing = atoi(str_skew); //
    
    /*----- end parameters for discrete systems -----*/

    
    
}



// reading system choice and neural network from file
void readfromfile_syschoice(const char * params_filename, char* sfx_filename, char* onnx_filename)
{
    const int LINESZ = 2048;
    char buff[LINESZ];
    char str_systype[LINESZ];
    
    int int_create_png;
    
    cout << "****** Reading system choice from file " <<  params_filename << " ******" << endl;
    FILE *params_file = fopen(params_filename,"r");
    if (params_file == NULL)
        cout << "Error reading " << params_filename << ": file not found" << endl;
    while (fgets(buff,LINESZ,params_file)) {
        sscanf(buff, "systype = %s\n", str_systype);
        sscanf(buff, "syschoice = %d\n", &syschoice);
        sscanf(buff, "nnfile-sfx = %s\n", sfx_filename);
        sscanf(buff, "nnfile-onnx = %s\n", onnx_filename);
        sscanf(buff, "nn-offset = %lf\n", &nn_offset);
        sscanf(buff, "nn-scaling = %lf\n", &nn_scaling_factor);
        sscanf(buff, "samples-per-dim = %d\n", &nb_sample_per_dim);
        sscanf(buff, "points-per-graph = %d\n", &points_per_graph);
        sscanf(buff, "create-png = %d\n", &int_create_png);
        sscanf(buff, "nb-initial-subdivisions = %d\n", &nb_subdiv_init);
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
    if (nn_scaling_factor)
        cout << "nn_scaling_factor =" << nn_scaling_factor;
    create_png = (int_create_png == 1);
    
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
   
}

/*
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
 
 */


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
    
    if (uncontrolled > 0) {
        out_approx << YAML::Key << "innerrobust";
        out_approx << YAML::Value << temp; //
        out_approx << YAML::Key << "outerrobust";
        out_approx << YAML::Value << temp; //
    }
    
    /*
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
    
    
if (sysdim >= 2)
{
    vector<double> temp2d(4);
    vector<double> temp3d(6);
    vector<double> tempskew(8);
    
    out_approx << YAML::Key << "outer2d";
    out_approx << YAML::Value << YAML::BeginSeq;
    

        
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            
            
            out_approx << YAML::BeginMap;
            
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << i;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << j;
            
            temp2d[0] = inf(x[i]); temp2d[1] = sup(x[i]); temp2d[2] = inf(x[j]); temp2d[3] = sup(x[j]);
            
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
            
            if (uncontrolled > 0) {
                out_approx << YAML::Key << "robskew";
                out_approx << YAML::Value << tempskew;
            }
            
            out_approx << YAML::EndMap;
        }
    }
    
    out_approx << YAML::EndSeq;
    
    
    out_approx << YAML::Key << "inner2d";
    out_approx << YAML::Value << YAML::BeginSeq;
    
    
        
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            
            
            out_approx << YAML::BeginMap;
            
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << i;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << j;
            
            temp2d[0] = inf(x[i]); temp2d[1] = sup(x[i]); temp2d[2] = inf(x[j]); temp2d[3] = sup(x[j]);
            
            out_approx << YAML::Key << "maxbox";
            out_approx << YAML::Value << temp2d;
            
        
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
            
            
            if (uncontrolled > 0) {
                out_approx << YAML::Key << "robbox";
                out_approx << YAML::Value << temp2d;
                
                out_approx << YAML::Key << "robskew";
                out_approx << YAML::Value << tempskew;
            }
            
            out_approx << YAML::EndMap;
        }
    }
    
    out_approx << YAML::EndSeq;
}
    
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
         //       if (uncontrolled > 0 || controlled > 0)
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
        Xinner_print[current_subdiv][0][i] = range_x;
        Xinner_robust_print[current_subdiv][0][i] = range_x;
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
    
    vector<double> temp2d(4);
    vector<double> temp3d(6);
    vector<double> tempskew(8);
    
    out_approx << YAML::Key << "outer2d";
    out_approx << YAML::Value << YAML::BeginSeq;
    

        
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            
            
            out_approx << YAML::BeginMap;
            
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << i;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << j;
            
            temp2d[0] = inf(x[i].convert_int()); temp2d[1] = sup(x[i].convert_int()); temp2d[2] = inf(x[j].convert_int()); temp2d[3] = sup(x[j].convert_int());
            
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
            
            if (uncontrolled > 0) {
                out_approx << YAML::Key << "robskew";
                out_approx << YAML::Value << tempskew;
            }
            
            out_approx << YAML::EndMap;
        }
    }
    
    out_approx << YAML::EndSeq;
            
            

    out_approx << YAML::Key << "inner2d";
    out_approx << YAML::Value << YAML::BeginSeq;
    

        
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=i+1 ; j < sysdim ; j++) {
            
            
            out_approx << YAML::BeginMap;
            
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << i;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << j;
            
            out_approx << YAML::Key << "maxbox";
            temp2d[0] = inf(x[i].convert_int()); temp2d[1] = sup(x[i].convert_int()); temp2d[2] = inf(x[j].convert_int()); temp2d[3] = sup(x[j].convert_int());
            out_approx << YAML::Value << temp2d;
            
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
    cout << "Visualize them with by : cd GUI; python3 Visu_output.py --interactive=0 --printvar=all; cd .." << endl;
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
    system("pwd");
    if (create_png)
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
        
        
    //    print_ErrorMeasures(current_iteration,d0);
        
        out_approx << YAML::EndMap;
    }
    
    if (no_hole)
        cout << "NO HOLE when joining the inner-approx tubes";
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
