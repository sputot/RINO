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
#include <cstring>

#include <stdlib.h>
#include <algorithm>

using namespace std;

#include "yaml-cpp/yaml.h"

//using namespace tinyxml2;

//using namespace fadbad;

// void read_system(const char * system_filename, char *sys_name);

//int systype; // 0 is ODE, 1 is DDE
//int syschoice; // choice of system to analyze

DiscreteFunc f; // for discrete systems - defined in discrete_system.h
OdeFunc obf;    // for odes: contains the differential system - defined in ode_def.h
DdeFunc bf;    // for ddes: contains the differential system - defined in ode_def.h
DdeJacFunc bbf; // for ddes

void read_system(int argc, char* argv[]);
vector<vector<interval>> estimate_reachset();

int main(int argc, char* argv[])
{
/*    DiscreteFunc f; // for discrete systems - defined in discrete_system.h
    OdeFunc obf;    // for odes: contains the differential system - defined in ode_def.h
    DdeFunc bf;    // for ddes: contains the differential system - defined in ode_def.h
    DdeJacFunc bbf; // for ddes
*/
    
    vector<vector<interval>> sampled_reachset;
    
    ReachSet RS;
   
    clock_t begin, end; //  = clock();
    double elapsed_secs_sampling;
    
    cout.unsetf(ios::floatfield);            // floatfield not set
    cout.precision(6);
    
    
/*    essai_sherlock();
    syschoice = 101;
    function_range(NULL);
    exit(0); */
    
    
    // reading system choices and parameters
    read_system(argc,argv);
   
    open_outputfiles();
    
    
    begin = clock();
    sampled_reachset = estimate_reachset();
    end = clock();
    elapsed_secs_sampling = double(end - begin) / CLOCKS_PER_SEC;
    
    begin = clock();
    
    /*******************************************************************************************/
   /************* Discrete Systems ************/
    if (systype == 2)
    {
        
        print_init_discrete(xinit,skewing);
        if (nb_steps == 1)
            RS = function_range(f,xinit,sampled_reachset);
        else {
            if (iter_method == 1)
                RS = discrete_dynamical(f,xinit,sampled_reachset,AEextension_order,skewing);
            else // (iter_method == 2)
                RS = discrete_dynamical_method2(f,xinit,sampled_reachset,skewing);
        }
       
    }
    else if (systype == 0 || systype == 1) // ODE or DDE
    {
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
            // resets initial_values, fullinputs and eps from saved values and current subdivisiion
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, initial_values_save, fullinputs_save, component_to_subdiv);
                        
            print_initstats(initial_values_aff,fullinputs_aff);  // print initial conditions and init XML
            
            if (systype == 0) // ODEs
            {
                HybridStep_ode cur_step = init_ode(obf,center_initial_values_aff,initial_values_aff,t_begin,tau,Taylor_order);
           
                int iter = 1;
                while (cur_step.tn < t_end-0.0001*t_end)
                {
                    // build Taylor model for Value and Jacobian and deduce guards for each active mode
                    cur_step.TM_build(params_aff,fullinputs_aff,center_fullinputs_aff);
                    RS = cur_step.TM_evalandprint_solutionstep(eps,cur_step.tn+tau,sampled_reachset[iter],current_subdiv);
                    cur_step.init_nextstep(params_aff,fullinputs_aff,tau);
                
                    if ((RS.Xouter[0] == interval::EMPTY()) || (RS.Xouter[0] == interval::ENTIRE())) {
                        printf("Terminated due to too large overestimation.\n");
                        break;
                    }
                    iter++;
                }
            }
            else if (systype == 1) // DDEs
            {
                HybridStep_dde prev_step = HybridStep_dde(bf,bbf,Taylor_order,t_begin,tau,delay,nb_subdiv_delay);
                prev_step.init_dde(sampled_reachset[0]);
                
                HybridStep_dde cur_step = prev_step.init_nextbigstep(tau);
                
                //******** Integration loop ************************
                int iter = nb_subdiv_delay+1;
                while (cur_step.tn+delay <= t_end)
                {
                    // build Taylor model for Value and Jacobian and deduce guards for each active mode
                    for (int j=0 ; j<nb_subdiv_delay ; j++)
                    {
                        cur_step.TM_build(j);  // j-th Taylor model valid on [tn+j*d0,tn+(j+1)*d0]
                        // eval at tn+tau
                        RS = cur_step.TM_evalandprint_solutionstep(j,eps,sampled_reachset[iter]); // cur_step.tn+(j+1)*tau/nb_subdiv);
                        
                        if ((RS.Xouter[0] == interval::EMPTY()) || (RS.Xouter[0] == interval::ENTIRE())) {
                            printf("Terminated due to too large overestimation.\n");
                            break;
                        }
                        
                        if (j < nb_subdiv_delay-1)
                            cur_step.init_nextsmallstep(j);
                        iter++;
                    }
                    cur_step = cur_step.init_nextbigstep(tau);
                }
            }
            // adding a white line separator between subdivisions in the output result (except for maximal outer approx which is computed by union of all subdivisions)
            // removed when passed to yaml but not tested since => check
        }
        
        if (nb_subdiv_init > 1)
            print_finalsolution(ceil((t_end-t_begin)/tau), delay);
    }

  
    
    // pb sur les fonctiosn trigos en formes affines a corriger a l'occasion ci-dessous:
//    vector<AAF> yp(sysdim);
//    yp[0] = -(5.0*cos(x[2].convert_int()) + fullinputs[0].convert_int());        // px' = v.cos(theta) + b1
//    yp[1] = -(5.0*sin(x[2].convert_int()) + fullinputs[1].convert_int());        // py' = v.sin(theta) + b2
//    yp[2] = -(param_inputs[3] + fullinputs[2]);    // theta' = a + b3
//    for (int i=0 ; i<sysdim ; i++)
//        cout << "yp[i]=" <<yp[i].convert_int() << " ";
//    cout << endl;
    
    
    
    out_approx << YAML::EndSeq;
    out_approx << YAML::EndMap;
    approxreachsetfile << out_approx.c_str();
    approxreachsetfile.close();
    
    end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "elapsed time (sec) =" << elapsed_secs << endl;
    
    
    ofstream summaryyamlfile;
    summaryyamlfile.open("output/sumup.yaml");
    YAML::Emitter out_summary;
    
    out_summary << YAML::BeginMap;
    
    out_summary << YAML::Key << "systype";
    out_summary << YAML::Value << systype;
    out_summary << YAML::Key << "sysdim";
    out_summary << YAML::Value << sysdim;
    
    out_summary << YAML::Key << "nb_sample_per_dim";
    out_summary << YAML::Value << nb_sample_per_dim;
    out_summary << YAML::Key << "nbmax_sampled_dim";
    out_summary << YAML::Value << nbmax_sampled_dim;
    
    // discrete systems
    if (systype == 2) {
        out_summary << YAML::Key << "nb_steps";
        out_summary << YAML::Value << nb_steps;
        out_summary << YAML::Key << "skew";
        out_summary << YAML::Value << skewing;
        out_summary << YAML::Key << "iter_method";
        out_summary << YAML::Value << iter_method;
    }
    
    if (systype == 0 || systype == 1) {
 //       out_summary << YAML::Key << "nb_timesteps";
 //       out_summary << YAML::Value << n_steps;
        out_summary << YAML::Key << "control_period";
        out_summary << YAML::Value << control_period;
        out_summary << YAML::Key << "dt";
        out_summary << YAML::Value << tau;
        out_summary << YAML::Key << "nb_subdiv_init";
        out_summary << YAML::Value << nb_subdiv_init;
    }
    
    out_summary << YAML::Key << "zouter";
    vector<double> aff_zouter(sysdim*2);
    for (int i=0; i<sysdim ; i++) {
        aff_zouter[2*i] = RS.Xouter[i].inf();
        aff_zouter[2*i+1] = RS.Xouter[i].sup();
    }
    out_summary << YAML::Value << aff_zouter;
    
    if (uncontrolled > 0) {
        out_summary << YAML::Key << "zouter_rob";
        for (int i=0; i<sysdim ; i++) {
            aff_zouter[2*i] = RS.Xouter_rob[i].inf();
            aff_zouter[2*i+1] = RS.Xouter_rob[i].sup();
        }
        out_summary << YAML::Value << aff_zouter;
    }
    
    out_summary << YAML::Key << "zinner";
    vector<double> aff_zinner(sysdim*2);
    for (int i=0; i<sysdim ; i++) {
        aff_zinner[2*i] = RS.Xinner[i].inf();
        aff_zinner[2*i+1] = RS.Xinner[i].sup();
    }
    out_summary << YAML::Value << aff_zinner;
    
    if (uncontrolled > 0) {
        out_summary << YAML::Key << "zinner_rob";
        for (int i=0; i<sysdim ; i++) {
            aff_zinner[2*i] = RS.Xinner_rob[i].inf();
            aff_zinner[2*i+1] = RS.Xinner_rob[i].sup();
        }
        out_summary << YAML::Value << aff_zinner;
    }
    
    
    out_summary << YAML::Key << "elapsed-secs";
    out_summary << YAML::Value << elapsed_secs;
    out_summary << YAML::Key << "elapsed-secs-sampling";
    out_summary << YAML::Value << elapsed_secs_sampling;
    //out_summary << YAML::Key << "elapsed-printing";
    //out_summary << YAML::Value << elapsed_printing;
    
    out_summary << YAML::EndMap;
    summaryyamlfile << out_summary.c_str();
    summaryyamlfile.close();
    
    
    
    
    
    std::string commandLineStr= "";
    for (int i=0;i<argc;i++) commandLineStr.append(std::string(argv[i]).append(" "));
    
    
    ofstream summaryfile;
    string sumupname = "output/sumup.txt";
    summaryfile.open(sumupname);
    summaryfile << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << std::endl
//                << "                           Summary   " << timestamp << std::endl
                << "                           Summary   " << endl
                << "                      Command-line   " << commandLineStr << std::endl
    << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << std::endl;
//                << "                          Dynamics │ " << system_network.dynamics_to_str() << std::endl
    if (systype == 0) {
        summaryfile << "                           Systype │ " << "ode"  << std::endl;
        summaryfile << "                          Dynamics │ " << syschoice << std::endl;
        summaryfile << "                Initial conditions │ " << initial_values << std::endl;
        summaryfile << "                            Inputs │ ";
        for (int i=0; i<inputsdim; i++)
            summaryfile <<"("<<inputs[i]<<","<<nb_inputs[i]<<") ";
        summaryfile << std::endl;
        summaryfile << "                     Time interval │ " << "[" << t_begin << "," << t_end << "]" << std::endl;
        summaryfile << "                      Control step │ " << control_period << std::endl;
        summaryfile << "                         Time step │ " << tau << std::endl;
        summaryfile << "                  System dimension │ " << sysdim << std::endl;
        summaryfile << "               Taylor Models order │ " << Taylor_order << std::endl;
        if (nb_subdiv_init >1)
            summaryfile << "     #subdivisions of input domain │ " << nb_subdiv_init << std::endl;
        if (sfx_filename)
            summaryfile << "               Neural network file │ " << sfx_filename << std::endl;
        else if (onnx_filename)
            summaryfile << "               Neural network file │ " << onnx_filename << std::endl;
        summaryfile << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;
        
        // le reste: a relire/calculer des fichiers YAML que je vais produire...
//        summaryfile << " Projection of under-approximation │ " << Xinner << std::endl;


    }
    else if (systype == 1) {
        summaryfile << "                           Systype │ " << "dde"  << std::endl;
        summaryfile << "                          Dynamics │ " << syschoice << std::endl;
        summaryfile << "                Initial conditions │ " << initial_values << std::endl;
        summaryfile << "                     Time interval │ " << "[" << t_begin << "," << t_end << "]" << std::endl;
        summaryfile << "                             Delay │ " << delay << std::endl;
        summaryfile << "                         Time step │ " << tau << std::endl;
        summaryfile << "                  System dimension │ " << sysdim << std::endl;
        summaryfile << "               Taylor Models order │ " << Taylor_order << std::endl;
        if (nb_subdiv_init >1)
            summaryfile << "     #subdivisions of input domain │ " << nb_subdiv_init << std::endl;
        summaryfile << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;

    }
    else if (systype == 2) {
        summaryfile << "                           Systype │ " << "discrete"  << std::endl;
        summaryfile << "                          Dynamics │ " << syschoice << std::endl;
        summaryfile << "                Initial conditions │ " << xinit << std::endl;
        summaryfile << "                            #steps │ " << nb_steps << std::endl;
        summaryfile << "                  System dimension │ " << sysdim << std::endl;
        summaryfile << "                  #samples per dim │ " << nb_sample_per_dim << std::endl;
        summaryfile << "                 AEextension order │ " << AEextension_order << std::endl;
        summaryfile << "                  Iteration method │ " << iter_method << std::endl;
        summaryfile << " Skewing/preconditiong joint range │ " << skewing << std::endl;
        if (sfx_filename)
            summaryfile << "               Neural network file │ " << sfx_filename << std::endl;
        else if (onnx_filename)
            summaryfile << "               Neural network file │ " << onnx_filename << std::endl;
        summaryfile << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;

    }
    
    summaryfile << "Sampled estimate of final reachset │ " << RS.Xsampled; // sampled_reachset[sampled_reachset.size()-2];  // or iter-1 would be better
    summaryfile << "                Over-approximation │ " << RS.Xouter;
    summaryfile << "               Under-approximation │ " << RS.Xinner;
    if (uncontrolled > 0) {
    summaryfile << "         Robust over-approximation │ " << RS.Xouter_rob;
    summaryfile << "        Robust under-approximation │ " << RS.Xinner_rob;
    }
    summaryfile << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;
    summaryfile << "       Elapsed analysis time (sec) │ " << elapsed_secs << endl;
    summaryfile << "       Elapsed sampling time (sec) │ " << elapsed_secs_sampling << endl;
  /*              << "                           Network │ " << argv[2] << std::endl
                << "                            Method │ " << method << std::endl
                << "                        # of steps │ " << n_steps << std::endl
                << "                  AAF approx. type │ " << aaftype << std::endl
                << "            # of ref. trajectories │ " << system_network.get_samples_amount() << std::endl
                << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;
    int once = 1;
    for(int j = 0; j < n_inputs; j++){
    summaryfile << "         x"<<j<<"  Projected UA coverage │ " << std::setw(3) << percent(stats["ua_proj_cov"][j][0]) << "/" << std::setw(3) << percent(stats["ua_proj_cov"][j][1]) << "/" << std::setw(3) << percent(stats["ua_proj_cov"][j][2]) << (once ? " (min/avg/max)" : "") << " faults : " << stats_proj_fault[j] << std::endl;
    if(once)
        once = 0;
    summaryfile << "                Robust UA coverage │ " << std::setw(3) << percent(stats["ua_rob_cov"][j][0]) << "/" << std::setw(3) << percent(stats["ua_rob_cov"][j][1]) << "/" << std::setw(3) << percent(stats["ua_rob_cov"][j][2]) << " faults : " << stats_rob_fault[j] << std::endl;
    summaryfile << "                      OA proximity │ " << std::setw(3) << percent(stats["oa_prox"][j][0]) << "/" << std::setw(3) << percent(stats["oa_prox"][j][1]) << "/" << std::setw(3) << percent(stats["oa_prox"][j][2]) << " faults : " << stats_oa_fault[j] << std::endl;
    }
    summaryfile << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;
    summaryfile << "                  Computation time │ " << std::chrono::duration_cast<std::chrono::seconds>(t_end_method - t_start).count() << "s" << std::endl
                << "                     Sampling time │ " << std::chrono::duration_cast<std::chrono::seconds>(t_end_sampling - t_end_method).count() << "s" << std::endl
                << "                             TOTAL │ " << std::chrono::duration_cast<std::chrono::seconds>(t_end_sampling - t_start).count() << "s" << std::endl;
    summaryfile << std::endl;
    for(int i = 0; i < argc; i ++){
        summaryfile << argv[i] << " ";
    } */
    summaryfile.close();
    
    //begin = clock();
    if (print_debug)
    {
      //  if (systype == 0 || systype == 1 || (systype == 2 && nb_steps>1))
            run_pythonscript_visualization();
     //   else if (systype == 2) // nb_steps = 1
     //       system("cd GUI; python3 Visu_function.py; cd ..");
    }
    //end = clock();
    //double elapsed_printing = double(end - begin) / CLOCKS_PER_SEC;
    
    // printing final summary
 //   print_finalstats(begin);
    
    
}








void read_system(int argc, char* argv[])
{
  //  char sfx_filename[LINESZ]={0};
  //  char onnx_filename[LINESZ]={0};
    
    systype = 1; // EDO = 0, DDE = 1, discrete systems = 2, DNN = 3
    syschoice = 1;
    nb_sample_per_dim = 20;
    nbmax_sampled_dim = 4;
    
    char * str_systype = getCmdOption(argv, argv + argc, "-systype");
    if (str_systype)
    {
        if (strcmp(str_systype,"ode")==0) {
            systype = 0;
            nb_sample_per_dim = 3;
            nbmax_sampled_dim = 3;
        }
        else if (strcmp(str_systype,"dde")==0) {
            systype = 1;
            nb_sample_per_dim = 3;
            nbmax_sampled_dim = 3;
        }
        else if (strcmp(str_systype,"discrete")==0) {
            systype = 2;
            nb_sample_per_dim = 20;
            nbmax_sampled_dim = 10;
        }
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

    if (systype == 2) // discrete systems
    {
        init_discrete_system(config_filename); // reading initial condition in xinit and parameters: config file erase hard-coded conditions
    }
    else if (systype == 0 || systype == 1) // odes or ddes
    {
        define_system_dim(); // defines value of sysdim: depends on syschoice --
        init_system(); // initializes param from hard-coded values
        if (config_filename) // called with configuration file: we overwrite the initialization of init_system
            read_parameters(config_filename);
        init_utils_inputs();

    }
    if (systype == 1)
        tau = delay/nb_subdiv_delay;
    
    // moving old output results directory and creating new one
    system("mv output output_sauv");
    system("rm -r output");
    system("mkdir output");
    
    
    if (systype == 0 || systype == 1)
    {
        // for subdivisions: saving the full ranges
        initial_values_save = vector<interval>(sysdim);
        fullinputs_save = vector<interval>(fullinputsdim);
        for (int i=0 ; i<sysdim; i++)
            initial_values_save[i] = initial_values[i];
        for (int i=0 ; i<fullinputsdim; i++)
            fullinputs_save[i] = fullinputs[i];
    
    }
    
}


vector<vector<interval>> estimate_reachset()
{
    cout << "Estimate reachset:" << endl;
    
    /* pour systemes discrets mais a factoriser ?
    int nbmax_sampled_dim = 4;
    int nb_sampled_dim = min(jacdim,nbmax_sampled_dim);
    vector<int> local_id(nb_sampled_dim);
    
    int nb_points = nb_sample_per_dim;
    // limit the number of sampled points
    for (int i=1; i < nb_sampled_dim ; i++)  // MODIF borne min : 1 => 0 (?)
        nb_points = nb_points * (nb_sample_per_dim);
    cout << "nb_points=" << nb_points << endl;
    */
    
    
    if (systype == 2) {
        return estimate_robust_reachset_discrete(f);
    }
    else if (systype == 0)
        return estimate_reachset(obf,nb_sample_per_dim);
    else if (systype == 1)
        return estimate_reachset_dde(bf,nb_sample_per_dim);
}








