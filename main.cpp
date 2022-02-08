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
    const int LINESZ = 2048;
    double tn;    // current time
    double tau;   // integration time step (fixed step for now)
    int order;    // order of Taylor expansion
    
    double d0; // = 1;   // delay in DDE
    int nb_subdiv; // = 10;   // number of Taylor models on [0,d0]
    char sfx_filename[LINESZ]={0};
    char onnx_filename[LINESZ]={0};
    //const char* onnx_filename;
    
    // for ODEs only
    vector<vector<interval>> sampled_reachset;
    
    // for discrete systems only
    int nb_steps = 1;
    int AEextension_order = 1;
    int iter_method = 1;
    bool skew = true; // compute skewboxes or regular boxes approx
    vector<interval> xinit;
    // end discrete systems only
    
    int nb_sample_per_dim; // for range estimation by sampling: # of samples per dimension
    
    
    
    ReachSet RS;
   
    double elapsed_secs_sampling;
    
/*    essai_sherlock();
    syschoice = 101;
    function_range(NULL);
    
    exit(0); */
    
    
    /********* DEFINING SYSTEM *******************/
    // default is running example of CAV 2018 paper
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
        readfromfile_syschoice(config_filename,sfx_filename,onnx_filename,nb_sample_per_dim);
    
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
    
    
    
    open_outputfiles();
    
    clock_t begin, end; //  = clock();
    begin = clock();
    
    /*******************************************************************************************/
   /************* Discrete Systems ************/
    if (systype == 2) {
        
        DiscreteFunc f;
       
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
            skew = atoi(str_skew); //
                
        xinit = init_discrete_system(); // reading hard-coded initial condition and reading config file
        
        // reading initial conditions from config file if any
        if (config_filename)
            read_parameters_discrete(config_filename,xinit,nb_steps,order,AEextension_order,iter_method,skew);

        sampled_reachset = estimate_reachset(f, nb_steps, xinit, nb_sample_per_dim);
        end = clock();
        
        elapsed_secs_sampling = double(end - begin) / CLOCKS_PER_SEC;
        
        begin = clock();

        print_init_discrete(xinit,skew);
        if (nb_steps == 1)
            function_range(f,xinit,sampled_reachset[1]);
        else {
            if (iter_method == 1)
                RS = discrete_dynamical(f,xinit,sampled_reachset,nb_steps,AEextension_order,skew);
            else // (iter_method == 2)
                RS = discrete_dynamical_method2(f,xinit,sampled_reachset,nb_steps,skew);
        }
       
    }
    
    else if (systype == 0 || systype == 1)
    {
    // *************** ODEs or DDEs *******************
    
    
    define_system_dim(); // defines value of sysdim: depends on syschoice -- reads from file if input at command-line
    if (config_filename) // called with configuration file: we overwrite the initialization of init_system
           readfromfile_nbsubdiv(config_filename, nb_subdiv_init);
    init_system(t_begin,t_end,tau,d0,nb_subdiv,order); // initializes param from hard-coded values
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
        
        cout << "params=" << params;
        cout << "Estimate reachset:" << endl;
        nb_sample_per_dim = 2;
        sampled_reachset = estimate_reachset_dde(bf,initial_values,t_begin,t_end,d0, nb_subdiv, nb_sample_per_dim);
        
        // ecrire ici un estimate_reachset_dde qui retourne soit sampled_reachset soit exact set quand il existe et qui afffiche les samples + la solution exacte en XML.
        
        end = clock();
        elapsed_secs_sampling = double(end - begin) / CLOCKS_PER_SEC;
        begin = clock();
        
        
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, initial_values_save, fullinputs_save, component_to_subdiv);
            
            // a changer de 0 en 1 comme pour ODE (l'iteration 0 servant a stocker l'instant initial non stocke sinon - c'est du moins le cas en ODE, a verifier pour DDE)
            current_iteration = 0;
            
            tn = t_begin;
            
            HybridStep_dde prev_step = HybridStep_dde(bf,bbf,order,tn,tau,d0,nb_subdiv);
            prev_step.init_dde(sampled_reachset[0]);
            
            HybridStep_dde cur_step = prev_step.init_nextbigstep(tau);
            
            //******** Integration loop ************************
            int iter = nb_subdiv+1;
            while (cur_step.tn+d0 <= t_end)
            {
                // build Taylor model for Value and Jacobian and deduce guards for each active mode
                for (int j=0 ; j<nb_subdiv ; j++)
                {
                    cur_step.TM_build(j);  // j-th Taylor model valid on [tn+j*d0,tn+(j+1)*d0]
                    // eval at tn+tau
                    RS = cur_step.TM_evalandprint_solutionstep(j,eps,sampled_reachset[iter]); // cur_step.tn+(j+1)*tau/nb_subdiv);
                    cout << endl;
                    
                    if (j < nb_subdiv-1)
                        cur_step.init_nextsmallstep(j);
                    iter++;
                }
                cur_step = cur_step.init_nextbigstep(tau);
                
            }
            
            // adding a white line separator between subdivisions in the output result (except for maximal outer approx which is computed by union of all subdivisions)
            // removed when passed to yaml but not tested since => check
        }
        
        if (nb_subdiv_init > 1)
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
        
        for (int j = 0; j < jacdim-sysdim ; j++)
             param_inputs[j] = fullinputs[j];
        
        cout << "params=" << params;
        cout << "Estimate reachset:" << endl;
        nb_sample_per_dim = 2;
        sampled_reachset = estimate_reachset(obf,initial_values,params,inputs,t_begin,t_end,tau, nb_sample_per_dim);
        
        end = clock();
        elapsed_secs_sampling = double(end - begin) / CLOCKS_PER_SEC;
        begin = clock();
        
        for (current_subdiv=1 ; current_subdiv<=nb_subdiv_init; current_subdiv++)
        {
           
            if (nb_subdiv_init > 1)
                init_subdiv(current_subdiv, initial_values_save, fullinputs_save, component_to_subdiv);
            
            set_initialconditions(param_inputs,param_inputs_center,x,xcenter,J);  //            setId(J0);
            
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
            
            vector<F<AAF>> temp(sysdim);
            for (int j=0 ; j<sysdim; j++)
                temp[j] = x[j];
           
            current_iteration = 1;
            HybridStep_ode cur_step = init_ode(obf,xcenter,x,J,tn,tau,order);
           
            int iter = 1;
            while (cur_step.tn < t_end-0.0001*t_end)
            {
                // build Taylor model for Value and Jacobian and deduce guards for each active mode
                cur_step.TM_build(params,param_inputs,param_inputs_center);
                RS = cur_step.TM_evalandprint_solutionstep(eps,cur_step.tn+tau,sampled_reachset[iter],current_subdiv);
                cur_step.init_nextstep(params,param_inputs,tau);
                
                if ((RS.Xouter[0] == interval::EMPTY()) || (RS.Xouter[0] == interval::ENTIRE())) {
                    printf("Terminated due to too large overestimation.\n");
                    break;
                }
                iter++;
            }
            // adding a white line separator between subdivisions in the output result (except for maximal outer approx which is computed by union of all subdivisions)
            // removed when passed to yaml but not tested since => check
            
        }
        
        if (nb_subdiv_init > 1)
            print_finalsolution(ceil((t_end-t_begin)/tau), d0);
        
    }
    }
    
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
    
   
    // discrete systems
    if (systype == 2) {
    
    out_summary << YAML::Key << "nb_steps";
    out_summary << YAML::Value << nb_steps;
    out_summary << YAML::Key << "nb_sample_per_dim";
    out_summary << YAML::Value << nb_sample_per_dim;
    out_summary << YAML::Key << "skew";
    out_summary << YAML::Value << skew;
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
    
    out_summary << YAML::Key << "zinner";
    vector<double> aff_zinner(sysdim*2);
    for (int i=0; i<sysdim ; i++) {
        aff_zinner[2*i] = RS.Xinner[i].inf();
        aff_zinner[2*i+1] = RS.Xinner[i].sup();
    }
    out_summary << YAML::Value << aff_zinner;
    
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
            summaryfile <<"("<<inputs[i].convert_int()<<","<<nb_inputs[i]<<") ";
        summaryfile << std::endl;
        summaryfile << "                     Time interval │ " << "[" << t_begin << "," << t_end << "]" << std::endl;
        summaryfile << "                      Control step │ " << control_period << std::endl;
        summaryfile << "                         Time step │ " << tau << std::endl;
        summaryfile << "                  System dimension │ " << sysdim << std::endl;
        summaryfile << "               Taylor Models order │ " << order << std::endl;
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
        summaryfile << "                             Delay │ " << d0 << std::endl;
        summaryfile << "                         Time step │ " << tau << std::endl;
        summaryfile << "                  System dimension │ " << sysdim << std::endl;
        summaryfile << "               Taylor Models order │ " << order << std::endl;
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
        summaryfile << " Skewing/preconditiong joint range │ " << skew << std::endl;
        if (sfx_filename)
            summaryfile << "               Neural network file │ " << sfx_filename << std::endl;
        else if (onnx_filename)
            summaryfile << "               Neural network file │ " << onnx_filename << std::endl;
        summaryfile << "───────────────────────────────────┼──────────────────────────────────────" << std::endl;

    }
    
    summaryfile << "Sampled estimate of final reachset │ " << RS.Xsampled; // sampled_reachset[sampled_reachset.size()-2];  // or iter-1 would be better
    summaryfile << "                Over-approximation │ " << RS.Xouter;
    summaryfile << "               Under-approximation │ " << RS.Xinner;
    
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
        if (systype == 0 || systype == 1 || (systype == 2 && nb_steps>1))
            run_pythonscript_visualization();
        else if (systype == 2) // nb_steps = 1
            system("cd GUI; python3 Visu_function.py; cd ..");
    }
    //end = clock();
    //double elapsed_printing = double(end - begin) / CLOCKS_PER_SEC;
    
    // printing final summary
 //   print_finalstats(begin);
    
    
}




















