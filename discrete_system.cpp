/* ============================================================================
 File   : discrete system.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file computes ranges for functions and discrete dynamical systems
 ============================================================================*/

#include "discrete_system.h"
#include "utils.h"
#include "matrix.h"
#include "inner.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <ctime>
#include <assert.h>
#include <cstring>
using namespace std;


vector<vector<vector<interval>>> constr_eps;  // constraints on noise symbols when partitioning 2D region in 4
vector<vector<vector<interval>>> eps_loc;     // consequence on [x]-x0  when partitioning 2D region in 4

vector<vector<vector<vector<interval>>>> constr_eps_discr;  // same but with additional discretization in each direction
vector<vector<vector<vector<interval>>>> eps_loc_discr;     // consequence on [x]-x0  when partitioning 2D region in 4
vector<vector<vector<vector<double>>>> extremity_eps_loc_discr;

int nb_discr, nb_discr1, nb_discr2;


vector<interval> init_discrete_system() // (char * config_filename)
{
    
    if (syschoice == 1) {
        jacdim = 1;  // number of input components
        sysdim = 1; // number of output components
    }
    if (syschoice == 111) {  // test sigmoid
        jacdim = 1;
        sysdim = 1;
    }
    else if (syschoice == 2) {
        jacdim = 2;
        sysdim = 1;
    }
    else if (syschoice == 3) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 4) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 5) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 6) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 7) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 8) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 9) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 10) {
        jacdim = 3;
        sysdim = 3;
    }
    else if (syschoice == 11) {
        jacdim = 3;
        sysdim = 3;
    }
    else if (syschoice == 12) {
        jacdim = 3;
        sysdim = 3;
    }
    else if (syschoice == 13) {
        jacdim = 3;  // number of input components
        sysdim = 3; // number of output components
    }
    else if (syschoice == 14) {
        jacdim = 1;  // number of input components
        sysdim = 1; // number of output components
    }
    else if (syschoice == 15) {  // test model parallelotope bundles HSCC 2016 p 303
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 16) {  // SIR epidemic model  - parallelotope bundles HSCC 2016 p 303
        jacdim = 3;
        sysdim = 3;
    }
    else if (syschoice == 17) {  // Honeybees model  - parallelotope bundles HSCC 2016 p 303-304
        jacdim = 5;
        sysdim = 5;
    }
    else if (syschoice == 18) {  // SIR epidemic model  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
        jacdim = 4;             // identical to 16 except parameters - gamma is now uncertain (x[3])
        sysdim = 3;
    }
    else if (syschoice == 19) {  // SIR epidemic model (first 2 dimensions) - parallelotope bundles HSCC 2016 p 303
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 20) {  // SIR epidemic model (first 2 dimensions)  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
        jacdim = 3;             // identical to 16 except parameters - gamma is now uncertain (x[3])
        sysdim = 2;
    }
    else if (syschoice == 21) {  // SIR epidemic model (first 2 dimensions)  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
        jacdim = 3;             // identical to 16 except parameters - gamma is now uncertain (x[3])
        sysdim = 3;
    }
    else if (syschoice == 22) {  // DNN CAV-DINO
        jacdim = 2;             //
        sysdim = 2;             // a transformer en sysdim = 1 si ca fonctionne ?
    }
    else if (syschoice == 23) {  // mountaincar (avec NN sfx format)
        jacdim = 2;             //
        sysdim = 2;
    }
    else if (syschoice == 231) {  // mountaincar (avec NN onnx format)
        jacdim = 2;             //
        sysdim = 2;
    }
    else if (syschoice == 100) { // neural network sfx format
        jacdim = NH.n_inputs;
        sysdim = NH.n_outputs;
    }
    else if (syschoice == 101) { // neural network onnx format
#if ONNX_active
        jacdim = CG.no_of_input_nodes;
        sysdim = CG.no_of_output_nodes;
#endif
    }
    
    vector<interval> res(jacdim);
    
    if (syschoice == 1) {
        res[0] = interval(2,3);
    }
    else if (syschoice == 111) { // test sigmoide
        res[0] = interval(-0.5,0.5);
    }
    else if (syschoice == 2) {
        res[0] = interval(2,3);
        res[1] = interval(2,3);
    }
    else if (syschoice == 3) { // example 3.5 Goldztejn
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
    }
    else if (syschoice == 4) { // example 3 CDC
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
    }
    else if (syschoice == 5) { // example 5.1 Goldstzjen
        res[0] = interval(0.9,1.1); // interval(0.99,1.01);
        res[1] = interval(0.9,1.1); // interval(0.99,1.01);
    }
    else if (syschoice == 6) { //
        res[0] = interval(0.9,1.5);
        res[1] = interval(0.9,1.5);
    }
    else if (syschoice == 7) { //
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
    }
    else if (syschoice == 8) { //
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
    }
    else if (syschoice == 9) { //
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
    }
    else if (syschoice == 10) {
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
        res[2] = interval(0.9,1.1);
    }
    else if (syschoice == 11) {
        res[0] = interval(0.9,1.1);
        res[1] = interval(0.9,1.1);
        res[2] = interval(0.9,1.1);
    }
    else if (syschoice == 12) {
        res[0] = interval(-1.0,1.0);
        res[1] = interval(-1.0,1.0);
        res[2] = interval(-1.0,1.0);
    }
    else if (syschoice == 13) {
        res[0] = interval(0,1.0);
        res[1] = interval(0.9,1.1);
        res[2] = interval(-0.25,0.25);
    }
    else if (syschoice == 14) {
        res[0] = interval(-0.25,0.25);
    }
    else if (syschoice == 15) {  // test model parallelotope bundles HSCC 2016 p 303
        
        res[0] = interval(0.05,0.1);
        res[1] = interval(0.99,1.0);
    //    nb_steps = 25;
    }
    else if (syschoice == 16) {  // SIR epidemic model  - parallelotope bundles HSCC 2016 p 303
        res[0] = interval(0.79,0.80);
        res[1] = interval(0.19,0.20);
        res[2] = interval(0.,0.1);
    //    nb_steps = 20; // 60;
    }
    else if (syschoice == 17) {  // Honeybees model  - parallelotope bundles HSCC 2016 p 303-304
        res[0] = interval(500.0,500.0); // interval(500.0,510.0);
        res[1] = interval(390.0,400.0);
        res[2] = interval(90.0,100.0);
       res[3] = interval(0.0,0.0);// res[3] = interval(-0.1,0.1);
       res[4] =  interval(0.0,0.0);// res[4] = interval(-0.1,0.1);
      //  nb_steps = 15;
    }
    else if (syschoice == 18) {         // SIR epidemic model  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
        res[0] = interval(0.79,0.80);   // identical to 16 except parameters - gamma is now uncertain (x[3])
        res[1] = interval(0.19,0.20);
        res[2] = interval(-0.,0.1);
        res[3] = interval(0.05,0.0675);  //interval(0.05,0.0675); // parameter gamma
      //  nb_steps = 30;
    }
    else if (syschoice == 19) {  // SIR epidemic model (first 2 dimensions) - parallelotope bundles HSCC 2016 p 303
        res[0] = interval(0.79,0.80);
        res[1] = interval(0.19,0.20);
 //       nb_steps = 20; // 60;
    }
    else if (syschoice == 20) {         // SIR epidemic model (first 2 dimensions)  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
        res[0] = interval(0.79,0.80);   // identical to 16 except parameters - gamma is now uncertain (x[3])
        res[1] = interval(0.19,0.20);
        res[2] = interval(0.05,0.0675);  //interval(0.05,0.0675); // parameter gamma
        //  nb_steps = 30;
    }
    else if (syschoice == 21) {         // SIR epidemic model (first 2 dimensions)  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
        res[0] = interval(1.,2.);  // interval(0.79,0.80);   // identical to 16 except parameters - gamma is now uncertain (x[3])
        res[1] = interval(3.,4.); // interval(0.19,0.20);
        res[2] = interval(0.,1.0); //interval(0.05,0.0675);  //interval(0.05,0.0675); // parameter gamma
        //  nb_steps = 30;
    }
    else if (syschoice == 22) {  // DNN CAV-DINO
        res[0] = interval(0.,0.1);
        res[1] = interval(0.,0.1);
    }
    else if (syschoice == 23) {  // mountain car (avec RNN sfx format)
        res[0] = interval(-0.5,-0.48); // interval(-0.5,-0.48);
       // res[1] = interval(0.,0.00);
        res[1] = interval(0.,0.001);
    }
    else if (syschoice == 231) {  // mountain car (avec RNN onnx format)
        res[0] = interval(-0.5,-0.48); // interval(-0.5,-0.48);
        // res[1] = interval(0.,0.00);
        res[1] = interval(0.,0.001);
    }
    else if (syschoice == 100) { // neural network sfx format
        // should be read from config file
   //     res[0] = interval(0,0.2);
    //    res[1] = interval(0,0.2);
    }
    else if (syschoice == 101) { // neural network onnx format
        // should be read from config file
        res[0] = interval(5,5.);
        res[1] = interval(-2.,-2.);
    }
    
    // reading initial conditions from config file
  /*  if (config_filename) {
        read_initialconditions(config_filename,res);
    } */
    return res;
}


// d0 and t_begin are for DDEs only, rest are common to ODE and DDE
void read_parameters_discrete(const char * params_filename, vector<interval> &res, int &nb_steps, int &order, int &AEextension_order, int &iter_method,bool &skew)
{
    const int LINESZ = 2048;
    char buff[LINESZ];
    char initialcondition[LINESZ];
    const char space[2] = " ";
    double a, b;
    int c;
    int skew_int;
    
    cout << "****** Reading initial conditions and parameters from file " <<  params_filename << " ******" << endl;
    FILE *params_file = fopen(params_filename,"r");
    if (params_file == NULL)
        cout << "Error reading " << params_filename << ": file not found" << endl;
    
    while (fgets(buff,LINESZ,params_file)) {
        sscanf(buff, "interactive-visualization = %d\n", &interactive_visualization);
        sscanf(buff, "order = %d\n", &order);
        sscanf(buff, "nbsteps = %d\n", &nb_steps);
        sscanf(buff, "AEextension-order = %d\n", &AEextension_order);
        sscanf(buff, "iter-method = %d\n", &iter_method);
        sscanf(buff, "skew = %d\n", &skew_int);
        sscanf(buff, "interactive-visualization = %d\n", &interactive_visualization);
        sscanf(buff, "refined-mean-value = %d\n", &refined_mean_value);
        if (sscanf(buff, "variables-to-display = %s\n", initialcondition) == 1)
        {
            for (int i=0; i< sysdim; i++)
                variables_to_display[i] = false;
            
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i;
            while( token != NULL ) {
                sscanf(token,"%d",&i);
                variables_to_display[i-1] = true;
             //   cout <<"input="<<inputs[i].convert_int()<<endl;
                token = strtok(NULL,space);
            }
        }
        //     sscanf(buff, "system = %s\n", sys_name);
        //      sscanf(buff, "initially = %[^\n]\n", initial_condition);   // tell separator is newline, otherwise by default it is white space
        if (sscanf(buff, "initial-values = %s\n", initialcondition) == 1)
        {
            char *token = strtok(buff,space);
            token = strtok(NULL,space);
            token = strtok(NULL,space);
            int i = 0;
            nb_subdiv_init = 1;
            component_to_subdiv = -1;
            component_to_subdiv2 = -1;
            // only one component can be subdivided for now (we keep the last if several...)
            while( token != NULL ) {
                if (sscanf(token,"([%lf,%lf],%d)",&a,&b,&c) == 3) {
                    if (component_to_subdiv > -1) // we already have one component to subdivide
                    {
                        component_to_subdiv2 = i;
                        cout << "component "<< i << " should also be subdivided " << endl;
                    }
                    else
                    {
                        nb_subdiv_init = c;
                        component_to_subdiv = i;
                        cout << "component "<< i << " should be dubdivided " << c << "times" << endl;
                    }
                }
                else
                    sscanf(token,"[%lf,%lf]",&a,&b);
                res[i] = interval(a,b);
                cout <<"initial_value="<<res[i]<<endl;
                i++;
                token = strtok(NULL,space);
            }
        }
    }
    skew = (skew_int == 1);
    
    // we don't want to save and print too many points
    printing_period = 1;
    if (nb_steps >= points_per_graph)
        printing_period = nb_steps / points_per_graph;
    
    fclose(params_file);
    
//    cout << "interactive-visualization=" << interactive_visualization << endl;
    
    // cout << "system name = " << sys_name << endl;
    cout << "****** End initial conditions reading ******" << endl << endl;
}







void print_outer_range(vector<F<AAF>> &z_o, vector<interval> &range) {
    interval interv_out;
    for (int i=0; i < sysdim ; i++) {
        interv_out = z_o[i].x().convert_int();
        cout << "z_o[" << i << "]=" << interv_out.inf() << ", " << interv_out.sup() << "] ";
        cout << " eta_o["<<i<<"]=" << (range[i].sup() - range[i].inf())/ (interv_out.sup() - interv_out.inf()) <<" "; // precision metric with respect to estimated range
    }
    cout << endl;
}

void print_outer_range(vector<interval> &z_o, vector<interval> &range) {
    for (int i=0; i < sysdim ; i++) {
        cout << "z_o[" << i << "]=" << z_o[i].inf() << ", " << z_o[i].sup() << "] ";
        cout << " eta_o["<<i<<"]=" << (range[i].sup() - range[i].inf())/ (z_o[i].sup() - z_o[i].inf()) <<" "; // precision metric with respect to estimated range
    }
    cout << endl;
}



// same as discrete_dynamical but with skew box for joint range
ReachSet discrete_dynamical(DiscreteFunc &f, vector<interval> &xinit, vector<vector<interval>> &estimated_range, int &nb_steps, int order, bool skew)
{
    vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    // preconditioning: A is center of Jacobian, C its inverse
    vector<vector<double>> A_o(sysdim,vector<double> (sysdim)), C_o(sysdim,vector<double> (sysdim));
    vector<vector<double>> A_i(sysdim,vector<double> (sysdim)), C_i(sysdim,vector<double> (sysdim));
    
    vector<interval> f0_o(sysdim), f0_i(sysdim);
    vector<vector<interval>> CJacf_o(sysdim, vector<interval>(jacdim));
    vector<vector<interval>> CJacf_i(sysdim, vector<interval>(jacdim));
    vector<interval> z0_o, z0_i;
    
    /*********************************** ORDER 2 Taylor model  ***********************************/
    vector<vector<interval>> dfdx0_i(sysdim,vector<interval>(jacdim)), dfdx0_o(sysdim,vector<interval>(jacdim));
    vector<vector<interval>> Cdfdx0_i(sysdim,vector<interval>(jacdim)), Cdfdx0_o(sysdim,vector<interval>(jacdim));
    vector<interval> zf0_o;
    vector<interval> zf0_i;
    vector<vector<vector<interval>>> Hessf_o(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector<vector<vector<interval>>> Hessf_i(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector<vector<vector<interval>>> CHessf_o(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector<vector<vector<interval>>> CHessf_i(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector <F<F<AAF>>> xff_o(jacdim), xff_i(jacdim), zff_o(sysdim), zff_i(sysdim);
    vector <F<interval>> x0ff_o(jacdim), x0ff_i(jacdim), z0ff_o(sysdim), z0ff_i(sysdim);
    /* end ORDER 2 */
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    vector<int> aux;
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    int varx = 0, vary = 1, varz = 2; // joint range we want to print/compute
    vector<vector<double>> output_skewedbox;
    
    // for now, no input, sysdim = jacdim
    z_inner = xinit;
    z_outer = xinit;
    
    
    for (int i=0 ; i<sysdim; i++) {
        A_i[i][i] = 1.0;
        A_o[i][i] = 1.0;
    }
    
    
    F<AAF> temp, temp1, temp2; // order 1
    F<interval> tempf, tempf1, tempf2; // order 2
    F<F<AAF>> tempaff, tempaff1, tempaff2;
    
    for (int i=0; i < jacdim ; i++) {
        x_i[i] = z_inner[i];
        x0_i[i] = mid(z_inner[i]);
        radx_i[i] = z_inner[i] - x0_i[i];
        x_o[i] = z_outer[i];
        x0_o[i] = mid(z_outer[i]);
        radx_o[i] = z_outer[i] - x0_o[i];
    }
    
    for (int step = 1; step <= nb_steps ; step++)
    {
        if (step % printing_period == 0) {
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "tn";
            out_approx << YAML::Value << step;
        }
        
        for (int i=0; i < jacdim ; i++)
            x_i[i].diff(i,jacdim);
        
        // x_i = A . x_i
        if (sysdim >=3 && skew) {
            temp1 = A_i[varx][varx]*x_i[varx] + A_i[varx][vary]*x_i[vary] + A_i[varx][varz]*x_i[varz];
            temp2 = A_i[vary][varx]*x_i[varx] + A_i[vary][vary]*x_i[vary] + A_i[vary][varz]*x_i[varz];
            x_i[varz] = A_i[varz][varx]*x_i[varx] + A_i[varz][vary]*x_i[vary] + A_i[varz][varz]*x_i[varz];
            x_i[varx] = temp1;
            x_i[vary] = temp2;
        }
        else if (sysdim >=2 && skew) {
            temp = A_i[varx][varx]*x_i[varx] + A_i[varx][vary]*x_i[vary];
            x_i[vary] = A_i[vary][vary]*x_i[vary] + A_i[vary][varx]*x_i[varx];
            x_i[varx] = temp;
        }
        
        
        z_i = f(x_i);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_i[i][j] = z_i[i].d(j).convert_int();// JacAff_i[i][j].convert_int();
        
        for (int i=0; i < jacdim ; i++)
            x_o[i].diff(i,jacdim);
        
        if (sysdim >=3 && skew) {
            temp1 = A_o[varx][varx]*x_o[varx] + A_o[varx][vary]*x_o[vary] + A_o[varx][varz]*x_o[varz];
            temp2 = A_o[vary][varx]*x_o[varx] + A_o[vary][vary]*x_o[vary] + A_o[vary][varz]*x_o[varz];
            x_o[varz] = A_o[varz][varx]*x_o[varx] + A_o[varz][vary]*x_o[vary] + A_o[varz][varz]*x_o[varz];
            x_o[varx] = temp1;
            x_o[vary] = temp2;
        }
        else if (sysdim >=2 && skew) {
            temp = A_o[varx][varx]*x_o[varx] + A_o[varx][vary]*x_o[vary];
            x_o[vary] = A_o[vary][vary]*x_o[vary] + A_o[vary][varx]*x_o[varx];
            x_o[varx] = temp;
        }
        
        z_o = f(x_o);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_o[i][j] = z_o[i].d(j).convert_int(); //JacAff_o[i][j].convert_int();
        
        
       
        
        if (order == 2)
        {
            /*********************************** ORDER 2 Taylor model  ***********************************/
            
          //  zf0_o = FFunc(dfdx0_o,x0_o);  // can't use it here because of preconditioning
   
            for (int i=0; i < jacdim ; i++) {
                x0ff_o[i] = x0_o[i];
                x0ff_o[i].diff(i,jacdim);
            }
            
            if (sysdim >=3 && skew) {
                tempf1 = A_o[varx][varx]*x0ff_o[varx] + A_o[varx][vary]*x0ff_o[vary] + A_o[varx][varz]*x0ff_o[varz];
                tempf2 = A_o[vary][varx]*x0ff_o[varx] + A_o[vary][vary]*x0ff_o[vary] + A_o[vary][varz]*x0ff_o[varz];
                x0ff_o[varz] = A_o[varz][varx]*x0ff_o[varx] + A_o[varz][vary]*x0ff_o[vary] + A_o[varz][varz]*x0ff_o[varz];
                x0ff_o[varx] = tempf1;
                x0ff_o[vary] = tempf2;
            }
            else if (sysdim >=2 && skew) {
                tempf = A_o[varx][varx]*x0ff_o[varx] + A_o[varx][vary]*x0ff_o[vary];
                x0ff_o[vary] = A_o[vary][vary]*x0ff_o[vary] + A_o[vary][varx]*x0ff_o[varx];
                x0ff_o[varx] = tempf;
            }
            
            z0ff_o = f(x0ff_o);
            
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    dfdx0_o[i][j]=z0ff_o[i].d(j);
            
            for (int i=0; i < jacdim ; i++) {
                x0ff_i[i] = x0_i[i];
                x0ff_i[i].diff(i,jacdim);
            }

            if (sysdim >=3 && skew) {
                tempf1 = A_i[varx][varx]*x0ff_i[varx] + A_i[varx][vary]*x0ff_i[vary] + A_i[varx][varz]*x0ff_i[varz];
                tempf2 = A_i[vary][varx]*x0ff_i[varx] + A_i[vary][vary]*x0ff_i[vary] + A_i[vary][varz]*x0ff_i[varz];
                x0ff_i[varz] = A_i[varz][varx]*x0ff_i[varx] + A_i[varz][vary]*x0ff_i[vary] + A_i[varz][varz]*x0ff_i[varz];
                x0ff_i[varx] = tempf1;
                x0ff_i[vary] = tempf2;
            }
            else if (sysdim >=2 && skew) {
                tempf = A_i[varx][varx]*x0ff_i[varx] + A_i[varx][vary]*x0ff_i[vary];
                x0ff_i[vary] = A_i[vary][vary]*x0ff_i[vary] + A_i[vary][varx]*x0ff_i[varx];
                x0ff_i[varx] = tempf;
            }
            
            z0ff_i = f(x0ff_i);
            
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    dfdx0_i[i][j]=z0ff_i[i].d(j);
            
         
            
            for (int i=0; i < jacdim ; i++) {
                xff_o[i] = z_outer[i];
                xff_o[i].diff(i,jacdim);          // first order
                xff_o[i].x().diff(i,jacdim);      // second order
            }
            
            if (sysdim >=3 && skew) {
                tempaff1 = A_o[varx][varx]*xff_o[varx] + A_o[varx][vary]*xff_o[vary] + A_o[varx][varz]*xff_o[varz];
                tempaff2 = A_o[vary][varx]*xff_o[varx] + A_o[vary][vary]*xff_o[vary] + A_o[vary][varz]*xff_o[varz];
                xff_o[varz] = A_o[varz][varx]*xff_o[varx] + A_o[varz][vary]*xff_o[vary] + A_o[varz][varz]*xff_o[varz];
                xff_o[varx] = tempaff1;
                xff_o[vary] = tempaff2;
            }
            else if (sysdim >=2 && skew) {
                tempaff = A_o[varx][varx]*xff_o[varx] + A_o[varx][vary]*xff_o[vary];
                xff_o[vary] = A_o[vary][vary]*xff_o[vary] + A_o[vary][varx]*xff_o[varx];
                xff_o[varx] = tempaff;
            }
            
            zff_o = f(xff_o);
            
            for (int i=0; i < jacdim ; i++) {
                xff_i[i] = z_inner[i];
                xff_i[i].diff(i,jacdim);          // first order
                xff_i[i].x().diff(i,jacdim);      // second order
            }
            
            if (sysdim >=3 && skew) {
                tempaff1 = A_i[varx][varx]*xff_i[varx] + A_i[varx][vary]*xff_i[vary] + A_i[varx][varz]*xff_i[varz];
                tempaff2 = A_i[vary][varx]*xff_i[varx] + A_i[vary][vary]*xff_i[vary] + A_i[vary][varz]*xff_i[varz];
                xff_i[varz] = A_i[varz][varx]*xff_i[varx] + A_i[varz][vary]*xff_i[vary] + A_i[varz][varz]*xff_i[varz];
                xff_i[varx] = tempaff1;
                xff_i[vary] = tempaff2;
            }
            else if (sysdim >=2 && skew) {
                tempaff = A_i[varx][varx]*xff_i[varx] + A_i[varx][vary]*xff_i[vary];
                xff_i[vary] = A_i[vary][vary]*xff_i[vary] + A_i[vary][varx]*xff_i[varx];
                xff_i[varx] = tempaff;
            }
            
            zff_i = f(xff_i);
            
            for (int i=0; i < sysdim ; i++) {
                for (int j=0; j < jacdim ; j++) {
                    for (int k=0; k < jacdim ; k++) {
                        Hessf_o[i][j][k] =zff_o[i].d(j).d(k).convert_int();
                        Hessf_i[i][j][k] =zff_i[i].d(j).d(k).convert_int();
                    }
                    //    cout << "JacAff[i][j]=" << JacAff[i][j] << JacAff[i][j].convert_int();
                }
            }
        }
        
        // x0 = A x0
        if (sysdim >=3 && skew) {
            interval tempi1 = A_i[varx][varx]*x0_i[varx] + A_i[varx][vary]*x0_i[vary] + A_i[varx][varz]*x0_i[varz];
            interval tempi2 = A_i[vary][varx]*x0_i[varx] + A_i[vary][vary]*x0_i[vary] + A_i[vary][varz]*x0_i[varz];
            x0_i[varz] = A_i[varz][varx]*x0_i[varx] + A_i[varz][vary]*x0_i[vary] + A_i[varz][varz]*x0_i[varz];
            x0_i[varx] = tempi1;
            x0_i[vary] = tempi2;
            
            tempi1 = A_o[varx][varx]*x0_o[varx] + A_o[varx][vary]*x0_o[vary] + A_o[varx][varz]*x0_o[varz];
            tempi2 = A_o[vary][varx]*x0_o[varx] + A_o[vary][vary]*x0_o[vary] + A_o[vary][varz]*x0_o[varz];
            x0_o[varz] = A_o[varz][varx]*x0_o[varx] + A_o[varz][vary]*x0_o[vary] + A_o[varz][varz]*x0_o[varz];
            x0_o[varx] = tempi1;
            x0_o[vary] = tempi2;
        }
        else if (sysdim >=2 && skew) {
            interval tempi = A_i[varx][varx]*x0_i[varx] + A_i[varx][vary]*x0_i[vary];
            x0_i[vary] = A_i[vary][vary]*x0_i[vary] + A_i[vary][varx]*x0_i[varx];
            x0_i[varx] = tempi;
        
            tempi = A_o[varx][varx]*x0_o[varx] + A_o[varx][vary]*x0_o[vary];
            x0_o[vary] = A_o[vary][vary]*x0_o[vary] + A_o[vary][varx]*x0_o[varx];
            x0_o[varx] = tempi;
        }
        
        z0_o = f(x0_o);
        z0_i = f(x0_i);
        
        /*********************************** Evaluation ***********************************/
        
        if (order == 1)
        {
            z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
            z_inner_proj = evaluate_innerrange(z0_i, radx_i, Jacf_i, true, aux);
            if (jacdim > sysdim)
                z_inner_proj_rob = evaluate_innerrange_robust(z0_i, radx_i, Jacf_i, true, aux);
        }
        else if (order == 2)
        {
            z_outer = evaluate_outerrange_order2(z0_o, radx_o, dfdx0_o, Hessf_o);
            z_inner_proj = evaluate_innerrange_order2(z0_i, radx_i, dfdx0_i, Hessf_i, true, aux);
            if (jacdim > sysdim)
                z_inner_proj_rob = evaluate_innerrange_order2_robust(z0_i, radx_i, dfdx0_i,  Hessf_i, true, aux);
        }
        else
            assert(false);
        
        for (int i=0; i < sysdim ; i++)
            z_outer[i] = intersect(z_outer[i],z_o[i].x().convert_int());
        
        if (step % printing_period == 0)
            print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step,estimated_range[step]);
        
        
        if (sysdim >= 2) {
            
            //   assert(sysdim <= 3); // higher dimension not treated yet
            for (int i=0 ; i<sysdim; i++)
                exist_quantified[i] = i;
            
            // fixed for now
            if (syschoice == 15) {
                exist_quantified[varx] = varx;
                exist_quantified[vary] = vary;
            }
            else if (syschoice == 16) {
                exist_quantified[varx] = varx;
                exist_quantified[vary] = vary;
                exist_quantified[2] = 2;
            }
            else if (syschoice == 17) {
                for (int j=0; j < jacdim ; j++)
                    exist_quantified[j] = j; // (j+2)%(jacdim);
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            else if (syschoice == 18) {
                for (int j=0; j < sysdim ; j++)
                    exist_quantified[j] = j; // (j+2)%(jacdim);
                exist_quantified[3] = 3; // because constant parameter
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            else if (syschoice == 20) {
                for (int j=0; j < sysdim ; j++)
                    exist_quantified[j] = j; // (j+2)%(jacdim);
                exist_quantified[2] = 1;
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            
            
            for (int i=0 ; i<sysdim; i++) {
                A_o[i][i] = 1.0;
                C_o[i][i] = 1.0;
                A_i[i][i] = 1.0;
                C_i[i][i] = 1.0;
            }
            
            if (sysdim >=3 && skew) {
                for (int i=0 ; i<sysdim; i++) {
                    A_o[i][i] = 1.0;
                    C_o[i][i] = 1.0;
                    A_i[i][i] = 1.0;
                    C_i[i][i] = 1.0;
                }
                // 3D preconditioner -  supposing for now that varx, vary, varz = 0, 1, 2
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++) {
                        if (order == 1) {
                            A_o[i][j] = mid(Jacf_o[i][j]);
                            A_i[i][j] = mid(Jacf_i[i][j]);
                        }
                        else {
                            A_o[i][j] = mid(dfdx0_o[i][j]);
                            A_i[i][j] = mid(dfdx0_i[i][j]);
                        }
                    }
                }
                
                // C is inverse of A
                // supposing for now that varx, vary, varz = 0, 1, 2
                double determinant = 0;
                for (int i = 0; i < 3; i++)
                    determinant = determinant + (A_o[0][i] * (A_o[1][(i+1)%3] * A_o[2][(i+2)%3] - A_o[1][(i+2)%3] * A_o[2][(i+1)%3]));
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++)
                        C_o[i][j] = ((A_o[(j+1)%3][(i+1)%3] * A_o[(j+2)%3][(i+2)%3]) - (A_o[(j+1)%3][(i+2)%3] * A_o[(j+2)%3][(i+1)%3]))/ determinant;
                }
                
                determinant = 0;
                for (int i = 0; i < 3; i++)
                    determinant = determinant + (A_i[0][i] * (A_i[1][(i+1)%3] * A_i[2][(i+2)%3] - A_i[1][(i+2)%3] * A_i[2][(i+1)%3]));
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++)
                        C_i[i][j] = ((A_i[(j+1)%3][(i+1)%3] * A_i[(j+2)%3][(i+2)%3]) - (A_i[(j+1)%3][(i+2)%3] * A_i[(j+2)%3][(i+1)%3]))/ determinant;
                }
            }
            else if (sysdim >=2 && skew)
            {
                // 2D preconditioner - just on components varx and vary
                if (order == 1) {
                    A_o[varx][varx] = mid(Jacf_o[varx][varx]);
                    A_o[vary][vary] = mid(Jacf_o[vary][vary]);
                    A_o[varx][vary] = mid(Jacf_o[varx][vary]);
                    A_o[vary][varx] = mid(Jacf_o[vary][varx]);
                    A_i[varx][varx] = mid(Jacf_i[varx][varx]);
                    A_i[vary][vary] = mid(Jacf_i[vary][vary]);
                    A_i[varx][vary] = mid(Jacf_i[varx][vary]);
                    A_i[vary][varx] = mid(Jacf_i[vary][varx]);
                }
                else if (order == 2) {
                    
                    A_o[varx][varx] = mid(dfdx0_o[varx][varx]);
                    A_o[vary][vary] = mid(dfdx0_o[vary][vary]);
                    A_o[varx][vary] = mid(dfdx0_o[varx][vary]);
                    A_o[vary][varx] = mid(dfdx0_o[vary][varx]);
                    A_i[varx][varx] = mid(dfdx0_i[varx][varx]);
                    A_i[vary][vary] = mid(dfdx0_i[vary][vary]);
                    A_i[varx][vary] = mid(dfdx0_i[varx][vary]);
                    A_i[vary][varx] = mid(dfdx0_i[vary][varx]);
                
                }
            
                // C is inverse of A
                double determinant = 1.0/(A_o[varx][varx]*A_o[vary][vary]-A_o[varx][vary]*A_o[vary][varx]);
                C_o[varx][varx] = determinant*A_o[vary][vary];
                C_o[varx][vary] = - determinant*A_o[varx][vary];
                C_o[vary][varx] = - determinant*A_o[vary][varx];
                C_o[vary][vary] = determinant*A_o[varx][varx];
                determinant = 1.0/(A_i[varx][varx]*A_i[vary][vary]-A_i[varx][vary]*A_i[vary][varx]);
                C_i[varx][varx] = determinant*A_i[vary][vary];
                C_i[varx][vary] = - determinant*A_i[varx][vary];
                C_i[vary][varx] = - determinant*A_i[vary][varx];
                C_i[vary][vary] = determinant*A_i[varx][varx];
            }
            
            // f0 = C * z0
            for (int i=0 ; i<sysdim; i++) {
                f0_o[i] = z0_o[i];
                f0_i[i] = z0_i[i];
            }
            
            if (sysdim >=3 && skew) {
                f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary] + C_o[varx][varz]*z0_o[varz];
                f0_o[vary] = C_o[vary][varx]*z0_o[varx] + C_o[vary][vary]*z0_o[vary] + C_o[vary][varz]*z0_o[varz];
                f0_o[varz] = C_o[varz][varx]*z0_o[varx] + C_o[varz][vary]*z0_o[vary] + C_o[varz][varz]*z0_o[varz];
                f0_i[varx] = C_i[varx][varx]*z0_i[varx] + C_i[varx][vary]*z0_i[vary] + C_i[varx][varz]*z0_i[varz];
                f0_i[vary] = C_i[vary][varx]*z0_i[varx] + C_i[vary][vary]*z0_i[vary] + C_i[vary][varz]*z0_i[varz];
                f0_i[varz] = C_i[varz][varx]*z0_i[varx] + C_i[varz][vary]*z0_i[vary] + C_i[varz][varz]*z0_i[varz];
            }
            else if (sysdim >=2 && skew)
            {
                f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary];
                f0_o[vary] = C_o[vary][vary]*z0_o[vary] + C_o[vary][varx]*z0_o[varx];
                f0_i[varx] = C_i[varx][varx]*z0_i[varx] + C_i[varx][vary]*z0_i[vary];
                f0_i[vary] = C_i[vary][vary]*z0_i[vary] + C_i[vary][varx]*z0_i[varx];
            }
            
            multMiMi(CJacf_o,C_o,Jacf_o);
            multMiMi(CJacf_i,C_i,Jacf_i);
            
            if (order == 2)
            {
                
                multMiMi(Cdfdx0_o,C_o,dfdx0_o);
                multMiMi(Cdfdx0_i,C_i,dfdx0_i);
                
                for (int j=0 ; j<jacdim; j++) {
                    for (int k=0 ; k<jacdim; k++) {
                        for (int i=0 ; i<sysdim; i++) {
                            CHessf_o[i][j][k] = Hessf_o[i][j][k];
                            CHessf_i[i][j][k] = Hessf_i[i][j][k];
                        }
                        if (sysdim >=3 && skew) {
                            CHessf_o[varx][j][k] = C_o[varx][varx]*Hessf_o[varx][j][k] + C_o[varx][vary]*Hessf_o[vary][j][k] + C_o[varx][varz]*Hessf_o[varz][j][k];
                            CHessf_i[varx][j][k] = C_i[varx][varx]*Hessf_i[varx][j][k] + C_i[varx][vary]*Hessf_i[vary][j][k] + C_i[varx][varz]*Hessf_i[varz][j][k];
                            CHessf_o[vary][j][k] = C_o[vary][varx]*Hessf_o[varx][j][k] + C_o[vary][vary]*Hessf_o[vary][j][k] + C_o[vary][varz]*Hessf_o[varz][j][k];
                            CHessf_i[vary][j][k] = C_i[vary][varx]*Hessf_i[varx][j][k] + C_i[vary][vary]*Hessf_i[vary][j][k] + C_i[vary][varz]*Hessf_i[varz][j][k];
                            CHessf_o[varz][j][k] = C_o[varz][varx]*Hessf_o[varx][j][k] + C_o[varz][vary]*Hessf_o[vary][j][k] + C_o[varz][varz]*Hessf_o[varz][j][k];
                            CHessf_i[varz][j][k] = C_i[varz][varx]*Hessf_i[varx][j][k] + C_i[varz][vary]*Hessf_i[vary][j][k] + C_i[varz][varz]*Hessf_i[varz][j][k];
                        }
                        else if (sysdim >=2 && skew) {
                            CHessf_o[varx][j][k] = C_o[varx][varx]*Hessf_o[varx][j][k] + C_o[varx][vary]*Hessf_o[vary][j][k];
                            CHessf_i[varx][j][k] = C_i[varx][varx]*Hessf_i[varx][j][k] + C_i[varx][vary]*Hessf_i[vary][j][k];
                            CHessf_o[vary][j][k] = C_o[vary][varx]*Hessf_o[varx][j][k] + C_o[vary][vary]*Hessf_o[vary][j][k];
                            CHessf_i[vary][j][k] = C_i[vary][varx]*Hessf_i[varx][j][k] + C_i[vary][vary]*Hessf_i[vary][j][k];
                        }
                    }
                }
                
            }
            
            // outer range
            vector<interval> temp_outer;
            if (order == 1)
                temp_outer = evaluate_outerrange(f0_o,radx_o,CJacf_o);
            else  if (order == 2)
                temp_outer = evaluate_outerrange_order2(f0_o, radx_o, Cdfdx0_o, CHessf_o);
            
            
            
            vector<interval> temp_inner;
            if (order == 1)
                temp_inner = evaluate_innerrange(f0_i,radx_i,CJacf_i,false,exist_quantified);
            else if (order == 2)
                temp_inner = evaluate_innerrange_order2(f0_i, radx_i, Cdfdx0_i, CHessf_i,false,exist_quantified);
            
            if (step % printing_period == 0) {
            if (skew)
            {
            out_approx << YAML::Key << "outer2d";
            out_approx << YAML::Value;
            out_approx << YAML::BeginSeq;
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << varx;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << vary;
            cout << "outer skewed box (mean value): ";
            output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A_o,  varx,  vary, step);
            out_approx << YAML::EndMap;
            out_approx << YAML::EndSeq;
            }
            
            
            out_approx << YAML::Key << "inner2d";
            out_approx << YAML::Value;
            out_approx << YAML::BeginSeq;
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << varx;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << vary;
            if (skew)
            {
                cout << "inner skewed box (mean value): ";
                print_pi(exist_quantified);
                output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A_i,  varx,  vary, step);
            }
            else
                print_innerbox(temp_inner, exist_quantified, varx, vary, step);
            out_approx << YAML::EndMap;
            out_approx << YAML::EndSeq;
            }
            
            z_inner = temp_inner;
            z_outer = temp_outer;
          
        }
        else
            z_inner = z_inner_proj;
        
        for (int i=0; i < sysdim ; i++) {
            x_i[i] = z_inner[i];
            x0_i[i] =  mid(z_inner[i]); // f0_i[i]; //
            radx_i[i] = z_inner[i] - x0_i[i];
            x_o[i] = z_outer[i];
            x0_o[i] =  mid(z_outer[i]); // f0_i[i]; //
            radx_o[i] = z_outer[i] - x0_o[i];
        }
        if (step % printing_period == 0)
            out_approx << YAML::EndMap;
    }
    
    ReachSet res = ReachSet(z_outer,z_inner);
    return res;
}



// computing at each step the sensitivity with respect to initial values
// same method as discrete_dynamical_method2 but using preconditioning uniquely for printing
ReachSet discrete_dynamical_method2(DiscreteFunc &f, vector<interval> &xinit, vector<vector<interval>> &estimated_range, int &nb_steps, bool skew)
{
    vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
 //   DiscreteFunc f;
    FDiff FFunc;
    
    
    
    //   estimate_range(f,xinit);
   // vector<vector<interval>> range(nb_steps+1,vector<interval>(sysdim));
  //  range =  estimate_reachset(f, nb_steps, xinit);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
 
    vector<vector<double>> output_skewedbox;
    
 
    
    // for now, no input, sysdim = jacdim
    z_outer = xinit;
    
    // preconditioning: A is center of Jacobian, C its inverse
    vector<vector<double>> A_o(sysdim,vector<double> (sysdim)), C_o(sysdim,vector<double> (sysdim));
    vector<interval> f0_o(sysdim);
    vector<vector<interval>> CJacf_o(sysdim, vector<interval>(jacdim));

    for (int i=0 ; i<sysdim; i++)
        A_o[i][i] = 1.0;
    
   
    
    for (int i=0; i < jacdim ; i++) {
        x_o[i] = z_outer[i];
        x0_o[i] = mid(z_outer[i]); //+(eps[i].sup()-mid(eps[i]))/2.0;
        radx_o[i] = z_outer[i] - x0_o[i];
    }
    
    for (int i=0; i < jacdim ; i++)
        x_o[i].diff(i,jacdim);
    
    vector<interval> z0_o = f(x0_o);
    for (int i=0; i < sysdim ; i++)
        x0_o[i] = z0_o[i];
    
    for (int step = 1; step <= nb_steps ; step++)
    {
        if (step % printing_period == 0)
        {
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "tn";
            out_approx << YAML::Value << step;
        }
        
        z_o = f(x_o);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++) {
                Jacf_o[i][j] = z_o[i].d(j).convert_int(); // JacAff_o[i][j].convert_int();
            }
        
        cout << "\n Step " << step << ": " << endl;
        cout << "Outer approx of f(x), direct evaluation:         ";
        print_outer_range(z_o,estimated_range[step]);
        
        
        /*********************************** Evaluation ***********************************/
        
        vector<int> aux;
        
        z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
        z_inner_proj = evaluate_innerrange(z0_o, radx_o, Jacf_o, true, aux);
        if (jacdim > sysdim)
            z_inner_proj_rob = evaluate_innerrange_robust(z0_o, radx_o, Jacf_o, true, aux);
        
        for (int i=0; i < sysdim ; i++)
            z_outer[i] = intersect(z_outer[i],z_o[i].x().convert_int());
        
        if (step % printing_period == 0)
            print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step,estimated_range[step]);
        
        
        if (sysdim >= 2) {
            // fixed for now
            //   assert(sysdim <= 3); // higher dimension not treated yet
            for (int i=0 ; i<sysdim; i++)
                exist_quantified[i] = i;
            
            if (syschoice == 18) {
                exist_quantified[3] = 1;  // has to be 3 ? Here I don't think so because we can evaluate at each step independently ? I don't understand !!!
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            else if (syschoice == 20) {
                exist_quantified[2] = 1;
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            
            
            for (int i=0 ; i<sysdim; i++) {
                A_o[i][i] = 1.0;
                C_o[i][i] = 1.0;
            }
            
            // f0 = C * z0
            for (int i=0 ; i<sysdim; i++)
                f0_o[i] = z0_o[i];
            
            int varx = 0, vary = 1, varz = 2; // joint range we want to print/compute
            
            if (step % printing_period == 0)
            {
            // SKEWED OUTER-APPROX
            if (skew)
            {
                
                out_approx << YAML::Key << "outer2d";
                out_approx << YAML::Value;
                out_approx << YAML::BeginSeq;
                
                
                for (varx=0; varx<sysdim; varx++)
                    for (vary=varx+1;vary<sysdim;vary++)
                    {
                        
                        for (int i=0 ; i<sysdim; i++) {
                            A_o[i][i] = 1.0;
                            C_o[i][i] = 1.0;
                        }
                    // 2D preconditioner - just on components varx and vary
                    A_o[varx][varx] = mid(Jacf_o[varx][varx]);
                    A_o[vary][vary] = mid(Jacf_o[vary][vary]);
                    A_o[varx][vary] = mid(Jacf_o[varx][vary]);
                    A_o[vary][varx] = mid(Jacf_o[vary][varx]);
           
                    // C is inverse of A
                    double determinant = 1.0/(A_o[varx][varx]*A_o[vary][vary]-A_o[varx][vary]*A_o[vary][varx]);
                    C_o[varx][varx] = determinant*A_o[vary][vary];
                    C_o[varx][vary] = - determinant*A_o[varx][vary];
                    C_o[vary][varx] = - determinant*A_o[vary][varx];
                    C_o[vary][vary] = determinant*A_o[varx][varx];
                    
                    for (int i=0 ; i<sysdim; i++)
                        f0_o[i] = z0_o[i];
                        
                    f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary];
                    f0_o[vary] = C_o[vary][vary]*z0_o[vary] + C_o[vary][varx]*z0_o[varx];
              
                
                // CJacf = C * Jacf
                multMiMi(CJacf_o,C_o,Jacf_o);
                
                
                // outer 2D range
                vector<interval> temp_outer = evaluate_outerrange(f0_o,radx_o,CJacf_o);
                
                out_approx << YAML::BeginMap;
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << varx;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << vary;
                cout << "outer skewed box (mean value) for x" << varx+1 <<" and x" << vary+1 <<": ";
                output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A_o,  varx,  vary, step);
                out_approx << YAML::EndMap;
                
                }
                out_approx << YAML::EndSeq;
            }
            
            // inner 2D range
            
            
            out_approx << YAML::Key << "inner2d";
            out_approx << YAML::Value;
            out_approx << YAML::BeginSeq;
            
            for (varx=0; varx<sysdim; varx++)
                for (vary=varx+1;vary<sysdim;vary++)
                {
                    
                        for (int i=0 ; i<sysdim; i++) {
                            A_o[i][i] = 1.0;
                            C_o[i][i] = 1.0;
                        }
             //   if (sysdim == 2) {
                    // 2D preconditioner - just on components varx and vary
                    A_o[varx][varx] = mid(Jacf_o[varx][varx]);
                    A_o[vary][vary] = mid(Jacf_o[vary][vary]);
                    A_o[varx][vary] = mid(Jacf_o[varx][vary]);
                    A_o[vary][varx] = mid(Jacf_o[vary][varx]);
           
                    // C is inverse of A
                    double determinant = 1.0/(A_o[varx][varx]*A_o[vary][vary]-A_o[varx][vary]*A_o[vary][varx]);
                    C_o[varx][varx] = determinant*A_o[vary][vary];
                    C_o[varx][vary] = - determinant*A_o[varx][vary];
                    C_o[vary][varx] = - determinant*A_o[vary][varx];
                    C_o[vary][vary] = determinant*A_o[varx][varx];
                    
                    for (int i=0 ; i<sysdim; i++)
                        f0_o[i] = z0_o[i];
                        
                    f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary];
                    f0_o[vary] = C_o[vary][vary]*z0_o[vary] + C_o[vary][varx]*z0_o[varx];
                
                // CJacf = C * Jacf
                multMiMi(CJacf_o,C_o,Jacf_o);
                        
                vector<interval> temp_inner = evaluate_innerrange(f0_o,radx_o,CJacf_o,false,exist_quantified);
                
                
                out_approx << YAML::BeginMap;
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << varx;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << vary;
                if (skew)
                {
                    cout << "inner skewed box (mean value): ";
                    print_pi(exist_quantified);
                    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A_o,  varx,  vary, step);
                }
                else
                    print_innerbox(temp_inner, exist_quantified, varx, vary, step);
                out_approx << YAML::EndMap;
                
            }
            out_approx << YAML::EndSeq;
        }
        
        } // end if (step % printing_period == 0)
        
        
        // initialize next iteration: z0 = f^n(x0)
        z0_o = f(x0_o);
        for (int i=0; i < sysdim ; i++)
            x0_o[i] = z0_o[i];
        
        for (int i=0; i < sysdim ; i++) {
            x_o[i] = z_o[i];
            //   x0_o[i] = mid(z_o[i].x().convert_int()); //+(eps[i].sup()-mid(eps[i]))/2.0;
            //       radx_o[i] = z_o[i].x().convert_int() - z0_o[i];
        }
        if (step % printing_period == 0)
                out_approx << YAML::EndMap;
    }
    
    ReachSet res = ReachSet(z_outer,z_inner);
    return res;
    
}





void function_range(DiscreteFunc &f, vector<interval> &xinit, vector<interval> &estimated_range) {
    
    int nb_steps;
   // vector<interval> xinit = init_discrete_system(nb_steps, config_filename); // initial condition
    
       
    nb_discr = 10;
    nb_discr1 = nb_discr;
    nb_discr2 = nb_discr;
    
    
    vector<F<AAF>> x(jacdim);
    vector<vector<interval>> Jacf(sysdim,vector<interval>(jacdim));
    vector<vector<AAF>> JacAff(sysdim,vector<AAF>(jacdim));
    vector<F<AAF>> z(sysdim), z1(sysdim), z2(sysdim);
    vector<interval> x0(jacdim), radx(jacdim);    // center and radius x-x0
    
    
//    DiscreteFunc f;
    
    //vector<interval> range(sysdim);
   // range = estimate_range(f,xinit);
   // range = estimate_reachset(f, 1, xinit);
    
    for (int i=0; i < jacdim ; i++) {
        x[i] = xinit[i];
        x0[i] = mid(xinit[i]); //+(eps[i].sup()-mid(eps[i]))/2.0;
        radx[i] = xinit[i] - x0[i];
    }
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    
    // vector<AAF> x_nn(jacdim);
    // vector<vector<AAF>> JacAff_nn(sysdim,vector<AAF>(jacdim));
  /*   vector<vector<F<AAF>>> JacAff_nn(sysdim,vector<F<AAF>>(jacdim));
     vector<vector<interval>> Jacf_nn(sysdim,vector<interval>(jacdim));
     vector<vector<F<AAF>>> net_outputs(NH.n_hidden_layers+2,vector<F<AAF>>(jacdim));
     
     // for (int i=0; i < jacdim ; i++)
     net_outputs[0] = x; // [i] = xinit[i];
     
     for (int k=0 ; k<NH.n_hidden_layers+1 ; k++ ) {
         net_outputs[k+1] = NH.L[k].eval_linear_layer(net_outputs[k]);
         for (int i=0; i < sysdim ; i++) {
             for (int j=0; j < jacdim ; j++) {  // ===> ATTENTION PAS JACDIM ICI MAIS LA DIMENSION INTERNE
                 JacAff_nn[i][j] = net_outputs[k+1][i].d(j); // x[i].d[j]
                 Jacf_nn[i][j] = JacAff_nn[i][j].x().convert_int();
             }
             net_outputs[k+1][i] = eval_activation(NH.L[k].activation,net_outputs[k+1][i]);
         }
         cout << "Jacf_nn evaluated on [x]" << endl;
         cout << Jacf_nn;
         for (int i=0; i < sysdim ; i++) {
             F<AAF> temp = net_outputs[k+1][i].x()*(1-net_outputs[k+1][i].x());
             for (int j=0; j < jacdim ; j++) {
                 JacAff_nn[i][j] =  temp*JacAff_nn[i][j]; // chain rule derivation (with s'(x) = s(x) (1 - s(x))) // net_outputs[k+1][i].d(j); //
                 Jacf_nn[i][j] = JacAff_nn[i][j].x().convert_int();
             }
         }
         cout << "Jacf_nn evaluated on [x]" << endl;
         cout << Jacf_nn;
     }
     vector<F<AAF>> z_nn = net_outputs[NH.n_hidden_layers+1];
         
     cout << "Outer approx of network, direct evaluation: " << endl;
     for (int i=0; i < sysdim ; i++)
         cout << z_nn[i].x().convert_int() << endl;
     
      for (int i=0; i < sysdim ; i++)
         for (int j=0; j < jacdim ; j++) {
             Jacf_nn[i][j] = JacAff_nn[i][j].x().convert_int();
         //    JacAff[i][j] = z[i].d(j);
          //           cout << "JacAff[i][j]=" << JacAff[i][j] << JacAff[i][j].convert_int();
     }
      cout << "Jacf_nn evaluated on [x]" << endl;
      cout << Jacf_nn;
     */
    
    
 /*   interval val = interval(-1,1);
    interval val2 = 1.0/(1.0 + exp(-val));
    cout << "interval image by sigmoid of [-1,1]=" << val2 << endl;
    
    AAF Aval = interval(-1,1);
    AAF Aval2 = 1.0/(1.0 + exp(-Aval));
    cout << "affine form range by sigmoid of [-1,1], CHEBYSCHEV mode=" << Aval2.convert_int()<< endl;
    cout << "affine form by sigmoid of [-1,1], CHEBYSCHEV mode=" << Aval2 << endl;
    
    val = interval(-2,2);
    val2 = 1.0/(1.0 + exp(-val));
    cout << "interval image by sigmoid of [-2,2]=" << val2 << endl;
    Aval = interval(-2,2);
    try
    {
    Aval2 = 1.0/(1.0 + exp(-Aval));
    cout << "affine form range by sigmoid of [-2,2], CHEBYSCHEV mode=" << Aval2.convert_int() << endl;
    cout << "affine form by sigmoid of [-2,2], CHEBYSHEV mode=" << Aval2 << endl;
    } catch(std::exception const& e)
    {
        cout << "ERREUR computing affine form by sigmoid of [-2,2] in CHEBYSHEV mode: " << e.what() << endl;
    }
    Aval2.setApproximationType(MINRANGE);
    Aval2 = 1.0/(1.0 + exp(-Aval));
    cout << "affine form range by sigmoid of [-2,2], MINMAX mode=" << Aval2.convert_int() << endl;
    cout << "affine form by sigmoid of [-2,2], MINMAX mode=" << Aval2 << endl;
    
    Aval2.setApproximationType(CHEBYSHEV);
    val = interval(-0.5,0.5);
    val2 = 1.0/(1.0 + exp(-val));
    cout << "interval image by sigmoid of [-0.5,0.5]=" << val2 << endl;
    Aval = interval(-0.5,0.5);
    Aval2 = 1.0/(1.0 + exp(-Aval));
    cout << "affine form range by sigmoid of [-0.5,0.5], CHEBYSCHEV mode=" << Aval2.convert_int() << endl;
    cout << "affine form by sigmoid of [-0.5,0.5], CHEBYSCHEV mode=" << Aval2 << endl;
    
    val = interval(0,1.);
    val2 = 1.0/(1.0 + exp(-val));
    cout << "interval image by sigmoid of [0,1]=" << val2 << endl;
    Aval = interval(0,1.);
    Aval2 = 1.0/(1.0 + exp(-Aval));
    cout << "affine form range by sigmoid of [0,1], CHEBYSCHEV mode=" << Aval2.convert_int() << endl;
    cout << "affine form by sigmoid of [0,1], CHEBYSCHEV mode=" << Aval2 << endl;
  */
    
    for (int i=0 ; i<sysdim; i++)
        z[i].x().setApproximationType(MINRANGE);
    z = f(x);
    
    //cout << "z=" <<z[0].x().convert_int();
   // cout << endl;
    
    cout << "Outer approx of f(x), direct evaluation: ";
    interval interv_out;
    for (int i=0; i < sysdim ; i++) {
        interv_out = z[i].x().convert_int();
        cout << "[" << interv_out.inf() << ", " << interv_out.sup() << "] ";
        cout << " eta["<<i<<"]=" << (estimated_range[i].sup() - estimated_range[i].inf())/ (interv_out.sup() - interv_out.inf()) <<" "; // precision metric with respect to estimated range
    }
    cout << endl;
    
    
    for (int i=0; i < sysdim ; i++)
        for (int j=0; j < jacdim ; j++) {
            Jacf[i][j] = z[i].d(j).convert_int();
            JacAff[i][j] = z[i].d(j);
            cout << "JacAff[i][j]=" << JacAff[i][j] << JacAff[i][j].convert_int();
        }
    cout << "Jacf evaluated on [x]: ";
    cout << Jacf;
    
    
   
 
    
    /*********************************** ORDER 2 Taylor model  ***********************************/
   
    FDiff FFunc;
  /*  vector<vector<AAF>> dfdx(sysdim,vector<AAF>(jacdim));
    vector<AAF> zf = FFunc(dfdx,initial_values);
    cout << "cf= " << zf[0].convert_int() << endl;
    cout << "dfdx= " << dfdx[0][0].convert_int() << endl; */
    
    vector<vector<interval>> dfdx0(sysdim,vector<interval>(jacdim));
    vector<interval> zf0 = FFunc(dfdx0,x0);
  //  cout << "x0=" << x0;
 //   cout << "cf= " << zf0;
  //  cout << "dfdx= " << dfdx0;
   // cout << "dfdx= " << dfdx0[1];
    
    
    vector<vector<vector<interval>>> Hessf(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector<vector<vector<AAF>>> HessAff(sysdim,vector<vector<AAF>>(jacdim,vector<AAF>(jacdim)));
    
    
    // (tres lointainement) inspiré de ExampleBADFAD1.cpp
    vector <F<F<AAF>>> xff(jacdim), zff(sysdim);
    
    for (int i=0; i < jacdim ; i++) {
        xff[i] = xinit[i];
        xff[i].diff(i,jacdim);          // first order
        xff[i].x().diff(i,jacdim);      // second order
    }
    zff = f(xff);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            for (int k=0; k < jacdim ; k++) {
                Hessf[i][j][k] =zff[i].d(j).d(k).convert_int();
                HessAff[i][j][k] =zff[i].d(j).d(k);
            }
            //    cout << "JacAff[i][j]=" << JacAff[i][j] << JacAff[i][j].convert_int();
        }
    }
    //    cout << "Jacf evaluated on all [x]" << endl;
    //    cout << Jacf;
 //   cout << "Hessian evaluated on all [x]" << endl;
 //   for (int i=0; i < sysdim ; i++)
  //      cout << Hessf[i];
    
    
    
    /*********************************** Evaluation ***********************************/
    
    vector<interval> z0 = f(x0);
    
 evaluate_projections(z0, radx, Jacf,estimated_range);
    //  evaluate_projections_subdiv(z0, JacAff);
    //  evaluate_projections_subdiv_discretize(z0, JacAff);
evaluate_projections_discretize_simultaneous(z0, radx, JacAff,estimated_range);
    // center and order 1 evaluated on x0, 2nd term on [x]
evaluate_projections_order2(z0, radx, dfdx0, Hessf, estimated_range);
    
    // here we actually want to use JacAff and not dfdx0 except for first subdivision
    // A LA FOIS FAUX ET IMPRECIS POUR LE MOMENT: A REPRENDRE
 evaluate_projections_order2_discretize_simultaneous(z0, radx, JacAff, dfdx0, HessAff, estimated_range);
    
    if (sysdim >= 2) {
 joint_ranges(z0,radx,Jacf,dfdx0,Hessf,0,1);
        //       twodim_discretization_by_quadrant();
        //      joint_ranges_subdiv(z0,JacAff,0,1);
        // joint_ranges_subdiv_discretize(z0,JacAff,0,1);
// joint_ranges_discretize_simultaneous(z0,radx,JacAff,0,1);
        preconditioned_joint_ranges(z0,radx,Jacf,dfdx0,Hessf,0,1);
        // preconditioned_joint_ranges_subdiv(z0,JacAff,0,1);
        //      preconditioned_joint_ranges_subdiv_discretize(z0,JacAff,0,1);
        preconditioned_joint_ranges_discretize_simultaneous(z0,radx,JacAff,dfdx0,HessAff,0,1);
    }
   
    if (systype == 3) {
        
        cout << "Direct range evaluation for neural network: ";
        
        if (NH.n_inputs > 0)
        { // network given as .sfx file
            // direct range by AAF
            vector<vector<AAF>> net_outputs(NH.n_hidden_layers+2); // outputs for each layer
            //vector<AAF> net_inputs(NH.n_inputs);
            net_outputs[0] = vector<AAF>(NH.n_inputs);
            //net_inputs[0] = interval(0,0.2);
            //net_inputs[1] = interval(0,0.2);
            
            for (int i=0 ; i<net_outputs[0].size(); i++)
                net_outputs[0][i] = xinit[i];
            for (int i=0 ; i<NH.n_hidden_layers+1 ; i++ ) {
                net_outputs[i+1] = NH.L[i].eval_layer(net_outputs[i]);
            }
            for (int i=-1 ; i<NH.n_hidden_layers+1 ; i++ ) {
                cout << "output layer " << i << " is " << net_outputs[i+1] << endl;
            }
        }
#if ONNX_active
        else if (CG.no_of_output_nodes > 0)
        {
            map<uint32_t, AAF> in, out;
            vector<AAF> net_outputs(CG.no_of_output_nodes);
            for (int i=1; i<=CG.no_of_input_nodes; i++)
                in.insert(make_pair(i, xinit[i-1]));
            computation_graph_evaluate_graph_abstract(CG,in,out);
            //    cout << "CG.no_of_input_nodes=" << CG.no_of_input_nodes << endl;
            //    cout << "CG.no_of_output_nodes=" << CG.no_of_output_nodes << endl;
            vector< uint32_t > in_nodes, out_nodes;
            CG.return_id_of_input_output_nodes(in_nodes , out_nodes );
            for (int i=0; i<CG.no_of_output_nodes; i++) {
                //  cout << "out_node=" << out_nodes[i] << endl;
                net_outputs[i] = out[out_nodes[i]];
                //cout << z[i];
            }
            cout << "output is " << net_outputs << endl;
        }
#endif
    }
    

    
//    system("cd GUI; python3 Visu_function.py; cd ..");
}



 
    
 






// range is the estimated range of system
void evaluate_projections(vector<interval> &z0, vector<interval> &radx,  vector<vector<interval>> &Jacf, vector<interval> &range) {
    interval inner_impro;
    vector<interval> z_inner, z_outer;
    vector<int> aux;
    
    z_outer = evaluate_outerrange(z0, radx, Jacf);
    z_inner = evaluate_innerrange(z0, radx, Jacf, true, aux);
    
  //  for (int i=0; i < sysdim ; i++) {
  //      z_outer[i] = evaluate_outerrange_x(z0, Jacf, i);
  //      z_inner[i] = evaluate_innerrange_x(z0, Jacf, true, aux, i);
  //  }
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value:                                           ";
  //  cout << z_outer;
   // cout << "gamma=";
    for (int i=0; i < sysdim ; i++)
        cout << " z_o[" << i << "]=" << z_outer[i] << " eta_o[" << i << "]=" << (range[i].sup() - range[i].inf())/ (z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to estimated range"
    cout << endl;
    cout << "Projection of inner range of f(x) by mean-value:                             ";
   // cout << z_inner;
    for (int i=0; i < sysdim ; i++)
        cout << " z_i[" << i << "]=" << z_inner[i] << " eta_i[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(range[i].sup() - range[i].inf()) << " gamma[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to
    cout << endl;
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}


void evaluate_projections_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf, vector<interval> &range) {
    interval inner_impro;
    vector<interval> z_inner, z_outer;
    vector<int> aux;
    
    z_outer = evaluate_outerrange_order2(z0, radx, Jacf, Hessf);
    z_inner = evaluate_innerrange_order2(z0, radx, Jacf, Hessf, true, aux);
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by order 2 Taylor model:                                 ";
    for (int i=0; i < sysdim ; i++)
        cout << " z_o[" << i << "]=" << z_outer[i] << " eta_o[" << i << "]=" << (range[i].sup() - range[i].inf())/ (z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to estimated range"
    cout << endl;
    //cout << z_outer;
    cout << "Projection of inner range of f(x) by order 2 Taylor model:                   ";
    for (int i=0; i < sysdim ; i++)
        cout << " z_i[" << i << "]=" << z_inner[i] << " eta_i[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(range[i].sup() - range[i].inf()) << " gamma[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to
    cout << endl;
//    cout << z_inner;
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}



// index1 is the first index that we subdivide, n1 the number of associated subdivisions
// evaluate projections on the 4 quadrant and join them
void evaluate_projections_subdiv(vector<interval> &z0, vector<interval> &radx,  vector<vector<AAF>> &JacAff) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = evaluate_outerrange_x_subdiv(z0, radx, JacAff, i);
        z_inner[i] = z0[i];
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++) {
                temp_inner = evaluate_innerrange_x_subdiv(z0, radx, JacAff, true, aux, i, i1, i2); // maximal inner approx on quadrant (i1,i2)
                // on maximal inner approx we can do interval hull
                z_inner[i] = interval_hull(z_inner[i],temp_inner);
            }
        }
    }
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value with 2 subdiv each input:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value with 2 subdiv each input: ";
    cout << z_inner;
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}

// index1 is the first index that we subdivide, n1 the number of associated subdivisions
// evaluate projections on the 4 quadrant and join them
void evaluate_projections_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = evaluate_outerrange_x_subdiv_discretize(z0, radx, JacAff, i);
        z_inner[i] = z0[i];
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++) {
                temp_inner = evaluate_innerrange_x_subdiv_discretize(z0, radx, JacAff, true, aux, i, i1, i2); // maximal inner approx on quadrant (i1,i2)
                // on maximal inner approx we can do interval hull
                z_inner[i] = interval_hull(z_inner[i],temp_inner);
            }
        }
    }
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value with 2 subdiv each input + discretization:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value with 2 subdiv each input + discretization: ";
    cout << z_inner;
    
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}


// evaluate projections
void evaluate_projections_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<interval> &range) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    z_outer = evaluate_outerrange_discretize_simultaneous(z0, radx, JacAff);
    z_inner = evaluate_innerrange_discretize_simultaneous(z0, radx, JacAff, true, aux);
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by centered mean-value with discretization:              ";
    for (int i=0; i < sysdim ; i++)
        cout << " z_o[" << i << "]=" << z_outer[i] << " eta_o[" << i << "]=" << (range[i].sup() - range[i].inf())/ (z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to estimated range"
    cout << endl;
   // cout << z_outer;
    cout << "Projection of inner range of f(x) by centered mean-value with discretization:";
    for (int i=0; i < sysdim ; i++)
        cout << " z_i[" << i << "]=" << z_inner[i] << " eta_i[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(range[i].sup() - range[i].inf()) << " gamma[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to
    cout << endl;
//    cout << z_inner;
    
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}


// evaluate projections
void evaluate_projections_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<interval>> &dfdx0, vector<vector<vector<AAF>>> &HessAff, vector<interval> &range) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    z_outer = evaluate_outerrange_order2_discretize_simultaneous(z0, radx, JacAff, dfdx0, HessAff);
    z_inner = evaluate_innerrange_order2_discretize_simultaneous(z0, radx, JacAff, dfdx0, HessAff, true, aux);
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by order 2 model with discretization:                    ";
    for (int i=0; i < sysdim ; i++)
        cout << " z_o[" << i << "]=" << z_outer[i] << " eta_o[" << i << "]=" << (range[i].sup() - range[i].inf())/ (z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to estimated range"
    cout << endl;
 //   cout << z_outer;
    cout << "Projection of inner range of f(x) by order 2 model with discretization:      ";
    for (int i=0; i < sysdim ; i++)
        cout << " z_i[" << i << "]=" << z_inner[i] << " eta_i[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(range[i].sup() - range[i].inf()) << " gamma[" << i << "]=" << (z_inner[i].sup() - z_inner[i].inf())/(z_outer[i].sup() - z_outer[i].inf()) <<" "; // precision metric with respect to
    cout << endl;
   // cout << z_inner;
    
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}



vector<interval> evaluate_outerrange(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf)
{
    vector<interval> z_outer(sysdim);
    
    // z_outer = z0 + Jacf*eps
    multMiVi(z_outer,Jacf,radx);
    addViVi(z_outer,z0);
    
    /*for (int i=0 ; i< sysdim ; i++)
    {
        z_outer[i] = z0[i];
        
        for (int j=0; j < jacdim ; j++) {
            z_outer[i] += Jacf[i][j]*eps[j];
        }
    } */
    return z_outer;
}


vector<interval> evaluate_outerrange_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf)
{
    
    vector<interval> z_outer(sysdim), z_temp(sysdim);
    
    // z0 + Jacf(x0) . x
    multMiVi(z_outer,Jacf,radx);
  //  cout << "z_outer" << z_outer;
    addViVi(z_outer,z0);
   // cout << "z_outer" << z_outer;
    
    // < x , Hessf([x]) . x >
    for (int i=0 ; i< sysdim ; i++)
    {
     //   multMiVi(z_temp,Hessf[i],eps);
     //   scalarproduct(z_temp,z_temp,eps);
        // replaced by below to exploit that eps^2 can only be positive
        z_temp[i] = 0;
        for (int j=0 ; j< jacdim ; j++) {
            for (int k=0 ; k< j ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
            z_temp[i] += Hessf[i][j][j]*interval(0,sup(radx[j])*sup(radx[j])) / 2.0;
            for (int k=j+1 ; k< jacdim ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
        }
    }
     addViVi(z_outer,z_temp);
  //  cout << "z_outer" << z_outer;
    return z_outer;
}


interval evaluate_outerrange_x_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int i) {
    interval z_temp, z_outer = z0[i];
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
            z_temp = z0[i];
            for (int j=0; j < jacdim ; j++) {
                z_temp += JacAff[i][j].convert_int(constr_eps[i1][i2],jacdim) * eps_loc[i1][i2][j];
            }
            
            z_outer = interval_hull(z_outer,z_temp);
        }
        
    }
    
    return z_outer;
}

// subdivise n times in the 2 inputs
// I think trouble is I progress in crabe (synchronously on the 2 components)
/*
interval evaluate_outerrange_x_subdiv_discretize_old(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, z_outer = z0[i];
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
            z_temp_partial = z0[i];
            // evaluated at the extremity on the n-1 first subdivisions
            for (int m=0; m < nb_discr-1 ; m++) {
                z_temp = z_temp_partial;
                for (int j=0; j < jacdim ; j++) {
                    if (m>0)
                    {
                        z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
                        z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (extremity_eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
                    }
                    else
                    {
                        z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * eps_loc_discr[i1][i2][m][j];
                        z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * extremity_eps_loc_discr[i1][i2][m][j];
                    }
                }
                // hull on all subdivisions
                z_outer = interval_hull(z_outer,z_temp);
            }
            // evaluated on interval on the last subdivision
            z_temp = z_temp_partial;
            for (int j=0; j < jacdim ; j++) {
                z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][nb_discr-1],jacdim) * (eps_loc_discr[i1][i2][nb_discr-1][j]-extremity_eps_loc_discr[i1][i2][nb_discr-2][j]);
            }
            
            z_outer = interval_hull(z_outer,z_temp);
        }
    }
    return z_outer;
}
 */

// the vector of constraints on eps
vector<interval> constraint_eps(int i1, int i2, int m1, int m2) {
    vector<interval> res(jacdim);
    
    for (int j=0; j < jacdim ; j++)
        res[j] = interval(-1,1);
    
    if (constr_eps[i1][i2][0].sup() > 0)
        res[0] = interval(((double)m1)/nb_discr1,(m1+1.0)/nb_discr1);
    else
        res[0] = -interval(((double)m1)/nb_discr1,(m1+1.0)/nb_discr1);
    
    if (constr_eps[i1][i2][1].sup() > 0)
        res[1] = interval(((double)m2)/nb_discr2,(m2+1.0)/nb_discr2);
    else
        res[1] = -interval(((double)m2)/nb_discr2,(m2+1.0)/nb_discr2);
    
    return res;
}

//interval concretize_eps(int i1, int i2, int m1, int m2, int k) {
//    return eps[k].mid() + eps[k].rad() * constr_eps_discr[i][j][m][k];
//}

double extremity_eps_max(vector<interval> &radx, int i1, int i2, int m1, int m2, int k) {
    if (constr_eps[i1][i2][k].sup() > 0) {
        if (k == 0)
            return radx[k].mid() + radx[k].rad() * (m1+1.0)/nb_discr1;
        else
            return radx[k].mid() + radx[k].rad() * (m2+1.0)/nb_discr2;
    }
    else {
        if (k == 0)
            return radx[k].mid() - radx[k].rad() * (m1+1.0)/nb_discr1;
        else
            return radx[k].mid() - radx[k].rad() * (m2+1.0)/nb_discr2;
    }
}

double extremity_eps_min(vector<interval> &radx, int i1, int i2, int m1, int m2, int k) {
    if (constr_eps[i1][i2][k].sup() > 0) {
        if (k == 0)
            return radx[k].mid() + radx[k].rad() * (double(m1))/nb_discr1;
        else
            return radx[k].mid() + radx[k].rad() * (double(m2))/nb_discr2;
    }
    else {
        if (k == 0)
            return radx[k].mid() - radx[k].rad() * (double(m1))/nb_discr1;
        else
            return radx[k].mid() - radx[k].rad() * (double(m2))/nb_discr2;
    }
}

/*
// subdivise n times in the 2 inputs
// I think trouble is I progress "in crab" (synchronously on the 2 components)
interval evaluate_outerrange_x_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, loc_eps, extremity_loc_eps_min, extremity_loc_eps_max, z_outer = z0[i];
    vector<interval> c_eps;
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
         //   z_temp_partial = z0[i];
            z_temp = z0[i]; // z_temp_partial;
            for (int m1=0; m1 < nb_discr1 ; m1++) {
                for (int m2=0; m2 < nb_discr2 ; m2++) {
                    c_eps = constraint_eps(i1,i2,m1,m2);
                    cout << "m1=" << m1 << " m2=" << m2 << " c_eps=" << c_eps;
                    for (int j=0; j < jacdim ; j++) {
                        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
                       // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                        extremity_loc_eps_min = extremity_eps_min(i1,i2,m1,m2,j);
                        z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
                        cout << "j=" << j << " z_temp=" << z_temp << " JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    }
                    // hull on all subdivisions
                    z_outer = interval_hull(z_outer,z_temp);
                    cout << "z_outer=" << z_outer << endl;
                }
            }
        }
    }
    return z_outer;
}
 */

// subdivise n times in the 2 inputs
interval evaluate_outerrange_x_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, loc_eps, extremity_loc_eps_min, extremity_loc_eps_max, z_outer = z0[i];
    vector<interval> c_eps;
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
            //   z_temp_partial = z0[i];
         //   if (i == 0)
          //      cout << "i1=" << i1 << " i2=" << i2 << endl;
            z_temp = z0[i]; // z_temp_partial;
            int m1 = 0; int m2 = 0;
            
            c_eps = constraint_eps(i1,i2,m1,m2);
            for (int j=0; j < jacdim ; j++) {
                loc_eps = radx[j].mid() + radx[j].rad() * c_eps[j];
                // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                extremity_loc_eps_min = extremity_eps_min(radx,i1,i2,m1,m2,j);
                z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
           /*     if (i == 0) {
                    cout << "j=" << j << " m1=" << m1 << " m2=" << m2 <<endl;
                    cout << "loc_eps=" << loc_eps << " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    cout << "JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " z_temp=" << z_temp << endl;
                }*/
            }
            // hull on all subdivisions
            //    z_outer = interval_hull(z_outer,z_temp);
            
            // first translate in x
            for (int m1=1; m1 < nb_discr1 ; m1++) {
                int j = 0;
                c_eps = constraint_eps(i1,i2,m1,m2);
                //for (int j=0; j < jacdim ; j++) {
                loc_eps = radx[j].mid() + radx[j].rad() * c_eps[j];
                // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                extremity_loc_eps_min = extremity_eps_min(radx,i1,i2,m1,m2,j);
                z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
         /*       if (i == 0) {
                    cout << "j=" << j << " m1=" << m1 << endl;
                    cout << "loc_eps=" << loc_eps<< " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    cout << "JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " z_temp=" << z_temp << endl;
                } */
            }
            
            // then translate in y but for whole [x]: Jacobian should be evaluated on all m1=0..nb_discr
            for (int m2=1; m2 < nb_discr2 ; m2++) {
                int j = 1;
                c_eps = constraint_eps(i1,i2,m1,m2);
                c_eps[0] = constr_eps[i1][i2][0];
                //for (int j=0; j < jacdim ; j++) {
                loc_eps = radx[j].mid() + radx[j].rad() * c_eps[j];
                // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                extremity_loc_eps_min = extremity_eps_min(radx,i1,i2,m1,m2,j);
                z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
           /*     if (i == 0) {
                    cout << "j=" << j << " m2=" << m2 << endl;
                    cout << "loc_eps=" << loc_eps<< " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    cout << "JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " z_temp=" << z_temp << endl;
                } */
            }
            
            // hull on all 4 subdivisions
            z_outer = interval_hull(z_outer,z_temp);
        //    cout << "z_outer=" << z_outer << endl;
        }
        
    }
    return z_outer;
}


// specializing JacAff given that we are between [x]^m and [x]^{m+1} in the subdivision
void constraint_eps(vector<vector<interval>> &Jac_m, vector<vector<AAF>> &JacAff, int m)
{
    vector<interval> c_eps(jacdim);
    
    for (int j=0; j < jacdim ; j++)
        c_eps[j] = interval(-1,1);
    
    if (m == 0)
    {
        c_eps[0] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        if (jacdim > 1)
            c_eps[1] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = JacAff[i][j].convert_int(c_eps,jacdim);
    }
    else
    {
        
        c_eps[0] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        if (jacdim > 1)
            c_eps[1] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = JacAff[i][j].convert_int(c_eps,jacdim);
        
        c_eps[0] = -interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
        
        if (jacdim > 1)
        {
            c_eps[0] = interval(-((double)m)/nb_discr,((double)m)/nb_discr);
            c_eps[1] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
            
            c_eps[1] = -interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
        }
    }
}

// specializing JacAff given that we are ON THE BORDER between [x]^m and [x]^{m+1} in the subdivision
void constraint_eps_border(vector<vector<interval>> &Jac_m, vector<vector<AAF>> &JacAff, int m)
{
    vector<interval> c_eps(jacdim);
    
    for (int j=0; j < jacdim ; j++)
        c_eps[j] = interval(-1,1);
    
    c_eps[0] = ((double)m)/nb_discr;
    c_eps[1] = interval(-((double)m)/nb_discr,((double)m)/nb_discr);
    for (int i=0 ; i<sysdim ; i++)
        for (int j=0 ; j<jacdim ; j++) {
            Jac_m[i][j] = JacAff[i][j].convert_int(c_eps,jacdim);
   //         cout << JacAff[i][j] << Jac_m[i][j];
        }
 //   cout << "Jac_m=" << Jac_m << endl;
    
    c_eps[0] = -((double)m)/nb_discr;
  //  c_eps[1] = interval(-((double)m)/nb_discr,((double)m)/nb_discr);
    for (int i=0 ; i<sysdim ; i++)
        for (int j=0 ; j<jacdim ; j++)
            Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
 //   cout << "Jac_m=" << Jac_m << endl;
    
    if (jacdim > 1)
    {
        c_eps[0] = interval(-((double)m)/nb_discr,((double)m)/nb_discr);
        c_eps[1] = ((double)m)/nb_discr;
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
 //       cout << "Jac_m=" << Jac_m << endl;
        
        c_eps[1] = -((double)m)/nb_discr;
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
 //       cout << "Jac_m=" << Jac_m << endl;
    }
    
}





// specializing HessAff given that we are between [x]^m and [x]^{m+1} in the subdivision
void constraint_eps(vector<vector<vector<interval>>> &Hessf, vector<vector<vector<AAF>>> &HessAff, int m)
{
    vector<interval> c_eps(jacdim);
    
    for (int j=0; j < jacdim ; j++)
        c_eps[j] = interval(-1,1);
    
    if (m == 0)
    {
        c_eps[0] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        if (jacdim > 1)
            c_eps[1] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                for (int k=0 ; k<jacdim ; k++)
                    Hessf[i][j][k] = HessAff[i][j][k].convert_int(c_eps,jacdim);
    }
    else
    {
        
        c_eps[0] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        if (jacdim > 1)
            c_eps[1] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                for (int k=0 ; k<jacdim ; k++)
                    Hessf[i][j][k] = HessAff[i][j][k].convert_int(c_eps,jacdim);
        
        c_eps[0] = -interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                for (int k=0 ; k<jacdim ; k++)
                    Hessf[i][j][k] = interval_hull(Hessf[i][j][k],HessAff[i][j][k].convert_int(c_eps,jacdim));
        
        if (jacdim > 1)
        {
            c_eps[0] = interval(-((double)m)/nb_discr,((double)m)/nb_discr);
            c_eps[1] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    for (int k=0 ; k<jacdim ; k++)
                        Hessf[i][j][k] = interval_hull(Hessf[i][j][k],HessAff[i][j][k].convert_int(c_eps,jacdim));
            
            c_eps[1] = -interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    for (int k=0 ; k<jacdim ; k++)
                        Hessf[i][j][k] = interval_hull(Hessf[i][j][k],HessAff[i][j][k].convert_int(c_eps,jacdim));
        }
    }
}



// subdivise n times in the 2 inputs with simultaneous progress in the 2 dimensions around the center
vector<interval> evaluate_outerrange_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff)
{
    vector<interval> c_eps;
    vector<interval> z_outer(sysdim);
    vector<interval> loc_eps(jacdim);
    
    vector<vector<interval>> Jac_m(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jac_m[i] = vector<interval>(jacdim);

    for (int j=0 ; j<jacdim; j++)
        loc_eps[j] = radx[j]/((double)nb_discr);
    
    for (int i=0 ; i<sysdim ; i++)
        z_outer[i] = z0[i];
    
    
    for (int m=0; m<nb_discr; m++)
    {
        // compute bounds for JacAff locally on the m-th subdivision and put them in Jac_m
        // on pourrait affiner en includant loc_eps dedans et en faisant le produit sur chaque sous-carre
        constraint_eps(Jac_m,JacAff,m);
        
        for (int i=0 ; i<sysdim ; i++)
        {
            for (int j=0; j < jacdim ; j++)
                z_outer[i] += Jac_m[i][j] * loc_eps[j];
        }
    }
    
    return z_outer;
}


// subdivise n times in the 2 inputs with simultaneous progress in the 2 dimensions around the center
vector<interval> evaluate_outerrange_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<interval>> &dfdx0, vector<vector<vector<AAF>>> &HessAff)
{
    vector<interval> c_eps;
    vector<interval> z_outer(sysdim), z_temp(sysdim);
    vector<interval> loc_eps(jacdim);
    
    vector<vector<interval>> Jac_m(sysdim,vector<interval>(jacdim));
    vector<vector<vector<interval>>> Hessf(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    
    for (int j=0 ; j<jacdim; j++)
        loc_eps[j] = radx[j]/((double)nb_discr);
    
    // first subdivision: z0 + <Jacf(x0) . dx>
    multMiVi(z_outer,dfdx0,loc_eps);
  //  cout << "z_outer=" << z_outer;
    addViVi(z_outer,z0);
  //  cout << "z_outer=" << z_outer;
    
    // sum of gradients on the bordier between the xi and xi+1
    for (int m=1; m<nb_discr; m++)
    {
        // compute bounds for JacAff on the bordier between the xi and xi+1 locally on the m-th subdivision and put them in Jac_m
        // A AFFINER !
        constraint_eps_border(Jac_m,JacAff,m);
  //      cout << "Jac_m=" << Jac_m;
        
        for (int i=0 ; i<sysdim ; i++)
        {
            for (int j=0; j < jacdim ; j++)
                z_outer[i] += Jac_m[i][j] * loc_eps[j];
        }
    }
//    cout << "z_outer=" << z_outer;
    
    // < x , Hessf([x]) . x >
    for (int m=0; m<nb_discr; m++)
    {
        constraint_eps(Hessf,HessAff,m);
        for (int i=0 ; i< sysdim ; i++)
        {
            //   multMiVi(z_temp,Hessf[i],eps);
            //   scalarproduct(z_temp,z_temp,eps);
            // replaced by below to exploit that eps^2 can only be positive
            z_temp[i] = 0;
            for (int j=0 ; j< jacdim ; j++) {
                for (int k=0 ; k< j ; k++)
                    z_temp[i] += Hessf[i][j][k]*loc_eps[j]*loc_eps[k];
                z_temp[i] += Hessf[i][j][j]*interval(0,sup(loc_eps[j])*sup(loc_eps[j])) / 2.0;
                for (int k=j+1 ; k< jacdim ; k++)
                    z_temp[i] += Hessf[i][j][k]*loc_eps[j]*loc_eps[k];
            }
        }
        addViVi(z_outer,z_temp);
    }
    
 //   cout << "z_outer=" << z_outer;
    return z_outer;
}



vector<interval> evaluate_innerrange_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified)
{
    vector<interval> c_eps;
    vector<interval> z_inner(sysdim);
    vector<interval> z_impro(sysdim);
    vector<interval> loc_eps(jacdim);
    
    vector<vector<interval>> Jac_m(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jac_m[i] = vector<interval>(jacdim);
    
    for (int j=0 ; j<jacdim; j++)
        loc_eps[j] = radx[j]/((double)nb_discr);
    
    for (int i=0 ; i<sysdim ; i++) {
        z_inner[i] = z0[i];
        z_impro[i] = 0;
    }
    
    
    for (int m=0; m<nb_discr; m++)
    {
        // compute bounds for JacAff locally on the m-th subdivision and put them in Jac_m
        constraint_eps(Jac_m,JacAff,m);
        
        for (int i=0 ; i<sysdim ; i++)
        {
            for (int j=0; j < jacdim ; j++) {
                if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
                    z_impro[i] +=  Kaucher_multeps(Jac_m[i][j],loc_eps[j]);
                else
                    z_inner[i] += Jac_m[i][j] * loc_eps[j];
            }
        }
    }
    for (int i=0 ; i<sysdim ; i++)
        z_inner[i] = Kaucher_add_pro_impro(z_inner[i],z_impro[i]);
    
    return z_inner;
    
}




vector<interval> evaluate_innerrange_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<interval>> &dfdx0, vector<vector<vector<AAF>>> &HessAff, bool maximal, vector<int> &exist_quantified)
{
    vector<interval> c_eps;
    vector<interval> z_inner(sysdim);
    vector<interval> z_impro(sysdim);
    vector<interval> z_pro(sysdim);
    vector<interval> loc_eps(jacdim);
    
    vector<vector<interval>> Jac_m(sysdim,vector<interval>(jacdim));
    vector<vector<vector<interval>>> Hessf(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));

    
    for (int j=0 ; j<jacdim; j++)
        loc_eps[j] = radx[j]/((double)nb_discr);
    
    // first subdivision order 1 : z0 + Jacf(x0) . x
    for (int i=0 ; i<sysdim ; i++)
    {
        z_impro[i] = 0;
        z_pro[i] = z0[i]; // + z_temp[i]; // 2nd order term is adversariam
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                z_impro[i] += Kaucher_multeps(dfdx0[i][j],loc_eps[j]);
            else
                z_pro[i] += dfdx0[i][j]*loc_eps[j];
        }
        
    }
    
 //   z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    
    
     // sum of gradients on the border between the xi and xi+1
    for (int m=1; m<nb_discr; m++)
    {
        // compute bounds for JacAff locally on the m-th subdivision and put them in Jac_m
        constraint_eps_border(Jac_m,JacAff,m);
        
        for (int i=0 ; i<sysdim ; i++)
        {
            for (int j=0; j < jacdim ; j++) {
                if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
                    z_impro[i] +=  Kaucher_multeps(Jac_m[i][j],loc_eps[j]);
                else
                    z_pro[i] += Jac_m[i][j] * loc_eps[j];
            }
        }
    }
    
    // order 2  < x , Hessf([x]) . x > is adversarial here
    for (int m=0; m<nb_discr; m++)
    {
        constraint_eps(Hessf,HessAff,m);
        for (int i=0 ; i< sysdim ; i++)
        {
            //   multMiVi(z_temp,Hessf[i],eps);
            //   scalarproduct(z_temp,z_temp,eps);
            // replaced by below to exploit that eps^2 can only be positive
            //    z_temp[i] = 0;
            for (int j=0 ; j< jacdim ; j++) {
                for (int k=0 ; k< j ; k++)
                    z_pro[i] += Hessf[i][j][k]*loc_eps[j]*loc_eps[k];
                z_pro[i] += Hessf[i][j][j]*interval(0,sup(loc_eps[j])*sup(loc_eps[j])) / 2.0;
                for (int k=j+1 ; k< jacdim ; k++)
                    z_pro[i] += Hessf[i][j][k]*loc_eps[j]*loc_eps[k];
            }
        }
    }
    
    for (int i=0 ; i<sysdim ; i++)
        z_inner[i] = Kaucher_add_pro_impro(z_pro[i],z_impro[i]);
    
    return z_inner;
    
    
    
    
    
    
    /*
    // order 2  < x , Hessf([x]) . x > is adversarial here
    for (int i=0 ; i< sysdim ; i++)
    {
        //   multMiVi(z_temp,Hessf[i],eps);
        //   scalarproduct(z_temp,z_temp,eps);
        // replaced by below to exploit that eps^2 can only be positive
        z_temp[i] = 0;
        for (int j=0 ; j< jacdim ; j++) {
            for (int k=0 ; k< j ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
            z_temp[i] += Hessf[i][j][j]*interval(0,sup(radx[j])*sup(radx[j])) / 2.0;
            for (int k=j+1 ; k< jacdim ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
        }
    }
    
    
    // order 1 : z0 + Jacf(x0) . x
    for (int i=0 ; i<sysdim ; i++)
    {
        inner_impro = 0;
        inner_pro = z0[i] + z_temp[i]; // 2nd order term is adversariam
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],radx[j]);
            else
                inner_pro += Jacf[i][j]*radx[j];
        }
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    }
    */
    
    
}



vector<interval> evaluate_innerrange(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, bool maximal, vector<int> &exist_quantified) {
    interval inner_impro, inner_pro;
    vector<interval> z_inner(sysdim);
    
    for (int i=0 ; i<sysdim ; i++)
    {
        inner_impro = 0;
        inner_pro = z0[i];
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],radx[j]);
            else
                inner_pro += Jacf[i][j]*radx[j];
        }
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    }
    
    return z_inner;
}

// same as evaluate_innerrange but robust to parameters
vector<interval> evaluate_innerrange_robust(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, bool maximal, vector<int> &exist_quantified) {
    interval inner_impro, inner_pro;
    vector<interval> z_inner(sysdim);
    
    for (int i=0 ; i<sysdim ; i++)
    {
        inner_impro = 0;
        inner_pro = z0[i];
        for (int j=0; j < sysdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],radx[j]);
            else
                inner_pro += Jacf[i][j]*radx[j];
        }
        for (int j=sysdim; j < jacdim ; j++)
            inner_pro += Jacf[i][j]*radx[j];
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    }
    
    return z_inner;
}


// same that evaluate_innerrange but for preconditioner C
vector<interval> evaluate_precond_innerrange(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<double>> C, int varx, int vary, bool maximal, vector<int> &exist_quantified) {
    interval inner_impro, inner_pro;
    vector<interval> z_inner(sysdim);
    
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    vector<vector<interval>> CJacf(sysdim, vector<interval>(jacdim));
    multMiMi(CJacf,C,Jacf);
    cout << "CJacf_i=" << endl;
    cout << CJacf;
    
    for (int i=0 ; i<sysdim ; i++)
    {
        inner_impro = 0;
        inner_pro = f0[i];
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(CJacf[i][j],radx[j]);
            else
                inner_pro += CJacf[i][j]*radx[j];
        }
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    }
    
    return z_inner;
}






// TODO A completer
vector<interval> evaluate_innerrange_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf, bool maximal, vector<int> &exist_quantified)
{
    interval inner_impro, inner_pro;
    vector<interval> z_inner(sysdim), z_temp(sysdim);
    
    // order 2  < x , Hessf([x]) . x > is adversarial here
    for (int i=0 ; i< sysdim ; i++)
    {
        //   multMiVi(z_temp,Hessf[i],eps);
        //   scalarproduct(z_temp,z_temp,eps);
        // replaced by below to exploit that eps^2 can only be positive
        z_temp[i] = 0;
        for (int j=0 ; j< jacdim ; j++) {
            for (int k=0 ; k< j ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
            z_temp[i] += Hessf[i][j][j]*interval(0,sup(radx[j])*sup(radx[j])) / 2.0;
            for (int k=j+1 ; k< jacdim ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
        }
    }
    
    
    // order 1 : z0 + Jacf(x0) . x
    for (int i=0 ; i<sysdim ; i++)
    {
        inner_impro = 0;
        inner_pro = z0[i] + z_temp[i]; // 2nd order term is adversariam
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],radx[j]);
            else
                inner_pro += Jacf[i][j]*radx[j];
        }
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    }
    
  //  addViVi(z_outer,z_temp);
    
    return z_inner;
}



// TODO A completer
vector<interval> evaluate_innerrange_order2_robust(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf, bool maximal, vector<int> &exist_quantified)
{
    interval inner_impro, inner_pro;
    vector<interval> z_inner(sysdim), z_temp(sysdim);
    
    // order 2  < x , Hessf([x]) . x > is adversarial here
    for (int i=0 ; i< sysdim ; i++)
    {
        //   multMiVi(z_temp,Hessf[i],eps);
        //   scalarproduct(z_temp,z_temp,eps);
        // replaced by below to exploit that eps^2 can only be positive
        z_temp[i] = 0;
        for (int j=0 ; j< jacdim ; j++) {
            for (int k=0 ; k< j ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
            z_temp[i] += Hessf[i][j][j]*interval(0,sup(radx[j])*sup(radx[j])) / 2.0;
            for (int k=j+1 ; k< jacdim ; k++)
                z_temp[i] += Hessf[i][j][k]*radx[j]*radx[k];
        }
    }
    
    
    // order 1 : z0 + Jacf(x0) . x
    for (int i=0 ; i<sysdim ; i++)
    {
        inner_impro = 0;
        inner_pro = z0[i] + z_temp[i]; // 2nd order term is adversariam
        for (int j=0; j < sysdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],radx[j]);
            else
                inner_pro += Jacf[i][j]*radx[j];
        }
        for (int j=sysdim; j < jacdim ; j++) {
            inner_pro += Jacf[i][j]*radx[j];
        }
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
    }
    
    //  addViVi(z_outer,z_temp);
    
    return z_inner;
}





// evaluation on one subdivision (defined by index1 and index2)
// defined per subdivision contrarily to over-approx because we cannot do the hull for inner approx
interval evaluate_innerrange_x_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2) {
    interval inner_impro, inner_pro;
    interval z_inner;
    
    inner_impro = 0;
    inner_pro = z0[i];
    for (int j=0; j < jacdim ; j++) {
        if (maximal || (exist_quantified[j] == i)) {// output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(constr_eps[index1][index2],jacdim),eps_loc[index1][index2][j]);
        }
        else {
            inner_pro += JacAff[i][j].convert_int(constr_eps[index1][index2],jacdim) * eps_loc[index1][index2][j];
        }
    }
    //   cout << endl;
    z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}

// evaluation on one subdivision (defined by index1 and index2)
// defined per subdivision contrarily to over-approx because we cannot do the hull for inner approx
/* interval evaluate_innerrange_x_subdiv_discretize_old(vector<interval> &z0,  vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2) {
    interval inner_impro, inner_pro;
    interval z_inner;
    
    inner_impro = 0;
    inner_pro = z0[i];
    for (int m=0 ; m<nb_discr; m++)
    {
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) {// output component i the one where j is existentially quantified {
                if (m>0)
                    inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim),(eps_loc_discr[index1][index2][m][j]-extremity_eps_loc_discr[index1][index2][m-1][j]));
                else
                    inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim),eps_loc_discr[index1][index2][m][j]);
            }
            else {
                if (m>0)
                    inner_pro += JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim) * (eps_loc_discr[index1][index2][m][j]-extremity_eps_loc_discr[index1][index2][m-1][j]);
                else
                    inner_pro += JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim) * eps_loc_discr[index1][index2][m][j];
            }
        }
    }
    //   cout << endl;
    z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}
*/
    
/*    z_temp_partial = z0[i];
    // evaluated at the extrimity on the n-1 first subdivisions
    for (int m=0; m < nb_discr-1 ; m++) {
        z_temp = z_temp_partial;
        for (int j=0; j < jacdim ; j++) {
            if (m>0)
            {
                z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
                z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (extremity_eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
            }
            else
            {
                z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * eps_loc_discr[i1][i2][m][j];
                z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * extremity_eps_loc_discr[i1][i2][m][j];
            }
        }
        // hull on all subdivisions
        z_outer = interval_hull(z_outer,z_temp);
        // shoudl be doing hull here I guess
    }
    // evaluated on interval on the last subdivision
    z_temp = z_temp_partial;
    for (int j=0; j < jacdim ; j++) {
        z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][nb_discr-1],jacdim) * (eps_loc_discr[i1][i2][nb_discr-1][j]-extremity_eps_loc_discr[i1][i2][nb_discr-2][j]);
    }
    
    z_outer = interval_hull(z_outer,z_temp);*/




// evaluation on one subdivision (defined by index1 and index2)
// defined per subdivision contrarily to over-approx because we cannot do the hull for inner approx
interval evaluate_innerrange_x_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2) {
    interval inner_impro, inner_pro, loc_eps;
    interval z_inner;
    vector<interval> c_eps(jacdim);
    double extremity_loc_eps_min;
    
    inner_impro = 0;
    inner_pro = z0[i];
    
    // mean-value on first block
    int m1 = 0, m2 = 0;
    c_eps = constraint_eps(index1,index2,m1,m2);
    for (int j=0; j < jacdim ; j++) {
        loc_eps = radx[j].mid() + radx[j].rad() * c_eps[j];
        extremity_loc_eps_min = extremity_eps_min(radx,index1,index2,m1,m2,j); // always 0 ?
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(c_eps,jacdim),loc_eps - extremity_loc_eps_min);
        else
            inner_pro += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    // first translate in first component x
    for (int m1=1; m1 < nb_discr1 ; m1++) {
        int j = 0;
        c_eps = constraint_eps(index1,index2,m1,m2);
        //for (int j=0; j < jacdim ; j++) {
        loc_eps = radx[j].mid() + radx[j].rad() * c_eps[j];
        // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
        extremity_loc_eps_min = extremity_eps_min(radx,index1,index2,m1,m2,j);
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(c_eps,jacdim),loc_eps - extremity_loc_eps_min);
        else
            inner_pro += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    // then translate in y but for whole [x]: Jacobian should be evaluated on all m1=0..nb_discr
    for (int m2=1; m2 < nb_discr2 ; m2++) {
        int j = 1;
        c_eps = constraint_eps(index1,index2,m1,m2);
        c_eps[0] = constr_eps[index1][index2][0];
        //for (int j=0; j < jacdim ; j++) {
        loc_eps = radx[j].mid() + radx[j].rad() * c_eps[j];
        // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
        extremity_loc_eps_min = extremity_eps_min(radx,index1,index2,m1,m2,j);
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(c_eps,jacdim),loc_eps - extremity_loc_eps_min);
        else
            inner_pro += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    //   cout << endl;
    z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}


    
    

void print_pi(vector<int> &exist_quantified) {
    cout << "(pi : ";
    for (int j = 0 ; j < jacdim ; j++)
        cout << j+1 << " -> " << exist_quantified[j]+1 << ", " ;
    cout << ") :  ";
}


// z0 = f(x0)
void joint_ranges(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<interval>> &Jacf0, vector<vector<vector<interval>>> &Hessf, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "inner_joint.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    // order 1
    vector<interval> inner = evaluate_innerrange(z0,radx,Jacf,false,exist_quantified);
    cout << "inner box mean-value ";
    print_pi(exist_quantified);
    cout << inner[varx] << "  " << inner[vary] << endl;
    outFile_jointinner << inf(inner[varx]) << "\t" << sup(inner[varx]) << "\t" << inf(inner[vary]) << "\t" << sup(inner[vary]) <<  endl;
    
    // order 2
    inner = evaluate_innerrange_order2(z0,radx,Jacf0,Hessf,false,exist_quantified);
    cout << "inner box order 2 ";
    print_pi(exist_quantified);
    cout << inner[varx] << "  " << inner[vary] << endl;
 //   outFile_jointinner << inf(inner[varx]) << "\t" << sup(inner[varx]) << "\t" << inf(inner[vary]) << "\t" << sup(inner[vary]) <<  endl;
    
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    // order 1
    inner = evaluate_innerrange(z0,radx,Jacf,false,exist_quantified);
    cout << "inner box mean-value ";
    print_pi(exist_quantified);
    cout << inner[varx] << "  " << inner[vary] << endl;
    outFile_jointinner << inf(inner[varx]) << "\t" << sup(inner[varx]) << "\t" << inf(inner[vary]) << "\t" << sup(inner[vary]) <<  endl;
    
    // order 2
    inner = evaluate_innerrange_order2(z0,radx,Jacf0,Hessf,false,exist_quantified);
    cout << "inner box order 2 ";
    print_pi(exist_quantified);
    cout << inner[varx] << "  " << inner[vary] << endl;
 //   outFile_jointinner << inf(inner[varx]) << "\t" << sup(inner[varx]) << "\t" << inf(inner[vary]) << "\t" << sup(inner[vary]) <<  endl;
    
    outFile_jointinner.close();
}




// z0 = f(x0)

void joint_ranges_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    interval temp_x, tempy, inner_x, inner_y;
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "inner_joint.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    
    cout << "inner box with subdiv"; print_pi(exist_quantified); cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv(z0,radx,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv(z0,radx,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    // TODO : en deriver une boite gloabel ?
   //  inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,varx);
   //  inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,vary);
    
   // cout << inner_x << "  " << inner_y << endl;
   // outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    cout << "inner box with subdiv ";
    print_pi(exist_quantified);
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv(z0,radx,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv(z0,radx,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    outFile_jointinner.close();
}


void joint_ranges_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    interval temp_x, tempy, inner_x, inner_y;
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "inner_joint.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    
    cout << "inner box "; print_pi(exist_quantified); cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv_discretize(z0,radx,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv_discretize(z0,radx,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    // TODO : en deriver une boite gloabel ?
    //  inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,varx);
    //  inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,vary);
    
    // cout << inner_x << "  " << inner_y << endl;
    // outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    cout << "inner box ";
    print_pi(exist_quantified);
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv_discretize(z0,radx,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv_discretize(z0,radx,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    outFile_jointinner.close();
}


void joint_ranges_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    vector<interval> z_inner;
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "inner_joint.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    
    cout << "inner box (mean value, 1d discretization): "; print_pi(exist_quantified);
    z_inner = evaluate_innerrange_discretize_simultaneous(z0,radx,JacAff,false,exist_quantified);
    cout << z_inner[0] << "  " << z_inner[1] << "       ";
    outFile_jointinner << inf(z_inner[0]) << "\t" << sup(z_inner[0]) << "\t" << inf(z_inner[1]) << "\t" << sup(z_inner[1]) <<  endl;
    cout << endl;
    // TODO : en deriver une boite gloabel ?
    //  inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,varx);
    //  inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,vary);
    
    // cout << inner_x << "  " << inner_y << endl;
    // outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    cout << "inner box (mean value, 1d discretization): "; print_pi(exist_quantified);
    z_inner = evaluate_innerrange_discretize_simultaneous(z0,radx,JacAff,false,exist_quantified);
    cout << z_inner[0] << "  " << z_inner[1] << "       ";
    outFile_jointinner << inf(z_inner[0]) << "\t" << sup(z_inner[0]) << "\t" << inf(z_inner[1]) << "\t" << sup(z_inner[1]) <<  endl;
    cout << endl;
    outFile_jointinner.close();
}



// f0 = f(x0),
// tests by Hansen Sengupta if a candidate box (fix,fixy) is in the joint range for f_varx, f_fary
/*void Hansen_Sengupta(vector<interval> &f0, vector<vector<interval>> &Jacf, interval &fix, interval &fi2, int varx, int vary)
{
    interval xtilde = x0 + gamma(Jacf)
}*/







void preconditioned_joint_ranges(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf,  vector<vector<interval>> &Jacf0, vector<vector<vector<interval>>> &Hessf, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    
    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim,vector<double> (sysdim)), C(sysdim,vector<double> (sysdim));
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(Jacf[varx][varx]);
    A[vary][vary] = mid(Jacf[vary][vary]);
    A[varx][vary] = mid(Jacf[varx][vary]);
    A[vary][varx] = mid(Jacf[vary][varx]);
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // CJacf = C * Jacf
    vector<vector<interval>> CJacf(sysdim, vector<interval>(jacdim));
   
    multMiMi(CJacf,C,Jacf);
    cout << "CJacf=" << endl;
    cout << CJacf;
    
    // CJacf = C * Jacf ===> TODO. plutot prendre un C directement a partir de Jacf0 ?
    vector<vector<interval>> CJacf0(sysdim, vector<interval>(jacdim));
    
    multMiMi(CJacf0,C,Jacf0);
    cout << "CJacf0=" << endl;
    cout << CJacf0;
  
    // CHessf = C * Hessf
    vector<vector<vector<interval>>> CHessf(sysdim, vector<vector<interval>>(jacdim, vector<interval>(jacdim)));
    
        for (int j=0 ; j<jacdim; j++) {
            for (int k=0 ; k<jacdim; k++) {
                for (int i=0 ; i<sysdim; i++)
                    CHessf[i][j][k] = Hessf[i][j][k];
                CHessf[varx][j][k] = C[varx][varx]*Hessf[varx][j][k] + C[varx][vary]*Hessf[vary][j][k];
                CHessf[vary][j][k] = C[vary][varx]*Hessf[varx][j][k] + C[vary][vary]*Hessf[vary][j][k];
        }
    }
    
    // outer range
    vector<interval> temp_outer = evaluate_outerrange(f0,radx,CJacf);
    out_approx << YAML::Key << "outer2d";
    out_approx << YAML::Value;
    out_approx << YAML::BeginSeq;
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "outer skewed box (mean value): ";
    output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    temp_outer = evaluate_outerrange_order2(f0,radx,CJacf0,CHessf);
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 2;
    cout << "outer skewed box (order 2): ";
    output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    out_approx << YAML::EndSeq;
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    vector<interval> temp_inner;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
     temp_inner = evaluate_innerrange(f0,radx,CJacf,false,exist_quantified);
    out_approx << YAML::Key << "inner2d";
    out_approx << YAML::Value;
    out_approx << YAML::BeginSeq;
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "inner skewed box (mean value): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 2;
    cout << "inner skewed box (order 2): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    temp_inner = evaluate_innerrange(f0,radx,CJacf,false,exist_quantified);
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "inner skewed box (mean value): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 2;
    cout << "inner skewed box (order 2): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        temp_inner = evaluate_innerrange(f0,radx,CJacf,false,exist_quantified);
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 1;
        cout << "inner skewed box (mean value): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
        temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 2;
        cout << "inner skewed box (order 2): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        temp_inner = evaluate_innerrange(f0,radx,CJacf,false,exist_quantified);
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 1;
        cout << "inner skewed box (mean value): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
        temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 2;
        cout << "inner skewed box (order 2): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
    }
    out_approx << YAML::EndSeq;
}



void preconditioned_joint_ranges_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;

    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim), C(sysdim);
    for (int i=0 ; i<sysdim; i++) {
        A[i] = vector<double> (sysdim);
        C[i] = vector<double> (sysdim);
    }
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(JacAff[varx][varx].convert_int());
    A[vary][vary] = mid(JacAff[vary][vary].convert_int());
    A[varx][vary] = mid(JacAff[varx][vary].convert_int());
    A[vary][varx] = mid(JacAff[vary][varx].convert_int());
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // CJacf = C * Jacf
    vector<vector<AAF>> CJacAff(sysdim);
    for (int i=0; i < sysdim ; i++)
        CJacAff[i] = vector<AAF>(jacdim);
    
    multMiMi(CJacAff,C,JacAff);
  //  cout << "CJacf=" << endl;
  //  cout << CJacf;
    
    
    // outer range
    interval temp_outer_x = evaluate_outerrange_x_subdiv(f0,radx,CJacAff,varx);
    interval temp_outer_y = evaluate_outerrange_x_subdiv(f0,radx,CJacAff,vary);
    out_approx << YAML::Key << "outer2d";
    out_approx << YAML::Value;
    out_approx << YAML::BeginSeq;
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    out_approx << YAML::EndSeq;
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    interval temp_inner_x, temp_inner_y;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            temp_inner_x = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
      //      cout << "f0= " << f0;
      //      cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            out_approx << YAML::Key << "inner2d";
            out_approx << YAML::Value;
            out_approx << YAML::BeginSeq;
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << varx;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << vary;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
            out_approx << YAML::EndMap;
        }
    }
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
    //        cout << "f0= " << f0 ;
            temp_inner_x = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
   //         cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << varx;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << vary;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
            out_approx << YAML::EndMap;
        }
    }
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
                out_approx << YAML::BeginMap;
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << varx;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << vary;
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
                out_approx << YAML::EndMap;
            }
        }
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
                out_approx << YAML::BeginMap;
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << varx;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << vary;
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
                out_approx << YAML::EndMap;
            }
        }
    }
    out_approx << YAML::EndSeq;
}



void preconditioned_joint_ranges_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim), C(sysdim);
    for (int i=0 ; i<sysdim; i++) {
        A[i] = vector<double> (sysdim);
        C[i] = vector<double> (sysdim);
    }
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(JacAff[varx][varx].convert_int());
    A[vary][vary] = mid(JacAff[vary][vary].convert_int());
    A[varx][vary] = mid(JacAff[varx][vary].convert_int());
    A[vary][varx] = mid(JacAff[vary][varx].convert_int());
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // CJacf = C * Jacf
    vector<vector<AAF>> CJacAff(sysdim);
    for (int i=0; i < sysdim ; i++)
        CJacAff[i] = vector<AAF>(jacdim);
    
    multMiMi(CJacAff,C,JacAff);
    //  cout << "CJacf=" << endl;
    //  cout << CJacf;
    
    
    // outer range
    interval temp_outer_x = evaluate_outerrange_x_subdiv_discretize(f0,radx,CJacAff,varx);
    interval temp_outer_y = evaluate_outerrange_x_subdiv_discretize(f0,radx,CJacAff,vary);
    out_approx << YAML::Key << "outer2d";
    out_approx << YAML::Value;
    out_approx << YAML::BeginSeq;
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    out_approx << YAML::EndSeq;
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    interval temp_inner_x, temp_inner_y;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
            //      cout << "f0= " << f0;
            //      cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            out_approx << YAML::Key << "inner2d";
            out_approx << YAML::Value;
            out_approx << YAML::BeginSeq;
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << varx;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << vary;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
            out_approx << YAML::EndMap;
        }
    }
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            //        cout << "f0= " << f0 ;
            temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
            //         cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            out_approx << YAML::BeginMap;
            out_approx << YAML::Key << "x1";
            out_approx << YAML::Value << varx;
            out_approx << YAML::Key << "x2";
            out_approx << YAML::Value << vary;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
            out_approx << YAML::EndMap;
        }
    }
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
                out_approx << YAML::BeginMap;
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << varx;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << vary;
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
                out_approx << YAML::EndMap;
            }
        }
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,radx,CJacAff,false,exist_quantified,vary,i1,i2);
                out_approx << YAML::BeginMap;
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << varx;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << vary;
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1);
                out_approx << YAML::EndMap;
            }
        }
    }
    out_approx << YAML::EndSeq;
    
}


void preconditioned_joint_ranges_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<interval>> &dfdx0, vector<vector<vector<AAF>>> &HessAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    vector<interval> temp_outer(sysdim);
    vector<interval> temp_inner(sysdim);
    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim), C(sysdim);
    for (int i=0 ; i<sysdim; i++) {
        A[i] = vector<double> (sysdim);
        C[i] = vector<double> (sysdim);
    }
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(JacAff[varx][varx].convert_int());
    A[vary][vary] = mid(JacAff[vary][vary].convert_int());
    A[varx][vary] = mid(JacAff[varx][vary].convert_int());
    A[vary][varx] = mid(JacAff[vary][varx].convert_int());
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // Cdfdx0 = C * dfdx0
    vector<vector<interval>> Cdfdx0(sysdim,vector<interval>(jacdim));
    multMiMi(Cdfdx0,C,dfdx0);
    
    // CJacf = C * Jacf
    vector<vector<AAF>> CJacAff(sysdim,vector<AAF>(jacdim));
    
    multMiMi(CJacAff,C,JacAff);
    //  cout << "CJacf=" << endl;
    //  cout << CJacf;
    
    // CHessf = C * Hessf
    vector<vector<vector<AAF>>> CHessAff(sysdim, vector<vector<AAF>>(jacdim, vector<AAF>(jacdim)));
    for (int j=0 ; j<jacdim; j++) {
        for (int k=0 ; k<jacdim; k++) {
            for (int i=0 ; i<sysdim; i++)
                CHessAff[i][j][k] = HessAff[i][j][k];
            CHessAff[varx][j][k] = C[varx][varx]*HessAff[varx][j][k] + C[varx][vary]*HessAff[vary][j][k];
            CHessAff[vary][j][k] = C[vary][varx]*HessAff[varx][j][k] + C[vary][vary]*HessAff[vary][j][k];
        }
    }
    
    // outer range
    temp_outer = evaluate_outerrange_discretize_simultaneous(f0,radx,CJacAff);
    out_approx << YAML::Key << "outer2d";
    out_approx << YAML::Value;
    out_approx << YAML::BeginSeq;
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "outer skewed box (mean-value, 1d discretization): ";
    output_skewedbox = print_skewbox(temp_outer[0], temp_outer[1], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    // outer range
    temp_outer = evaluate_outerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff);
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 2;
    cout << "outer skewed box (order 2, 1d discretization): ";
    output_skewedbox = print_skewbox(temp_outer[0], temp_outer[1], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    out_approx << YAML::EndSeq;
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
   // interval temp_inner_x, temp_inner_y;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    
    out_approx << YAML::Key << "inner2d";
    out_approx << YAML::Value;
    out_approx << YAML::BeginSeq;
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "inner skewed box (mean-value, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 2;
    cout << "inner skewed box (order 2, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 1;
    cout << "inner skewed box (mean-value, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    out_approx << YAML::BeginMap;
    out_approx << YAML::Key << "x1";
    out_approx << YAML::Value << varx;
    out_approx << YAML::Key << "x2";
    out_approx << YAML::Value << vary;
    out_approx << YAML::Key << "order";
    out_approx << YAML::Value << 2;
    cout << "inner skewed box (order 2, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
    out_approx << YAML::EndMap;
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 1;
        cout << "inner skewed box (mean-value, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 2;
        cout << "inner skewed box (order 2, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 1;
        cout << "inner skewed box (mean-value, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
        
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "x1";
        out_approx << YAML::Value << varx;
        out_approx << YAML::Key << "x2";
        out_approx << YAML::Value << vary;
        out_approx << YAML::Key << "order";
        out_approx << YAML::Value << 2;
        cout << "inner skewed box (order 2, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1);
        out_approx << YAML::EndMap;
    }
    
    out_approx << YAML::EndSeq;
}



void twodim_discretization_by_quadrant(vector<interval> &radx) {
    
    // begin for subdivisions
    constr_eps = vector<vector<vector<interval>>>(jacdim);
    eps_loc = vector<vector<vector<interval>>>(jacdim);
    for (int i=0; i < jacdim ; i++) {
        constr_eps[i] = vector<vector<interval>>(jacdim);
        eps_loc[i] = vector<vector<interval>>(jacdim);
        for (int j=0; j < jacdim ; j++) {
            constr_eps[i][j] = vector<interval>(jacdim);
            eps_loc[i][j] = vector<interval>(jacdim);
            for (int k=0; k < jacdim ; k++) {
                constr_eps[i][j][k] = interval(-1,1);
                eps_loc[i][j][k] = radx[k];
            }
        }
    }
    
    constr_eps[0][0][0] = interval(-1,0);
    constr_eps[0][0][1] = interval(-1,0);
    constr_eps[0][1][0] = interval(-1,0);
    constr_eps[0][1][1] = interval(0,1);
    constr_eps[1][0][0] = interval(0,1);
    constr_eps[1][0][1] = interval(-1,0);
    constr_eps[1][1][0] = interval(0,1);
    constr_eps[1][1][1] = interval(0,1);
    
    for (int i=0; i < jacdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            for (int k=0; k < jacdim ; k++)
                eps_loc[i][j][k] = radx[k].mid() + radx[k].rad() * constr_eps[i][j][k];
        }
    }
    // end subdiv
    
    // additional discretization of each component
    constr_eps_discr = vector<vector<vector<vector<interval>>>>(jacdim);
    eps_loc_discr = vector<vector<vector<vector<interval>>>>(jacdim);
    extremity_eps_loc_discr = vector<vector<vector<vector<double>>>>(jacdim);
    for (int i=0; i < jacdim ; i++) {
        constr_eps_discr[i] = vector<vector<vector<interval>>>(jacdim);
        eps_loc_discr[i] = vector<vector<vector<interval>>>(jacdim);
        extremity_eps_loc_discr[i] = vector<vector<vector<double>>>(jacdim);
        for (int j=0; j < jacdim ; j++) {
            constr_eps_discr[i][j] = vector<vector<interval>>(nb_discr);
            eps_loc_discr[i][j] = vector<vector<interval>>(nb_discr);
            extremity_eps_loc_discr[i][j] = vector<vector<double>>(nb_discr);
            for (int m=0; m < nb_discr ; m++) {
                constr_eps_discr[i][j][m] = vector<interval>(jacdim);
                eps_loc_discr[i][j][m] = vector<interval>(jacdim);
                extremity_eps_loc_discr[i][j][m] = vector<double>(jacdim);
                for (int k=0; k < jacdim ; k++) {
                    constr_eps_discr[i][j][m][k] = interval(-1,1);
                    eps_loc_discr[i][j][m][k] = radx[k];
                    // completer pour extremity_eps_loc_discr ?
                }
            }
        }
    }
    
    for (int m=0; m < nb_discr ; m++) {
        constr_eps_discr[1][1][m][0] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        constr_eps_discr[1][1][m][1] = constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][0][m][0] = -constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][0][m][1] = -constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][1][m][0] = -constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][1][m][1] = constr_eps_discr[1][1][m][0];
        constr_eps_discr[1][0][m][0] = constr_eps_discr[1][1][m][0];
        constr_eps_discr[1][0][m][1] = -constr_eps_discr[1][1][m][0];
        extremity_eps_loc_discr[1][1][m][0] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[1][1][m][1] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][0][m][0] = -(m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][0][m][1] = -(m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][1][m][0] = -(m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][1][m][1] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[1][0][m][0] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[1][0][m][1] = -(m+1.0)/nb_discr;
    }
    
    
    for (int i=0; i < jacdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            for (int m=0; m < nb_discr ; m++) {
                for (int k=0; k < jacdim ; k++) {
                    eps_loc_discr[i][j][m][k] = radx[k].mid() + radx[k].rad() * constr_eps_discr[i][j][m][k];
                    extremity_eps_loc_discr[i][j][m][k] = radx[k].mid() + radx[k].rad() * extremity_eps_loc_discr[i][j][m][k];
                    //      cout << "constr_eps_discr[i][j][m][k]=" << constr_eps_discr[i][j][m][k] << endl;
                    //     cout << "eps_loc_discr[i][j][m][k]=" << eps_loc_discr[i][j][m][k] << endl;
                }
            }
        }
    }
}


vector<interval> estimate_range(DiscreteFunc &f, vector<interval> &xinit, int discr)
{
    
    
    // estimation of the range of f
    vector<double> input(jacdim);
    vector<double> output(sysdim);
    vector<double> max_output(sysdim);
    vector<double> min_output(sysdim);
    vector<interval> range(sysdim);
    
    cout << "Initial condition x: " << xinit;
    cout << endl;
    
  //  int discr = 10;
    for (int i=0; i < jacdim ; i++)
        input[i] = xinit[i].mid();
    for (int i=0; i < sysdim ; i++) {
        max_output = f(input);
        min_output = f(input);
    }
    
    ofstream samplesreachsetfile;
    samplesreachsetfile.open("output/samplesreachset.yaml");
    
    YAML::Emitter out_samples;
    out_samples << YAML::BeginMap;
    
    out_samples << YAML::Key << "samples";
    out_samples << YAML::Value << YAML::BeginSeq;
    

    for (int i1=0; i1 <= discr ; i1++) {
        input[0] = xinit[0].inf() + (2.0*i1*xinit[0].rad())/discr;
        if (jacdim > 1) {
            for (int i2=0; i2 <= discr ; i2++) {
                input[1] = xinit[1].inf() + (2.0*i2*xinit[1].rad())/discr;
                if (jacdim > 2) {
                    for (int i3=0; i3 <= discr ; i3++) {
                        input[2] = xinit[2].inf() + (2.0*i3*xinit[2].rad())/discr;
                        if (jacdim > 3) {
                            for (int i4=0; i4 <= discr ; i4++) {
                                input[3] = xinit[3].inf() + (2.0*i4*xinit[3].rad())/discr;
                                output = f(input);
                                out_samples << YAML::BeginMap;
                                out_samples << YAML::Key << "tn";
                                out_samples << YAML::Value << 1;
                                out_samples << YAML::Key << "sample";
                                out_samples << YAML::Value << output;
                                out_samples << YAML::EndMap;
                                for (int i=0; i < sysdim ; i++) {
                                    if (output[i] < min_output[i])
                                        min_output[i] = output[i];
                                    else if (output[i] > max_output[i])
                                        max_output[i] = output[i];
                                }
                            }
                        }
                        else
                        {
                        output = f(input);
                            out_samples << YAML::BeginMap;
                            out_samples << YAML::Key << "tn";
                            out_samples << YAML::Value << 1;
                            out_samples << YAML::Key << "sample";
                            out_samples << YAML::Value << output;
                            out_samples << YAML::EndMap;
                        for (int i=0; i < sysdim ; i++) {
                            if (output[i] < min_output[i])
                                min_output[i] = output[i];
                            else if (output[i] > max_output[i])
                                max_output[i] = output[i];
                        }
                        }
                    }
                }
                else
                {
                    output = f(input);
                    out_samples << YAML::BeginMap;
                    out_samples << YAML::Key << "tn";
                    out_samples << YAML::Value << 1;
                    out_samples << YAML::Key << "sample";
                    out_samples << YAML::Value << output;
                    out_samples << YAML::EndMap;
                    for (int i=0; i < sysdim ; i++) {
                        if (output[i] < min_output[i])
                            min_output[i] = output[i];
                        else if (output[i] > max_output[i])
                            max_output[i] = output[i];
                    }
                    
                }
            }
        }
        else {
            output = f(input);
            out_samples << YAML::BeginMap;
            out_samples << YAML::Key << "tn";
            out_samples << YAML::Value << 1;
            out_samples << YAML::Key << "sample";
            out_samples << YAML::Value << output;
            out_samples << YAML::EndMap;
            
            for (int i=0; i < sysdim ; i++) {
                if (output[i] < min_output[i])
                    min_output[i] = output[i];
                else if (output[i] > max_output[i])
                    max_output[i] = output[i];
            }
        }
    }
    cout << "Estimated range f(x): ";
    for (int i=0; i < sysdim ; i++) {
        cout << "z["<<i << "]=[" << min_output[i] << ", " << max_output[i] <<"]  ";
        range[i] = interval(min_output[i],max_output[i]);
    }
    cout << endl;
    
    return range;
}

// estimate the range of the n iterates f(x) ... f^n(xn)
vector<vector<interval>> estimate_reachset(DiscreteFunc &f, int n, vector<interval> &xinit, int discr) {
    
 
    
    // int discr = 20;
    int nb_points = discr+1;
    
    // limit the number of sampled points
    for (int i=1; i < min(jacdim,4) ; i++)  // MODIF borne min : 1 => 0 (?)
        nb_points = nb_points * (discr+1);
    
    // estimation of the range of f
    vector<vector<double>> input(nb_points,vector<double>(jacdim));  //  the iterates f^n(x_j)
    vector<vector<double>> output(nb_points,vector<double>(sysdim));
    
    vector<vector<double>> max_output(n+1,vector<double>(sysdim));  // store the min and max for each iterate
    vector<vector<double>> min_output(n+1,vector<double>(sysdim));
    vector<vector<interval>> range(n+1,vector<interval>(sysdim));
    
    
    cout << "Initial condition x: " << xinit;
    cout << endl;
    
    // choosing the sampling points in the initial box
    int cur_point = 0;
    for (int i1=0; i1 <= discr ; i1++)
    {
        input[cur_point][0] = xinit[0].inf() + (2.0*i1*xinit[0].rad())/discr;
        if (jacdim > 1)
        {
            for (int i2=0; i2 <= discr ; i2++)
            {
                input[cur_point][0] = xinit[0].inf() + (2.0*i1*xinit[0].rad())/discr;
                input[cur_point][1] = xinit[1].inf() + (2.0*i2*xinit[1].rad())/discr;
                if (jacdim > 2)
                {
                    for (int i3=0; i3 <= discr ; i3++)
                    {
                        input[cur_point][0] = xinit[0].inf() + (2.0*i1*xinit[0].rad())/discr;
                        input[cur_point][1] = xinit[1].inf() + (2.0*i2*xinit[1].rad())/discr;
                        input[cur_point][2] = xinit[2].inf() + (2.0*i3*xinit[2].rad())/discr;
                        // to limit the number of sampled points
                        if (jacdim > 3) {
                            for (int i4=0; i4 <= discr ; i4++)
                            {
                                input[cur_point][0] = xinit[0].inf() + (2.0*i1*xinit[0].rad())/discr;
                                input[cur_point][1] = xinit[1].inf() + (2.0*i2*xinit[1].rad())/discr;
                                input[cur_point][2] = xinit[2].inf() + (2.0*i3*xinit[2].rad())/discr;
                                input[cur_point][3] = xinit[3].inf() + (2.0*i4*xinit[3].rad())/discr;
                                
                                // input[cur_point][3] = (xinit[3].inf()+xinit[3].sup())/2.0;
                                if (jacdim > 4) {
                                    if (xinit[4].inf() != xinit[4].sup())
                                        printf("warning, case not fully implemented");
                                    input[cur_point][4] = (xinit[4].inf()+xinit[4].sup())/2.0;
                                }
                        //        cur_point++; // AJOUT
                            }
                            
                        }
                        cur_point++;
                    }
                }
                else
                    cur_point++;
            }
        }
        else
            cur_point++;
    }
  //  cout << "nb_points=" << nb_points << " cur_point=" << cur_point << endl;
    
    
    
    
    ofstream samplesreachsetfile;
    samplesreachsetfile.open("output/samplesreachset.yaml");
    
    YAML::Emitter out_samples;
    out_samples << YAML::BeginMap;
    
    out_samples << YAML::Key << "samples";
    out_samples << YAML::Value << YAML::BeginSeq;
    
    
    for (int iter=1 ; iter <=n ; iter++)
    {
        for (int i=0; i < sysdim ; i++) {
            max_output[iter] = f(input[0]);
            min_output[iter] = f(input[0]);
        }
        
        for (cur_point=0 ; cur_point<nb_points; cur_point++)
        {
            output[cur_point] = f(input[cur_point]);
            
            if (iter % printing_period == 0)
            {
                out_samples << YAML::BeginMap;
                out_samples << YAML::Key << "tn";
                out_samples << YAML::Value << iter;
                out_samples << YAML::Key << "sample";
                out_samples << YAML::Value << output[cur_point];
                out_samples << YAML::EndMap;
            }
        
            for (int i=0; i < sysdim ; i++) {
                if (output[cur_point][i] < min_output[iter][i])
                    min_output[iter][i] = output[cur_point][i];
                if (output[cur_point][i] > max_output[iter][i])
                    max_output[iter][i] = output[cur_point][i];
            }
            // initializing next step (iter)
            for (int i=0; i < sysdim ; i++)
                input[cur_point][i] = output[cur_point][i];
        }
        
        if (iter % printing_period == 0)
        {
            cout << "Estimated reachable set f^n(x) at step " << iter << " is: ";
            for (int i=0; i < sysdim ; i++) {
                cout << "z["<<i << "]=[" << min_output[iter][i] << ", " << max_output[iter][i] <<"]  ";
            }
            cout << endl;
        }
        for (int i=0; i < sysdim ; i++)
            range[iter][i] = interval(min_output[iter][i],max_output[iter][i]);
    }
    
    out_samples << YAML::EndSeq;
    out_samples << YAML::EndMap;
    samplesreachsetfile << out_samples.c_str();
    samplesreachsetfile.close();
    
    
    return range;
}




// for discrete dynamical systems
void print_projections(vector<interval> &z_inner, vector<interval> &z_inner_rob, vector<interval> &z_outer, int step, vector<interval> &range)
{
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    
    
    vector<double> temp(2*sysdim);
    
    out_approx << YAML::Key << "outer";
    for (int i=0 ; i<sysdim ; i++) {
        temp[2*i] = inf(z_outer[i]);
        temp[2*i+1] = sup(z_outer[i]);
    }
    out_approx << YAML::Value << temp; // Xouter does not work because of interval type (I guess I could solve this but...)
    
    out_approx << YAML::Key << "inner";
    for (int i=0 ; i<sysdim ; i++) {
        temp[2*i] = inf(z_inner[i]);
        temp[2*i+1] = sup(z_inner[i]);
    }
    out_approx << YAML::Value << temp;
    
    if (jacdim > sysdim) {
        out_approx << YAML::Key << "innerrobust";
        for (int i=0; i < sysdim ; i++) {
            temp[2*i] = inf(z_inner_rob[i]);
            temp[2*i+1] = sup(z_inner_rob[i]);
        }
        out_approx << YAML::Value << temp;
    }
    
    
    // error measures
    vector<double> eta(sysdim);
    vector<double> gamma(sysdim);
    cout << "Outer range of f(x) by mean-value:               ";
    for (int i=0; i < sysdim ; i++) {
        cout << "z_o[" << i << "]=" << z_outer[i].inf() << ", " << z_outer[i].sup() << "] ";
        eta[i] = (range[i].sup() - range[i].inf())/ (z_outer[i].sup() - z_outer[i].inf()); // precision metric with respect to estimated range
        cout << " eta_o["<<i<<"]=" << eta[i] <<" ";
    }
    cout << endl;
    out_approx << YAML::Key << "etaouter";
    out_approx << YAML::Value << eta;
    
    cout << "Projection of inner range of f(x) by mean-value:";
    for (int i=0; i < sysdim ; i++) {
        cout << " z_i[" << i << "]=" << z_inner[i];
        eta[i] = (z_inner[i].sup() - z_inner[i].inf())/(range[i].sup() - range[i].inf());
        gamma[i] = (z_inner[i].sup() - z_inner[i].inf())/(z_outer[i].sup() - z_outer[i].inf());
        cout << " eta_i[" << i << "]=" << eta[i];
        cout << " gamma[" << i << "]=" << gamma[i] <<" ";
    }
    cout << endl;
    out_approx << YAML::Key << "etainner";
    out_approx << YAML::Value << eta;
    out_approx << YAML::Key << "gamma";
    out_approx << YAML::Value << gamma;
    
}

// for discrete dynamical systems
void print_innerbox(vector<interval> &inner, vector<int> &exist_quantified, int varx, int vary, int step)
{
    vector<double> temp(4);
    
    out_approx << YAML::Key << "maxbox";
    temp[0] = inf(inner[varx]); temp[1] = sup(inner[varx]); temp[2] = inf(inner[vary]); temp[3] = sup(inner[vary]) ;
    out_approx << YAML::Value << temp;
    
    print_pi(exist_quantified);
    cout << inner[varx] << "  " << inner[vary] << endl;
    //outFile_jointinner << step << "\t" << inf(inner[varx]) << "\t" << sup(inner[varx]) << "\t" << inf(inner[vary]) << "\t" << sup(inner[vary]) <<  endl;
    
    //outFile_jointinner.close();
}

vector<vector<double>> print_skewbox(interval &temp_inner_x, interval &temp_inner_y, vector<vector<double>> &A, int varx, int vary, int step) {
    // resulting quadrilatere A * inner is an inner approximation of the range of f
    vector<vector<double>> output_skewedbox;
    
    output_skewedbox = compute_skewbox(temp_inner_x,temp_inner_y,A,varx,vary);
    
    
    //  cout << "inner skewed box: (pi : 1 -> " << exist_quantified[0]+1 << ", 2 -> " << exist_quantified[1]+1 << ")" << endl;
    for (int i=0; i<4; i++)
        cout << "(" << output_skewedbox[i][0] <<", " << output_skewedbox[i][1] << ") ";
    cout << endl;
    
    vector<double> temp2(8);
    for (int p=0; p<=3; p++) {
        temp2[2*p] = output_skewedbox[p][0];
        temp2[2*p+1] = output_skewedbox[p][1];
    }
    
    out_approx << YAML::Key << "maxskew";
    out_approx << YAML::Value << temp2;
    
    return output_skewedbox;
}
