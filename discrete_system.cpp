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
using namespace std;


vector<vector<vector<interval>>> constr_eps;  // constraints on noise symbols when partitioning 2D region in 4
vector<vector<vector<interval>>> eps_loc;     // consequence on [x]-x0  when partitioning 2D region in 4

vector<vector<vector<vector<interval>>>> constr_eps_discr;  // same but with additional discretization in each direction
vector<vector<vector<vector<interval>>>> eps_loc_discr;     // consequence on [x]-x0  when partitioning 2D region in 4
vector<vector<vector<vector<double>>>> extremity_eps_loc_discr;

int nb_discr, nb_discr1, nb_discr2;

vector<interval> init_discrete_system(int &nb_steps)
{
    
    if (syschoice == 1) {
        jacdim = 1;  // number of input components
        sysdim = 1; // number of output components
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
    
    vector<interval> res(jacdim);
    
    if (syschoice == 1) {
        res[0] = interval(2,3);
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
        res[2] = interval(0.,0.0);
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
    return res;
}




// iterate function range by computing joint range as input of next iterate
// order 1 (Mean-value) or order 2 (2nd order Taylor model)
void discrete_dynamical(int &nb_steps, int order) {
    
    vector<interval> res(jacdim);
    
    //  int nb_steps;
    
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    
    //  vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
    vector<AAF> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    /*********************************** ORDER 2 Taylor model  ***********************************/
    vector<vector<interval>> dfdx0_i(sysdim,vector<interval>(jacdim)), dfdx0_o(sysdim,vector<interval>(jacdim));
    vector<interval> zf0_o;
    vector<interval> zf0_i;
    vector<vector<vector<interval>>> Hessf_o(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector<vector<vector<interval>>> Hessf_i(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector <F<F<AAF>>> xff_o(jacdim), xff_i(jacdim), zff_o(sysdim), zff_i(sysdim);
    /* end ORDER 2 */
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    vector<int> exist_quantified(jacdim);
    int varx = 0, vary = 1;
    
    DiscreteFunc f;
    FDiff FFunc;
    
    
    open_outputfiles();
    
    //   estimate_range(f,xinit);
     estimate_reachset(f, nb_steps, xinit);
    
    clock_t begin = clock();
    cout << "begin clock" << endl;
    
    // for now, no input, sysdim = jacdim
    z_inner = xinit;
    z_outer = xinit;
    
    // initial state
    print_projections(z_inner,z_inner,z_outer,0);
    print_innerbox(z_inner, exist_quantified, varx, vary, 0);
    
    for (int i=0; i < jacdim ; i++) {
        x_o[i] = z_outer[i];
        x_i[i] = z_inner[i];
        x0_o[i] = mid(z_outer[i]); //+(eps[i].sup()-mid(eps[i]))/2.0;
        x0_i[i] = mid(z_inner[i]);
        radx_o[i] = z_outer[i] - x0_o[i];
        radx_i[i] = z_inner[i] - x0_i[i];
    }
    
    for (int step = 1; step <= nb_steps ; step++)
    {
        // ESTIMER RANGE DE f^n pâr estimate_range(f,xinit,n) ?
        
        z_o = FFunc(JacAff_o,x_o);
        z_i = FFunc(JacAff_i,x_i);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++) {
                Jacf_o[i][j] = JacAff_o[i][j].convert_int();
                Jacf_i[i][j] = JacAff_i[i][j].convert_int();
            }
        //   cout << "Jacf evaluated on all [x]" << endl;
        //    cout << Jacf_o;
        //    cout << Jacf_i;
        
        cout << "Outer approx of f(x), direct evaluation, step "<< step << ": ";
        for (int i=0; i < sysdim ; i++)
            cout << z_o[i].convert_int() << " ";
        cout << endl;
        
        vector<interval> z0_o = f(x0_o);
        vector<interval> z0_i = f(x0_i);
        
        if (order == 2)
        {
            /*********************************** ORDER 2 Taylor model  ***********************************/
            zf0_o = FFunc(dfdx0_o,x0_o);
            zf0_i = FFunc(dfdx0_i,x0_i);
            
            for (int i=0; i < jacdim ; i++) {
                xff_o[i] = z_outer[i];
                xff_o[i].diff(i,jacdim);          // first order
                xff_o[i].x().diff(i,jacdim);      // second order
            }
            zff_o = f(xff_o);
            for (int i=0; i < jacdim ; i++) {
                xff_i[i] = z_inner[i];
                xff_i[i].diff(i,jacdim);          // first order
                xff_i[i].x().diff(i,jacdim);      // second order
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
        
        /*********************************** Evaluation ***********************************/
        
        vector<int> aux;
        // for each input, index of the output in which it is existentially quantified
        
        
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
            z_outer[i] = intersect(z_outer[i],z_o[i].convert_int());
        
        print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step);
        
        
        if (sysdim >= 2) {
            
            //   assert(sysdim <= 3); // higher dimension not treated yet
            
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
                    exist_quantified[j] = j;
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            
            // order 1 - maximal reachability with respect to parameters
            if (order == 1)
                z_inner = evaluate_innerrange(z0_i, radx_i,Jacf_i,false,exist_quantified);
            else if (order == 2)
                z_inner = evaluate_innerrange_order2(z0_i, radx_i, dfdx0_i, Hessf_i, false, exist_quantified);
            cout << "inner box mean-value ";
            print_innerbox(z_inner, exist_quantified, varx, vary, step);
        }
        else
            z_inner = z_inner_proj;
        
        for (int i=0; i < sysdim ; i++) {
            x_o[i] = z_outer[i];
            x_i[i] = z_inner[i];
            x0_o[i] = mid(z_outer[i]); //+(eps[i].sup()-mid(eps[i]))/2.0;
            x0_i[i] = mid(z_inner[i]);
            radx_o[i] = z_outer[i] - x0_o[i];
            radx_i[i] = z_inner[i] - x0_i[i];
        }
        
        
    }
    
    print_finalstats(begin);
    
    
    system("cd GUI; python3 Visu_discrete.py; cd ..");
    
}



// computing at each step the sensitivity with respect to initial values
void discrete_dynamical_method2(int &nb_steps) {
    
    vector<interval> res(jacdim);
    
    //  int nb_steps;
    
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    
    //  vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
    vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    
    DiscreteFunc f;
    FDiff FFunc;
    
    
    open_outputfiles();
    
    //   estimate_range(f,xinit);
   //  estimate_reachset(f, nb_steps, xinit);
    
    clock_t begin = clock();
    cout << "begin clock" << endl;
    
    // initial state
    print_projections(xinit,xinit,xinit,0);
    
    // for now, no input, sysdim = jacdim
    z_outer = xinit;
    
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
        // ESTIMER RANGE DE f^n pâr estimate_range(f,xinit,n) ?
    /*    cout << "x_o=";
        for (int i=0; i < sysdim ; i++)
           cout << x_o[i].x(); */
        z_o = f(x_o);
        
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++) {
                Jacf_o[i][j] = z_o[i].d(j).convert_int(); // JacAff_o[i][j].convert_int();
            }
        //  cout << "Jacf evaluated on all [x]" << endl;
        //   cout << Jacf_o;
        
        cout << "Outer approx of f(x), direct evaluation, step "<< step << ": ";
        for (int i=0; i < sysdim ; i++)
            cout << z_o[i].x().convert_int() << " ";
        cout << endl;
        
        /*********************************** Evaluation ***********************************/
        
        vector<int> aux;
        
        z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
        z_inner_proj = evaluate_innerrange(z0_o, radx_o, Jacf_o, true, aux);
        if (jacdim > sysdim)
            z_inner_proj_rob = evaluate_innerrange_robust(z0_o, radx_o, Jacf_o, true, aux);
        
        for (int i=0; i < sysdim ; i++)
            z_outer[i] = intersect(z_outer[i],z_o[i].x().convert_int());
        
        print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step);
        
        
        // initialize next iteration: z0 = f^n(x0)
        z0_o = f(x0_o);
        for (int i=0; i < sysdim ; i++)
            x0_o[i] = z0_o[i];
        
        for (int i=0; i < sysdim ; i++) {
            x_o[i] = z_o[i];
         //   x0_o[i] = mid(z_o[i].x().convert_int()); //+(eps[i].sup()-mid(eps[i]))/2.0;
     //       radx_o[i] = z_o[i].x().convert_int() - z0_o[i];
        }
        
        
    }
    
    print_finalstats(begin);
    system("cd GUI; python3 Visu_discrete.py; cd ..");
}



// computing at each step the sensitivity with respect to initial values
// same method as discrete_dynamical_method2 but using preconditioning uniquely for printing
void discrete_dynamical_method2_preconditioned(int &nb_steps) {
    
    vector<interval> res(jacdim);
    
    //  int nb_steps;
    
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    
    //  vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
    vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    DiscreteFunc f;
    FDiff FFunc;
    
    
    open_outputfiles();
    
    //   estimate_range(f,xinit);
  //  estimate_reachset(f, nb_steps, xinit);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    int varx = 0, vary = 1; // joint range we want to print/compute
 /*   if (syschoice == 18) {
        varx = 1;
        vary = 2;
    } */
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skinner_joint.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skouter_joint.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
   
    
    clock_t begin = clock();
    cout << "begin clock" << endl;
    
    // for now, no input, sysdim = jacdim
    z_outer = xinit;
    
    // preconditioning: A is center of Jacobian, C its inverse
    vector<vector<double>> A_o(sysdim,vector<double> (sysdim)), C_o(sysdim,vector<double> (sysdim));
    vector<interval> f0_o(sysdim);
    vector<vector<interval>> CJacf_o(sysdim, vector<interval>(jacdim));

    for (int i=0 ; i<sysdim; i++)
        A_o[i][i] = 1.0;
    
    // initial state
    print_projections(xinit,xinit,xinit,0);
    
    output_skewedbox = print_skewbox(xinit[varx], xinit[vary], A_o,  varx,  vary, 0, outFile_skewedouter);
    output_skewedbox = print_skewbox(xinit[varx], xinit[vary], A_o,  varx,  vary, 0, outFile_skewedinner);
    
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
        // ESTIMER RANGE DE f^n pâr estimate_range(f,xinit,n) ?
        /*    cout << "x_o=";
         for (int i=0; i < sysdim ; i++)
         cout << x_o[i].x(); */
        z_o = f(x_o);
        
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++) {
                Jacf_o[i][j] = z_o[i].d(j).convert_int(); // JacAff_o[i][j].convert_int();
            }
        //  cout << "Jacf evaluated on all [x]" << endl;
        //   cout << Jacf_o;
        
        cout << "Outer approx of f(x), direct evaluation, step "<< step << ": ";
        for (int i=0; i < sysdim ; i++)
            cout << z_o[i].x().convert_int() << " ";
        cout << endl;
        
        /*********************************** Evaluation ***********************************/
        
        vector<int> aux;
        
        z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
        z_inner_proj = evaluate_innerrange(z0_o, radx_o, Jacf_o, true, aux);
        if (jacdim > sysdim)
            z_inner_proj_rob = evaluate_innerrange_robust(z0_o, radx_o, Jacf_o, true, aux);
        
        print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step);
        
        
        
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
                    exist_quantified[j] = j;
                // if varx = 1, vary = 2
               // exist_quantified[2] = 0;
                exist_quantified[3] = 1;
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
          
            
            // f0 = C * z0
            for (int i=0 ; i<sysdim; i++)
                f0_o[i] = z0_o[i];
            f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary];
            f0_o[vary] = C_o[vary][vary]*z0_o[vary] + C_o[vary][varx]*z0_o[varx];
            
          
            // CJacf = C * Jacf
            multMiMi(CJacf_o,C_o,Jacf_o);
            
            // outer range
            vector<interval> temp_outer = evaluate_outerrange(f0_o,radx_o,CJacf_o);
            cout << "outer skewed box (mean value): ";
            output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A_o,  varx,  vary, step, outFile_skewedouter);
            
            vector<interval> temp_inner = evaluate_innerrange(f0_o,radx_o,CJacf_o,false,exist_quantified);
            cout << "inner skewed box (mean value): ";
            print_pi(exist_quantified);
            output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A_o,  varx,  vary, step, outFile_skewedinner);
          
        }
        
        
        
        // initialize next iteration: z0 = f^n(x0)
        z0_o = f(x0_o);
        for (int i=0; i < sysdim ; i++)
            x0_o[i] = z0_o[i];
        
        for (int i=0; i < sysdim ; i++) {
            x_o[i] = z_o[i];
            //   x0_o[i] = mid(z_o[i].x().convert_int()); //+(eps[i].sup()-mid(eps[i]))/2.0;
            //       radx_o[i] = z_o[i].x().convert_int() - z0_o[i];
        }
        
        
    }
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
    
    print_finalstats(begin);
    system("cd GUI; python3 Visu_discrete.py; cd ..");
}


// same as discrete_dynamical but with skew box for joint range
void discrete_dynamical_preconditioned(int &nb_steps, int order) {
    
    vector<interval> res(jacdim);
    
//    int nb_steps;
    
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    
      vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
    //vector<AAF> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    /*********************************** ORDER 2 Taylor model  ***********************************/
    vector<vector<interval>> dfdx0_i(sysdim,vector<interval>(jacdim)), dfdx0_o(sysdim,vector<interval>(jacdim));
    vector<interval> zf0_o;
    vector<interval> zf0_i;
    vector<vector<vector<interval>>> Hessf_o(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector<vector<vector<interval>>> Hessf_i(sysdim,vector<vector<interval>>(jacdim,vector<interval>(jacdim)));
    vector <F<F<AAF>>> xff_o(jacdim), xff_i(jacdim), zff_o(sysdim), zff_i(sysdim);
    /* end ORDER 2 */
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    DiscreteFunc f;
    FDiff FFunc, FFunc2;
    
    open_outputfiles();
    
    //   estimate_range(f,xinit);
   // estimate_reachset(f, nb_steps, xinit);
    
    clock_t begin = clock();
    cout << "begin clock" << endl;
    
    vector<int> aux;
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    int varx = 0, vary = 1; // joint range we want to print/compute
    vector<vector<double>> output_skewedbox;
    
   
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skinner_joint.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skouter_joint.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
    
    
    // for now, no input, sysdim = jacdim
    z_inner = xinit;
    z_outer = xinit;
    
    
    
    // preconditioning: A is center of Jacobian, C its inverse
    vector<vector<double>> A_o(sysdim,vector<double> (sysdim)), C_o(sysdim,vector<double> (sysdim));
    vector<vector<double>> A_i(sysdim,vector<double> (sysdim)), C_i(sysdim,vector<double> (sysdim));
    
    vector<interval> f0_o(sysdim), f0_i(sysdim);
    vector<vector<interval>> CJacf_o(sysdim, vector<interval>(jacdim));
    vector<vector<interval>> CJacf_i(sysdim, vector<interval>(jacdim));
    vector<interval> z0_o, z0_i;
    
    for (int i=0 ; i<sysdim; i++) {
        A_i[i][i] = 1.0;
        A_o[i][i] = 1.0;
    }
    
    // initial state
    print_projections(z_inner,z_inner,z_outer,0);
    //   print_innerbox(z_inner, exist_quantified, varx, vary, 0);
    output_skewedbox = print_skewbox(xinit[varx], xinit[vary], A_o,  varx,  vary, 0, outFile_skewedouter);
    output_skewedbox = print_skewbox(xinit[varx], xinit[vary], A_o,  varx,  vary, 0, outFile_skewedinner);
    
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
        for (int i=0; i < jacdim ; i++)
            x_i[i].diff(i,jacdim);
        
        // x_i = A . x_i
        F<AAF> temp = A_i[varx][varx]*x_i[varx] + A_i[varx][vary]*x_i[vary];
        x_i[vary] = A_i[vary][vary]*x_i[vary] + A_i[vary][varx]*x_i[varx];
        x_i[varx] = temp;
        
        z_i = f(x_i);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_i[i][j] = z_i[i].d(j).convert_int();// JacAff_i[i][j].convert_int();
        
        for (int i=0; i < jacdim ; i++)
            x_o[i].diff(i,jacdim);
        
        temp = A_o[varx][varx]*x_o[varx] + A_o[varx][vary]*x_o[vary];
        x_o[vary] = A_o[vary][vary]*x_o[vary] + A_o[vary][varx]*x_o[varx];
        x_o[varx] = temp;
        
        z_o = f(x_o);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_o[i][j] = z_o[i].d(j).convert_int(); //JacAff_o[i][j].convert_int();
        
        
        // x0 = A x0
        interval tempi = A_i[varx][varx]*x0_i[varx] + A_i[varx][vary]*x0_i[vary];
        x0_i[vary] = A_i[vary][vary]*x0_i[vary] + A_i[vary][varx]*x0_i[varx];
        x0_i[varx] = tempi;
        
        tempi = A_o[varx][varx]*x0_o[varx] + A_o[varx][vary]*x0_o[vary];
        x0_o[vary] = A_o[vary][vary]*x0_o[vary] + A_o[vary][varx]*x0_o[varx];
        x0_o[varx] = tempi;
        
        z0_o = f(x0_o);
        z0_i = f(x0_i);
        
        if (order == 2)
        {
            /*********************************** ORDER 2 Taylor model  ***********************************/
            zf0_o = FFunc(dfdx0_o,x0_o);
            zf0_i = FFunc(dfdx0_i,x0_i);
            
            for (int i=0; i < jacdim ; i++) {
                xff_o[i] = z_outer[i];
                xff_o[i].diff(i,jacdim);          // first order
                xff_o[i].x().diff(i,jacdim);      // second order
            }
            zff_o = f(xff_o);
            for (int i=0; i < jacdim ; i++) {
                xff_i[i] = z_inner[i];
                xff_i[i].diff(i,jacdim);          // first order
                xff_i[i].x().diff(i,jacdim);      // second order
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
        
        
        /*********************************** Evaluation ***********************************/
        
        z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
        z_inner_proj = evaluate_innerrange(z0_i, radx_i, Jacf_i, true, aux);
        if (jacdim > sysdim)
            z_inner_proj_rob = evaluate_innerrange_robust(z0_i, radx_i, Jacf_i, true, aux);
        
        print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step);
        
        
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
                exist_quantified[3] = 1;
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
            // 2D preconditioner - just on components varx and vary
            A_o[varx][varx] = mid(Jacf_o[varx][varx]);
            A_o[vary][vary] = mid(Jacf_o[vary][vary]);
            A_o[varx][vary] = mid(Jacf_o[varx][vary]);
            A_o[vary][varx] = mid(Jacf_o[vary][varx]);
            A_i[varx][varx] = mid(Jacf_i[varx][varx]);
            A_i[vary][vary] = mid(Jacf_i[vary][vary]);
            A_i[varx][vary] = mid(Jacf_i[varx][vary]);
            A_i[vary][varx] = mid(Jacf_i[vary][varx]);
            
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
            
            // f0 = C * z0
            for (int i=0 ; i<sysdim; i++)
                f0_o[i] = z0_o[i];
            f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary];
            f0_o[vary] = C_o[vary][vary]*z0_o[vary] + C_o[vary][varx]*z0_o[varx];
            
            for (int i=0 ; i<sysdim; i++)
                f0_i[i] = z0_i[i];
            f0_i[varx] = C_i[varx][varx]*z0_i[varx] + C_i[varx][vary]*z0_i[vary];
            f0_i[vary] = C_i[vary][vary]*z0_i[vary] + C_i[vary][varx]*z0_i[varx];
            
            // CJacf = C * Jacf
            multMiMi(CJacf_o,C_o,Jacf_o);
            multMiMi(CJacf_i,C_i,Jacf_i);
            
            // outer range
            vector<interval> temp_outer = evaluate_outerrange(f0_o,radx_o,CJacf_o);
            cout << "outer skewed box (mean value): ";
            output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A_o,  varx,  vary, step, outFile_skewedouter);
            
            vector<interval> temp_inner = evaluate_innerrange(f0_i,radx_i,CJacf_i,false,exist_quantified);
            cout << "inner skewed box (mean value): ";
            print_pi(exist_quantified);
            output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A_i,  varx,  vary, step, outFile_skewedinner);
            
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
        
    }
    
    
    print_finalstats(begin);
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
    
    system("cd GUI; python3 Visu_discrete.py; cd ..");
}


// old version - to remove
// same as discrete_dynamical but with skew box for joint range
void discrete_dynamical_preconditioned(int &nb_steps) {
    
    vector<interval> res(jacdim);
    
    //    int nb_steps;
    
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    
    vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
    //vector<AAF> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    DiscreteFunc f;
    FDiff FFunc, FFunc2;
    
    open_outputfiles();
    
    //   estimate_range(f,xinit);
    // estimate_reachset(f, nb_steps, xinit);
    
    clock_t begin = clock();
    cout << "begin clock" << endl;
    
    vector<int> aux;
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    int varx = 0, vary = 1; // joint range we want to print/compute
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skinner_joint.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skouter_joint.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
    
    
    // for now, no input, sysdim = jacdim
    z_inner = xinit;
    z_outer = xinit;
    
    // preconditioning: A is center of Jacobian, C its inverse
    vector<vector<double>> A_o(sysdim,vector<double> (sysdim)), C_o(sysdim,vector<double> (sysdim));
    vector<vector<double>> A_i(sysdim,vector<double> (sysdim)), C_i(sysdim,vector<double> (sysdim));
    
    vector<interval> f0_o(sysdim), f0_i(sysdim);
    vector<vector<interval>> CJacf_o(sysdim, vector<interval>(jacdim));
    vector<vector<interval>> CJacf_i(sysdim, vector<interval>(jacdim));
    vector<interval> z0_o, z0_i;
    
    for (int i=0 ; i<sysdim; i++) {
        A_i[i][i] = 1.0;
        A_o[i][i] = 1.0;
    }
    
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
        for (int i=0; i < jacdim ; i++)
            x_i[i].diff(i,jacdim);
        
        // x_i = A . x_i
        F<AAF> temp = A_i[varx][varx]*x_i[varx] + A_i[varx][vary]*x_i[vary];
        x_i[vary] = A_i[vary][vary]*x_i[vary] + A_i[vary][varx]*x_i[varx];
        x_i[varx] = temp;
        
        z_i = f(x_i);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_i[i][j] = z_i[i].d(j).convert_int();// JacAff_i[i][j].convert_int();
        
        for (int i=0; i < jacdim ; i++)
            x_o[i].diff(i,jacdim);
        
        temp = A_o[varx][varx]*x_o[varx] + A_o[varx][vary]*x_o[vary];
        x_o[vary] = A_o[vary][vary]*x_o[vary] + A_o[vary][varx]*x_o[varx];
        x_o[varx] = temp;
        
        z_o = f(x_o);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_o[i][j] = z_o[i].d(j).convert_int(); //JacAff_o[i][j].convert_int();
        
        
        // x0 = A x0
        interval tempi = A_i[varx][varx]*x0_i[varx] + A_i[varx][vary]*x0_i[vary];
        x0_i[vary] = A_i[vary][vary]*x0_i[vary] + A_i[vary][varx]*x0_i[varx];
        x0_i[varx] = tempi;
        
        tempi = A_o[varx][varx]*x0_o[varx] + A_o[varx][vary]*x0_o[vary];
        x0_o[vary] = A_o[vary][vary]*x0_o[vary] + A_o[vary][varx]*x0_o[varx];
        x0_o[varx] = tempi;
        
        z0_o = f(x0_o);
        z0_i = f(x0_i);
        
        
        /*********************************** Evaluation ***********************************/
        
        z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
        z_inner_proj = evaluate_innerrange(z0_i, radx_i, Jacf_i, true, aux);
        if (jacdim > sysdim)
            z_inner_proj_rob = evaluate_innerrange_robust(z0_i, radx_i, Jacf_i, true, aux);
        
        print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step);
        
        
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
                exist_quantified[3] = 1;
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
            // 2D preconditioner - just on components varx and vary
            A_o[varx][varx] = mid(Jacf_o[varx][varx]);
            A_o[vary][vary] = mid(Jacf_o[vary][vary]);
            A_o[varx][vary] = mid(Jacf_o[varx][vary]);
            A_o[vary][varx] = mid(Jacf_o[vary][varx]);
            A_i[varx][varx] = mid(Jacf_i[varx][varx]);
            A_i[vary][vary] = mid(Jacf_i[vary][vary]);
            A_i[varx][vary] = mid(Jacf_i[varx][vary]);
            A_i[vary][varx] = mid(Jacf_i[vary][varx]);
            
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
            
            // f0 = C * z0
            for (int i=0 ; i<sysdim; i++)
                f0_o[i] = z0_o[i];
            f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary];
            f0_o[vary] = C_o[vary][vary]*z0_o[vary] + C_o[vary][varx]*z0_o[varx];
            
            for (int i=0 ; i<sysdim; i++)
                f0_i[i] = z0_i[i];
            f0_i[varx] = C_i[varx][varx]*z0_i[varx] + C_i[varx][vary]*z0_i[vary];
            f0_i[vary] = C_i[vary][vary]*z0_i[vary] + C_i[vary][varx]*z0_i[varx];
            
            // CJacf = C * Jacf
            multMiMi(CJacf_o,C_o,Jacf_o);
            multMiMi(CJacf_i,C_i,Jacf_i);
            
            // outer range
            vector<interval> temp_outer = evaluate_outerrange(f0_o,radx_o,CJacf_o);
            cout << "outer skewed box (mean value): ";
            output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A_o,  varx,  vary, step, outFile_skewedouter);
            
            vector<interval> temp_inner = evaluate_innerrange(f0_i,radx_i,CJacf_i,false,exist_quantified);
            cout << "inner skewed box (mean value): ";
            print_pi(exist_quantified);
            output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A_i,  varx,  vary, step, outFile_skewedinner);
            
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
        
    }
    
    
    print_finalstats(begin);
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
    
    system("cd GUI; python3 Visu_discrete.py; cd ..");
}


void discrete_dynamical_preconditioned_3d(int &nb_steps) {
    
    vector<interval> res(jacdim);
    
   // int nb_steps;
    
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    vector<F<AAF>> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);  // o for outer/over, i for inner/under
  //  vector<AAF> x_o(jacdim), z_o(sysdim), x_i(jacdim), z_i(sysdim);
    vector<vector<AAF>> JacAff_o(sysdim,vector<AAF>(jacdim)), JacAff_i(sysdim,vector<AAF>(jacdim));
    vector<vector<interval>> Jacf_o(sysdim,vector<interval>(jacdim)), Jacf_i(sysdim,vector<interval>(jacdim));
    vector<interval> x0_o(jacdim), radx_o(jacdim), x0_i(jacdim), radx_i(jacdim);    // center and radius x-x0;
    
    
    vector<interval> z_inner, z_inner_proj, z_inner_proj_rob, z_outer;
    
    DiscreteFunc f;
 //   FDiff FFunc, FFunc2;
    
    open_outputfiles();
    
    //   estimate_range(f,xinit);
    estimate_reachset(f, nb_steps, xinit);
    
    clock_t begin = clock();
    cout << "begin clock" << endl;
    
    vector<int> aux;
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    int varx = 0, vary = 1, varz = 2; // joint range we want to print/compute
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skinner_joint.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skouter_joint.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
    
    
    // for now, no input, sysdim = jacdim
    z_inner = xinit;
    z_outer = xinit;
    
    // preconditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A_o(sysdim,vector<double> (sysdim)), C_o(sysdim,vector<double> (sysdim));
    vector<vector<double>> A_i(sysdim,vector<double> (sysdim)), C_i(sysdim,vector<double> (sysdim));
    
    vector<interval> f0_o(sysdim), f0_i(sysdim);
    vector<vector<interval>> CJacf_o(sysdim, vector<interval>(jacdim));
    vector<vector<interval>> CJacf_i(sysdim, vector<interval>(jacdim));
    vector<interval> z0_o, z0_i;
    
    for (int i=0 ; i<sysdim; i++) {
        A_i[i][i] = 1.0;
        A_o[i][i] = 1.0;
    }
    
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
        
        for (int i=0; i < jacdim ; i++)
            x_i[i].diff(i,jacdim);
        
        // x_i = A . x_i
        F<AAF> temp1 = A_i[varx][varx]*x_i[varx] + A_i[varx][vary]*x_i[vary] + A_i[varx][varz]*x_i[varz];
        F<AAF> temp2 = A_i[vary][varx]*x_i[varx] + A_i[vary][vary]*x_i[vary] + A_i[vary][varz]*x_i[varz];
        x_i[varz] = A_i[varz][varx]*x_i[varx] + A_i[varz][vary]*x_i[vary] + A_i[varz][varz]*x_i[varz];
        x_i[varx] = temp1;
        x_i[vary] = temp2;
        
        z_i = f(x_i);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_i[i][j] = z_i[i].d(j).convert_int();// JacAff_i[i][j].convert_int();
        
        
        for (int i=0; i < jacdim ; i++)
            x_o[i].diff(i,jacdim);
        
        temp1 = A_o[varx][varx]*x_o[varx] + A_o[varx][vary]*x_o[vary] + A_o[varx][varz]*x_o[varz];
        temp2 = A_o[vary][varx]*x_o[varx] + A_o[vary][vary]*x_o[vary] + A_o[vary][varz]*x_o[varz];
        x_o[varz] = A_o[varz][varx]*x_o[varx] + A_o[varz][vary]*x_o[vary] + A_o[varz][varz]*x_o[varz];
        x_o[varx] = temp1;
        x_o[vary] = temp2;
        
        z_o = f(x_o);
        
        for (int i=0; i < sysdim ; i++)
            for (int j=0; j < jacdim ; j++)
                Jacf_o[i][j] = z_o[i].d(j).convert_int();// JacAff_i[i][j].convert_int();
        
        
        
       
        
        // x0 = A x0
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
        
        z0_o = f(x0_o);
        z0_i = f(x0_i);
        
        
        /*********************************** Evaluation ***********************************/
        
        z_outer = evaluate_outerrange(z0_o, radx_o, Jacf_o);
        z_inner_proj = evaluate_innerrange(z0_i, radx_i, Jacf_i, true, aux);
        
        if (jacdim > sysdim)
            z_inner_proj_rob = evaluate_innerrange_robust(z0_i, radx_i, Jacf_i, true, aux);
        
        print_projections(z_inner_proj,z_inner_proj_rob,z_outer,step);
        
        
        if (sysdim >= 2) {
            
            //   assert(sysdim <= 3); // higher dimension not treated yet
            
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
                    exist_quantified[j] = j;
                //       cout << "inner box mean-value "; print_innerbox(z_inner, exist_quantified, 1, 2);
            }
            
            
            for (int i=0 ; i<sysdim; i++) {
                A_o[i][i] = 1.0;
                C_o[i][i] = 1.0;
                A_i[i][i] = 1.0;
                C_i[i][i] = 1.0;
            }
            // 2D preconditioner - just on components varx and vary
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++) {
                    A_o[i][j] = mid(Jacf_o[i][j]);
                    A_i[i][j] = mid(Jacf_i[i][j]);
                }
            }
        /*    A_o[varx][varx] = mid(Jacf_o[varx][varx]);
            A_o[vary][vary] = mid(Jacf_o[vary][vary]);
            A_o[varz][varz] = mid(Jacf_o[varz][varz]);
            A_o[varx][vary] = mid(Jacf_o[varx][vary]);
            A_o[varx][varz] = mid(Jacf_o[varx][varz]);
            A_o[vary][varx] = mid(Jacf_o[vary][varx]);
            A_o[vary][varz] = mid(Jacf_o[vary][varz]);
            A_o[varz][varx] = mid(Jacf_o[varz][varx]);
            A_o[varz][vary] = mid(Jacf_o[varz][vary]);
            A_i[varx][varx] = mid(Jacf_i[varx][varx]);
            A_i[vary][vary] = mid(Jacf_i[vary][vary]);
            A_i[varz][varz] = mid(Jacf_i[varz][varz]);
            A_i[varx][vary] = mid(Jacf_i[varx][vary]);
            A_i[varx][varz] = mid(Jacf_i[varx][varz]);
            A_i[vary][varx] = mid(Jacf_i[vary][varx]);
            A_i[vary][varz] = mid(Jacf_i[vary][varz]);
            A_i[varz][varx] = mid(Jacf_i[varz][varx]);
            A_i[varz][vary] = mid(Jacf_i[varz][vary]); */
            
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
            
            
            // f0 = C * z0
            for (int i=0 ; i<sysdim; i++)
                f0_o[i] = z0_o[i];
            f0_o[varx] = C_o[varx][varx]*z0_o[varx] + C_o[varx][vary]*z0_o[vary] + C_o[varx][varz]*z0_o[varz];
            f0_o[vary] = C_o[vary][varx]*z0_o[varx] + C_o[vary][vary]*z0_o[vary] + C_o[vary][varz]*z0_o[varz];
            f0_o[varz] = C_o[varz][varx]*z0_o[varx] + C_o[varz][vary]*z0_o[vary] + C_o[varz][varz]*z0_o[varz];
            
            for (int i=0 ; i<sysdim; i++)
                f0_i[i] = z0_i[i];
            f0_i[varx] = C_i[varx][varx]*z0_i[varx] + C_i[varx][vary]*z0_i[vary] + C_i[varx][varz]*z0_i[varz];
            f0_i[vary] = C_i[vary][varx]*z0_i[varx] + C_i[vary][vary]*z0_i[vary] + C_i[vary][varz]*z0_i[varz];
            f0_i[varz] = C_i[varz][varx]*z0_i[varx] + C_i[varz][vary]*z0_i[vary] + C_i[varz][varz]*z0_i[varz];
            
            // CJacf = C * Jacf
            multMiMi(CJacf_o,C_o,Jacf_o);
            multMiMi(CJacf_i,C_i,Jacf_i);
            
            // outer range
            vector<interval> temp_outer = evaluate_outerrange(f0_o,radx_o,CJacf_o);
            cout << "outer skewed box (mean value): ";
            output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A_o,  varx,  vary, step, outFile_skewedouter);
            
            vector<interval> temp_inner = evaluate_innerrange(f0_i,radx_i,CJacf_i,false,exist_quantified);
            cout << "inner skewed box (mean value): ";
            print_pi(exist_quantified);
            output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A_i,  varx,  vary, step, outFile_skewedinner);
            
            z_inner = temp_inner;
            z_outer = temp_outer;
            
            
            for (int i=0; i < sysdim ; i++) {
                x_i[i] = temp_inner[i];
                x0_i[i] =  mid(temp_inner[i]); // f0_i[i]; //
                radx_i[i] = temp_inner[i] - x0_i[i];
                x_o[i] = temp_outer[i];
                x0_o[i] =  mid(temp_outer[i]); // f0_i[i]; //
                radx_o[i] = temp_outer[i] - x0_o[i];
            }
            
        }
        else
            z_inner = z_inner_proj;
        
        
    }
    
    
    print_finalstats(begin);
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
    
    system("cd GUI; python3 Visu_discrete.py; cd ..");
}







void function_range(void) {
    
    int nb_steps;
    vector<interval> xinit = init_discrete_system(nb_steps); // initial condition
    
    nb_discr = 10;
    nb_discr1 = nb_discr;
    nb_discr2 = nb_discr;
    
    
    vector<F<AAF>> x(jacdim);
    vector<vector<interval>> Jacf(sysdim,vector<interval>(jacdim));
    vector<vector<AAF>> JacAff(sysdim,vector<AAF>(jacdim));
    vector<F<AAF>> z(sysdim), z1(sysdim), z2(sysdim);
    vector<interval> x0(jacdim), radx(jacdim);    // center and radius x-x0
    
    
    DiscreteFunc f;
    
   
    estimate_range(f,xinit);
    
    for (int i=0; i < jacdim ; i++) {
        x[i] = xinit[i];
        x0[i] = mid(xinit[i]); //+(eps[i].sup()-mid(eps[i]))/2.0;
        radx[i] = xinit[i] - x0[i];
    }
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    
    z = f(x);
    

    
    cout << "Outer approx of f(x), direct evaluation: ";
    for (int i=0; i < sysdim ; i++)
        cout << z[i].x().convert_int() << " ";
    cout << endl;
    
    
    for (int i=0; i < sysdim ; i++)
        for (int j=0; j < jacdim ; j++) {
            Jacf[i][j] = z[i].d(j).convert_int();
            JacAff[i][j] = z[i].d(j);
        //    cout << "JacAff[i][j]=" << JacAff[i][j] << JacAff[i][j].convert_int();
        }
    cout << "Jacf evaluated on all [x]" << endl;
    cout << Jacf;
    
    
   
 
    
    /*********************************** ORDER 2 Taylor model  ***********************************/
   
    FDiff FFunc;
  /*  vector<vector<AAF>> dfdx(sysdim,vector<AAF>(jacdim));
    vector<AAF> zf = FFunc(dfdx,initial_values);
    cout << "cf= " << zf[0].convert_int() << endl;
    cout << "dfdx= " << dfdx[0][0].convert_int() << endl; */
    
    vector<vector<interval>> dfdx0(sysdim,vector<interval>(jacdim));
    vector<interval> zf0 = FFunc(dfdx0,x0);
    cout << "x0=" << x0;
 //   cout << "cf= " << zf0;
    cout << "dfdx= " << dfdx0[0];
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
    
 evaluate_projections(z0, radx, Jacf);
    //  evaluate_projections_subdiv(z0, JacAff);
    //  evaluate_projections_subdiv_discretize(z0, JacAff);
evaluate_projections_discretize_simultaneous(z0, radx, JacAff);
    // center and order 1 evaluated on x0, 2nd term on [x]
evaluate_projections_order2(z0, radx, dfdx0, Hessf);
    
    // here we actually want to use JacAff and not dfdx0 except for first subdivision
    // A LA FOIS FAUX ET IMPRECIS POUR LE MOMENT: A REPRENDRE
 evaluate_projections_order2_discretize_simultaneous(z0, radx, JacAff, dfdx0, HessAff);
    
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
   
    
    
    
    
    
    // F<AAF> zdiff = z[0].d(0);
    // zdiff.diff(0,jacdim);
    // cout << "zdiff=" << zdiff.d(0).convert_int();
    
    
   // vector<F<AAF>> x(jacdim);
  // vector<vector<interval>> Jacf(sysdim,vector<interval>(jacdim));
/*    vector<vector<F<AAF>>> JacFAff(sysdim,vector<F<AAF>>(jacdim));
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            JacFAff[i][j] = z[i].d(j);
            cout << "JacFAff[i][j] =" << JacFAff[i][j].x().convert_int() << " ";
            for (int k=0; k < jacdim ; k++) {
                Laplf[i][j][k] = JacFAff[i][j].d(k).convert_int();
                cout << "Lapl[i][j][k]=" << Laplf[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
 */
    
    /*********************************/
/*
    cout << endl;
   // cout << "Improved mean-value:" << endl;
    vector<vector<interval>> Jacf_refined(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jacf_refined[i] = vector<interval>(jacdim);
    
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    x[jacdim-1] = mid(initial_values[jacdim-1].convert_int());
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    z1 = f(x);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim-1 ; j++) { // first column evaluated on ([x1],x2_0)
            Jacf_refined[i][j] = z1[i].d(j).convert_int();
            JacAff[i][j] = z1[i].d(j);
        }
        for (int j=jacdim-1; j < jacdim ; j++) { // rest evaluated on box
            Jacf_refined[i][j] = z[i].d(j).convert_int();
            JacAff[i][j] = z[i].d(j);
        }
    }
    cout << "Jacf: 1st column evaluated on ([x1],x2_0), 2nd column on [x]" << endl;
    cout << Jacf_refined;

    evaluate_projections(z0, Jacf_refined);
    evaluate_projections_subdiv(z0, JacAff);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
 */
 
    
 /*
    cout << endl;
   // cout << "Improved mean-value:" << endl;
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    x[0] = mid(initial_values[0].convert_int());
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    z2 = f(x);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim-1 ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z2[i].d(j).convert_int();
    }
    cout << "Jacf: 1st column evaluated on [x], 2nd column on (x0_0,[x1])" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
  */

/*
    cout << endl;
 //   cout << "Improved mean-value:" << endl;
    
    for (int i=0; i < sysdim-1 ; i++) {
        for (int j=0; j < jacdim-1 ; j++)
            Jacf_refined[i][j] = z1[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
    }
    for (int i=sysdim-1; i < sysdim ; i++) {
        for (int j=0; j < jacdim ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z2[i].d(j).convert_int();
    }
    cout << "Jacf: f1: 1st column evaluated on ([x0],x1_0), 2nd column on [x], f2: 1st column evaluated on [x], 2nd column on (x0_0,[x1])" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
  */
    
   /*
    cout << endl;
    //   cout << "Improved mean-value:" << endl;
    
    for (int i=0; i < sysdim-1 ; i++) {
        for (int j=0; j < jacdim-1 ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z2[i].d(j).convert_int();
    }
    for (int i=sysdim-1; i < sysdim ; i++) {
        for (int j=0; j < jacdim ; j++)
            Jacf_refined[i][j] = z1[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
    }
    cout << "Jacf: f1: 1st column evaluated on [x], 2nd column on (x0_0,[x1]), f2: 1st column evaluated on ([x0],x1_0), 2nd column on [x]" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
    */
    
    
    system("cd GUI; python3 Visu_function.py; cd ..");
}




void evaluate_projections(vector<interval> &z0, vector<interval> &radx,  vector<vector<interval>> &Jacf) {
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
    
    cout << "Outer range of f(x) by mean-value:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value: ";
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


void evaluate_projections_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf) {
    interval inner_impro;
    vector<interval> z_inner, z_outer;
    vector<int> aux;
    
    z_outer = evaluate_outerrange_order2(z0, radx, Jacf, Hessf);
    z_inner = evaluate_innerrange_order2(z0, radx, Jacf, Hessf, true, aux);
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by order 2 Taylor model:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by order 2 Taylor model: ";
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
void evaluate_projections_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    z_outer = evaluate_outerrange_discretize_simultaneous(z0, radx, JacAff);
    z_inner = evaluate_innerrange_discretize_simultaneous(z0, radx, JacAff, true, aux);
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by centered mean-value with discretization:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by centered mean-value with discretization: ";
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
void evaluate_projections_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<interval>> &dfdx0, vector<vector<vector<AAF>>> &HessAff) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    z_outer = evaluate_outerrange_order2_discretize_simultaneous(z0, radx, JacAff, dfdx0, HessAff);
    z_inner = evaluate_innerrange_order2_discretize_simultaneous(z0, radx, JacAff, dfdx0, HessAff, true, aux);
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by order 2 model with discretization:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by order 2 model with discretization: ";
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
    cout << "z_outer" << z_outer;
    addViVi(z_outer,z0);
    cout << "z_outer" << z_outer;
    
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
    cout << "z_outer" << z_outer;
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
// I think trouble is I progress in crabe (synchronously on the 2 components)
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
    cout << "z_outer=" << z_outer;
    addViVi(z_outer,z0);
    cout << "z_outer=" << z_outer;
    
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
    cout << "z_outer=" << z_outer;
    
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
    
    cout << "z_outer=" << z_outer;
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
                if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
                if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
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
        if (maximal || (exist_quantified[j] == i)) {// output component i the one where j is existentially quantified {
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
            if (maximal || (exist_quantified[j] == i)) {// output component i the one where j is existentially quantified {
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
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
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
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
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
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
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
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
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
    cout << "outer skewed box (mean value): ";
    output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A,  varx,  vary, -1, outFile_skewedouter);
    
    temp_outer = evaluate_outerrange_order2(f0,radx,CJacf0,CHessf);
    cout << "outer skewed box (order 2): ";
    output_skewedbox = print_skewbox(temp_outer[varx], temp_outer[vary], A,  varx,  vary, -1, outFile_skewedouter);
    
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
    cout << "inner skewed box (mean value): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
    
    temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
    cout << "inner skewed box (order 2): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
    
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
    cout << "inner skewed box (mean value): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
    
    temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
    cout << "inner skewed box (order 2): ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
    
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
        cout << "inner skewed box (mean value): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
        
        temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
        cout << "inner skewed box (order 2): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
        
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
        cout << "inner skewed box (mean value): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
        
        temp_inner = evaluate_innerrange_order2(f0,radx,CJacf0,CHessf,false,exist_quantified);
        cout << "inner skewed box (order 2): ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[varx], temp_inner[vary], A,  varx,  vary, -1, outFile_skewedinner);
    }
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}



void preconditioned_joint_ranges_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
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
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, -1, outFile_skewedouter);
    
    
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
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
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
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
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
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
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
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
            }
        }
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}



void preconditioned_joint_ranges_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
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
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, -1, outFile_skewedouter);
    
    
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
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
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
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
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
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
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
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, -1, outFile_skewedinner);
            }
        }
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}


void preconditioned_joint_ranges_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<interval>> &dfdx0, vector<vector<vector<AAF>>> &HessAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    vector<interval> temp_outer(sysdim);
    vector<interval> temp_inner(sysdim);
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
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
    cout << "outer skewed box (mean-value, 1d discretization): ";
    output_skewedbox = print_skewbox(temp_outer[0], temp_outer[1], A,  varx,  vary, -1, outFile_skewedouter);
    
    // outer range
    temp_outer = evaluate_outerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff);
    cout << "outer skewed box (order 2, 1d discretization): ";
    output_skewedbox = print_skewbox(temp_outer[0], temp_outer[1], A,  varx,  vary, -1, outFile_skewedouter);
    
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
    
    cout << "inner skewed box (mean-value, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    
    cout << "inner skewed box (order 2, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    
    
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
    cout << "inner skewed box (mean-value, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    
    cout << "inner skewed box (order 2, 1d discretization): ";
    print_pi(exist_quantified);
    temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    
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
        
        cout << "inner skewed box (mean-value, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
        
        cout << "inner skewed box (order 2, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
        
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
        cout << "inner skewed box (mean-value, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
        
        cout << "inner skewed box (order 2, 1d discretization): ";
        print_pi(exist_quantified);
        temp_inner = evaluate_innerrange_order2_discretize_simultaneous(f0,radx,CJacAff,Cdfdx0,CHessAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}

void preconditioned_joint_ranges_discretize_simultaneous_sauv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    vector<interval> temp_outer(sysdim);
    vector<interval> temp_inner(sysdim);
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
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
    temp_outer = evaluate_outerrange_discretize_simultaneous(f0,radx,CJacAff);
    cout << "outer skewed box (mean-value, 1d discretization): ";
    output_skewedbox = print_skewbox(temp_outer[0], temp_outer[1], A,  varx,  vary, -1, outFile_skewedouter);
    
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
    
    cout << "inner skewed box (mean-value, 1d discretization): ";
    print_pi(exist_quantified);
    cout << endl;
    temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    
    
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
    cout << "inner skewed box (mean-value, 1d discretization): ";
    print_pi(exist_quantified);
    cout << endl;
    temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
    output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    
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
        
        cout << "inner skewed box (mean-value, 1d discretization): ";
        print_pi(exist_quantified);
        
        temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
        
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
        cout << "inner skewed box (mean-value, 1d discretization): ";
        print_pi(exist_quantified);
        
        temp_inner = evaluate_innerrange_discretize_simultaneous(f0,radx,CJacAff,false,exist_quantified);
        output_skewedbox = print_skewbox(temp_inner[0], temp_inner[1], A,  varx,  vary, -1, outFile_skewedinner);
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
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


void estimate_range(DiscreteFunc &f, vector<interval> &xinit) {
    
    system("rm -r output");
    system("mkdir output");
    
    // estimation of the range of f
    vector<double> input(jacdim);
    vector<double> output(sysdim);
    vector<double> max_output(sysdim);
    vector<double> min_output(sysdim);
    
    
    cout << "Initial condition x: " << xinit;
    cout << endl;
    
    int discr = 100;
    for (int i=0; i < jacdim ; i++)
        input[i] = xinit[i].mid();
    for (int i=0; i < sysdim ; i++) {
        max_output = f(input);
        min_output = f(input);
    }
    
    ofstream outFile_xi;
    stringstream file_name;
    file_name.str("");
    file_name << "output/xi.out";
    outFile_xi.open(file_name.str().c_str());
    for (int i1=0; i1 <= discr ; i1++) {
        input[0] = xinit[0].inf() + (2.0*i1*xinit[0].rad())/discr;
        if (jacdim > 1) {
            for (int i2=0; i2 <= discr ; i2++) {
                input[1] = xinit[1].inf() + (2.0*i2*xinit[1].rad())/discr;
                if (jacdim > 2) {
                    for (int i3=0; i3 <= discr ; i3++) {
                        input[2] = xinit[2].inf() + (2.0*i3*xinit[2].rad())/discr;
                        
                        output = f(input);
                        for (int i=0; i < sysdim ; i++)
                            outFile_xi << output[i] <<  "\t";
                        outFile_xi << endl;
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
                    for (int i=0; i < sysdim ; i++)
                        outFile_xi << output[i] <<  "\t";
                    outFile_xi << endl;
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
            for (int i=0; i < sysdim ; i++)
                outFile_xi << output[i] <<  "\t";
            outFile_xi << endl;
            for (int i=0; i < sysdim ; i++) {
                if (output[i] < min_output[i])
                    min_output[i] = output[i];
                else if (output[i] > max_output[i])
                    max_output[i] = output[i];
            }
        }
    }
    outFile_xi.close();
    cout << "Estimated range f(x): ";
    for (int i=0; i < sysdim ; i++)
        cout << "z["<<i << "]=[" << min_output[i] << ", " << max_output[i] <<"]  ";
    cout << endl;
    
}

// estimate the range of the n iterates f(x) ... f^n(xn)
void estimate_reachset(DiscreteFunc &f, int n, vector<interval> &xinit) {
    
    system("rm -r output");
    system("mkdir output");
    
    int discr = 10;
    int nb_points = discr+1;
    
    // limit the number of sampled points
    for (int i=1; i < min(jacdim,3) ; i++)
        nb_points = nb_points * (discr+1);
    
    // estimation of the range of f
    vector<vector<double>> input(nb_points,vector<double>(jacdim));  //  the iterates f^n(x_j)
    vector<vector<double>> output(nb_points,vector<double>(sysdim));
    
    vector<vector<double>> max_output(n+1,vector<double>(sysdim));  // store the min and max for each iterate
    vector<vector<double>> min_output(n+1,vector<double>(sysdim));
    
    
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
                           if (xinit[3].inf() != xinit[3].sup())
                               printf("warning, case not fully implemented");
                            input[cur_point][3] = (xinit[3].inf()+xinit[3].sup())/2.0;
                            if (jacdim > 4) {
                                if (xinit[4].inf() != xinit[4].sup())
                                    printf("warning, case not fully implemented");
                                input[cur_point][4] = (xinit[4].inf()+xinit[4].sup())/2.0;
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
    
    
    ofstream outFile_xi;
    stringstream file_name;
    file_name.str("");
    file_name << "output/xi.out";
    outFile_xi.open(file_name.str().c_str());
    
    for (int iter=1 ; iter <=n ; iter++)
    {
        for (int i=0; i < sysdim ; i++) {
            max_output[iter] = f(input[0]);
            min_output[iter] = f(input[0]);
        }
        
        for (cur_point=0 ; cur_point<nb_points; cur_point++)
        {
            output[cur_point] = f(input[cur_point]);
            outFile_xi << iter <<  "\t";
            for (int i=0; i < sysdim ; i++)
                outFile_xi << output[cur_point][i] <<  "\t";
            outFile_xi << endl;
            for (int i=0; i < sysdim ; i++) {
                if (output[cur_point][i] < min_output[iter][i])
                    min_output[iter][i] = output[cur_point][i];
                else if (output[cur_point][i] > max_output[iter][i])
                    max_output[iter][i] = output[cur_point][i];
            }
            // initializing next step (iter)
            for (int i=0; i < sysdim ; i++)
                input[cur_point][i] = output[cur_point][i];
        }
        
        cout << "Estimated reachable set f^n(x) at step " << iter << " is: ";
        for (int i=0; i < sysdim ; i++)
            cout << "z["<<i << "]=[" << min_output[iter][i] << ", " << max_output[iter][i] <<"]  ";
        cout << endl;
    }
    outFile_xi.close();
    
    
}




// for discrete dynamical systems
void print_projections(vector<interval> &z_inner, vector<interval> &z_inner_rob, vector<interval> &z_outer, int step)
{
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value: ";
    cout << z_inner;
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << step << "\t"  << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << step << "\t" << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
    
    if (jacdim > sysdim) {
        for (int i=0; i < sysdim ; i++) {
            file_name.str("");
            file_name << "output/x" << i+1 << "inner_robust.out";
            outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
            outFile_outer_mean_value[i] << step << "\t" << inf(z_inner_rob[i]) << "\t" << sup(z_inner_rob[i]) << endl;
            outFile_outer_mean_value[i].close();
        }
    }
    
}

// for discrete dynamical systems
void print_innerbox(vector<interval> &inner, vector<int> &exist_quantified, int varx, int vary, int step)
{
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "inner_joint.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    print_pi(exist_quantified);
    cout << inner[varx] << "  " << inner[vary] << endl;
    outFile_jointinner << step << "\t" << inf(inner[varx]) << "\t" << sup(inner[varx]) << "\t" << inf(inner[vary]) << "\t" << sup(inner[vary]) <<  endl;
    
    outFile_jointinner.close();
}

vector<vector<double>> print_skewbox(interval &temp_inner_x, interval &temp_inner_y, vector<vector<double>> &A, int varx, int vary, int step, ofstream &outFile) {
    // resulting quadrilatere A * inner is an inner approximation of the range of f
    vector<vector<double>> output_skewedbox;
    
    output_skewedbox = compute_skewbox(temp_inner_x,temp_inner_y,A,varx,vary);
    
    /*   output_skewedbox[0][0] = inf(temp_inner_x)*A[varx][varx] + inf(temp_inner_y)*A[varx][vary];
     output_skewedbox[0][1] = inf(temp_inner_x)*A[vary][varx] + inf(temp_inner_y)*A[vary][vary];
     output_skewedbox[1][0] = inf(temp_inner_x)*A[varx][varx] + sup(temp_inner_y)*A[varx][vary];
     output_skewedbox[1][1] = inf(temp_inner_x)*A[vary][varx] + sup(temp_inner_y)*A[vary][vary];
     output_skewedbox[2][0] = sup(temp_inner_x)*A[varx][varx] + sup(temp_inner_y)*A[varx][vary];
     output_skewedbox[2][1] = sup(temp_inner_x)*A[vary][varx] + sup(temp_inner_y)*A[vary][vary];
     output_skewedbox[3][0] = sup(temp_inner_x)*A[varx][varx] + inf(temp_inner_y)*A[varx][vary];
     output_skewedbox[3][1] = sup(temp_inner_x)*A[vary][varx] + inf(temp_inner_y)*A[vary][vary];*/
    
    //  cout << "inner skewed box: (pi : 1 -> " << exist_quantified[0]+1 << ", 2 -> " << exist_quantified[1]+1 << ")" << endl;
    for (int i=0; i<4; i++)
        cout << "(" << output_skewedbox[i][0] <<", " << output_skewedbox[i][1] << ") ";
    cout << endl;
    
    if (step == -1)
    {
        for (int i=0; i<3; i++)
            outFile << output_skewedbox[i][0] << "\t" << output_skewedbox[i][1] << "\t" ;
        outFile << output_skewedbox[3][0] << "\t" << output_skewedbox[3][1] << endl ;
    }
    else
    {
        outFile << step << "\t";
        for (int i=0; i<3; i++)
            outFile << output_skewedbox[i][0] << "\t" << output_skewedbox[i][1] << "\t" ;
        outFile << output_skewedbox[3][0] << "\t" << output_skewedbox[3][1] << endl ;
    }
    
    return output_skewedbox;
}
