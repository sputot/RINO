/* ============================================================================
 File   : discrete system.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file computes ranges for functions and discrete dynamical systems
 ============================================================================*/


#ifndef DISCRETE_SYSTEM_H
#define DISCRETE_SYSTEM_H

#include "filib_interval.h"
#include "tadiff.h"
#include "badiff.h"
#include "fadiff.h"
#include "fadbad_aa.h"
#include "ode_def.h"

class DiscreteFunc {
public:
    template <class C>
    vector<C> operator()(vector<C> x) {
        vector<C> z(sysdim);
        
        if (syschoice == 1)
            z[0] = x[0]*x[0] - x[0];
        else if (syschoice == 2)
            z[0] = x[1]*x[1] - 2.0*x[0];
        else if (syschoice == 3) // example 3.5 Goldztejn
        {
            z[0] = 81.0*x[0]*x[0] + x[1]*x[1] + 18.0*x[0]*x[1] - 100.0;
            z[1] = x[0]*x[0] + 81.0*x[1]*x[1] + 18.0*x[0]*x[1] - 100.0;
        }
        else if (syschoice == 4) // example 3 CDC
        {
            z[0] = 5.0*x[0]*x[0] + x[1]*x[1]  - 2.0*x[0]*x[1] - 4.0;
            z[1] = x[0]*x[0] + 5.0*x[1]*x[1] - 2.0*x[0]*x[1] - 4.0;
        }
        else if (syschoice == 5) { // example 5.1 Goldstzjen - surprising that square joint inner-approx is empty !?
            z[0] = x[0]*x[0]*x[0]*x[0]*x[0]*x[0] + x[1]*x[1]*x[1]*x[1]*x[1]*x[1] + x[0]*x[1] - 3.0;
            z[1] = x[0]*x[0]*x[0]*x[0]*x[0]*x[0] - x[1]*x[1]*x[1]*x[1]*x[1]*x[1] - x[0]*x[1] + 1.0;
        }
        else if (syschoice == 6) { //
            z[0] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1] + 2.0*x[0]*x[1] - 4.0;
            z[1] = x[0]*x[0]*x[0] - x[1]*x[1]*x[1] - 2.0*x[0]*x[1] + 2.0;
        }
        else if (syschoice == 7) { // skewed inner-approx is not very good - try to see with zonotope ?
            z[0] = 2.0*x[0]*x[0] + 2.0*x[1]*x[1] - 2.0*x[0]*x[1] - 2.0;
            z[1] = x[0]*x[0] - x[1]*x[1] + 4.0*x[0]*x[1] - 4.0;
        }
        else if (syschoice == 8) { // skewed inner-approx is empty - try to see with zonotope ?
            z[0] = 2.0*x[0]*x[0] + 1.0*x[1]*x[1] - 2.0*x[0]*x[1] - 1.0;
            z[1] = x[0]*x[0] - x[1]*x[1] + 3.0*x[0]*x[1] - 3.0;
        }
        else if (syschoice == 9) { //
            z[0] = 2.0*x[0]*x[0] - x[0]*x[1] - 1.0;
            z[1] = x[0]*x[0] + x[1]*x[1]  - 2.0;
        }
        else if (syschoice == 10) { // (z0,z1) empty ???
            z[0] = 2.0*x[0]*x[0] - x[0]*x[1] + x[0]*x[2] + x[2]*x[2] - 3.0;
            z[1] = x[0]*x[0] + x[1]*x[1] - x[2]*x[2] - 1.0;
            z[2] = x[0] + x[1] - 2.0;
        }
        else if (syschoice == 11) { // ok
            z[0] = 2.0*x[0]*x[0] - x[0]*x[1] - 1.0;
            z[1] = x[0]*x[0] + x[1]*x[1]  - 2.0;
            z[2] = x[2] - 1.0;
        }
        else if (syschoice == 12) { // ok
            z[0] = x[0];
            z[1] = 0.707*x[1] + 0.707*x[1]*x[1]  - 2.0;
            z[2] = x[2] - 1.0;
        }
        else if (syschoice == 13) {
            z[0] = 4.0*(1.0-x[0]+x[0]*x[2])-2.0+4.0*x[0]*(x[1]-1.0)-2.0*(x[2]-1.0);
            z[1] = -4.0*(x[1]-1.0)+4.0*x[1]*(x[1]-1.0);
         //   z[2] = 4.0*(1.0-x[0]+x[0]*x[2])-2.0*x[1]+4.0*x[0]*(x[1]-1.0);
            z[2] =x[2]*x[2]*x[2] + x[2]*x[2] + x[2] + 1.0;
        }
        else if (syschoice == 14) {
            z[0] = x[0]*x[0]*x[0] + x[0]*x[0] + x[0] + 1.0;
        }
        else if (syschoice == 15) {  // test model - parallelotope bundles HSCC 2016 p 303
            double Delta = 0.01; // 0.01;
            z[0] = x[0] + (0.5*x[0]*x[0] - 0.5*x[1]*x[1])*Delta;
            z[1] = x[1] + 2.0*x[0]*x[1]*Delta;
        }
        else if (syschoice == 16) {  // SIR epidemic model  - parallelotope bundles HSCC 2016 p 303
            double beta = 0.34;
            double Delta = 0.5;
            double gamma = 0.05;
            z[0] = x[0] * (1.0 - beta*x[1]*Delta);
            z[1] = x[1] * (1.0 + (beta*x[0]-gamma)*Delta);
            z[2] = x[2] + gamma*x[1]*Delta;
        }
        else if (syschoice == 17) {  // Honeybees model  - parallelotope bundles HSCC 2016 p 303-304
            double beta1 = 0.001;
            double beta2 = 0.001;
            double gamma = 0.3;
            double delta = 0.5;
            double alpha = 0.7;
            double Delta = 0.01;
            z[0] = x[0] * (1.0 - (beta1*x[1] + beta2*x[2]) * Delta);
            z[1] = x[1] * (1.0 + (beta1*x[0] - gamma + delta*beta1*x[3] + alpha*beta1*x[4]) * Delta);
            z[2] = x[2] * (1.0 + (beta2*x[0] - gamma + delta*beta2*x[4] + alpha*beta2*x[3]) * Delta);
            z[3] = x[3] * (1.0 - (delta*beta1*x[1] + alpha*beta2*x[2]) * Delta) + gamma*x[1]*Delta;
            z[4] = x[4] * (1.0 - (delta*beta2*x[2] + alpha*beta1*x[1]) * Delta) + gamma*x[2]*Delta;
        }
        else if (syschoice == 18) {  // SIR epidemic model  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
            double beta = 0.34;       // identical to 16 except parameters - gamma is now uncertain (x[3])
            double h = 1.0;
            //            double gamma : x[3] = interval(0.05,0.0675);
            z[0] = x[0] * (1.0 - beta*x[1]*h);
            z[1] = x[1] * (1.0 + (beta*x[0]-x[3])*h);
            z[2] = x[2] + x[3]*x[1]*h;
        }
        else if (syschoice == 19) {  // SIR epidemic model  - first 2 dimensions - parallelotope bundles HSCC 2016 p 303
            double beta = 0.34;
            double Delta = 0.5;
            double gamma = 0.05;
            z[0] = x[0] * (1.0 - beta*x[1]*Delta);
            z[1] = x[1] * (1.0 + (beta*x[0]-gamma)*Delta);
        }
        else if (syschoice == 20) {  // SIR epidemic model  - first 2 dimensions  - Parameter Synthesis for Polynomial Biological Models HSCC 2014 p 239
            double beta = 0.34;       // identical to 16 except parameters - gamma is now uncertain (x[3])
            double h = 1.0;
            //            double gamma : x[2] = interval(0.05,0.0675);
            z[0] = x[0] * (1.0 - beta*x[1]*h);
            z[1] = x[1] * (1.0 + (beta*x[0]-x[2])*h); // (beta*x[0]-x[2])*h);
        }
        else if (syschoice == 21) {  // essai avec param
            double beta = 0.34;       //
            double h = 1.0;
            //            double gamma : x[2] = interval(0.05,0.0675);
            z[0] = x[0]; // * (1.0 - beta*x[1]*h);
            z[1] = x[1] + x[0]+ x[2]; // (beta*x[0]-x[2])*h);
            z[2] = x[2];
        }
        return z;
    }
};


//template <class C>
struct FDiff
{
    template <class C>
    vector<C> operator()(vector<vector<C>> &o_dfdx, vector<C> x) {
        vector<F<C>> loc_x(jacdim);  // Initialize arguments
        vector<C> res_f(sysdim);
        
        for (int i=0 ; i<jacdim ; i++) {
            loc_x[i] = x[i];
            loc_x[i].diff(i,jacdim);            // Differentiate wrt. x
        }
        DiscreteFunc func;             // Instantiate functor
        vector<F<C>> f(func(loc_x));  // Evaluate function and record DAG
        
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                o_dfdx[i][j]=f[i].d(j);      // Value of df/dx
     //   o_dfdy=f.d(1);      // Value of df/dy */
        for (int i=0 ; i<sysdim ; i++)
            res_f[i] = f[i].x();
        return res_f;
    }
};

extern vector<vector<vector<interval>>> constr_eps;  // constraints on noise symbols
extern vector<vector<vector<interval>>> eps_loc;     // consequence on eps=[x]-x0

extern vector<vector<vector<vector<interval>>>> constr_eps_discr;  // same but with additional discretization in each direction
extern vector<vector<vector<vector<interval>>>> eps_loc_discr;     // consequence on [x]-x0  when partitioning 2D region in 4
extern vector<vector<vector<vector<double>>>> extremity_eps_loc_discr;

extern int nb_discr, nb_discr1, nb_discr2;

vector<interval> init_discrete_system(int &nb_steps);
void discrete_dynamical(int &nb_steps, int order);
void discrete_dynamical_method2(int &nb_steps);
void discrete_dynamical_method2_preconditioned(int &nb_steps);
void discrete_dynamical_preconditioned(int &nb_steps, int order);
void discrete_dynamical_preconditioned_3d(int &nb_steps);
void function_range(void);

void constraint_eps(vector<vector<interval>> &Jac_m, vector<vector<AAF>> &JacAff, int m);

void evaluate_projections(vector<interval> &z0, vector<interval> &radx,  vector<vector<interval>> &Jacf);
void evaluate_projections_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf);
void evaluate_projections_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff);
void evaluate_projections_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff);
void evaluate_projections_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff);
void evaluate_projections_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<vector<AAF>>> &HessAff);

vector<interval> evaluate_outerrange(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf);
vector<interval> evaluate_outerrange_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf);
interval evaluate_outerrange_x_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int i);
//interval evaluate_outerrange_x_subdiv_discretize_old(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i);
interval evaluate_outerrange_x_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int i);
vector<interval> evaluate_outerrange_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff);
vector<interval> evaluate_outerrange_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<vector<AAF>>> &HessAff);

// for only one component: used for joint range
interval evaluate_innerrange_x_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2);
//interval evaluate_innerrange_x_subdiv_discretize_old(vector<interval> &z0,  vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2);
interval evaluate_innerrange_x_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2);
vector<interval> evaluate_innerrange(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, bool maximal, vector<int> &exist_quantified);
vector<interval> evaluate_innerrange_robust(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, bool maximal, vector<int> &exist_quantified);
vector<interval> evaluate_precond_innerrange(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<double>> C, int varx, int vary, bool maximal, vector<int> &exist_quantified);
vector<interval> evaluate_innerrange_order2(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf, bool maximal, vector<int> &exist_quantified);
vector<interval> evaluate_innerrange_order2_robust(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<vector<interval>>> &Hessf, bool maximal, vector<int> &exist_quantified);
vector<interval> evaluate_innerrange_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified);
vector<interval> evaluate_innerrange_order2_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, vector<vector<vector<AAF>>> &HessAfff, bool maximal, vector<int> &exist_quantified);

void joint_ranges(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf,  vector<vector<interval>> &Jacf0, vector<vector<vector<interval>>> &Hessf, int varx, int vary);
void joint_ranges_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary);
void joint_ranges_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary);
void joint_ranges_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary);

void print_pi(vector<int> &exist_quantified);
vector<vector<double>> print_skewbox(interval &temp_inner_x, interval &temp_inner_y, vector<vector<double>> &A, int varx, int vary, int step, ofstream &outFile);

void preconditioned_joint_ranges(vector<interval> &z0, vector<interval> &radx, vector<vector<interval>> &Jacf, vector<vector<interval>> &Jacf0, vector<vector<vector<interval>>> &Hessf, int varx, int vary);
void preconditioned_joint_ranges_subdiv(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAfff, int varx, int vary);
void preconditioned_joint_ranges_subdiv_discretize(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary);
void preconditioned_joint_ranges_discretize_simultaneous(vector<interval> &z0, vector<interval> &radx, vector<vector<AAF>> &JacAff, int varx, int vary);

void twodim_discretization_by_quadrant(vector<interval> &radx);

// estimation of exact image by sampling
void estimate_range(DiscreteFunc &f, vector<interval> &xinit);
void estimate_reachset(DiscreteFunc &f, int n, vector<interval> &xinit);


// for discrete-time dynamical systems
void print_projections(vector<interval> &z_inner, vector<interval> &z_inner_rob, vector<interval> &z_outer, int step);
void print_innerbox(vector<interval> &inner, vector<int> &exist_quantified, int varx, int vary, int step);

#endif
