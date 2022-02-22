/* ============================================================================
 File   : ode_integr.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file builds the Taylor models and defines the main reachability analysis functions for ODEs
 - compute Taylor model coefficients by automatic differentiation from the ODE systems
 - compute a priori rough enclosures of the flow and its Jacobian (fixpoint functions)
 - compute refined tubes of the flow and its jacobian using the rough enclosure and the Taylor models
 ============================================================================ */

#ifndef ODE_INTEGR_H
#define ODE_INTEGR_H

#include "filib_interval.h"
#include "tadiff.h"
#include "fadiff.h"

#include "fadbad_aa.h"
#include "utils.h"
//#include "hybrid.h"
//#include "taylor.h"

#include <iostream>
#include <ostream>
#include <fstream>


using namespace std;



// class used for computation of Taylor coefficients and Jacobian with FADBAD++
class OdeVar
{
public:
    vector<T<F<AAF>>> cst_params;
    vector<T<F<AAF>>> param_inputs; // params of the ODE that do appear in the Jacobian
    vector<T<F<AAF>>> control_inputs; // params of the ODE that don't appear in the Jacobian such as control given by NN output
    vector<T<F<AAF> >> x; // Independent variables
    vector<T<F<AAF> >> xp;  // Dependent variables

 
    OdeVar() : cst_params(paramsdim), param_inputs(inputsdim), control_inputs(nncontroldim), x(sysdim), xp(sysdim)
    {
    }
    
    OdeVar(OdeFunc f)  : cst_params(paramsdim), param_inputs(inputsdim), control_inputs(nncontroldim), x(sysdim), xp(sysdim)
    {
        for (int i=0; i<nncontroldim; i++)
            control_inputs[i] = nncontrol[i];
        for (int i=0; i<paramsdim; i++)
            cst_params[i] = params[i];
        f(xp,cst_params,param_inputs,control_inputs,x); // record DAG at construction:
    }
    
    OdeVar(OdeFunc f, vector<F<AAF>> control_in)  : cst_params(paramsdim), param_inputs(inputsdim), control_inputs(nncontroldim), x(sysdim), xp(sysdim)
    {
        for (int i=0; i<nncontroldim; i++)
            control_inputs[i] = control_in[i].x();
        for (int i=0; i<paramsdim; i++)
            cst_params[i] = params[i];
        f(xp,cst_params,param_inputs,control_inputs,x); // record DAG at construction:
    }
    
 
    
    void reset()
    {
        for (int i=0 ; i<sysdim ; i++)
            xp[i].reset();
    }
    
    void print(int min, int max);
    void printTM(int order);    // TM with interval coeefficients format
    void printAAFTM(int order);  // TM with AAF coeefficients format
};

// class used for computation of Taylor coefficients with FADBAD++
class Ode
{
public:
    //  T<AAF> x[sysdim];
    vector<T<AAF>> cst_params;
    vector<T<AAF>> param_inputs;
    vector<T<AAF>> control_inputs;
    vector< T<AAF>> x; // Independent variables
    vector< T<AAF>> xp;  // Dependent variables
    
    Ode()  : cst_params(paramsdim), param_inputs(inputsdim), control_inputs(nncontroldim), x(sysdim), xp(sysdim)
    {
        
    }
    
    Ode(OdeFunc f)  : cst_params(paramsdim), param_inputs(inputsdim), control_inputs(nncontroldim), x(sysdim), xp(sysdim)
    {
        // record DAG at construction:
        for (int i=0; i<nncontroldim; i++)
            control_inputs[i] = nncontrol[i];
        for (int i=0; i<paramsdim; i++)
            cst_params[i] = params[i];
        f(xp,cst_params,param_inputs,control_inputs,x);
    }
    
    Ode(OdeFunc f, vector<AAF> control_in)  : cst_params(paramsdim), param_inputs(inputsdim), control_inputs(nncontroldim), x(sysdim), xp(sysdim)
    {
        // record DAG at construction:
        for (int i=0; i<nncontroldim; i++)
            control_inputs[i] = control_in[i];
        for (int i=0; i<paramsdim; i++)
            cst_params[i] = params[i];
        f(xp,cst_params,param_inputs,control_inputs,x);
    }
    
    void reset()
    {
        for (int i=0 ; i<sysdim ; i++)
            xp[i].reset();
    }
    
    void print(int min, int max);
    void printTM(int order);
    
};




// TM valid on time interval [tn,tn+tau]
class TM_val
{
public:
    
  //  OdeFunc bf;   // def of switched system - really useful ?
    Ode ode_x;      // Taylor coeff 0 to k-1
    Ode ode_g;      // k-th Taylor coeff
    int order;
    vector<AAF> x;  // value at time tn
    double tn;
    double tau;
    vector<AAF> xp1;  // value at time tn+tau
    
    TM_val(OdeFunc _bf, vector<AAF> &param_inputs, int o, double t, double h) :  ode_x(_bf), ode_g(_bf), order(o), x(sysdim), tn(t), tau(h), xp1(sysdim){}
    
    TM_val(Ode &o_x, Ode &o_g, int o, double t, double h) : ode_x(o_x), ode_g(o_g),  order(o), x(sysdim),  tn(t), tau(h), xp1(sysdim) {}
    
    TM_val(Ode &o_x, Ode &o_g, int o, vector<AAF> &_x, double t, double h) :  ode_x(o_x), ode_g(o_g), order(o), x(_x),  tn(t), tau(h), xp1(sysdim)  {}
    
    // builds ode_x, ode_g
    void build(OdeFunc _bf, vector<AAF> &params, vector<AAF> &param_inputs /*DiscreteTrans &prev_trans, bool active_discrete_trans*/);
  
    // eval TM at time tn+h and return value in res
    void eval(vector<AAF> &res, double h);
    
    // eval TM at time tn+tau and store value in xp1
    void eval();
    
    // eval range on time interval [tn,tn+tau] and store in xrange
   // void eval_range();
    // eval range on time interval [tn+tau1,tn+tau2] and return value in res
   // void eval_range(vector<AAF> &res, double tau1, double tau2);
    
  //  void change_stepsize(double _tau); // tau := _tau
    void init_nextstep(double _tau); // set x := xp1, tau = _tau
 //   void reinit_nextstep(vector<AAF> &_x, double _tn, double _tau); //
};

// TM valid on time interval [tn,tn+tau]
class TM_Jac
{
public:
    
//    OdeFunc bf;   // def of switched system - really useful ?
    vector<OdeVar> odeVAR_x;      // computed, Taylor coeff 0 to k-1 -- will be evaluated on [x1]...[xi],x0i+1,...x0n
    OdeVar odeVAR_g;      // computed, for k-th Taylor coeff
    vector<vector<AAF>> J_rough; // computed, for k-th Taylor coeff - will be evaluated on all [x]
    int order;
    vector<AAF> x;          // input, value at time tn
    vector<AAF> x0;          // center value at time tn, coming from TM_val
    vector<vector<AAF>> J;  // input, jacobian at time tn: \partial z / \partial z0 (dimension jacdim \times jacdim)
    double tn;
    double tau;
    vector<AAF> xp1;  // computed, value at time tn+tau
    vector<vector<AAF>> Jp1; // computed, Jac at time tn+tau
 //   vector<AAF> xrange; // value on time range [tn,tn+tau]
 //   vector<vector<AAF>> AAF Jrange[sysdim][sysdim]; // computed, Jac on time range [tn,tn+tau]
    
   /* TM_Jac(OdeFunc _bf, vector<AAF> &param_inputs, int o, double t, double h) :  odeVAR_x(_bf), odeVAR_g(_bf), J_rough(sysdim,vector<AAF>(jacdim)), order(o), x(sysdim), J(sysdim,vector<AAF>(jacdim)), tn(t), tau(h), xp1(sysdim) , Jp1(sysdim,vector<AAF>(jacdim)){}
    
    TM_Jac(OdeVar &o_x, OdeVar &o_g, int o, double t, double h) : odeVAR_x(o_x), odeVAR_g(o_g), J_rough(sysdim,vector<AAF>(jacdim)),  order(o), x(sysdim), J(sysdim,vector<AAF>(jacdim)),  tn(t), tau(h), xp1(sysdim) , Jp1(sysdim,vector<AAF>(jacdim)){} */
    
    TM_Jac(vector<OdeVar> &o_x, OdeVar &o_g, int o, vector<AAF> &_x, vector<AAF> &_x0, vector<vector<AAF>> &_J, double t, double h) :  odeVAR_x(o_x), odeVAR_g(o_g), J_rough(sysdim,vector<AAF>(jacdim)), order(o),  x(_x), x0(_x0), J(_J), tn(t), tau(h), xp1(sysdim) , Jp1(sysdim,vector<AAF>(jacdim)) {
       // assign(J,_J);
    }
    
    void build(OdeFunc _bf, vector<AAF> &params, vector<AAF> &param_inputs, vector<AAF> &param_inputs_center);
   
    // eval TM at time tn+h and return value
    void eval_val(vector<AAF> &res, double h);
    // eval range on time interval [tn+tau1,tn+tau2] and return value in res
//    void eval_val_range(vector<AAF> &res, double tau1, double tau2);
    
    void eval_Jac(vector<vector<AAF>> &J_res , double h);
//    void eval_Jac_range(AAF J_res[sysdim][sysdim], double tau1, double tau2);
    
    void eval();
    // eval range value and Jacobian on time interval [tn,tn+tau] and store in xrange and Jrange
 //   void eval_range();
    
  //  void change_stepsize(double _tau); // tau := _tau
    void init_nextstep(double _tau, vector<AAF> &_x0); // set x := xp1 ; J := Jp1
  //  void reinit_nextstep(vector<AAF> &_x, AAF _J[sysdim][sysdim], double _tn, double _tau);
};


// printing the Taylor coefficients and Jacobians of orders between min and max

// compute an a priori rough enclosure of flow (affine forms)
vector<AAF> fixpoint(OdeFunc bf, vector<AAF> &params, vector<AAF> &param_inputs, vector<AAF> &x0, double tau);

// compute an a priori  rough enclosure of Jacobian of flow
//iMatrix fixpoint(iMatrix &Jac1_g_rough, iMatrix &J0, double tau);
void fixpoint(vector<vector<AAF>> &J_rough , vector<vector<AAF>> &Jac1_g_rough, vector<vector<AAF>> &J0, double tau);









/* building Taylor models */

// TM of the solution of ODE starting from the center of IC
// outputs are ode_x0 (Taylor coefficients from 0 to order-1) and ode_g0 (for the remainder term,
// evaluated using a priori enclosure on [tn,tn+tau])
// if active_discrete_trans is true, join range (for last coeff) with that coming from previous tranistion
void build_TM_center(Ode &ode_x0, Ode &ode_g0, OdeFunc bf, vector<AAF> &x0,  double tn, double tau, int order);


// TM of the Jacobian with respect to IC of the solution of the ODE for all the range of IC
// outputs are odeVAR_x (Taylor coefficients from 0 to order-1), odeVAR_g (for the remainder term,
// evaluated using a priori enclosure on [tn,tn+tau]), and J_rough, also used for the remainder term
// if active_discrete_trans is true, join range (for last coeff) with that coming from previous tranistion
void build_TM_Jacobian(OdeVar &odeVAR_x, OdeVar &odeVAR_g, vector<vector<AAF>> &J_rough , OdeFunc bf, vector<AAF> &x, vector<vector<AAF>> &J0, double tn, double tau, int order);


/* end building Taylor models */



class HybridStep_ode
{
public:
    OdeFunc bf;
    
    TM_val TMcenter;
    TM_Jac TMJac;
    
    int order;
    double tn;
    double tau;
    
    
  /*  HybridStep_ode(OdeFunc _bf, vector<AAF> &param_inputs, double _tn, double _tau, int _order) : bf(_bf), TMcenter(_bf, param_inputs, _order, _tn, _tau), TMJac(jacdim), tn(_tn), tau(_tau), order(_order)
    {
        for (int j=0 ; j<jacdim; j++)
            TMJac[j] = TMJac(_bf, param_inputs, _order, _tn, _tau);
    } */
    
    HybridStep_ode(OdeFunc _bf, TM_val &_TMcenter, TM_Jac &_TMJac, double _tn, double _tau, int _order) : bf(_bf), TMcenter(_TMcenter), TMJac(_TMJac), tn(_tn), tau(_tau), order(_order)
    {
        
    }
    
   // void TM_buildval();
    //void TM_buildJac();
    // build from TMcenter and TMJac Taylor Model for outer-approximation of guard expression from mode to dest on [tn;tn+1] and store in gamma[dest]
    void TM_build(vector<AAF> &params, vector<AAF> &param_inputs,vector<AAF> &param_inputs_center);
   // void TM_build();
    
    void TMval_eval(vector<AAF> &res, double h);
    void TM_evalval();
    void TMJac_eval(vector<AAF> &res, double h);
    void TMJac_eval(vector<vector<AAF>> &J_res, double h);
    void TM_evalJac();
    void TM_eval();
    
    ReachSet TM_evalandprint_solutionstep(vector<interval> &eps, double tnp1, vector<interval> &sampled_reachset, int current_subdiv);
    void init_nextstep(vector<AAF> &params, vector<AAF> &param_inputs, double _tau);
    
    
    void print_solutionstep(vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xcenter, vector<interval> &sampled_reachset, int current_subdiv);
    
    
    void eval_valandJacobian_nn(vector<AAF> x, vector<AAF> &param_inputs, double tn, double tau, vector<vector<AAF>> J);

};

HybridStep_ode init_ode(OdeFunc _bf, vector<AAF> &x0,  vector<AAF> &x, double _tn, double _tau, int _order);

// estimate the range of the n iterates (same stepsize as for reachability analysis)
vector<vector<interval>> estimate_reachset(OdeFunc &obf, int discr);
vector<double> RK(OdeFunc &obf, vector<double> &yn, vector<double> &params, vector<double> &param_inputs, vector<double> &control_inputs, double h);

#endif
