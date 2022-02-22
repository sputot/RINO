/* ============================================================================
 File   : dde_integr.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 Defines the classes and functions used to actually perform reachability analysis function for DDEs
 ============================================================================ */
#ifndef DDE_H
#define DDE_H

#include "filib_interval.h"
#include "tadiff.h"
#include "fadiff.h"

#include "fadbad_aa.h"
//#include "hybrid.h"
#include "matrix.h"
#include "ode_def.h"
//#include "taylor.h"

#include <iostream>
#include <ostream>
#include <fstream>

using namespace std;

vector<vector<interval>> estimate_reachset_dde(DdeFunc &bf, int discr);

// class used for computation of Taylor coefficients and Jacobian with FADBAD++
class DdeVar
{
public:
    vector<T<F<AAF> >> x; // Independent variables
    vector<T<F<AAF> >> xp;  // Dependent variables
    
    DdeVar() : x(sysdim), xp(sysdim)
    {
        
    }
    
    DdeVar(DdeFunc f, vector< T<F<AAF>>> &_x_prev, vector<AAF> &param_inputs)  : x(sysdim), xp(sysdim)
    {
        f(xp,x,_x_prev,param_inputs); // record DAG at construction:
    }
    void reset()
    {
        for (int i=0 ; i<sysdim ; i++)
            xp[i].reset();
    }
    
    void print(int min, int max);
    void printTM(int order);    // TM with interval coeefficients format
};


// class used for computation of Taylor coefficients with FADBAD++
class Dde
{
public:
    vector< T<AAF>> x; // Independent variables
    vector< T<AAF>> xp;  // Dependent variables
    
    Dde()  : x(sysdim), xp(sysdim)
    {
        
    }
    
    Dde(DdeFunc f, vector< T<AAF>> &_x_prev, vector<AAF> &param_inputs)  : x(sysdim), xp(sysdim)
    {
        // record DAG at construction:
        f(xp,x,_x_prev,param_inputs);
    }
    void reset()
    {
        for (int i=0 ; i<sysdim ; i++)
            xp[i].reset();
    }
    
    
    void print(int min, int max);
    void printTM(int order);
    
};

class Dde_Jac
{
public:
    vector<vector< T<AAF>>> J; // Independent variables
    vector<vector< T<AAF>>> Jp;  // Dependent variables
    
    Dde_Jac() : J(sysdim,vector<T<AAF>>(jacdim)), Jp(sysdim,vector<T<AAF>>(jacdim)) {}
    
    Dde_Jac(DdeJacFunc f, vector<vector< T<AAF>>> &J_prev, vector< T<AAF>> &x, vector< T<AAF>> &x_prev) : J(sysdim,vector<T<AAF>>(jacdim)), Jp(sysdim,vector<T<AAF>>(jacdim))
    {
        f(Jp,J,J_prev,x,x_prev);
    }
    void reset()
    {
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jp[i][j].reset();
    }
};

//void dde_printinit();

// compute an a priori rough enclosure of flow (affine forms)
vector<AAF> fixpoint(DdeFunc bf, vector<AAF> &x0, vector<AAF> &x0_prec, vector<AAF> &param_inputs, double tau);

// compute an a priori  rough enclosure of Jacobian of flow
//iMatrix fixpoint(iMatrix &Jac1_g_rough, iMatrix &J0, double tau);
//void fixpoint(AAF J_rough[sysdim][sysdim], AAF Jac1_g_rough[sysdim][sysdim], AAF J0[sysdim][sysdim], double tau);

// compute an a priori rough enclosure of Jacobian (affine forms) - for variational version
vector<vector<AAF>> fixpoint(DdeJacFunc bf, vector<vector<AAF>> &J0, vector<vector<AAF>> &J0_prec, vector<AAF> &x0, vector<AAF> &x0_prec, double tau);


// enclosure by Taylor model of order order of the flow g_i(X0) at time tau
void gTaylor(vector<AAF> &g, Dde ode_x0, Dde ode_g0, double tau, int order);


// enclosure by Taylor model of order order of the Jacobian of the flow g_i(X0)
//void JTaylor(AAF J_res[sysdim][sysdim], DdeVar odeVAR_x, DdeVar odeVAR_g, AAF J0[sysdim][sysdim], AAF J_rough[sysdim][sysdim], double tau, int order);


class Dde_TM_val
{
public:
    vector<Dde> prev_dde_x;      // at each sub-time step Taylor coeff 0 to k-1
    vector<Dde> prev_dde_g;      // at each sub-time step k-th Taylor coeff
    vector<Dde> dde_x;      // at each sub-time step Taylor coeff 0 to k-1
    vector<Dde> dde_g;      // at each sub-time step k-th Taylor coeff
    vector<vector<AAF>> x;  // value at time tn+j*tau
    vector<vector<AAF>> xp1;  // value at time tn+d0
    int p;       // nb of Taylor models stored on [tn,tn+d0]
    
    Dde_TM_val(DdeFunc bf, vector<Dde> _prev_dde_x, vector<Dde> _prev_dde_g, int _p) : prev_dde_x(_prev_dde_x), prev_dde_g(_prev_dde_g),
    dde_x(vector<Dde>(_p+1,Dde())), dde_g(vector<Dde>(_p+1,Dde())), x(_p+1,vector<AAF>(sysdim)),
    xp1(_p+1,vector<AAF>(sysdim)), p(_p)
    {
        // Taylor Models
        for (int i=0; i<p+1 ; i++)
        {
            if (innerapprox == 0)
            {
                dde_x[i] = Dde(bf, prev_dde_x[i].x,fullinputs_aff);
                dde_g[i] = Dde(bf, prev_dde_g[i].x,fullinputs_aff);
            }
            else
            {
                dde_x[i] = Dde(bf, prev_dde_x[i].x,center_fullinputs_aff);
                dde_g[i] = Dde(bf, prev_dde_g[i].x,center_fullinputs_aff);
            }
        }
    }
    
    Dde_TM_val(DdeFunc bf, int _p) : prev_dde_x(vector<Dde>(_p+1,Dde())),prev_dde_g(vector<Dde>(_p+1,Dde())),
    dde_x(vector<Dde>(_p+1,Dde())), dde_g(vector<Dde>(_p+1,Dde())), x(_p+1,vector<AAF>(sysdim)),
    xp1(_p+1,vector<AAF>(sysdim)), p(_p)
    {
        // Taylor Models
        for (int i=0; i<p+1 ; i++)
        {
            dde_x[i] = Dde();
            dde_g[i] = Dde();
        }
    }
    
    // compute Taylor models corresponding to initial conditions
    void init_dde(double t0, vector<AAF> &ix,  vector<AAF> &inputs, double tau, int order);
    void build(DdeFunc bf, vector<AAF> &param_inputs, int s, double d0, double tau, int order);
    void eval(vector<AAF> &res, int s, double h, int order);
    void init_nextsmallstep(int s);
};

//class Dde_TM_Jac;

class Dde_TM_Jac
{
public:
   // Dde_TM_Jac *prev;
    vector<Dde> prev_ddeVAR_x;      // at each sub-time step Taylor coeff 0 to k-1
    vector<Dde> prev_ddeVAR_g;      // at each sub-time step k-th Taylor coeff
    vector<Dde> ddeVAR_x;      // at each sub-time step Taylor coeff 0 to k-1
    vector<Dde> ddeVAR_g;      // at each sub-time step k-th Taylor coeff
    vector<Dde_Jac> prev_ddeJAC_x;
    vector<Dde_Jac> prev_ddeJAC_g;
    vector<Dde_Jac> ddeJAC_x;
    vector<Dde_Jac> ddeJAC_g;
  //  vector<vector<vector<AAF>>> J_rough;
    vector<vector<vector<AAF>>> J;
   // vector< float **> J;
    vector<vector<AAF>> x;  // value at time tn+j*tau
    vector<vector<AAF>> xp1;  // value at time tn+d0
    vector<vector<vector<AAF>>> Jp1;
    int p;       // nb of Taylor models stored on [tn,tn+d0]
    
    Dde_TM_Jac(DdeFunc _bf, DdeJacFunc _bjf, vector<Dde> &_prev_ddeVAR_x, vector<Dde> &_prev_ddeVAR_g, vector<Dde_Jac> &_prev_ddeJAC_x,
               vector<Dde_Jac> &_prev_ddeJAC_g, int _p) : prev_ddeVAR_x(_prev_ddeVAR_x), prev_ddeVAR_g(_prev_ddeVAR_g),
    ddeVAR_x(vector<Dde>(_p+1,Dde())), ddeVAR_g(vector<Dde>(_p+1,Dde())), prev_ddeJAC_x(_prev_ddeJAC_x), prev_ddeJAC_g(_prev_ddeJAC_g),
    ddeJAC_x(vector<Dde_Jac>(_p+1,Dde_Jac())), ddeJAC_g(vector<Dde_Jac>(_p+1,Dde_Jac())),
    J(_p+1,vector<vector<AAF>>(sysdim,vector<AAF>(jacdim)))
 , x(_p+1,vector<AAF>(sysdim)), xp1(_p+1,vector<AAF>(sysdim)), Jp1(_p+1,vector<vector<AAF>>(sysdim,vector<AAF>(jacdim))), p(_p)
    {
        // Taylor Models
        for (int i=0; i<p+1 ; i++)
        {
            ddeVAR_x[i] = Dde(_bf, prev_ddeVAR_x[i].x,inputs_aff);
            ddeVAR_g[i] = Dde(_bf, prev_ddeVAR_g[i].x,inputs_aff);
            ddeJAC_x[i] = Dde_Jac(_bjf, prev_ddeJAC_x[i].J, ddeVAR_x[i].x, prev_ddeVAR_x[i].x);
            ddeJAC_g[i] = Dde_Jac(_bjf, prev_ddeJAC_g[i].J, ddeVAR_g[i].x, prev_ddeVAR_g[i].x);
        }
    }
    
    Dde_TM_Jac(DdeFunc _bf, DdeJacFunc _bjf, int _p) : prev_ddeVAR_x(vector<Dde>(_p+1,Dde())),prev_ddeVAR_g(vector<Dde>(_p+1,Dde())), /*prev(NULL),*/ ddeVAR_x(vector<Dde>(_p+1,Dde())), ddeVAR_g(vector<Dde>(_p+1,Dde())), prev_ddeJAC_x(vector<Dde_Jac>(_p+1,Dde_Jac())),
    prev_ddeJAC_g(vector<Dde_Jac>(_p+1,Dde_Jac())), ddeJAC_x(vector<Dde_Jac>(_p+1,Dde_Jac())), ddeJAC_g(vector<Dde_Jac>(_p+1,Dde_Jac())),
     J(_p+1,vector<vector<AAF>>(sysdim,vector<AAF>(jacdim)))
    , x(_p+1,vector<AAF>(sysdim)), xp1(_p+1,vector<AAF>(sysdim)), Jp1(_p+1,vector<vector<AAF>>(sysdim,vector<AAF>(jacdim))), p(_p)
    {
        // Taylor Models
        for (int i=0; i<p+1 ; i++)
        {
            ddeVAR_x[i] = Dde();
            ddeVAR_g[i] = Dde();
            ddeJAC_x[i] = Dde_Jac();
            ddeJAC_g[i] = Dde_Jac();
        }
    }
    
    // compute Taylor models corresponding to initial conditions on [-d0,0] or [0,d0] depending on the value of t0 (-d0 or 0)
    void init_dde(double t0, vector<AAF> &x, vector<AAF> &inputs, double tau, int order);
    
    void build(DdeFunc bf, DdeJacFunc bbf, int s, double d0, double tau, int order);
    // eval TM at tn+tau and store in xp1
    // void eval_val();
    // eval TM at time tn+h and return value
    void eval_val(vector<AAF> &res, int s, double h, int order);
    
    // void eval_Jac();
    void eval_Jac(vector<vector<AAF>> &J_res, int s, double h, int order);
    
    void eval(int s, double h, int order);
    
    void init_nextsmallstep(int s); // set x := xp1 ; J := Jp1
};



// TM valid on time interval [tn,tn+d0], with discretization with step tau
class HybridStep_dde
{
public:
    DdeFunc bf;        // def of switched system - really useful ?
    DdeJacFunc bbf;    // eq satisfied by the Jacobian
    
    Dde_TM_val TMcenter;
    Dde_TM_Jac TMJac;

    int order;
    double tn;   // current time
    double tau;  // elementary step
    double d0;   // delay occurring in DDE
    int p;       // nb of Taylor models stored on [tn,tn+d0]
    
    HybridStep_dde(DdeFunc _bf, DdeJacFunc _bbf, vector<Dde> _prev_dde_x, vector<Dde> _prev_dde_g, vector<Dde> _prev_ddeVAR_x, vector<Dde> _prev_ddeVAR_g, vector<Dde_Jac> &_prev_ddeJAC_x, vector<Dde_Jac> &_prev_ddeJAC_g,  int _order, double _tn, double _tau, double _d0, int _p) : bf(_bf), bbf(_bbf), TMcenter(_bf,_prev_dde_x,_prev_dde_g,_p), TMJac(_bf,_bbf,_prev_ddeVAR_x,_prev_ddeVAR_g,_prev_ddeJAC_x,_prev_ddeJAC_g,_p), order(_order), tn(_tn), tau(_tau), d0(_d0), p(_p)
    {   }
   
    HybridStep_dde(DdeFunc _bf, DdeJacFunc _bbf, int _order, double _tn, double _tau, double _d0, int _p) : bf(_bf), bbf(_bbf), TMcenter(_bf,_p), TMJac(_bf,_bbf,_p),
    order(_order), tn(_tn), tau(_tau), d0(_d0), p(_p)
    {   }
    
 
    
    // compute Taylor models corresponding to initial conditions
    void init_dde(vector<interval> &sampled_reachset);
    
    // builds ode_x, ode_g
    void TM_build(int s);
    // eval TM at tn+tau and store in xp1
    void TM_eval(int s);
  
    ReachSet TM_evalandprint_solutionstep(int s, vector<interval> &eps, vector<interval> &sampled_reachset);
    
    void init_nextsmallstep(int s);
    HybridStep_dde init_nextbigstep(double tau);
    
    void print_solutionstep(int s, vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xcenter, vector<interval> &sampled_reachset);
    
};


#endif
