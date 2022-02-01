/* ============================================================================
 File   : dde_integr.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 Defines the classes and functions used to actually perform reachability analysis function for DDEs (Delay Differenial Equations)
 ============================================================================ */

#include "filib_interval.h"

#include "tadiff.h"
#include "fadiff.h"
#include "fadbad_aa.h"
#include "utils.h"
#include "matrix.h"
#include "inner.h"
#include "dde_integr.h"
#include <math.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
using namespace std;
//using namespace fadbad;




// printing the Taylor coefficients  of orders between min and max
void Dde::print(int min, int max) {
    for(int i=min;i<=max;i++) {
        for (int j=0 ; j<sysdim ; j++)
            cout << "dde.x[" << j << "][" << i << "]=" << x[j][i] << "\t" ;
        cout << endl;
    }
}

void Dde::printTM(int order) {
    //  cout << "printTM" << endl;
    for (int j=0 ; j<sysdim ; j++) {
        cout << "x[" << j << "]=" << x[j][0].convert_int() << " + ";
        
        for (int i=1;i<order;i++) {
            cout << x[j][i].convert_int() << " t^" << i << " + " ;
        }
        cout << x[j][order].convert_int() << " t^" << order << endl ;
    }
}

vector<AAF> fixpoint(DdeFunc bf, vector<AAF> &x0, vector<AAF> &x0_prec, vector<AAF> &param_inputs, double tau)
{
// calcul de x satisfaisant x0 + [0,tau][f](x,x-tau) \subseteq x
    // valeur de x-tau = ici evaluation sur le range [0,tau]
    vector<AAF> y0(sysdim);
    vector<AAF> y1(sysdim);
    vector<AAF> fx0(sysdim);
    int iter;
    interval widen, coeff;
    
  //  cout << "Entering fixpoint with x0[0]=" << x0[0].convert_int() << "x0_prec[0]=" << x0_prec[0].convert_int() << endl;
    
    bf(fx0,x0,x0_prec,param_inputs);
    
    for (int i=0; i<sysdim ; i++)
        y1[i] = x0[i] + fx0[i].convert_int();
    
    y0 = x0;

    iter = 1;
    widen = interval(-1,1);

    while (iter <=1 || !subseteq(y1,y0))
    {
        if (iter > 25)
            coeff = 1;
        else if (iter > 20)
            coeff = 0.1;
        else if (iter > 15)
            coeff = 0.01;
        else if (iter > 10)
            coeff = 0.001;
    else if (iter > 5)
        coeff = 0.0001;
    else if (iter > 2)
        coeff = 0.00001;
        
        if (iter > 2)
        {
            for (int i=0; i<sysdim ; i++)
                y1[i] = y1[i] + coeff*widen*y1[i].convert_int();
        }
        y0 = y1;
    
        bf(fx0,y0,x0_prec,param_inputs);
   // cout << "fx[0]=" << fx0[0] << "\t"  << fx0[1] << endl;
 //   cout << "y1=" << y1[0] << "\t"  << y1[1] << endl;
    for (int i=0; i<sysdim ; i++)
        y1[i] = x0[i] + interval(0,tau)*fx0[i].convert_int();
      
    iter = iter+1;
 }
//cout << "y0 (end of fixpoint)=" << y0[0].convert_int() << " iter=" << iter << endl;
return y0;
}



// compute an a priori rough enclosure of Jacobian (affine forms) - for variational version
vector<vector<AAF>> fixpoint(DdeJacFunc bf, vector<vector<AAF>> &J0, vector<vector<AAF>> &J0_prec, vector<AAF> &x0, vector<AAF> &x0_prec, double tau)
{
    // calcul de x satisfaisant x0 + [0,tau][f](x,x-tau) \subseteq x
    // valeur de x-tau = ici evaluation sur le range [0,tau]
    vector<vector<AAF>> Jy0(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> Jy1(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> fJ0(sysdim, vector<AAF>(jacdim));
    int iter;
    interval widen, coeff;
    
   // cout << "Entering fixpoint with x0[0]=" << x0[0].convert_int() << "x0_prec[0]=" << x0_prec[0].convert_int() << endl;
    
    bf(fJ0,J0,J0_prec,x0,x0_prec);
    
    for (int i=0; i<sysdim ; i++)
        for (int j=0; j<jacdim ; j++)
        Jy1[i][j] = J0[i][j] + fJ0[i][j].convert_int();
    
    Jy0 = J0;
    
    iter = 1;
    widen = interval(-1,1);
    
    while (iter <=1 || !subseteq(Jy1,Jy0))
    {
        if (iter > 25)
            coeff = 1;
        else if (iter > 20)
            coeff = 0.1;
        else if (iter > 15)
            coeff = 0.01;
        else if (iter > 10)
            coeff = 0.001;
        else if (iter > 5)
            coeff = 0.0001;
        else if (iter > 2)
            coeff = 0.00001;
        
        if (iter > 2)
        {
            for (int i=0; i<sysdim ; i++)
                for (int j=0; j<jacdim ; j++)
                Jy1[i][j] = Jy1[i][j] + coeff*widen*Jy1[i][j].convert_int();
        }
        Jy0 = Jy1;
        
        //bf(fx0,y0,x0_prec);
        bf(fJ0,Jy0,J0_prec,x0,x0_prec);
        // cout << "fx[0]=" << fx0[0] << "\t"  << fx0[1] << endl;
        //   cout << "y1=" << y1[0] << "\t"  << y1[1] << endl;
        for (int i=0; i<sysdim ; i++)
            for (int j=0; j<jacdim ; j++)
            Jy1[i][j] = J0[i][j] + interval(0,tau)*fJ0[i][j].convert_int();
        
        iter = iter+1;
    }
    //cout << "y0 (end of fixpoint)=" << y0[0].convert_int() << " iter=" << iter << endl;
    return Jy0;
}


//FadbadVarODE(int n,TFfunction f,void*param= 0);

// enclosure by Taylor model of order order of the flow g_i(center(X0)) at time tau
void gTaylor(vector<AAF> &g, Dde dde_x, Dde dde_g, double tau, int order)
{
    double taui = tau;
    for (int j=0 ; j<sysdim; j++)
        g[j] = (dde_x.x[j])[0];
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            g[j] += dde_x.x[j][i]*taui;
        taui *= tau;
    }
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        g[j] += dde_g.x[j][order]*taui;
}




// sets initial conditions and parameters for the ODE system
void  Dde_TM_val::init_dde(double t0, vector<AAF> &ix, vector<AAF> &inputs, double tau, int order)
{
    // when solution is given explicitely on [-d0,tau]
    
    // init by computing in x Taylor model of solution on [-d0,0]
    T<AAF> t;

 //   cout << "Start initialisation TM_val with ix0=" << beta[0].convert_int() << endl;
    t = t0; // -d0; // -d0 + k.tau
    for (int j=0; j<p+1 ; j++) // extended all vectors to size p+1 to be able to store these initial conditions - think again
    {
        t[1] = 1; // Taylor-expand wrt. t (dt/dt=1)
        dde_x[j].x = Initfunc(t,ix,inputs);    // InitByFunc(dde_x[j].x,t); // marche po :(
      
        for (int i=0 ; i<sysdim ; i++)
            dde_x[j].x[i].eval(order);
        
        for (int i= 0 ; i<sysdim; i++) {
            x[j][i] = dde_x[j].x[i][0];
             cout << "x["<<j<<"]["<<i<<"]="<<x[j][i].convert_int() << endl;
        }
    
        t = t + tau;
        if (j>0)
            xp1[j-1] = x[j];
    }
//    xp1[p-1] = x[p];
    
    t = interval(t0,t0+tau); // interval(-d0,-d0+tau); // -d0 + k.tau
    for (int j=0; j<p ; j++) // extended all vectors to size p+1 to be able to store these initial conditions - think again
    {
        t[1] = 1; // Taylor-expand wrt. t (dt/dt=1)
        dde_g[j].x = Initfunc(t,ix,inputs);
        for (int i=0 ; i<sysdim ; i++)
            dde_g[j].x[i].eval(order);
        t = t + tau;
    }
    
    cout << "End initialisation TM_val" << endl<< endl;
}


// A completer. sets initial conditions and parameters for the ODE system
void  Dde_TM_Jac::init_dde(double t0, vector<AAF> &ix, vector<AAF> &inputs, double tau, int order)
{
    vector<T<F<AAF>>> beta_initial(sysdim);
    vector<T<F<AAF>>> beta_inputs(jacdim-sysdim);
    for (int i= 0 ; i<sysdim; i++)
        beta_initial[i] = ix[i];
    for (int i= 0 ; i<jacdim-sysdim; i++)
        beta_inputs[i] = inputs[i];
    
    vector<DdeVar> x_init(p+1);
    vector<DdeVar> g_init(p+1);
    
    // initialize jacobian of the flow (on interval x)
  //  setId(J[0]);
    
    // when solution is given explicitely on [-d0,tau]
    
    // init by computing in x Taylor model of solution on [-d0,0]
    T<F<AAF>> t;
    
   // cout << "Start initialisation TM_Jac with ix=" << beta[0][0].x().convert_int() << endl;
    t = t0; // -d0; // -d0 + k.tau
    for (int j=0; j<p+1 ; j++) // extended all vectors to size p+1 to be able to store these initial conditions - think again
    {
        t[1] = 1; // Taylor-expand wrt. t (dt/dt=1)
 
         // specify the variables to differentiate
        for (int i=0 ; i<sysdim ; i++)
            beta_initial[i][0].diff(i,jacdim);   // ddeVAR_x[j].x[i][0].diff(i,sysdim);
        for (int i= 0 ; i<jacdim-sysdim; i++)
            beta_inputs[i][0].diff(sysdim+i,jacdim);
        
        x_init[j].x  = Initfunc(t,beta_initial,beta_inputs); // ddeVAR_x[j].x = Initfunc(t,beta);
        
        
        // InitByFunc(dde_x[j].x,t); // marche po :(
        for (int i=0 ; i<sysdim ; i++)
            x_init[j].x[i].eval(order); // ddeVAR_x[j].x[i].eval(order);
        
        
        for (int i= 0 ; i<sysdim; i++) {
            x[j][i] = x_init[j].x[i][0].x(); // ddeVAR_x[j].x[i][0].x();
            ddeVAR_x[j].x[i] = x_init[j].x[i][0].x();
      //      cout << "x["<<j<<"]["<<i<<"]="<<x[j][i].convert_int() << endl;
        }
        for (int i=0 ; i<sysdim ; i++)
            for (int k=0 ; k<jacdim ; k++) {
                J[j][i][k] = x_init[j].x[i][0].d(k); // ddeVAR_x[j].x[i][0].d(k);
                ddeJAC_x[j].J[i][k] = x_init[j].x[i][0].d(k);
          //      cout << "J["<<j<<"]["<<i<<"]["<<k<<"]="<<J[j][i][k].convert_int() << endl;
            }
        
        
        for(int l=1;l<=order;l++) {
            for (int i=0 ; i<sysdim ; i++)
            {
                ddeVAR_x[j].x[i][l] = x_init[j].x[i][l].x();
                for (int k=0 ; k<jacdim ; k++) {
                    ddeJAC_x[j].J[i][k][l] = x_init[j].x[i][l].d(k);
                }
            }
         //    cout << "(1/k!)*(d^"<<l<<"f/dx^" << l << ")=" <<  x_init[j].x[0][l].d(0).convert_int() /*ddeVAR_x[j].x[0][i].d(0).convert_int()*/ << endl; // The i'th taylor coefficient
        }
        
        t = t + tau;
        if (j>0)
        {
            xp1[j-1] = x[j];
            Jp1[j-1] = J[j];
        }
    }
    
    
    t = interval(t0,t0+tau); // interval(-d0,-d0+tau); // -d0 + k.tau
    for (int j=0; j<p ; j++) // extended all vectors to size p+1 to be able to store these initial conditions - think again
    {
        t[1] = 1; // Taylor-expand wrt. t (dt/dt=1)
    //    for (int i=0 ; i<jacdim ; i++)
    //        beta[i][0].diff(i,jacdim);
         g_init[j].x  = Initfunc(t,beta_initial,beta_inputs); //  ddeVAR_g[j].x = Initfunc(t,beta);
        
        for (int i=0 ; i<sysdim ; i++)
            g_init[j].x[i].eval(order); //  ddeVAR_g[j].x[i].eval(order);
        
        for (int i= 0 ; i<sysdim; i++) {
            ddeVAR_g[j].x[i] = g_init[j].x[i][0].x(); // ddeVAR_g[j].x[i][0].x();
        }
        
        for (int i=0 ; i<sysdim ; i++)
            for (int k=0 ; k<jacdim ; k++) {
                 ddeJAC_g[j].J[i][k] = g_init[j].x[i][0].d(k); // ddeVAR_g[j].x[i][0].d(k);
              //  cout << "J_rough["<<j<<"]["<<i<<"]["<<k<<"]="<<J_rough[j][i][k].convert_int() << endl;
            }
        
        for(int l=1;l<=order;l++) {
            for (int i=0 ; i<sysdim ; i++)
            {
                ddeVAR_g[j].x[i][l] = g_init[j].x[i][l].x();
                for (int k=0 ; k<jacdim ; k++) {
                    ddeJAC_g[j].J[i][k][l] = g_init[j].x[i][l].d(k);
                }
            }
        }

        
        t = t + tau;
        
        // TM_evalandprint_solutionstep(j);
    }
    
    cout << "End initialisation TM_Jac" << endl<< endl;
}


// sets initial conditions and parameters for the ODE system
void  HybridStep_dde::init_dde(/*vector<AAF> &x, vector<AAF> &x0, vector<interval> &eps*/)
{
    vector<interval> Xouter(sysdim),Xouter_robust(sysdim),Xouter_minimal(sysdim),Xinner(sysdim),Xinner_robust(sysdim),Xinner_minimal(sysdim),Xcenter(sysdim);
    
    
    
    if (innerapprox == 0) {
        TMcenter.init_dde(tn,initial_values,fullinputs,tau,order);
        
        for (int i= 0 ; i<sysdim; i++)
            Xouter[i] = TMcenter.x[0][i].convert_int();
        print_solutionstep(-1,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter); // initial solution t=-d0
        
        // print exact range
        for (int j=1; j<=p ; j++) {
            for (int i= 0 ; i<sysdim; i++)
                Xouter[i] = TMcenter.x[j][i].convert_int();
            print_solutionstep(j-1,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter); // initial solution t=-d0
        }
    }
    else
    {
   //     cout << "eps=" << eps[0] << endl;
        TMcenter.init_dde(tn,center_initial_values,center_fullinputs,tau,order);
        TMJac.init_dde(tn,initial_values,fullinputs,tau,order);
        
        //
        for (int i= 0 ; i<sysdim; i++) {
            Xouter[i] = TMJac.x[0][i].convert_int();
            Xcenter[i] = TMcenter.x[0][i].convert_int();
        }
        //
        out_approx << YAML::BeginMap;
        out_approx << YAML::Key << "tn";
        out_approx << YAML::Value << tn;
        
        InnerOuter(Xinner,Xinner_robust,Xinner_minimal,Xouter,Xouter_robust,Xouter_minimal,TMcenter.x[0],TMJac.J[0],eps,tn+tau); //x0p1,Jp1,eps);
        intersectViVi(Xouter,TMJac.x[0]);
        
        print_solutionstep(-1,Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter); // initial solution t=-d0
        
        out_approx << YAML::EndMap;
        
        for (int j=0; j<p ; j++)
            TM_evalandprint_solutionstep(j,eps);
    }
    
  //  out_approx << YAML::EndMap;
 //   cout << " end print, p=" << p << endl;
}

// building the Taylor model for the solution of ODE starting from the center of initial conditions
// outputs are ode_x0 (Taylor coefficients from 0 to order-1) and ode_g0 (giving the remainder Taylor coefficient,
// computed using an a priori enclosure of the solution on [tn,tn+tau])
// if active_discrete_trans is true, then we join with range of previous transition for the enclosure on the range (last coeff of TM)
// tn here is the final time on which this TM should be valid
// buyilding s-th Taylor model on interval [tn,tn+d0], ie on [tn+s*tau,tn+(s+1)*d0/p]
void Dde_TM_val::build(DdeFunc bf, vector<AAF> &param_inputs, int s, double d0, double tau, int order)
{
    vector<AAF> g_rough(sysdim); // rough estimation of solution of current mode on [tn,tn+tau]
    int i, j;
    
    // compute Taylor coefficients / Lie derivatives starting from center of initial conditions
  //  dde_x[s].reset();
    // record DDE
    bf(dde_x[s].xp,dde_x[s].x,prev_dde_x[s].x,param_inputs);
    
    for (j=0 ; j<sysdim ; j++)
        dde_x[s].x[j][0]=x[s][j]; // initialize;
    
 //   cout << "building AD : dde_x[s].x=" << dde_x[s].x[0][0].convert_int() <<  endl;
   
    // computing a priori enclosure on [tn,tn+tau] of solution starting from center of initial conditions
    // x-tau should eb the enclosure on the interval step
    // evaluate the Lie derivatives on this enclosure: we will use the last coefficients as remainder terms of Taylor model
  //  dde_g[s].reset();
    
    vector<AAF> x_prev_range(sysdim);
    for (j=0 ; j<sysdim ; j++)
        x_prev_range[j] = prev_dde_g[s].x[j][0]; // dde_g[s].x_prev[j][0];
  
    g_rough = fixpoint(bf,x[s],x_prev_range,param_inputs,tau); // a verifier !
    for (int i=0 ; i<sysdim ; i++)
        g_rough[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term
   
    // record DDE
    bf(dde_g[s].xp,dde_g[s].x,prev_dde_g[s].x,param_inputs);
    
    for (j=0 ; j<sysdim ; j++)
        dde_g[s].x[j][0] = g_rough[j]; // initialize with enclosure of the flow on the time step
    
    
    for(i=0;i<order;i++)
    {
        // Evaluate the i'th Taylor coefficient of the r.h.s. of the ODE:
        for (j=0 ; j<sysdim ; j++) {
            dde_x[s].xp[j].eval(i);
            dde_g[s].xp[j].eval(i);
            
            // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
            dde_x[s].x[j][i+1]=dde_x[s].xp[j][i]/double(i+1);
            dde_g[s].x[j][i+1]=dde_g[s].xp[j][i]/double(i+1);
            // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
        }
    }
}

//
void Dde_TM_Jac::build(DdeFunc bf, DdeJacFunc bbf, int s, double d0, double tau, int order)
{
    int i,j, k;
    vector<AAF> g_rough(sysdim);
    vector<vector<AAF>> Jac1_g_rough(sysdim, vector<AAF>(jacdim));
    
//    cout << "entering Dde_TM_Jac::build s=" << s << endl;
    
 //   ddeVAR_x[s].reset();
 //   ddeJAC_x[s].reset();
  //  bf(ddeVAR_x[s].xp,ddeVAR_x[s].x,prev_ddeVAR_x[s].x);
  
   //  beta[0].diff(0,sysdim);   // necessaire ou une fois pour toutes dans la fonction init ou autre ?
    
    for (j=0 ; j<sysdim ; j++)
        ddeVAR_x[s].x[j][0]=x[s][j]; // initialize;
    
    for (i=0 ; i<sysdim ; i++)
        for (j=0 ; j<jacdim ; j++)
            ddeJAC_x[s].J[i][j] = J[s][i][j];
    
    bf(ddeVAR_x[s].xp,ddeVAR_x[s].x,prev_ddeVAR_x[s].x,fullinputs);
    
    // ne sert plus a rien?
    bbf(ddeJAC_x[s].Jp,ddeJAC_x[s].J,prev_ddeJAC_x[s].J,ddeVAR_x[s].x,prev_ddeVAR_x[s].x);
    
  //    cout << "ddeVAR_x[s].x = " << ddeVAR_x[s].x[0][0].convert_int() << " prev_ddeVAR_x[s].x =" << prev_ddeVAR_x[s].x[0][0].convert_int() << endl;
  //   cout << "ddeVAR_x[s].J = " << ddeJAC_x[s].J[0][0][0].convert_int() << " prev_ddeVAR_x[s].J =" << prev_ddeJAC_x[s].J[0][0][0].convert_int() << endl;
  //   cout << "ddeVAR_x[s].J = " << ddeJAC_x[s].J[0][1][0].convert_int() << " prev_ddeVAR_x[s].J =" << prev_ddeJAC_x[s].J[0][1][0].convert_int() << endl;
  //   cout << "ddeVAR_x[s].J = " << ddeJAC_x[s].J[1][0][0].convert_int() << " prev_ddeVAR_x[s].J =" << prev_ddeJAC_x[s].J[1][0][0].convert_int() << endl;
  //   cout << "ddeVAR_x[s].J = " << ddeJAC_x[s].J[1][1][0].convert_int() << " prev_ddeVAR_x[s].J =" << prev_ddeJAC_x[s].J[1][1][0].convert_int() << endl;
    
    vector<AAF> x_prev_range(sysdim);
    for (j=0 ; j<sysdim ; j++)
        x_prev_range[j] = prev_ddeVAR_g[s].x[j][0]; // dde_g[s].x_prev[j][0];
    g_rough = fixpoint(bf,x[s],x_prev_range,fullinputs,tau);
    
    for (int i=0 ; i<sysdim ; i++)
        g_rough[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term
    
    // record DDE
   // ddeVAR_g[s].reset();
   // ddeJAC_g[s].reset();

    // ne sert plus a rien ?
    for (j=0 ; j<sysdim ; j++)
        ddeVAR_g[s].x[j][0] = g_rough[j]; // initialize with enclosure of the flow on the time step
     bf(ddeVAR_g[s].xp,ddeVAR_g[s].x,prev_ddeVAR_g[s].x,fullinputs);
    
    for (i=0 ; i<sysdim ; i++)
        for (j=0 ; j<jacdim ; j++)
            Jac1_g_rough[i][j] = prev_ddeJAC_g[s].J[i][j][0];
    vector<vector<AAF>> J_rough =  fixpoint(bbf, J[s], Jac1_g_rough , g_rough, x_prev_range , tau);
    
    for (i=0 ; i<sysdim ; i++)
        for (j=0 ; j<jacdim ; j++)
            ddeJAC_g[s].J[i][j] = J_rough[i][j];
    bbf(ddeJAC_g[s].Jp,ddeJAC_g[s].J,prev_ddeJAC_g[s].J,ddeVAR_g[s].x,prev_ddeVAR_g[s].x);
  //  bbf(ddeJAC_g[s].Jp,ddeJAC_g[s].J,prev_ddeJAC_g[s].J,ddeVAR_g[s].x,prev_ddeVAR_g[s].x);
    
    // specify the variables to differentiate
  //  for (int i=0 ; i<sysdim ; i++)
    
   //  beta[0].diff(0,sysdim);
  //  ddeVAR_x[s].x[0][0].diff(0,sysdim);
  //  ddeVAR_g[s].x[0][0].diff(0,sysdim);
    
    for(i=0;i<order;i++)
    {
        // Evaluate the i'th Taylor coefficient of the r.h.s. of the ODE:
        for (j=0 ; j<sysdim ; j++) {
            ddeVAR_x[s].xp[j].eval(i);
            ddeVAR_g[s].xp[j].eval(i);
        
            // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
            ddeVAR_x[s].x[j][i+1]=ddeVAR_x[s].xp[j][i]/double(i+1);
            ddeVAR_g[s].x[j][i+1]=ddeVAR_g[s].xp[j][i]/double(i+1);
            // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
            
         /*   if (i == 0) {
                cout << "Jaci[0][0][0] = "<< ddeVAR_x[s].x[j][i].d(0).convert_int();
                cout << "Jaci[1][0][0] = "<< ddeVAR_x[s].x[j][i+1].d(0).convert_int();
            }*/
        }
        
        for (j=0 ; j<sysdim ; j++) {
            for (k=0 ; k<jacdim ; k++) {
                ddeJAC_x[s].Jp[j][k].eval(i);
                ddeJAC_g[s].Jp[j][k].eval(i);
            }
        }
        
        for (j=0 ; j<sysdim ; j++) {
            for (k=0 ; k<jacdim ; k++) {
                // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
                ddeJAC_x[s].J[j][k][i+1]=ddeJAC_x[s].Jp[j][k][i]/double(i+1);
                ddeJAC_g[s].J[j][k][i+1]=ddeJAC_g[s].Jp[j][k][i]/double(i+1);
                // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
            }
        }
    }
  
//    cout << "exiting Dde_TM_Jac::build" << endl;
}

// building the Taylor model for the solution of ODE starting from the center of initial conditions
// outputs are ode_x0 (Taylor coefficients from 0 to order-1) and ode_g0 (giving the remainder Taylor coefficient,
// computed using an a priori enclosure of the solution on [tn,tn+tau])
// if active_discrete_trans is true, then we join with range of previous transition for the enclosure on the range (last coeff of TM)
// tn here is the final time on which this TM should be valid
void HybridStep_dde::TM_build(int s)
{
//    cout << "Building TM valid on [" <<tn+s*tau << "," << tn+(s+1)*tau << "]" << endl;
    if (innerapprox == 0)
        TMcenter.build(bf,fullinputs,s,d0,tau,order);
    else
        TMcenter.build(bf,center_fullinputs,s,d0,tau,order);

    if (innerapprox == 1)
        TMJac.build(bf,bbf,s,d0,tau,order);
    
}


// eval s-th Taylor model
void Dde_TM_val::eval(vector<AAF> &res, int s, double h, int order)
{
    double taui = h;
    
    for (int j=0 ; j<sysdim; j++)
        res[j] = (dde_x[s].x[j])[0];
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            res[j] += dde_x[s].x[j][i]*taui;
        taui *= h;
    }
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        res[j] += dde_g[s].x[j][order]*taui;
}


// eval s-th Taylor model
void Dde_TM_Jac::eval_val(vector<AAF> &res, int s, double h, int order)
{
    double taui = h;
    
    for (int j=0 ; j<sysdim; j++)
        res[j] = ddeVAR_x[s].x[j][0];
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            res[j] += ddeVAR_x[s].x[j][i]*taui;
        taui *= h;
    }
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        res[j] += ddeVAR_g[s].x[j][order]*taui;
 //    cout << "Value is x["<<s<<"][0][0]="<<res[0].convert_int() << endl;
}


// eval s-th Taylor model
void Dde_TM_Jac::eval_Jac(vector<vector<AAF>> &J_res, int s, double h, int order) // vector<AAF> &res, int s, double h, int order)
{
    double taui = h;
    vector<vector<AAF>> aux(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> Jaci(sysdim, vector<AAF>(jacdim));
    
    J_res = J[s];
    //for (int j=0 ; j<sysdim; j++)
    //    for (int k=0 ; k<sysdim; k++)
    //        J_res[j][k] = ddeVAR_x[s].x[j][0].d(k); // J[s];
//    cout << "Beginning of eval_jac at local time " << s*taui << endl;
//    cout << "J["<<s<<"][0][0]="<<J_res[0][0].convert_int() << endl;
    
    // pour les premiers termes, Jac^i(x)*J*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            for (int k=0 ; k<jacdim; k++)
                Jaci[j][k] = ddeJAC_x[s].J[j][k][i];

        
            aux = Jaci;
            scaleM(aux,taui); // scaleM(aux,taui);
            addMiMi(J_res,aux);
        
    
        taui *= h;
    }
    // pour le dernier terme, Jac^i(g_rough)*J_rough*tau^i
    for (int j=0 ; j<sysdim; j++)
        for (int k=0 ; k<jacdim; k++)
            Jaci[j][k] = ddeJAC_g[s].J[j][k][order];
  
    {
        aux = Jaci;
        scaleM(aux,taui);
        addMiMi(J_res,aux);
    }
  
    
 //   J_res = J_res + aux;
//    cout << "Final sum is J["<<s<<"][0][0]="<<J_res[0][0].convert_int() << endl << endl;
    
}

// A completer.
void Dde_TM_Jac::eval(int s, double h, int order) // vector<AAF> &res, int s, double h, int order)
{
    eval_val(xp1[s],s,h,order);
    eval_Jac(Jp1[s],s,h,order);
}

// eval s-th Taylor model and initialize s+1-th
void HybridStep_dde::TM_eval(int s)
{
    TMcenter.eval(TMcenter.xp1[s],s,tau,order);
    if (innerapprox == 1)
        TMJac.eval(s,tau,order);
}



void HybridStep_dde::print_solutionstep(int s, vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xouter_minimal, vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xinner_minimal, vector<interval> &Xcenter)
{
    
    cout << "print_solutionstep at t=" << tn+(s+1)*tau << ": " << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "Xouter_maximal[" << i <<"]=" << Xouter[i] << "\t";
    for (int i=0 ; i<sysdim ; i++)
        cout << "Xinner_maximal[" << i <<"]=" << Xinner[i] << "\t";
    if (uncontrolled > 0) {
        for (int i=0 ; i<sysdim ; i++)
            cout << "Xouter_robust[" << i <<"]=" << Xouter_robust[i] << "\t";
        for (int i=0 ; i<sysdim ; i++)
            cout << "Xinner_robust[" << i <<"]=" << Xinner_robust[i] << "\t";
    }
    if (controlled > 0 || uncontrolled > 0) {
        for (int i=0 ; i<sysdim ; i++)
            cout << "Xouter_minimal[" << i <<"]=" << Xouter_minimal[i] << "\t";
        for (int i=0 ; i<sysdim ; i++)
            cout << "Xinner_minimal[" << i <<"]=" << Xinner_minimal[i] << "\t";
    }
    
    //out_approx << YAML::BeginMap;
    //out_approx << YAML::Key << "tn";
    //out_approx << YAML::Value << tn+(s+1)*tau;
    
    vector<double> temp(2*sysdim);
    
    if (nb_subdiv_init ==1) {
        out_approx << YAML::Key << "outer";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xouter[i].inf();
            temp[2*i+1] = Xouter[i].sup();
        }
        out_approx << YAML::Value << temp; // Xouter does not work because of interval type (I guess I could solve this but...)
        out_approx << YAML::Key << "center";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xcenter[i].inf();
            temp[2*i+1] = Xcenter[i].sup();
        }
        out_approx << YAML::Value << temp; // Xcenter;
    }
    if (innerapprox == 1)
    {
        out_approx << YAML::Key << "inner";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xinner[i].inf();
            temp[2*i+1] = Xinner[i].sup();
        }
        out_approx << YAML::Value << temp; // Xinner;
        if (uncontrolled > 0) {
        out_approx << YAML::Key << "outerrobust";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xouter_robust[i].inf();
            temp[2*i+1] = Xouter_robust[i].sup();
        }
        out_approx << YAML::Value << temp; // Xouter_robust;
        out_approx << YAML::Key << "innerrobust";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xinner_robust[i].inf();
            temp[2*i+1] = Xinner_robust[i].sup();
        }
        out_approx << YAML::Value << temp; // Xinner_robust;
    }
    if (controlled > 0 || uncontrolled > 0) {
        out_approx << YAML::Key << "outerminimal";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xouter_minimal[i].inf();
            temp[2*i+1] = Xouter_minimal[i].sup();
        }
        out_approx << YAML::Value << temp; // Xouter_minimal;
        out_approx << YAML::Key << "innerminimal";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xinner_minimal[i].inf();
            temp[2*i+1] = Xinner_minimal[i].sup();
        }
        out_approx << YAML::Value << temp; // Xinner_minimal;
    }
    }
   // out_approx << YAML::EndMap;
    
    for (int i=0 ; i<sysdim ; i++) {
        // saving result
        Xouter_print[current_subdiv][current_iteration][i] = Xouter[i];
        Xouter_robust_print[current_subdiv][current_iteration][i] = Xouter_robust[i];
        Xouter_minimal_print[current_subdiv][current_iteration][i] = Xouter_minimal[i];
        Xinner_print[current_subdiv][current_iteration][i] = Xinner[i];
        Xinner_robust_print[current_subdiv][current_iteration][i] = Xinner_robust[i];
        Xinner_minimal_print[current_subdiv][current_iteration][i] = Xinner_minimal[i];
        t_print[current_iteration] = tn+(s+1)*tau;
    }
     current_iteration++;
}


// print solution at time cur_step.tn+(j+1)*tau/nb_subdiv
//  A completer pour la sous-approx
ReachSet HybridStep_dde::TM_evalandprint_solutionstep(int s, vector<interval> &eps)
{
    vector<interval> Xouter(sysdim),Xouter_robust(sysdim),Xouter_minimal(sysdim),Xinner(sysdim),Xinner_robust(sysdim),Xinner_minimal(sysdim),Xcenter(sysdim);
    ReachSet res;
    
    TM_eval(s);
    out_approx << YAML::BeginMap;
    
    out_approx << YAML::Key << "tn";
    out_approx << YAML::Value << tn+(s+1)*tau;
    
    if (innerapprox == 0)
    {
        for (int i = 0 ; i<sysdim ; i++) {
            TMcenter.xp1[s][i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
            TMcenter.xp1[s][i].sumup(tol_noise); // group small terms
            Xouter[i] = TMcenter.xp1[s][i].convert_int();
        }
        
        print_solutionstep(s,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter);
    }
    else // A modifier
    {
        for (int i = 0 ; i<sysdim ; i++) {
            TMcenter.xp1[s][i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
            TMcenter.xp1[s][i].sumup(tol_noise); // group small terms
            Xcenter[i] = TMcenter.xp1[s][i].convert_int();
        }
       
        InnerOuter(Xinner,Xinner_robust,Xinner_minimal,Xouter,Xouter_robust,Xouter_minimal,TMcenter.xp1[s],TMJac.Jp1[s],eps,tn+tau); //x0p1,Jp1,eps);
        cout << "without quadrature: ";
        cout << "Xouter=" << Xouter;
        cout << "Xinner=" << Xinner;
     //   print_solutionstep(s,Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter);
        
        InnerOuter_discretize(Xinner,Xinner_robust,Xinner_minimal,Xouter,Xouter_robust,Xouter_minimal,TMcenter.xp1[s],TMJac.Jp1[s],eps,tn+tau);
        cout << "with quadrature: ";
        cout << "Xouter=" << Xouter;
        cout << "Xinner=" << Xinner;
       //  print_solutionstep(s,Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter);
        
        
     //   cout << "Xouter[0]" << Xouter[0] << endl;
        intersectViVi(Xouter,TMJac.xp1[s]);
       // for (int i = 0 ; i<sysdim ; i++)
       //     Xouter[i] = TMJac.xp1[s][i].convert_int();
        cout << "with intersection with direct solution: ";
        print_solutionstep(s,Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter);
    }
    out_approx << YAML::EndMap;
    
    vector<interval> Xsampled(sysdim);
    for (int i=0; i<sysdim; i++)
        Xsampled[i] = interval::EMPTY();
    res = ReachSet(Xsampled,Xouter,Xinner);
    return res;
//    print_solutionstep(Xouter,Xinner,Xcenter);
}

void Dde_TM_val::init_nextsmallstep(int s)
{
    x[s+1] = xp1[s];
}


void Dde_TM_Jac::init_nextsmallstep(int s)
{
    x[s+1]= xp1[s];
    J[s+1] = Jp1[s];
 //   for (int i=0 ; i<sysdim ; i++)
  //      for (int j=0 ; j<sysdim ; j++)
   //         J[s+1][i][j] = Jp1[s][i][j];
}



void HybridStep_dde::init_nextsmallstep(int s)
{
    TMcenter.init_nextsmallstep(s);
    if (innerapprox == 1)
        TMJac.init_nextsmallstep(s);
   // TMcenter.x[s+1] = TMcenter.xp1[s];
}


//
HybridStep_dde HybridStep_dde::init_nextbigstep(double tau)
{
    HybridStep_dde res(bf,bbf,TMcenter.dde_x,TMcenter.dde_g,TMJac.ddeVAR_x,TMJac.ddeVAR_g,TMJac.ddeJAC_x, TMJac.ddeJAC_g,order,tn+d0,tau,d0,p);
    
    for (int i=0 ; i<sysdim ; i++) {
        TMcenter.xp1[p-1][i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
        TMcenter.xp1[p-1][i].sumup(tol_noise); // group small terms
    }
    res.TMcenter.x[0] = TMcenter.xp1[p-1];
    
    //
    if (innerapprox == 1)
    {
        for (int i=0 ; i<sysdim ; i++) {
            TMJac.xp1[p-1][i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
            TMJac.xp1[p-1][i].sumup(tol_noise); // group small terms
        }
        for (int i=0 ; i<sysdim ; i++)
            for (int k=0 ; k<jacdim ; k++) {
                TMJac.Jp1[p-1][i][k].compact();
                TMJac.Jp1[p-1][i][k].sumup(tol_noise);
            }
        res.TMJac.x[0] = TMJac.xp1[p-1]; // matrix assignment (ok or should copy each coeff ?)
        
         res.TMJac.J[0] = TMJac.Jp1[p-1];
       // for (int i=0 ; i<sysdim ; i++)
       //     for (int j=0 ; j<sysdim ; j++)
       // res.TMJac.J[0][i][j] = TMJac.Jp1[p-1][i][j];
        
//        cout << "res.TMJac.JddeVAR_x[s].J = " << res.TMJac.J[0][0][0].convert_int()  << endl;
//        cout << "res.TMJac.JddeVAR_x[s].J = " << res.TMJac.J[0][0][1].convert_int()  << endl;
//        cout << "res.TMJac.JddeVAR_x[s].J = " << res.TMJac.J[0][1][0].convert_int()  << endl;
//        cout << "res.TMJac.JddeVAR_x[s].J = " << res.TMJac.J[0][1][1].convert_int() << endl;
        
    }
    return res;
}




