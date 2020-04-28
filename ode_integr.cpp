/* ============================================================================
 File   : ode_integr.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file builds the Taylor models and defines the main reachability analysis functions for ODEs
 - compute Taylor model coefficients by automatic differentiation from the ODE systems
 - compute a priori rough enclosures of the flow and its Jacobian (fixpoint functions)
 - compute refined tubes of the flow and its jacobian using the rough enclosure and the Taylor models
 ============================================================================ */

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"
#include "matrix.h"
#include "inner.h"
#include "ode_integr.h"
#include "ode_def.h"


#include <iostream>
#include <ostream>
#include <fstream>
#include <ctime>
#include <assert.h>
using namespace std;





// printing the Taylor coefficients and Jacobians of orders between min and max
void OdeVar::print(int min, int max) {
    for(int i=min;i<=max;i++) {
        for (int j=0 ; j<sysdim ; j++)
            cout << "odeVAR.x[" << j << "][" << i << "]=" << x[j][i].x() << "\t";
        cout << endl;
        for (int j=0 ; j<sysdim ; j++) {
            for (int k=0 ; k<sysdim ; k++)
                cout << "odeVAR.Jacx[" << j << "]x[" << k << "][" << i << "]=" << x[j][i].d(k) << "\t";
            cout << endl;
            cout << endl;
        }
        cout << endl;
    }
}

void OdeVar::printTM(int order) {
    //    cout << "printTM" << endl;
    for (int j=0 ; j<sysdim ; j++) {
        cout << "x[" << j << "]=" << x[j][0].x().convert_int() << " + ";
        
        for (int i=1;i<order;i++) {
            cout << x[j][i].x().convert_int()  << " t^" << i << " + " ;
        }
        cout << x[j][order].x().convert_int()  << " t^" << order << endl ;
    }
}

void OdeVar::printAAFTM(int order) {
    unsigned int new_index;
    //    cout << "printTM" << endl;
    for (int j=0 ; j<sysdim ; j++) {
        new_index = (x[j][0].x()).sumup((unsigned)sysdim); // attention modifie x ...
        cout << "x[" << j << "]=" << x[j][0].x() << " + ";
        
        for (int i=1;i<order;i++) {
            new_index = (x[j][i].x()).sumup((unsigned)sysdim); // attention modifie x ...
            cout << x[j][i].x() << " t^" << i << " + " ;
        }
        new_index = (x[j][order].x()).sumup((unsigned)sysdim); // attention modifie x ...
        cout << x[j][order].x() << " t^" << order << endl ;
    }
}

// printing the Taylor coefficients  of orders between min and max
void Ode::print(int min, int max) {
    for(int i=min;i<=max;i++) {
        for (int j=0 ; j<sysdim ; j++)
            cout << "ode.x[" << j << "][" << i << "]=" << x[j][i] << "\t" ;
        cout << endl;
    }
}

void Ode::printTM(int order) {
    //  cout << "printTM" << endl;
    for (int j=0 ; j<sysdim ; j++) {
        cout << "x[" << j << "]=" << x[j][0].convert_int() << " + ";
        
        for (int i=1;i<order;i++) {
            cout << x[j][i].convert_int() << " t^" << i << " + " ;
        }
        cout << x[j][order].convert_int() << " t^" << order << endl ;
    }
}


// building the Taylor model for the solution of ODE starting from the center of initial conditions
// outputs are ode_x0 (Taylor coefficients from 0 to order-1) and ode_g0 (giving the remainder Taylor coefficient,
// computed using an a priori enclosure of the solution on [tn,tn+tau])
// if active_discrete_trans is true, then we join with range of previous transition for the enclosure on the range (last coeff of TM)
// tn here is the final time on which this TM should be valid
void TM_val::build(OdeFunc _bf, vector<AAF> &param_inputs) {
    vector<AAF> g_rough(sysdim); // rough estimation of solution of current mode on [tn,tn+tau]
    int i, j;
    double offset = (tn-t_begin)/t_end;
    
    // compute Taylor coefficients / Lie derivatives starting from center of initial conditions
    ode_x.reset();
    for (j=0 ; j<sysdim ; j++)
        ode_x.x[j][0]=x[j]; // initialize with center;
    for (j=0 ; j<inputsdim ; j++)
        ode_x.param_inputs[j][0]=param_inputs[index_param_inv[j]+floor(offset*nb_inputs[j])];
    
    // computing a priori enclosure on [tn,tn+tau] of solution starting from center of initial conditions
    g_rough = fixpoint(_bf,param_inputs,x,tau);
    
    
        for (int i=0 ; i<sysdim ; i++)
            g_rough[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term
    
    
    // evaluate the Lie derivatives on this enclosure: we will use the last coefficients as remainder terms of Taylor model
    ode_g.reset();
    for (j=0 ; j<sysdim ; j++)
        ode_g.x[j][0] = g_rough[j]; // initialize with enclosure of the flow on the time step
    for (j=0 ; j<inputsdim ; j++)
        ode_g.param_inputs[j][0]=param_inputs[index_param_inv[j]+floor(offset*nb_inputs[j])];
   
    
    for(i=0;i<order;i++)
    {
        // Evaluate the i'th Taylor coefficient of the r.h.s. of the ODE:
        for (j=0 ; j<sysdim ; j++) {
            ode_x.xp[j].eval(i);
            ode_g.xp[j].eval(i);
            
            // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
            ode_x.x[j][i+1]=ode_x.xp[j][i]/double(i+1);
            ode_g.x[j][i+1]=ode_g.xp[j][i]/double(i+1);
            // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
        }
        for (j=0 ; j<inputsdim ; j++) {
            ode_x.param_inputs[j][i+1]=0;
            ode_g.param_inputs[j][i+1]=0;
        }
    }
}


// eval TM at time tn+h and return value
void TM_val::eval(vector<AAF> &res, double h)
{
    double taui = h;
    
    for (int j=0 ; j<sysdim; j++)
        res[j] = (ode_x.x[j])[0];
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            res[j] += ode_x.x[j][i]*taui;
        taui *= h;
    }
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        res[j] += ode_g.x[j][order]*taui;
    
 // A REVOIR
/*    for (int j=0 ; j<sysdim; j++)
    {
        if (is_variable[j]) // variable param with value always in the same range
        res[j] = ode_x.x[j][0];
    } */
}

// eval TM at tn+tau and store in xp1 and J
void TM_val::eval()
{
    eval(xp1,tau);
}



void TM_val::init_nextstep(double _tau)
{
    tn = tn + tau;
    tau = _tau;
    
    for (int i=0 ; i<sysdim ; i++) {
        xp1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
        xp1[i].sumup(tol_noise); // group small terms
   
        // A REVOIR
//        if (is_variable[i]) // removing dependencies to previous occurence in time dependent input
//            xp1[i] = xp1[i].convert_int();
    }
    x = xp1;
}






// TM of the Jacobian with respect to IC of the solution of the ODE for all the range of IC
// outputs are odeVAR_x (Taylor coefficients from 0 to order-1), odeVAR_g (for the remainder term,
// evaluated using a priori enclosure on [tn,tn+tau]), and J_rough, also used for the remainder term
void TM_Jac::build(OdeFunc _bf, vector<AAF> &param_inputs, vector<AAF> &param_inputs_center) {
    vector<AAF> g_rough(sysdim);
    vector<vector<AAF>> Jac1_g_rough(jacdim, vector<AAF>(jacdim)); // \partial f_i / \partial (z_j,beta_k) => only the sysdim first lines are relevant
    // vector<T<F<AAF>>> tf_param_inputs(jacdim-sysdim);
    int i, j, k;
    double offset = (tn-t_begin)/t_end;
    
    // compute Lie derivatives and Jacobian on CI: these will give coefficients 0 to order-1 of TM for Jacobian
    for (k=0 ; k<sysdim+inputsdim; k++)
    {
        odeVAR_x[k].reset();
        for (j=0 ; j<sysdim ; j++) {
            if ((k < j) && (refined_mean_value))
                odeVAR_x[k].x[j][0] = x0[j];
            else{
                odeVAR_x[k].x[j][0] = x[j];
     //           cout << "k=" << k << " j=" << j <<" odeVAR_x[k].x[j][0] =" << x[j].convert_int() << endl;
            }
        }
        
        for (j=0 ; j<inputsdim ; j++) {
            if ((k < j+sysdim) && (refined_mean_value))
                odeVAR_x[k].param_inputs[j][0]=param_inputs_center[index_param_inv[j]+floor(offset*nb_inputs[j])];
            else {
                odeVAR_x[k].param_inputs[j][0]=param_inputs[index_param_inv[j]+floor(offset*nb_inputs[j])];
          //      cout << "k=" << k << " j=" << j <<" odeVAR_x[k].param_inputs[j][0] =" << param_inputs[index_param_inv[j]+floor(offset*nb_inputs[j])].convert_int() << endl;
            }
            
           // cout << "index of param_inputs=" << index_param_inv[0]+floor(offset*nb_inputs[0]) << endl;
           // cout << "param inputs=" << param_inputs[index_param_inv[0]+floor(offset*nb_inputs[0])] << endl;
        }
    }
    
    // compute an a priori enclosure on [tn,tn+tau] of solution
    g_rough = fixpoint(_bf,param_inputs,x,tau);
    
    for (i=0 ; i<sysdim ; i++) {
        g_rough[i].compact();
        g_rough[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term (for efficiency)
    }
    
    // evaluate the Lie derivatives and Jacobian on this a priori enclosure on [tn,tn+tau]: we will use the last coefficients as remainder terms of Taylor model
    odeVAR_g.reset();
    for (j=0 ; j<sysdim ; j++)
        odeVAR_g.x[j][0] = g_rough[j];
    for (j=0 ; j<inputsdim ; j++)
        odeVAR_g.param_inputs[j][0]=param_inputs[index_param_inv[j]+floor(offset*nb_inputs[j])];
    
    // specify the variables to differentiate - we can probably differentiate only on the component we will be interested in - see later !
    for (k=0 ; k<sysdim; k++)
        odeVAR_x[k].x[k][0].diff(0,1);
    for (k=0 ; k<inputsdim ; k++)
        odeVAR_x[k+sysdim].param_inputs[k][0].diff(0,1); // index_param: correspondance for variable parameter
   //     odeVAR_x[k+sysdim].param_inputs[index_param[k]][0].diff(0,1); // index_param: correspondance for variable parameter
    
    
    for (j=0 ; j<sysdim ; j++)
        odeVAR_g.x[j][0].diff(j,sysdim+inputsdim);
    for (j=0 ; j<inputsdim ; j++)
        odeVAR_g.param_inputs[j][0].diff(sysdim+j,sysdim+inputsdim);
    
    
    
    for(i=0;i<order;i++)
    {
        // Evaluate the i'th Taylor coefficient of the r.h.s. of the ODE:
        for (k=0 ; k<sysdim+inputsdim; k++)
            for (j=0 ; j<sysdim ; j++) {
                odeVAR_x[k].xp[j].eval(i);
                // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
                odeVAR_x[k].x[j][i+1]=odeVAR_x[k].xp[j][i]/double(i+1);
                // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
            }
        for (j=0 ; j<sysdim ; j++) {
            odeVAR_g.xp[j].eval(i);
            // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
            odeVAR_g.x[j][i+1]=odeVAR_g.xp[j][i]/double(i+1);
            // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
        }
        // zero derivative: important, does not work otherwise...
        for (k=0 ; k<sysdim+inputsdim; k++)
            for (j=0 ; j<inputsdim ; j++)
                odeVAR_x[k].param_inputs[j][i+1]=0;
        for (j=0 ; j<inputsdim ; j++)
            odeVAR_g.param_inputs[j][i+1]=0;
        
    }
    
    
//    for (int k=0 ; k<inputsdim ; k++) {
//        cout << "k+sysdim+floor(offset*nb_inputs[j])" << k+sysdim+floor(offset*nb_inputs[j]) << endl;
 //       Jaci[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])] = odeVAR_g.x[j][order].d(k);
    
   
    // first order coefficient Jac1_g_rough of the  Taylor model for the Jacobian, evaluated on a priori enclosure of the flow (for remainder term)
    for (j=0 ; j<sysdim ; j++) {
        for (int k=0 ; k<sysdim ; k++)
        {
            Jac1_g_rough[j][k] = odeVAR_g.x[j][1].d(k);
            Jac1_g_rough[j][k].compact();
            Jac1_g_rough[j][k].sumup(tol_noise);
       //     cout << "Jac1_g_rough[i][k] in build: " << Jac1_g_rough[j][k].convert_int() << endl;
       //     cout << "Jac1_g_rough[i][k] in build: " << odeVAR_x.x[j][1].d(k).convert_int() << endl;
        }
        for (int k=0 ; k<inputsdim ; k++)
        {
            Jac1_g_rough[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])] = odeVAR_g.x[j][1].d(k);
            Jac1_g_rough[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])].compact();
            Jac1_g_rough[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])].sumup(tol_noise);
            //     cout << "Jac1_g_rough[i][k] in build: " << Jac1_g_rough[j][k].convert_int() << endl;
            //     cout << "Jac1_g_rough[i][k] in build: " << odeVAR_x.x[j][1].d(k).convert_int() << endl;
        }
    }
    
    // deduce a priori enclosure J_rough of the Jacobian (used to compute the remainder)
    fixpoint(J_rough, Jac1_g_rough,J,tau);
    
  
        for (i=0 ; i<sysdim ; i++)
            for (int k=0 ; k<jacdim ; k++)
            {
                J_rough[i][k].compact();
                J_rough[i][k].sumup(tol_noise); // group terms smaller than eps*radius in a new error term (for efficiency)
         //       cout << "J_rough[i][k] in build: " << J_rough[i][k].convert_int() << endl;
            }
}



// eval TM at time tn+h and return value in res
// the full value is in odeVar_x[jacdim-1], the others are evaluated partially in the center (improved mean-value theorem)
void TM_Jac::eval_val(vector<AAF> &res, double h)
{
    double taui = h;
    
    // only the last differentiation is computed with full intervals (improved mean value)
    for (int j=0 ; j<sysdim; j++)
        res[j] = odeVAR_x[sysdim+inputsdim-1].x[j][0].x();
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            res[j] += odeVAR_x[sysdim+inputsdim-1].x[j][i].x()*taui;
        taui *= h;
    }
    
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        res[j] += odeVAR_g.x[j][order].x()*taui;
    
    // A REVOIR
  /*  for (int j=0 ; j<sysdim; j++)
    {
    if (is_variable[j]) // variable param with value always in the same range
        res[j] = odeVAR_x[sysdim+inputsdim-1].x[j][0].x();
    } */
 //   for (int j=0 ; j<sysdim; j++)
  //      cout << "res[j]=" << res[j].convert_int() << endl;
}




// eval TM at tn+h and return in J_res
void TM_Jac::eval_Jac(vector<vector<AAF>> &J_res, double h)
{
    double taui = h;
    double offset = (tn-t_begin)/t_end;
    vector<vector<AAF>> aux(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> Jaci(sysdim, vector<AAF>(jacdim));
    
    J_res = J;
    
    // pour les premiers termes, Jac^i(x)*J*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++) {
            for (int k=0 ; k<sysdim; k++) {
                Jaci[j][k] = odeVAR_x[k].x[j][i].d(0);  //  on each structure differentiate to only one variable
       //         cout << "order=" << i << " Jaci["<< j <<"]["<< k << "]=" << Jaci[j][k] << endl;
            }
            for (int k=0 ; k<inputsdim ; k++) {
                Jaci[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])] = odeVAR_x[sysdim+k].x[j][i].d(0);
     //           cout << "order=" << i << " Jaci["<< j <<"]["<< sysdim+index_param_inv[k]+floor(offset*nb_inputs[k]) << "]=" << Jaci[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])] << endl;
            }
        }
        
        multJacfzJaczz0(aux,Jaci,J); // Jac^i(x)*J
        scaleJacfz(aux,taui);
        addJacfzJacfz(J_res,aux);
        taui *= h;
    }
    
    // pour le dernier terme, Jac^i(g_rough)*J_rough*tau^i
    for (int j=0 ; j<sysdim; j++) {
        for (int k=0 ; k<sysdim; k++) {
            Jaci[j][k] = odeVAR_g.x[j][order].d(k);
        }
        for (int k=0 ; k<inputsdim ; k++) {
     //       cout << "k+sysdim+floor(offset*nb_inputs[k])" << k+sysdim+floor(offset*nb_inputs[k]) << endl;
            Jaci[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])] = odeVAR_g.x[j][order].d(k);
        }
    }
    
    //  print_interv("Jaci",Jaci);
    multJacfzJaczz0(aux,Jaci,J_rough); // Jac^i(g_rough)*J_rough
    
    //    print_interv("J_rough",J_rough);
    scaleJacfz(aux,taui);
    addJacfzJacfz(J_res,aux);
     /*   for (int j=0 ; j<sysdim; j++)
            for (int k=0 ; k<jacdim; k++)
                cout << "J_res[j][k]=" << J_res[j][k].convert_int() << endl; */
}



// eval TM at tn+tau and store in xp1 and J
void TM_Jac::eval()
{
    eval_val(xp1,tau);
    eval_Jac(Jp1,tau);
}

// get _x0 computed by TM_Val
void TM_Jac::init_nextstep(double _tau, vector<AAF> &_x0)
{
    tn = tn + tau;
    tau = _tau;
    
    x0 = _x0;
    
    for (int i=0 ; i<sysdim ; i++) {
        xp1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
        xp1[i].sumup(tol_noise); // group small terms
        
        // A REVOIR
    //    if (is_variable[i]) // removing dependencies to previous occurence in time dependent input
     //       xp1[i] = xp1[i].convert_int();
    }
    x = xp1;
    
    for (int i=0 ; i<sysdim ; i++)
        for (int k=0 ; k<jacdim ; k++) {
            Jp1[i][k].compact();
            Jp1[i][k].sumup(tol_noise);
            
            // AJOUTER QQCHOSE POUR ENLEVER DEPENDENCES AUSSI ?
        }
    
    J = Jp1; // J = Jp1
}



vector<AAF> fixpoint(OdeFunc bf, vector<AAF> &param_inputs, vector<AAF> &x0, double tau)
{
    // calcul de x satisfaisant x0 + [0,tau][f](x) \subseteq x
    vector<AAF> y0(sysdim);
    vector<AAF> y1(sysdim);
    vector<AAF> fx0(sysdim);
    int iter;
    interval widen, coeff;
    
    bf(fx0,param_inputs,x0);
    
    for (int i=0; i<sysdim ; i++)
        y1[i] = x0[i] + interval(0,tau)*fx0[i].convert_int();  // modif (*tau)
    
    
    
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
        
      //  cout << "iter=" << iter << "fx0[1]=" << fx0[1].convert_int() << "\t" "y0[1]=" << y0[1].convert_int() << endl;
        
        bf(fx0,param_inputs,y0);
        // cout << "fx[0]=" << fx0[0] << "\t"  << fx0[1] << endl;
        //   cout << "y1=" << y1[0] << "\t"  << y1[1] << endl;
        for (int i=0; i<sysdim ; i++)
            y1[i] = x0[i] + interval(0,tau)*fx0[i].convert_int();
        
        
       
        iter = iter+1;
    }
    //cout << "y0 (end of fixpoint)=" << y0[0] << "\t"  << y0[1] << " iter=" << iter << endl;
    return y0;
}


// computing a priori enclosure y0 for Jacobian, such that J0 + [0,tau]*Jac1_g_rough*y0 \subseteq y0
void fixpoint(vector<vector<AAF>> &y0, vector<vector<AAF>> &Jac1_g_rough, vector<vector<AAF>> &J0, double tau)
{
    // calcul de x satisfaisant J0 + [0,tau]Jac1(g)*x \subseteq x
    vector<vector<AAF>> J1(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> fJ0(sysdim, vector<AAF>(jacdim));
    int iter;
    interval widen, coeff;
    
    multJacfzJaczz0(fJ0,Jac1_g_rough,J0);
   
    
    for (int i=0; i<sysdim ; i++)
        for (int j=0; j<jacdim ; j++)
            J1[i][j] = J0[i][j] + interval(0,tau)*fJ0[i][j].convert_int();
    
    y0 = J0;
    
    iter = 1;
    widen = interval(-1,1);
    
    while (iter <=1 || !subseteq(J1,y0)) {
        
        //     cout << "J1=";
        //     print(J1);
        //    cout << "y0=";
        //    print(y0);
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
                for (int j=0; j<jacdim ; j++) {
                    J1[i][j] = J1[i][j] + coeff*widen*J1[i][j].convert_int();
                }
            //   print_interv("fixpoint J1",J1);
        }
        y0 = J1;
    
        
        multJacfzJaczz0(fJ0,Jac1_g_rough,y0); // fx0 = f(y0)
        
        for (int i=0; i<sysdim ; i++)
            for (int j=0; j<jacdim ; j++) {
                J1[i][j] = J0[i][j] + interval(0,tau)*fJ0[i][j].convert_int();
                J1[i][j].compact();
                J1[i][j].sumup(tol_noise);
            }
        
        iter = iter+1;
    }
  //  cout << "nb_iter = " << iter << endl;

}







HybridStep_ode init_ode(OdeFunc bf, vector<AAF> &x0, vector<AAF> &x, vector<vector<AAF>> &J0, double tn, double tau, int order)
{
    vector<OdeVar> odeVAR_x(sysdim+inputsdim);
    OdeVar odeVAR_g = OdeVar(bf);
    Ode ode_x0 = Ode(bf);
    Ode ode_g0 = Ode(bf);
    vector<AAF> init_x;
    
    if (innerapprox == 1)
       init_x = x0;
    else
        init_x = x;
    
    
    
    TM_val TMcenter = TM_val(ode_x0, ode_g0, order, init_x, tn, tau);
    //if (innerapprox == 1)
    
    if (innerapprox == 1) {
        for (int k = 0; k<sysdim+inputsdim; k++)
            odeVAR_x[k] = OdeVar(bf);
    }
    TM_Jac TMJac = TM_Jac(odeVAR_x, odeVAR_g, order, x, init_x, J0, tn, tau);
    
    HybridStep_ode res = HybridStep_ode(bf,TMcenter,TMJac,tn,tau,order);
    return res;
}




void HybridStep_ode::init_nextstep(double _tau)
{
    TMcenter.init_nextstep(_tau);
    // pass the center of TMcenter to TMJac
    if (innerapprox == 1)
        TMJac.init_nextstep(_tau,TMcenter.x);
    
    tn = tn + tau;
    tau = _tau;
}


void HybridStep_ode::TM_build(vector<AAF> &param_inputs,vector<AAF> &param_inputs_center)
{
    TMcenter.build(bf,param_inputs_center);
    if (innerapprox == 1)
        TMJac.build(bf,param_inputs,param_inputs_center);
}

// eval s-th Taylor model and initialize s+1-th
void HybridStep_ode::TM_eval()
{
    TMcenter.eval();
    if (innerapprox == 1)
        TMJac.eval();
}




void HybridStep_ode::print_solutionstep(vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xouter_minimal, vector<interval> &Xinner,  vector<interval> &Xinner_robust, vector<interval> &Xinner_minimal, vector<interval> &Xcenter)
{
    double aux, minwidth_ratio = 1.0;  // min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
    double mean_dist; // mean value on the xi of max distance between inner and outer approximations
    double tnp1 = tn + tau;
    
    
    if (print_debug)
    {
        cout << "t=" << tnp1 << endl;
        
        for (int i=0 ; i<sysdim ; i++) {
            cout << "Xouter_maximal[" << i <<"]=" << Xouter[i] << "\t";
            cout << "Xinner_maximal[" << i <<"]=" << Xinner[i] << "\t";
        }
        cout << endl;
        if (uncontrolled > 0)
        {
            for (int i=0 ; i<sysdim ; i++) {
                cout << "Xouter_robust[" << i <<"]=" << Xouter_robust[i] << "\t";
                cout << "Xinner_robust[" << i <<"]=" << Xinner_robust[i] << "\t";
            }
            cout << endl;
        }
        if (controlled > 0 || uncontrolled > 0)
        {
            for (int i=0 ; i<sysdim ; i++) {
                cout << "Xouter_minimal[" << i <<"]=" << Xouter_minimal[i] << "\t";
                cout << "Xinner_minimal[" << i <<"]=" << Xinner_minimal[i] << "\t";
            }
            cout << endl;
        }
    }
  
    
    for (int i=0 ; i<sysdim ; i++) {
        if (print_debug)
        {
            if (nb_subdiv_init ==1)
            {
                outFile_outer[i] << tnp1 << "\t" << inf(Xouter[i]) << "\t" << sup(Xouter[i]) << endl;
                outFile_center[i] << tnp1 << "\t" << inf(Xcenter[i]) << "\t" << sup(Xcenter[i]) << endl;
            }
            if (uncontrolled > 0) {
                outFile_outer_robust[i] << tnp1 << "\t" << inf(Xouter_robust[i]) << "\t" << sup(Xouter_robust[i]) << endl;
                outFile_inner_robust[i] << tnp1 << "\t" << inf(Xinner_robust[i]) << "\t" << sup(Xinner_robust[i]) << endl;
            }
            if (controlled > 0 || uncontrolled > 0) {
                outFile_outer_minimal[i] << tnp1 << "\t" << inf(Xouter_minimal[i]) << "\t" << sup(Xouter_minimal[i]) << endl;
                outFile_inner_minimal[i] << tnp1 << "\t" << inf(Xinner_minimal[i]) << "\t" << sup(Xinner_minimal[i]) << endl;
            }
            
            outFile_inner[i] << tnp1 << "\t" << inf(Xinner[i]) << "\t" << sup(Xinner[i]) << endl;
            //    outFile_inner_joint[i] << tnp1 << "\t" << inf(Xinner_joint[i]) << "\t" << sup(Xinner_joint[i]) << endl;
        }
        
        // saving result
        Xouter_print[current_subdiv][current_iteration][i] = Xouter[i];
        Xouter_robust_print[current_subdiv][current_iteration][i] = Xouter_robust[i];
        Xouter_minimal_print[current_subdiv][current_iteration][i] = Xouter_minimal[i];
        Xinner_print[current_subdiv][current_iteration][i] = Xinner[i];
        Xinner_robust_print[current_subdiv][current_iteration][i] = Xinner_robust[i];
        Xinner_minimal_print[current_subdiv][current_iteration][i] = Xinner_robust[i];
        t_print[current_iteration] = tnp1;
    }
    current_iteration++;
    
    if (print_debug)
    {
        minwidth_ratio = (sup(Xinner[0])-inf(Xinner[0]))/(sup(Xouter[0])-inf(Xouter[0]));
        for (int i=1 ; i<sysdim ; i++) {
            aux = (sup(Xinner[i])-inf(Xinner[i]))/(sup(Xouter[i])-inf(Xouter[i]));
            if (minwidth_ratio > aux)
                minwidth_ratio = aux;
        }
        if (tnp1 != 0)
            outFile_width_ratio << tnp1 << "\t" << minwidth_ratio << endl;
    }
}



void HybridStep_ode::TM_evalandprint_solutionstep(vector<interval> &eps, double tnp1)
{
    assert (tn <= tnp1);
    assert (tnp1 <= tn+tau);
    
    vector<interval> Xouter(sysdim), Xouter_robust(sysdim), Xouter_minimal(sysdim), Xinner(sysdim), Xinner_robust(sysdim), Xinner_minimal(sysdim), Xcenter(sysdim);
    
    // eval and store at time tnp1
    TM_eval();
    
    if (innerapprox == 0)
    {
        for (int i = 0 ; i<sysdim ; i++) {
            TMcenter.xp1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
            TMcenter.xp1[i].sumup(tol_noise); // group small terms
            Xouter[i] = TMcenter.xp1[i].convert_int();
        }
        print_solutionstep(Xouter,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter);
    }
    else {
    // deduce inner-approx by mean-value thm
        for (int i = 0 ; i<sysdim ; i++) {
            TMcenter.xp1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
            TMcenter.xp1[i].sumup(tol_noise); // group small terms
            Xcenter[i] = TMcenter.xp1[i].convert_int();
        }
        InnerOuter(Xinner,Xinner_robust,Xinner_minimal,Xouter,Xouter_robust,Xouter_minimal,TMcenter.xp1,TMJac.Jp1,eps,tn+tau);
        cout << "without quadrature: ";
        cout << "Xouter=" << Xouter;
        cout << "Xinner=" << Xinner;
        
    /*    InnerOuter_discretize(Xinner,Xinner_robust,Xinner_minimal,Xouter,Xouter_robust,Xouter_minimal,TMcenter.xp1,TMJac.Jp1,eps,tn+tau);
        cout << "with quadrature: ";
        cout << "Xouter=" << Xouter;
        cout << "Xinner=" << Xinner; */
       
      //  InnerOuter(Xinner,Xouter,TMcenter.xp1,TMJac.Jp1,eps);
     //   for (int i = 0 ; i<sysdim ; i++)
     //   cout << "before intersect: Xouter_maximal[" << i <<"]=" << Xouter[i] << "\t";
        
        
        intersectViVi(Xouter,TMJac.xp1);
        
        cout << "with intersection with direct solution: ";
        print_solutionstep(Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter);
        
    // print_solutionstep_ode(Xouter,Xinner,Xcenter,tnp1);
    }
}





