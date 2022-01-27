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
#include <float.h>
using namespace std;

#include "yaml-cpp/yaml.h"



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
    
    
 /*   vector<F<AAF>> temp(sysdim);
    taui = h;
    
    // only the last differentiation is computed with full intervals (improved mean value)
    for (int j=0 ; j<sysdim; j++)
        temp[j] = odeVAR_x[sysdim+inputsdim-1].x[j][0].x();
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            temp[j] += odeVAR_x[sysdim+inputsdim-1].x[j][i].x()*taui;
        taui *= h;
    }
    
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        temp[j] += odeVAR_g.x[j][order].x()*taui;
    
    // vector<F<AAF>>
    nn_outputs = NH.eval_network(temp); */
    
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
    vector<vector<AAF>> aux2(inputsdim, vector<AAF>(sysdim));
    vector<vector<AAF>> Jaci(sysdim, vector<AAF>(jacdim));
    
    J_res = J;
    
    // partial u(z) / partial z0
    for (int i=0 ; i<inputsdim ; i++) {
        for (int j=0 ; j<sysdim; j++) {
            aux2[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                aux2[i][j] += AAF(Jac_param_inputs[i][k])*J[k][j];
  //          cout << "aux2["<< i <<"]["<< j << "]=" << aux2[i][j].convert_int() << endl;
        }
    }
    
    //if (syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 492 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388)
     //   addJacfzJacfz(J_res,Jac_params);
    
    // pour les premiers termes, Jac^i(x)*J*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++) {
            for (int k=0 ; k<sysdim; k++) {
                Jaci[j][k] = odeVAR_x[k].x[j][i].d(0);  //  on each structure differentiate to only one variable
                
          //      cout << "order=" << i << " Jaci["<< j <<"]["<< k << "]=" << Jaci[j][k].convert_int() << endl;
            }
            for (int k=0 ; k<inputsdim ; k++) {
                Jaci[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])] = odeVAR_x[sysdim+k].x[j][i].d(0);
            //    cout << "order=" << i << " Jaci["<< j <<"]["<< sysdim+index_param_inv[k]+floor(offset*nb_inputs[k]) << "]=" << Jaci[j][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])].convert_int() << endl;
            }
        }
        
        
        
        
    //    for (int k=0 ; k<sysdim+inputsdim; k++)
    //        for (int j=0 ; j<inputsdim ; j++)
     //           cout << "order=" << i << " deriv param inputs"<<  odeVAR_x[k].param_inputs[j][i].d(0).convert_int() << " deriv param inputs"<<  odeVAR_x[k].param_inputs[j][i].d(1).convert_int() << endl;
        
       //  multJacfzJaczz0(aux,Jaci,J); // Jac^i(x)*J
        multJacfzuJaczz0Jacuz0(aux,Jaci,aux2,J,offset); // Jac^i(x)*J  - includes the case where u can be a control functio,
        
    //    for (int i=0 ; i<Jaci.size(); i++)
    //        for (int j=0 ; j<Jaci[i].size(); j++)
    //            cout << "Jaci["<<i<<"]["<<j<<"]="<<Jaci[i][j].convert_int() << endl;
    /*    for (int j=0 ; j<sysdim; j++)
            for (int k=0 ; k<jacdim; k++)
                cout << "aux["<<j<<"]["<<k<<"]=" << aux[j][k].convert_int() << endl; */
        
        // TODO. Adding term of the jacobian due to NN control.
        if ((i == 1) && ((syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388|| syschoice == 1111)))
            addMiMi(aux,Jac_params);
         //   for (int j=0 ; j<sysdim; j++)
        //            for (int k=0 ; k<jacdim; k++)
//                        cout << "Jac_params["<<j<<"]["<<k<<"]=" << Jac_params[j][k].convert_int() << endl;
        
        if ((i == 2) && ((syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388|| syschoice == 1111)))
                addMiMi(aux,Jac_params_order2);
    
        // In fact we should compute this for all orders...
        
  //      for (int j=0 ; j<sysdim; j++)
   //             for (int k=0 ; k<jacdim; k++)
    //                cout << "i=" << i << " aux["<<j<<"]["<<k<<"]=" << aux[j][k].convert_int() << endl;
        
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
  //      for (int j=0 ; j<sysdim; j++)
  //          for (int k=0 ; k<jacdim; k++)
  //              cout << "J_res["<<j<<"]["<<k<<"]=" << J_res[j][k].convert_int() << endl;
    

    
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
    
    bf(fx0,param_inputs,params,x0);
    
    for (int i=0; i<sysdim ; i++)
        y1[i] = x0[i] + interval(0,tau)*fx0[i].convert_int();  // modif (*tau)
    
 //   bool save_recompute = recompute_control;
 //   recompute_control = false;
    
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
        
        bf(fx0,param_inputs,params,y0);
        // cout << "fx[0]=" << fx0[0] << "\t"  << fx0[1] << endl;
        //   cout << "y1=" << y1[0] << "\t"  << y1[1] << endl;
        for (int i=0; i<sysdim ; i++)
            y1[i] = x0[i] + interval(0,tau)*fx0[i].convert_int();
        
        
       
        iter = iter+1;
    }
    
//    recompute_control = save_recompute;
    
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
    
    if (syschoice == 491)
        params = NH.eval_network(syst_to_nn(x));
    else    if (syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 492 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388 || syschoice == 1111|| syschoice == 1112)
    {
        params = NH.eval_network(x);
        cout << "params=" << params[0].convert_int() << " x=" << x[0].convert_int() << endl;
        
        
    }
    
    vector<F<AAF>> xf(sysdim);
    for (int i=0 ; i<sysdim; i++)
        xf[i] = x[i];
    for (int i=0; i < sysdim ; i++)
        xf[i].diff(i,sysdim);    // differentiate to x   // TODO. we still need to add the inputs (inputsdim) in the dimensions that we are differentiatiing to // something like for (int k=0 ; k<inputsdim ; k++) { Jac[j][sysdim+...]
    vector<F<AAF>> paramsf;
    if (syschoice == 491) {
        vector<F<AAF>> resf = vector<F<AAF>>(NH.n_inputs);
        resf[0] = 30.0;
        resf[1] = 1.4;
        resf[2] = xf[4]; // v_ego
        resf[3] = xf[0]-xf[3]; // x_lead-x_ego
        resf[4] = xf[1]-xf[4];
        paramsf = NH.eval_network(resf);
    }
    else
        paramsf = NH.eval_network(xf);
    
    OdeVar odeVAR_g = OdeVar(bf,paramsf);
    Ode ode_x0 = Ode(bf,params);
    Ode ode_g0 = Ode(bf,params);
    vector<AAF> init_x;
    
    if (innerapprox == 1)
       init_x = x0;
    else
        init_x = x;
    
    
    
    TM_val TMcenter = TM_val(ode_x0, ode_g0, order, init_x, tn, tau);
    //if (innerapprox == 1)
    
    if (innerapprox == 1) {
        for (int k = 0; k<sysdim+inputsdim; k++)
            odeVAR_x[k] = OdeVar(bf,paramsf);
    }
    TM_Jac TMJac = TM_Jac(odeVAR_x, odeVAR_g, order, x, init_x, J0, tn, tau);
    
    
    
    HybridStep_ode res = HybridStep_ode(bf,TMcenter,TMJac,tn,tau,order);
    
    if (syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 492 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388 || syschoice == 1111|| syschoice == 1112)
    {
        vector<vector<AAF>> J(sysdim, vector<AAF>(sysdim+inputsdim));  // should be jacdim but not yet defined ?
        for (int i=0; i<sysdim; i++)
            J[i][i] = 1.0;
        res.eval_valandJacobian_nn(x,inputs,0,tau,J);
    }
    
    return res;
}




void HybridStep_ode::init_nextstep(vector<AAF> &param_inputs, double _tau)
{
   
    
    TMcenter.init_nextstep(_tau);
    // pass the center of TMcenter to TMJac
    if (innerapprox == 1)
        TMJac.init_nextstep(_tau,TMcenter.x);
    
    tn = tn + tau;
    tau = _tau;
    
    if (syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 492 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388|| syschoice == 1111|| syschoice == 1112) {
       
        eval_valandJacobian_nn(TMJac.xp1,param_inputs,tn,tau,TMJac.J); // J_rough ou J ???
       /* if (control_period > 0) // control update rate is not the same as time step
        {
            if ((tn/control_period >= (int)(tn/control_period)) && ((tn-tau)/control_period < (int)(tn/control_period))) {
                
                vector<F<AAF>> xf(sysdim);
                for (int i=0 ; i<sysdim; i++)
                    xf[i] = TMJac.x[i];
                vector<F<AAF>> paramsf = NH.eval_network(syst_to_nn(xf));
                for (int i=0 ; i<paramsf.size(); i++)
                    params[i] = paramsf[i].x();
                // params = NH.eval_network(syst_to_nn(TMJac.x));
                cout << "recomputing params at time " << tn << endl;
            }
        }
        else {
            vector<F<AAF>> xf(sysdim);
            for (int i=0 ; i<sysdim; i++)
                xf[i] = TMJac.x[i];
            vector<F<AAF>> paramsf = NH.eval_network(syst_to_nn(xf));
            for (int i=0 ; i<paramsf.size(); i++)
                params[i] = paramsf[i].x();
           // params = NH.eval_network(syst_to_nn(TMJac.x));
            cout << "recomputing params at time " << tn << endl;
        }
        cout << "params=" << params[0].convert_int() << endl;
        */
    
        // no need to do this if we do not modify params ?
  /*      cout << "setting new bf" << endl;
        bf = OdeFunc();  // to take into account new params if they have changed
        // resetting
        TMcenter.ode_x = Ode(bf);
        TMcenter.ode_g = Ode(bf);
        TMJac.odeVAR_x = vector<OdeVar>(sysdim+inputsdim);
        for (int k = 0; k<sysdim+inputsdim; k++)
            TMJac.odeVAR_x[k] = OdeVar(bf);
        TMJac.odeVAR_g = OdeVar(bf); */
    }
    else if ((control_period == 0) || (tn == 0) || ((tn/control_period >= (int)(tn/control_period)) && ((tn-tau)/control_period < (int)(tn/control_period))))  // control update rate is not the same as time step
    {
        vector<AAF> y_temp(sysdim);
        vector<AAF> control_inputs;
        recompute_control = true;
        bf(y_temp,param_inputs,control_inputs,TMJac.xp1);  // to evaluate new control value
        recompute_control = false;
    }
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


// sets params and Jac_params (Jac params should contain only df/du . du/dx and not all df/dx hence the decomposition implemented below) - separate control and the rest of dependencies because control period is different from integration period
// Jac_params_order2 contains d/dt (df/du . du/dx) = d/dx(Jac_order_1).(dx/dt)
void HybridStep_ode::eval_valandJacobian_nn(vector<AAF> x, vector<AAF> &param_inputs, double tn, double tau, vector<vector<AAF>> J)
{
    if (syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 ||syschoice == 491 || syschoice == 492 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388|| syschoice == 1111|| syschoice == 1112)
    {
        
        if ((control_period == 0) || (tn == 0) || ((tn/control_period >= (int)(tn/control_period)) && ((tn-tau)/control_period < (int)(tn/control_period))))  // control update rate is not the same as time step
        {
                    
            // evaluate (\partial u / partial x) separately
            
            vector<F<AAF>> xf(sysdim);
            for (int i=0 ; i<sysdim; i++)
                xf[i] = x[i];
            for (int i=0; i < sysdim ; i++)
                xf[i].diff(i,sysdim);    // differentiate to x
            // TODO. we still need to add the inputs (inputsdim) in the dimensions that we are differentiatiing to
            // something like for (int k=0 ; k<inputsdim ; k++) { Jac[j][sysdim+...]
                    
            vector<F<AAF>> paramsf;
            if (syschoice == 491) {
                vector<F<AAF>> resf = vector<F<AAF>>(NH.n_inputs);
                resf[0] = 30.0;
                resf[1] = 1.4;
                resf[2] = xf[4]; // v_ego
                resf[3] = xf[0]-xf[3]; // x_lead-x_ego
                resf[4] = xf[1]-xf[4];
                paramsf = NH.eval_network(resf);
            }
            else
                paramsf = NH.eval_network(xf);
            for (int i=0 ; i<sysdim_params; i++) {
                params[i] = paramsf[i].x();
                cout << "params[i]=" << params[i].convert_int() << endl;
            }
                        
            vector<vector<AAF>> auxu = vector<vector<AAF>>(sysdim, vector<AAF>(sysdim));
            for (int i=0; i < sysdim_params ; i++)
                for (int j=0; j < sysdim ; j++)
                    auxu[i][j] = paramsf[i].d(j);
                  
            
            // evaluate (\partial f) / (\partial u)
            OdeFunc f = OdeFunc();
            vector<vector<AAF>> auxf = vector<vector<AAF>>(sysdim, vector<AAF>(sysdim_params));
            
            
            // new experiment for order 2
            vector<F<F<AAF>>> ffu(sysdim_params);
            vector<F<F<AAF>>> ffxff(sysdim);
            vector<F<F<AAF>>> fftemp(sysdim);
            vector<F<F<AAF>>> ffparam_inputs(inputsdim);
            
            for (int i=0 ; i<sysdim; i++)
                ffxff[i] = x[i];
            for (int i=0; i < inputsdim ; i++)
                ffparam_inputs[i] = param_inputs[i];
            for (int i=0; i < sysdim ; i++) {
                ffxff[i].diff(i,sysdim+sysdim_params);          // differentiate to x, first order
                ffxff[i].x().diff(i,sysdim+sysdim_params);      // second order
            }
            
            for (int i=0; i < sysdim_params ; i++) {
                ffu[i] = params[i];
                ffu[i].diff(sysdim+i,sysdim+sysdim_params);    // differentiate to u, first order
                ffu[i].x().diff(sysdim+i,sysdim+sysdim_params);  // second order
            }
            
            f(fftemp,ffparam_inputs,ffu,ffxff);
            
            for (int i=0; i < sysdim ; i++)
                for (int j=0; j < sysdim_params ; j++)
                    auxf[i][j] = fftemp[i].d(sysdim+j).x();   // \partial f . partial u
           
        //    if (syschoice == 471)
        //        auxf[1][0] = 4 * x[1] * x[1];
         //   if (syschoice == 1112)
               // auxf[1][0] = 1;
          //      auxf[1][0] = x[1];
            
            // compute (\partial f) / (\partial u) . (\partial u / partial x)
            vector<vector<AAF>> aux = vector<vector<AAF>>(sysdim, vector<AAF>(sysdim));
            for (int i=0 ; i<sysdim; i++)
                for (int j=0; j < sysdim ; j++) {
                    aux[i][j] = 0;
                    for (int k=0; k < sysdim_params ; k++)
                        aux[i][j] += auxf[i][k]*auxu[k][j];
                   // cout << "partial derivatives: " << aux[i][j].convert_int() << endl;
                }
            
            // compte aux . J and store the result in Jac_params
             multMiMi(Jac_params,aux,J);
            
        
            
            // EVALUATE ORDER 2 (Lie derivative of Jacobian)
            vector<vector<AAF>> aux_order2 = vector<vector<AAF>>(sysdim, vector<AAF>(sysdim));
            vector<vector<F<AAF>>> f_order2 = vector<vector<F<AAF>>>(sysdim, vector<F<AAF>>(sysdim));
            // if J is the 1st order jacobian, evaluate for all i,j in sysdim, \sum_k (partial J_ij / partial z_k) . (dz_k.dt) where dz_k.dt=f_k
               for (int i=0 ; i<sysdim; i++)
                    for (int j=0; j < sysdim ; j++) {
                        f_order2[i][j] = 0;
                        for (int k=0; k < sysdim_params ; k++)
                            f_order2[i][j] += fftemp[i].d(sysdim+k)*auxu[k][j];
                    }
            
            for (int i=0 ; i<sysdim; i++)
                for (int j=0; j < sysdim ; j++) {
                    aux_order2[i][j] = 0;
                    for (int k=0; k < sysdim ; k++)
                        aux_order2[i][j] += f_order2[i][j].d(k)*fftemp[k].x().x();
               //     cout << "partial derivatives order2 " << i << " " << j << " new new : " << aux_order2[i][j].convert_int() << endl;
                }
       
            // compte aux . J and store the result in Jac_params
            multMiMi(Jac_params_order2,aux_order2,J);
            
            // WE COULD/SHOULD PROCEED WITH HIGHER ORDERS....
            
     /*       for (int i=0 ; i<sysdim; i++)
                for (int j=0 ; j<Jac_params[i].size(); j++)
                    cout << "Jac_params["<<i<<"]["<<j<<"]=" << Jac_params[i][j].convert_int() << endl;
            
            for (int i=0 ; i<sysdim; i++)
                for (int j=0 ; j<Jac_params_order2[i].size(); j++)
                    cout << "Jac_params_order2["<<i<<"]["<<j<<"]=" << Jac_params_order2[i][j].convert_int() << endl;
         */
  
        
        cout << "setting new bf at time " << tn << endl;
              bf = OdeFunc();  // to take into account new params if they have changed
              // resetting
              TMcenter.ode_x = Ode(bf,params);
              TMcenter.ode_g = Ode(bf,params);
              TMJac.odeVAR_x = vector<OdeVar>(sysdim+inputsdim);
              for (int k = 0; k<sysdim+inputsdim; k++)
                  TMJac.odeVAR_x[k] = OdeVar(bf,paramsf);
              TMJac.odeVAR_g = OdeVar(bf,paramsf);
            
    }
    }
        
}




void HybridStep_ode::print_solutionstep(vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xouter_minimal, vector<interval> &Xinner,  vector<interval> &Xinner_robust, vector<interval> &Xinner_minimal, vector<interval> &Xcenter, vector<interval> &sampled_reachset, int current_subdiv)
{
    //double mean_dist; // mean value on the xi of max distance between inner and outer approximations
    double tnp1 = tn + tau;
    
    
    if (print_debug)
    {
        cout << "t=" << tnp1 << endl;
    
        if (innerapprox == 0)
        {
            for (int i=0 ; i<sysdim ; i++)
                cout << "Xouter_maximal[" << i <<"]=" << Xouter[i] << "\t";
        }
        else
        {
            for (int i=0 ; i<sysdim ; i++) {
                cout.precision(6);
                cout << "Xouter_maximal[" << i <<"]=[" << Xouter[i].inf() << ", " << Xouter[i].sup() << "] \t";
                //printf("%.6f\t", Xouter[i]);
                cout << "Xinner_maximal[" << i <<"]=[" << Xinner[i].inf() << ", " << Xinner[i].sup() << "] \t"; //<<"]=" << Xinner[i] << "\t";
                cout << "Sampled estim.[" << i <<"]=[" << sampled_reachset[i].inf() << ", " << sampled_reachset[i].sup() << "] \t"; //<<"]=" << sampled_reachset[i] << "\t";
                cout << " eta_o["<<i<<"]=" << (sampled_reachset[i].sup() - sampled_reachset[i].inf())/ (Xouter[i].sup() - Xouter[i].inf()) << "\t";
                cout << " eta_i[" << i << "]=" << (Xinner[i].sup() - Xinner[i].inf())/(sampled_reachset[i].sup() - sampled_reachset[i].inf()) << "\t";
                cout << " gamma[" << i << "]=" << (Xinner[i].sup() - Xinner[i].inf())/(Xouter[i].sup() - Xouter[i].inf());
                if (! subseteq(sampled_reachset[i],Xouter[i]))
                    cout << "WARNING: problem with outer-estimated approximation!\n";
                if (! subseteq(Xinner[i],sampled_reachset[i]))
                    cout << "WARNING: problem with inner-estimated approximation!\n";
                cout << endl;
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
    }
  
   // out_approx << YAML::BeginMap;
    
    
    
    vector<double> temp(2*sysdim);
    
    if (nb_subdiv_init ==1) {
        out_approx << YAML::Key << "outer";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xouter[i].inf();
            temp[2*i+1] = Xouter[i].sup();
        }
        out_approx << YAML::Value << temp; // Xouter does not work because of interval type (I guess I could solve this but...)
    }
    
    if (innerapprox == 1)
    {
        out_approx << YAML::Key << "inner";
        for (int i=0 ; i<sysdim ; i++) {
            temp[2*i] = Xinner[i].inf();
            temp[2*i+1] = Xinner[i].sup();
        }
        out_approx << YAML::Value << temp; // Xinner;
        if (nb_subdiv_init ==1)
        {
            out_approx << YAML::Key << "center";
            for (int i=0 ; i<sysdim ; i++) {
                temp[2*i] = Xcenter[i].inf();
                temp[2*i+1] = Xcenter[i].sup();
            }
            out_approx << YAML::Value << temp; // Xcenter;
        }
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
        
        // error measures
        vector<double> temp2(sysdim);
        for (int i=0 ; i<sysdim ; i++)
            temp2[i] = (sampled_reachset[i].sup() - sampled_reachset[i].inf())/ (Xouter[i].sup() - Xouter[i].inf());
        out_approx << YAML::Key << "etaouter";
        out_approx << YAML::Value << temp2;
        for (int i=0 ; i<sysdim ; i++)
            temp2[i] = (Xinner[i].sup() - Xinner[i].inf())/(sampled_reachset[i].sup() - sampled_reachset[i].inf());
        out_approx << YAML::Key << "etainner";
        out_approx << YAML::Value << temp2;
        for (int i=0 ; i<sysdim ; i++)
            temp2[i] = (Xinner[i].sup() - Xinner[i].inf())/(Xouter[i].sup() - Xouter[i].inf());
        out_approx << YAML::Key << "gamma";
        out_approx << YAML::Value << temp2;
        
        
    }
 //   out_approx << YAML::EndMap;
    
    
    for (int i=0 ; i<sysdim ; i++) {
        // saving result
        Xouter_print[current_subdiv][current_iteration][i] = Xouter[i];
        if (innerapprox == 1)
        {
            Xouter_robust_print[current_subdiv][current_iteration][i] = Xouter_robust[i];
            Xouter_minimal_print[current_subdiv][current_iteration][i] = Xouter_minimal[i];
            Xinner_print[current_subdiv][current_iteration][i] = Xinner[i];
            Xinner_robust_print[current_subdiv][current_iteration][i] = Xinner_robust[i];
            Xinner_minimal_print[current_subdiv][current_iteration][i] = Xinner_robust[i];
        }
        t_print[current_iteration] = tnp1;
    }
    current_iteration++;
    
    
}



vector<interval> HybridStep_ode::TM_evalandprint_solutionstep(vector<interval> &eps, double tnp1, vector<interval> &sampled_reachset, int current_subdiv)
{
    assert (tn <= tnp1);
    assert (tnp1 <= tn+tau);
    
    vector<interval> Xouter(sysdim), Xouter_robust(sysdim), Xouter_minimal(sysdim), Xinner(sysdim), Xinner_robust(sysdim), Xinner_minimal(sysdim), Xcenter(sysdim);
    
    // eval and store at time tnp1
    TM_eval();
    
    out_approx << YAML::BeginMap;
    
    out_approx << YAML::Key << "tn";
    out_approx << YAML::Value << tnp1;
    if (nb_subdiv_init > 1)
    {
        out_approx << YAML::Key << "currentsubdiv";
        out_approx << YAML::Value << current_subdiv;
    }
    
    if (innerapprox == 0)
    {
        for (int i = 0 ; i<sysdim ; i++) {
            TMcenter.xp1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
            TMcenter.xp1[i].sumup(tol_noise); // group small terms
            Xouter[i] = TMcenter.xp1[i].convert_int();
        }
        print_solutionstep(Xouter,Xouter,Xouter,Xouter,Xouter,Xouter,Xouter,sampled_reachset,current_subdiv);
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
        print_solutionstep(Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter,sampled_reachset,current_subdiv);
        
    // print_solutionstep_ode(Xouter,Xinner,Xcenter,tnp1);
    }
    out_approx << YAML::EndMap;
    
    return Xouter;
}


// setting the k-th control_values in fullinputs and param_inputs
// used when we need to set control values dynamically and not statically before execution as done until now
// EN FAIT j'ai l'impression que cette fonction soit n'a pas ete terminee soit ne fait pas ce qu'elle est prevue pour
// TODO. A vori si je corrige pour utiliser le control_input ou pas au final (dans les champs)
void HybridStep_ode::set_controlinput(vector<AAF> &param_inputs, vector<AAF> &param_inputs_center, const vector<AAF> &control_input, int k)
{
    // the part which is statically done in init_system
    interval temp;
    for (int i=0; i<inputsdim; i++) {
        cout << index_param_inv[i]+k << endl;
        fullinputs[index_param_inv[i]+k] = control_input[i];
        temp = control_input[i].convert_int();
        center_fullinputs[index_param_inv[i]+k] = mid(temp);
        eps[sysdim+index_param_inv[i]+k] = temp - mid(temp);
    }
    
    // the part which is statically done in set_initialconditions
    for (int i=0; i<inputsdim; i++) {
        param_inputs[index_param_inv[i]+k] = fullinputs[index_param_inv[i]+k];
        if (innerapprox == 1)
            param_inputs_center[index_param_inv[i]+k] = center_fullinputs[index_param_inv[i]+k];
        else
            param_inputs_center[index_param_inv[i]+k] = fullinputs[index_param_inv[i]+k];
    }
}

// setting the k-th control_values in fullinputs and param_inputs
// used when we need to set control values dynamically and not statically before execution as done until now
void HybridStep_ode::set_controlinput_regression(vector<AAF> &param_inputs, vector<AAF> &param_inputs_center, const vector<AAF> &control_input, const vector<interval> &control_input_uncertainty, vector<vector<interval>> Jac_control_input, int k)
{
    // the part which is statically done in init_system
    interval temp;
    for (int i=0; i<inputsdim; i++) {
        cout << index_param_inv[i]+k << endl;
        fullinputs[index_param_inv[i]+k] = control_input_uncertainty[i];
        temp = control_input_uncertainty[i];
        center_fullinputs[index_param_inv[i]+k] = mid(temp);
        eps[sysdim+index_param_inv[i]+k] = temp - mid(temp);
    }
    
    // the part which is statically done in set_initialconditions
    for (int i=0; i<inputsdim; i++) {
        param_inputs[index_param_inv[i]+k] = control_input[i] + control_input_uncertainty[i];
        for (int j=0 ; j<sysdim; j++)
           // Jac_param_inputs[index_param_inv[i]+k][j] = Jac_control_input[i][j];
            Jac_param_inputs[i][j] = Jac_control_input[i][j];
        if (innerapprox == 1)
            param_inputs_center[index_param_inv[i]+k] = mid(control_input[i].convert_int()) + center_fullinputs[index_param_inv[i]+k];
        else
            param_inputs_center[index_param_inv[i]+k] = control_input[i] + control_input_uncertainty[i];
    }
}


// estimate the range of the n iterates (same stepsize as for reachability analysis)
vector<vector<interval>> estimate_reachset(OdeFunc &obf, vector<AAF> &initial_values, vector<AAF> &param_inputs, double t_begin, double t_end, double tau, int discr)
{
    int n = (t_end - t_begin)/tau + 1;
   // int discr = 2;
    int nb_points = discr+1;
    
    
    
    vector<interval> xinit(sysdim);
    for (int i=0 ; i<sysdim ; i++)
        xinit[i] = initial_values[i].convert_int();
    
    // limit the number of sampled points
    for (int i=1; i < min(jacdim,4) ; i++)  // MODIF borne min : 1 => 0 (?)
        nb_points = nb_points * (discr+1);
    
    
    vector<vector<double>> input(nb_points,vector<double>(sysdim));  //  the iterates f^n(x_j)
    vector<vector<double>> output(nb_points,vector<double>(sysdim));
    
    vector<vector<double>> max_output(n+1,vector<double>(sysdim));  // store the min and max for each iterate
    vector<vector<double>> min_output(n+1,vector<double>(sysdim));
    vector<vector<interval>> range(n+1,vector<interval>(sysdim));

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

    
    ofstream samplesreachsetfile;
    samplesreachsetfile.open("output/samplesreachset.yaml");
    
    YAML::Emitter out_samples;
    out_samples << YAML::BeginMap;
    
    out_samples << YAML::Key << "samples";
    out_samples << YAML::Value << YAML::BeginSeq;
    
    double tn;
    
    
    vector<AAF> control_inputs(sysdim_params);
    vector<AAF> input_dummy(sysdim), output_dummy(sysdim);
    vector<double> control(sysdim_params);
    
    recompute_control = false;
    
    int iter = 1;
    for (tn=t_begin ; tn <=t_end ; tn = tn+tau)
    {
       
        for (int i=0; i < sysdim ; i++) {
            max_output[iter][i] = -DBL_MAX; // RK(obf,input[0],param_inputs,control_inputs,tau);
            min_output[iter][i] = DBL_MAX; // RK(obf,input[0],param_inputs,control_inputs,tau);
        }
        iter++;
    }
    
    for (cur_point=0 ; cur_point<nb_points; cur_point++)
    {
        iter = 1;
            
        for (tn=t_begin ; tn <t_end-0.00001*t_end ; tn = tn+tau)
        {
                
            if ((control_period == 0) || (tn == 0) || ((tn/control_period >= (int)(tn/control_period)) && ((tn-tau)/control_period < (int)(tn/control_period))))  // control update rate is not the same as time step
            {
                if (syschoice == 461 ||syschoice == 451 || syschoice == 471 || syschoice == 481 || syschoice == 482 || syschoice == 483 || syschoice == 484 || syschoice == 492 || syschoice == 493 || syschoice == 381 || syschoice == 382|| syschoice == 383 || syschoice == 384|| syschoice == 385|| syschoice == 386 || syschoice == 387|| syschoice == 388|| syschoice == 1111|| syschoice == 1112)
                {
                    control = NH.eval_network(input[cur_point]);
                    for (int i=0; i < control.size() ; i++)
                        params[i] = control[i];
                 }
                else if (syschoice == 491)
                {
                    vector<double> res = vector<double>(NH.n_inputs);
                    res[0] = 30.0;
                    res[1] = 1.4;
                    res[2] = input[cur_point][4]; // v_ego
                    res[3] = input[cur_point][0]-input[cur_point][3]; // x_lead-x_ego
                    res[4] = input[cur_point][1]-input[cur_point][4];
                    control = NH.eval_network(res);
                    for (int i=0; i < control.size() ; i++)
                        params[i] = control[i];
                }
            }
            
            // param_inputs,control_inputs unused for now ?
            output[cur_point] = RK(obf,input[cur_point],param_inputs,params,tau);
            //outFile_xi << iter <<  "\t";
            //for (int i=0; i < sysdim ; i++)
            //    outFile_xi << output[cur_point][i] <<  "\t";
            //outFile_xi << endl;
            
            out_samples << YAML::BeginMap;
            out_samples << YAML::Key << "tn";
            out_samples << YAML::Value << tn+tau;
            out_samples << YAML::Key << "sample";
            out_samples << YAML::Value << output[cur_point];
            out_samples << YAML::EndMap;
           // out_samples << iter << output[cur_point];
            
            for (int i=0; i < sysdim ; i++) {
                if (output[cur_point][i] < min_output[iter][i])
                    min_output[iter][i] = output[cur_point][i];
                if (output[cur_point][i] > max_output[iter][i])
                    max_output[iter][i] = output[cur_point][i];
            }
            // initializing next step (iter)
            for (int i=0; i < sysdim ; i++)
                input[cur_point][i] = output[cur_point][i];
                
            iter++;
        }
        
    }
    
    iter = 1;
    for (tn=t_begin ; tn <=t_end ; tn = tn+tau)
    {
        cout << "Estimated reachable set at tn = " << tn << " is: ";
        for (int i=0; i < sysdim ; i++) {
            cout << "z["<<i << "]=[" << min_output[iter][i] << ", " << max_output[iter][i] <<"]  ";
            
        }
        cout << endl;
        for (int i=0; i < sysdim ; i++)
            range[iter][i] = interval(min_output[iter][i],max_output[iter][i]);
        
        iter++;
    }
    
    // resetting params to initial condition

    if (syschoice == 491)
        params = NH.eval_network(syst_to_nn(initial_values));
    else
        params = NH.eval_network(initial_values);
    
    out_samples << YAML::EndSeq;
    out_samples << YAML::EndMap;
    samplesreachsetfile << out_samples.c_str();
    samplesreachsetfile.close();
    
    // outFile_xi.close();
    
    return range;
}

// ponctual values but AAF necessary to comply with obf...
vector<double> RK(OdeFunc &obf, vector<double> &yn, vector<AAF> &param_inputs, vector<AAF> &control_inputs, double h)
{
    vector<AAF> k1(sysdim), kb(sysdim), k2(sysdim), k3(sysdim), k4(sysdim), ynloc(sysdim);
    vector<double> ynp1(sysdim);
    
    for (int i=0 ; i<sysdim; i++)
        ynloc[i] = yn[i];
    
    //    k1 = h * f(yn);
    obf(k1,param_inputs,control_inputs,ynloc);
    for (int i=0 ; i<sysdim; i++)
        k1[i] = h*k1[i];
    
    // k2 = h * (f((x0+h/2), (y0+k1/2)));
    for (int i=0 ; i<sysdim; i++)
        kb[i] = yn[i] + k1[i]/2.0;
    obf(k2,param_inputs,control_inputs,kb);
    for (int i=0 ; i<sysdim; i++)
        k2[i] = h*k2[i];
    
    // k3 = h * (f((x0+h/2), (y0+k2/2)));
    for (int i=0 ; i<sysdim; i++)
        kb[i] = yn[i] + k2[i]/2.0;
    obf(k3,param_inputs,control_inputs,kb);
    for (int i=0 ; i<sysdim; i++)
        k3[i] = h*k3[i];
    
    // k4 = h * (f((x0+h), (y0+k3)));
    for (int i=0 ; i<sysdim; i++)
        kb[i] = yn[i] + k3[i];
    obf(k4,param_inputs,control_inputs,kb);
    for (int i=0 ; i<sysdim; i++)
        k4[i] = h*k4[i];
    
 /*   k = (k1+2*k2+2*k3+k4)/6;
      yn = y0 + k; */
    for (int i=0 ; i<sysdim; i++)
        ynp1[i] = yn[i] + (k1[i].convert_int().inf()+2.0*k2[i].convert_int().inf()+2.0*k3[i].convert_int().inf()+k4[i].convert_int().inf())/6.0;
    return ynp1;
}
