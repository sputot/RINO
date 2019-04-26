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
void TM_val::build(OdeFunc _bf) {
    vector<AAF> g_rough(sysdim); // rough estimation of solution of current mode on [tn,tn+tau]
    int i, j;
    
    // compute Taylor coefficients / Lie derivatives starting from center of initial conditions
    ode_x.reset();
    for (j=0 ; j<sysdim ; j++)
        ode_x.x[j][0]=x[j]; // initialize with center;
    
    // computing a priori enclosure on [tn,tn+tau] of solution starting from center of initial conditions
    g_rough = fixpoint(_bf,x,tau);
    
    
        for (int i=0 ; i<sysdim ; i++)
            g_rough[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term
    
    
    // evaluate the Lie derivatives on this enclosure: we will use the last coefficients as remainder terms of Taylor model
    ode_g.reset();
    for (j=0 ; j<sysdim ; j++)
        ode_g.x[j][0] = g_rough[j]; // initialize with enclosure of the flow on the time step
    
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
    
    for (int j=0 ; j<sysdim; j++)
    {
        if (is_variable[j]) // variable param with value always in the same range
        res[j] = ode_x.x[j][0];
    }
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
        
        if (is_variable[i]) // removing dependencies to previous occurence in time dependent input
            xp1[i] = xp1[i].convert_int();
    }
    x = xp1;
}






// TM of the Jacobian with respect to IC of the solution of the ODE for all the range of IC
// outputs are odeVAR_x (Taylor coefficients from 0 to order-1), odeVAR_g (for the remainder term,
// evaluated using a priori enclosure on [tn,tn+tau]), and J_rough, also used for the remainder term
void TM_Jac::build(OdeFunc _bf) {
    vector<AAF> g_rough(sysdim);
    vector<vector<AAF>> Jac1_g_rough(sysdim, vector<AAF>(sysdim));
    int i, j;
    
    // compute Lie derivatives and Jacobian on CI: these will give coefficients 0 to order-1 of TM for Jacobian
    odeVAR_x.reset();
    for (j=0 ; j<sysdim ; j++)
        odeVAR_x.x[j][0] = x[j];
    
    // compute an a priori enclosure on [tn,tn+tau] of solution
    g_rough = fixpoint(_bf,x,tau);
    
 
        for (i=0 ; i<sysdim ; i++) {
            g_rough[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term (for efficiency)
        }
    
    // evaluate the Lie derivatives and Jacobian on this a priori enclosure on [tn,tn+tau]: we will use the last coefficients as remainder terms of Taylor model
    odeVAR_g.reset();
    for (j=0 ; j<sysdim ; j++)
        odeVAR_g.x[j][0] = g_rough[j];
    
    // specify the variables to differentiate
    for (j=0 ; j<jacdim ; j++) {
        odeVAR_x.x[j][0].diff(j,jacdim);
        odeVAR_g.x[j][0].diff(j,jacdim);
    }

    for(i=0;i<order;i++)
    {
        // Evaluate the i'th Taylor coefficient of the r.h.s. of the ODE:
        for (j=0 ; j<sysdim ; j++) {
            odeVAR_x.xp[j].eval(i);
            odeVAR_g.xp[j].eval(i);
            
            // Since d(x,y,z,p)/dt=f(x,y,z,p) we have
            odeVAR_x.x[j][i+1]=odeVAR_x.xp[j][i]/double(i+1);
            odeVAR_g.x[j][i+1]=odeVAR_g.xp[j][i]/double(i+1);
            // ode.x[0]...ode.x[10] now contains the Taylor-coefficients of the solution of the ODE.
        }
    }
   
    // first order coefficient Jac1_g_rough of the  Taylor model for the Jacobian, evaluated on a priori enclosure of the flow (for remainder term)
    for (j=0 ; j<sysdim ; j++)
        for (int k=0 ; k<sysdim ; k++)
        {
            Jac1_g_rough[j][k] = odeVAR_g.x[j][1].d(k);
      //      cout << "Jac1_g_rough[i][k] in build: " << Jac1_g_rough[j][k].convert_int() << endl;
       //     cout << "Jac1_g_rough[i][k] in build: " << odeVAR_x.x[j][1].d(k).convert_int() << endl;
        }
   
   // for (j=0 ; j<sysdim ; j++)
   //     for (int k=sysdim ; k<jacdim ; k++)
   //  cout << "differentiation in build: " << odeVAR_x.x[j][1].d(k).convert_int() << endl;
    
    // deduce a priori enclosure J_rough of the Jacobian (used to compute the remainder)
    fixpoint(J_rough, Jac1_g_rough,J,tau);
    
  
        for (i=0 ; i<sysdim ; i++)
            for (int k=0 ; k<jacdim ; k++)
            {
                J_rough[i][k].sumup(tol_noise); // group terms smaller than eps*radius in a new error term (for efficiency)
         //       cout << "J_rough[i][k] in build: " << J_rough[i][k].convert_int() << endl;
            }
    
}



// eval TM at time tn+h and return value in res
void TM_Jac::eval_val(vector<AAF> &res, double h)
{
    double taui = h;
    
    for (int j=0 ; j<sysdim; j++)
        res[j] = odeVAR_x.x[j][0].x();
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            res[j] += odeVAR_x.x[j][i].x()*taui;
        taui *= h;
    }
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        res[j] += odeVAR_g.x[j][order].x()*taui;
    
    for (int j=0 ; j<sysdim; j++)
    {
    if (is_variable[j]) // variable param with value always in the same range
        res[j] = odeVAR_x.x[j][0].x();
    }
}




// eval TM at tn+h and return in J_res
void TM_Jac::eval_Jac(vector<vector<AAF>> &J_res, double h)
{
    double taui = h;
    vector<vector<AAF>> aux(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> Jaci(sysdim, vector<AAF>(sysdim));
    
    J_res = J;
    
    // pour les premiers termes, Jac^i(x)*J*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            for (int k=0 ; k<sysdim; k++)
                Jaci[j][k] = odeVAR_x.x[j][i].d(k);
        multMiMi(aux,Jaci,J); // Jac^i(x)*J
        scaleM(aux,taui);
        addMiMi(J_res,aux);
        taui *= h;
    }
    // pour le dernier terme, Jac^i(g_rough)*J_rough*tau^i
    for (int j=0 ; j<sysdim; j++)
        for (int k=0 ; k<sysdim; k++)
            Jaci[j][k] = odeVAR_g.x[j][order].d(k);
    //  print_interv("Jaci",Jaci);
    multMiMi(aux,Jaci,J_rough); // Jac^i(g_rough)*J_rough
    //    print_interv("J_rough",J_rough);
    scaleM(aux,taui);
    addMiMi(J_res,aux);
}



// eval TM at tn+tau and store in xp1 and J
void TM_Jac::eval()
{
    eval_val(xp1,tau);
    eval_Jac(Jp1,tau);
}


void TM_Jac::init_nextstep(double _tau)
{
    tn = tn + tau;
    tau = _tau;
    
    for (int i=0 ; i<sysdim ; i++) {
        xp1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
        xp1[i].sumup(tol_noise); // group small terms
        
        if (is_variable[i]) // removing dependencies to previous occurence in time dependent input
            xp1[i] = xp1[i].convert_int();
    }
    x = xp1;
    
    for (int i=0 ; i<sysdim ; i++)
        for (int k=0 ; k<sysdim ; k++) {
            Jp1[i][k].compact();
            Jp1[i][k].sumup(tol_noise);
            
            // AJOUTER QQCHOSE POUR ENLEVER DEPENDENCES AUSSI ?
        }
    
    J = Jp1; // J = Jp1
}



vector<AAF> fixpoint(OdeFunc bf, vector<AAF> &x0, double tau)
{
    // calcul de x satisfaisant x0 + [0,tau][f](x) \subseteq x
    vector<AAF> y0(sysdim);
    vector<AAF> y1(sysdim);
    vector<AAF> fx0(sysdim);
    int iter;
    interval widen, coeff;
    
    bf(fx0,x0);
    
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
        
        bf(fx0,y0);
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
    
    multMiMi(fJ0,Jac1_g_rough,J0);
    
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
                for (int j=0; j<jacdim ; j++)
                    J1[i][j] = J1[i][j] + coeff*widen*J1[i][j].convert_int();
            //   print_interv("fixpoint J1",J1);
        }
        y0 = J1;
        
        multMiMi(fJ0,Jac1_g_rough,y0); // fx0 = f(y0)
        
        for (int i=0; i<sysdim ; i++)
            for (int j=0; j<jacdim ; j++)
                J1[i][j] = J0[i][j] + interval(0,tau)*fJ0[i][j].convert_int();
        
        iter = iter+1;
    }
}





//FadbadVarODE(int n,TFfunction f,void*param= 0);

// enclosure by Taylor model of order order of the flow g_i(center(X0)) at time tau
void gTaylor(vector<AAF> &g, Ode ode_x0, Ode ode_g0, double tau, int order)
{
    double taui = tau;
    for (int j=0 ; j<sysdim; j++)
        g[j] = (ode_x0.x[j])[0];
    
    // pour les premiers termes, on utilise dL^i(x0)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            g[j] += ode_x0.x[j][i]*taui;
        taui *= tau;
    }
    // pour le dernier terme, on utilise le reste dL^i(g0_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        g[j] += ode_g0.x[j][order]*taui;
}


// enclosure by Taylor model of order order of the flow g_i(X0) at time tau
void gTaylor(vector<AAF> &g, OdeVar odeVAR_x, OdeVar odeVAR_g, double tau, int order)
{
    double taui = tau;
    //  AAF* temp_g = new AAF[2];
    for (int j=0 ; j<sysdim; j++)
        g[j] = odeVAR_x.x[j][0].x();
    // pour les premiers termes, on utilise dL^i(x)*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            g[j] += odeVAR_x.x[j][i].x()*taui;
        taui *= tau;
    }
    // pour le dernier terme, on utilise le reste dL^i(g_rough)*tau^i
    for (int j=0 ; j<sysdim; j++)
        g[j] += odeVAR_g.x[j][order].x()*taui;
}




void JTaylor(vector<vector<AAF>> &J_res, OdeVar odeVAR_x, OdeVar odeVAR_g, vector<vector<AAF>> &J0, vector<vector<AAF>> &J_rough, double tau, int order)
{
    double taui = tau;
    vector<vector<AAF>> aux(sysdim, vector<AAF>(jacdim));
    vector<vector<AAF>> Jaci(sysdim, vector<AAF>(sysdim));
    //  iMatrix aux, Jaci;
    //  sizeM(aux,sysdim);
    //  sizeM(Jaci,sysdim);
    
    J_res = J0;
    
    // pour les premiers termes, Jac^i(x)*J0*tau^i
    for (int i=1; i<order ; i++)
    {
        for (int j=0 ; j<sysdim; j++)
            for (int k=0 ; k<sysdim; k++)
                Jaci[j][k] = odeVAR_x.x[j][i].d(k);
        multMiMi(aux,Jaci,J0); // Jac^i(x)*J0
        scaleM(aux,taui);
        addMiMi(J_res,aux);
        taui *= tau;
    }
    // pour le dernier terme, Jac^i(g_rough)*J_rough*tau^i
    for (int j=0 ; j<sysdim; j++)
        for (int k=0 ; k<sysdim; k++)
            Jaci[j][k] = odeVAR_g.x[j][order].d(k);
    //  print_interv("Jaci",Jaci);
    multMiMi(aux,Jaci,J_rough); // Jac^i(g_rough)*J_rough
    //    print_interv("J_rough",J_rough);
    scaleM(aux,taui);
    addMiMi(J_res,aux);
    for (int j=0 ; j<sysdim; j++)
        for (int k=0 ; k<sysdim; k++)
            if (isInfinite(J_res[j][k].convert_int()))
                cout << "Jac_res is infinite!" << endl;
}





/* init next time step of integration:
//     x0, x0p1 are the solution of the ODE starting from the center of initial conditions
//     x, xp1 are the solution of the ODE starting from the full range of initial conditions
//     J0, Jtau are the Jacobian of the solution with respect to initial conditions 
 */
void init_nextstep_ode(vector<AAF> &x0, vector<AAF> &x0p1, vector<AAF> &x, vector<AAF> &xp1, vector<vector<AAF>> &J0, vector<vector<AAF>> &Jtau)
{
    for (int i=0 ; i<sysdim ; i++)
        x0p1[i].compact();  // compact the affine form: remove zero coefficients (there are some ...)
    for (int i=0 ; i<sysdim ; i++) {
        x0p1[i].sumup(10*tol_noise); // group small terms : we keep for now more terms than for x and J, to be adapted
        if (is_variable[i]) // removing dependencies to previous occurence in time dependent input
            x0p1[i] = x0p1[i].convert_int();
    }
    x0 = x0p1;
 //   cout << "x0=" << x0[0].convert_int() << "\t"  << x0[1].convert_int() << endl;
    
    for (int i=0 ; i<sysdim ; i++) {
        xp1[i].sumup(tol_noise); // group terms smaller than eps*radius in a new error term (for efficiency)
        if (is_variable[i]) // removing dependencies to previous occurence in time dependent input
            xp1[i] = xp1[i].convert_int();
    }
    x = xp1;
    
    for (int i=0 ; i<sysdim ; i++)
        for (int k=0 ; k<sysdim ; k++)
            Jtau[i][k].sumup(tol_noise);
    
    J0 = Jtau; // J0 = Jtau
    cout << endl << endl;
    
}






// printing the inner and outer approx at each integration step (more precisely at time tnp1)
void print_solutionstep_ode(vector<interval> &Xouter, vector<interval> &Xinner, vector<interval> &Xcenter, double tnp1)
{
    double aux, minwidth_ratio = 1.0;  // min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
    double mean_dist; // mean value on the xi of max distance between inner and outer approximations
    
    cout << "t=" << tnp1 << endl;
    
    cout << "Xouter=" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "Xouter[" << i <<"]=" << Xouter[i] << "\t";
    cout << endl;
    
    cout << "Xinner=" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "Xinner[" << i <<"]=" << Xinner[i] << "\t";
    cout << endl;
    
    for (int i=0 ; i<sysdim ; i++) {
        outFile_outer[i] << tnp1 << "\t" << inf(Xouter[i]) << "\t" << sup(Xouter[i]) << endl;
        outFile_center[i] << tnp1 << "\t" << inf(Xcenter[i]) << "\t" << sup(Xcenter[i]) << endl;
        outFile_inner[i] << tnp1 << "\t" << inf(Xinner[i]) << "\t" << sup(Xinner[i]) << endl;
   //     outFile_inner_robust[i] << tnp1 << "\t" << inf(Xinner_robust[i]) << "\t" << sup(Xinner_robust[i]) << endl;
    }
    minwidth_ratio = (sup(Xinner[0])-inf(Xinner[0]))/(sup(Xouter[0])-inf(Xouter[0]));
    for (int i=1 ; i<sysdim ; i++) {
        aux = (sup(Xinner[i])-inf(Xinner[i]))/(sup(Xouter[i])-inf(Xouter[i]));
        if (minwidth_ratio > aux)
            minwidth_ratio = aux;
    }
    if (tnp1 != 0)
    outFile_width_ratio << tnp1 << "\t" << minwidth_ratio << endl;
    
  /*  mean_dist = 0.0;
    for (int i=0 ; i<sysdim ; i++) {
        if (! isEmpty(Xinner[i])) {
            if (sup(Xouter[i])-sup(Xinner[i]) > inf(Xinner[i])-inf(Xouter[i]))
                aux = sup(Xouter[i])-sup(Xinner[i]);
            else
                aux = inf(Xinner[i])-inf(Xouter[i]);
            mean_dist += aux;
            cout << "distance between inner and outer approx on component " << i <<" is " << aux << endl;
        }
    }
    
    mean_dist = mean_dist / sysdim;
    if (mean_dist != 0)
        cout << "mean distance between inner and outer approx is " << mean_dist << endl; */
}




HybridStep_ode init_ode(OdeFunc bf, vector<AAF> &x0, vector<AAF> &x, vector<vector<AAF>> &J0, double tn, double tau, int order)
{
    OdeVar odeVAR_x = OdeVar(bf);
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
    TM_Jac TMJac = TM_Jac(odeVAR_x, odeVAR_g, order, x, J0, tn, tau);
    
    HybridStep_ode res = HybridStep_ode(bf,TMcenter,TMJac,tn,tau,order);
    return res;
}




void HybridStep_ode::init_nextstep(double _tau)
{
    
    TMcenter.init_nextstep(_tau);
    if (innerapprox == 1)
        TMJac.init_nextstep(_tau);
    
    tn = tn + tau;
    tau = _tau;
    

}


void HybridStep_ode::TM_build()
{
    TMcenter.build(bf);
    if (innerapprox == 1)
        TMJac.build(bf);
}

// eval s-th Taylor model and initialize s+1-th
void HybridStep_ode::TM_eval()
{
    TMcenter.eval();
    if (innerapprox == 1)
        TMJac.eval();
}




void HybridStep_ode::print_solutionstep(vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xouter_minimal, vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xinner_minimal, vector<interval> &Xcenter)
{
    double aux, minwidth_ratio = 1.0;  // min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
    double mean_dist; // mean value on the xi of max distance between inner and outer approximations
    double tnp1 = tn + tau;
    
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
    
  
    
    for (int i=0 ; i<sysdim ; i++) {
        if (nb_subdiv_init ==1)
        {
        outFile_outer[i] << tnp1 << "\t" << inf(Xouter[i]) << "\t" << sup(Xouter[i]) << endl;
        outFile_center[i] << tnp1 << "\t" << inf(Xcenter[i]) << "\t" << sup(Xcenter[i]) << endl;
        }
        outFile_outer_robust[i] << tnp1 << "\t" << inf(Xouter_robust[i]) << "\t" << sup(Xouter_robust[i]) << endl;
        outFile_outer_minimal[i] << tnp1 << "\t" << inf(Xouter_minimal[i]) << "\t" << sup(Xouter_minimal[i]) << endl;
        outFile_inner[i] << tnp1 << "\t" << inf(Xinner[i]) << "\t" << sup(Xinner[i]) << endl;
        outFile_inner_robust[i] << tnp1 << "\t" << inf(Xinner_robust[i]) << "\t" << sup(Xinner_robust[i]) << endl;
        outFile_inner_minimal[i] << tnp1 << "\t" << inf(Xinner_minimal[i]) << "\t" << sup(Xinner_minimal[i]) << endl;
        
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
    
    minwidth_ratio = (sup(Xinner[0])-inf(Xinner[0]))/(sup(Xouter[0])-inf(Xouter[0]));
    for (int i=1 ; i<sysdim ; i++) {
        aux = (sup(Xinner[i])-inf(Xinner[i]))/(sup(Xouter[i])-inf(Xouter[i]));
        if (minwidth_ratio > aux)
            minwidth_ratio = aux;
    }
    if (tnp1 != 0)
    outFile_width_ratio << tnp1 << "\t" << minwidth_ratio << endl;
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
        InnerOuter(Xinner,Xinner_robust,Xinner_minimal,Xouter,Xouter_robust,Xouter_minimal,TMcenter.xp1,TMJac.Jp1,eps);
      //  InnerOuter(Xinner,Xouter,TMcenter.xp1,TMJac.Jp1,eps);
        intersectViVi(Xouter,TMJac.xp1);
    
    print_solutionstep(Xouter,Xouter_robust,Xouter_minimal,Xinner,Xinner_robust,Xinner_minimal,Xcenter);
    // print_solutionstep_ode(Xouter,Xinner,Xcenter,tnp1);
    }
}





