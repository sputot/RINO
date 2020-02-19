/* ============================================================================
 File   : discrete system.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file computes ranges for functions and discrete dynamical systems
 ============================================================================*/

#include "discrete_system.h"
#include "ode_def.h"
#include "matrix.h"
#include "inner.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <ctime>
#include <assert.h>
using namespace std;


void range_discrete_system(void) {
    sysdim = 1;
    vector<F<AAF>> x(sysdim);
    vector<vector<interval>> Jacf(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jacf[i] = vector<interval>(sysdim);
    vector<F<AAF>> z(sysdim);
    vector<interval> x0(sysdim);    // center
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    DiscreteFunc f;
    
    eps = vector<interval>(sysdim);
    initial_values = vector<AAF>(sysdim);
    
    initial_values[0] = interval(2,3);
    for (int i=0; i < sysdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    
    for (int i=0; i < sysdim ; i++) {
        x[i].diff(i,sysdim);    // differentiate to x
    }
    
    z = f(x);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < sysdim ; j++) {
            Jacf[i][j] = z[i].d(j).convert_int();
        }
    }
    
    vector<interval> z0 = f(x0);
    interval inner_impro;
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = z0[i];
        inner_impro = 0;
        for (int j=0; j < sysdim ; j++) {
            z_outer[i] += Jacf[i][j]*eps[j];
            inner_impro += Kaucher_multeps(Jacf[i][j],eps[j]);
        }
        z_inner[i] = Kaucher_add_pro_impro(z0[i],inner_impro);
    }
    
    for (int i=0; i < sysdim ; i++) {
        cout << "outer range of f direct evaluation:" << z[i].x().convert_int() << endl;
        cout << "outer range of f mean-value:" << z_outer[i] << endl;
        cout << "inner range of f mean-value:" << z_inner[i] << endl;
    }
    // cout << "gradient of f:" << z[0].d(0).convert_int() << endl;
}
