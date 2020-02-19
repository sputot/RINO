/* ============================================================================
 File   : discrete system.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file computes ranges for functions and discrete dynamical systems
 ============================================================================*/

#include "discrete_system.h"
#include "matrix.h"
#include "inner.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <ctime>
#include <assert.h>
using namespace std;


void range_discrete_system(void) {
    
    if (syschoice == 1) {
        jacdim = 1;
        sysdim = 1;
    }
    else if (syschoice == 2) {
        jacdim = 2;
        sysdim = 1;
    }
    initial_values = vector<AAF>(jacdim);
    
    vector<bool> is_existential(jacdim);
    for (int i=0; i < jacdim ; i++)
        is_existential[i] = true;
        
    if (syschoice == 1) {
        initial_values[0] = interval(2,3);
    }
    else if (syschoice == 2) {
        initial_values[0] = interval(2,3);
        initial_values[1] = interval(2,3);
    }
    
    vector<F<AAF>> x(jacdim);
    vector<vector<interval>> Jacf(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jacf[i] = vector<interval>(jacdim);
    vector<F<AAF>> z(sysdim);
    vector<interval> x0(jacdim);    // center
    
    DiscreteFunc f;
    
    eps = vector<interval>(jacdim);
    
    
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    
    z = f(x);
    
    for (int i=0; i < jacdim ; i++)
        cout << "outer range of f(x) direct evaluation:" << z[i].x().convert_int() << endl;
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            Jacf[i][j] = z[i].d(j).convert_int();
        }
    }
    
    vector<interval> z0 = f(x0);
    
    evaluate_ranges(z0,Jacf,is_existential);
    is_existential[0] = false;
    evaluate_ranges(z0,Jacf,is_existential);
    
    // cout << "gradient of f:" << z[0].d(0).convert_int() << endl;
}


// A FINIR : AJOUTER LE CAS OU DES VAR SONT UNIVERSELLES
void evaluate_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, vector<bool> &is_existential) {
    interval inner_impro;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = z0[i];
        inner_impro = 0;
        for (int j=0; j < jacdim ; j++) {
            z_outer[i] += Jacf[i][j]*eps[j];
            inner_impro += Kaucher_multeps(Jacf[i][j],eps[j]);
        }
        z_inner[i] = Kaucher_add_pro_impro(z0[i],inner_impro);
    }
    
    string str = "";
    for (int i=0; i < jacdim ; i++)
        if (!is_existential[i])
            str = str+to_string(i);
        
    for (int i=0; i < sysdim ; i++) {
        cout << "outer range of f(x"<< str <<") by mean-value:" << z_outer[i] << endl;
        cout << "inner range of f(x) by mean-value:" << z_inner[i] << endl;
    }
}
