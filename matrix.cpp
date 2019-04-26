/* ============================================================================
 File   : matrix.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file defines rather standard vector/matrix operations over intervals and affine forms (AAF)
 ============================================================================ */

#include "filib_interval.h"
#include "fadbad_aa.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include "matrix.h"
#include "ode_def.h"
using namespace std;




bool subseteq (vector<vector<AAF>> &M1, vector<vector<AAF>> &M2)
{
    for (int i=0; i<sysdim ; i++) {
        for (int j=0; j<jacdim ; j++) {
            if (! subseteq(M1[i][j].convert_int(),M2[i][j].convert_int()))
                return false;
        }
    }
    return true;
}


bool subseteq(vector<AAF> &x, vector<AAF> &y) {
    for (int i=0; i<sysdim ; i++) {
        if (! subseteq(x[i].convert_int(),y[i].convert_int()))
            return false;
    }
    return true;
}

bool subseteq(vector<interval> &x, vector<interval> &y) {
    for (int i=0; i<sysdim ; i++) {
        if (! subseteq(x[i],y[i]))
            return false;
    }
    return true;
}

// beware that arguments are modified

vector<AAF> hull(vector<AAF> &x, vector<AAF> &y) {
    vector<AAF> res(sysdim);
    for (int i=0; i<sysdim ; i++) {
     //   x[i].sumup(tol_noise);
      //  y[i].sumup(tol_noise);
        res[i] = hull(x[i],y[i]);
    }
    return res;
}


void hull(vector<interval> &res, vector<interval> &x, vector<interval> &y)
{
    for (int i=0; i<sysdim ; i++)
        res[i] = hull(x[i],y[i]);
}

// beware that arguments are modified

void hull(vector<AAF> &res, vector<AAF> &x, vector<AAF> &y) {
    //cout << "hull" << endl;
    for (int i=0; i<sysdim ; i++) {
     //   x[i].sumup(tol_noise);
     //   y[i].sumup(tol_noise);
    //    cout << "x[i]" << x[i].tnt() << endl;
       // cout << "x[i]" << x[i] << endl;
     //   cout << "y[i]" << y[i].convert_int() << endl;
       // cout << "y[i]" << y[i] << endl;
        res[i] = hull(x[i],y[i]);
      //  cout << "res[i]" << res[i].convert_int() << endl;
        //cout << "res[i]" << res[i] << endl;
    }
}






void multMiVi(vector<interval> &y, vector<vector<interval>> &A, vector<interval> &x) {
    for (int i=0 ; i<sysdim ; i++) {
        y[i] = 0.0;
        for (int j=0 ; j<jacdim; j++)
            y[i]=y[i]+A[i][j]*x[j];
    }
}

// produit d'une matrice (n,n) par une (n,p) -> res (n,p)
void multMiMi(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<jacdim; j++) {
            z[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
    }
}



void scaleM(vector<vector<AAF>> &x, double d) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<jacdim; j++) {
            x[i][j]=x[i][j]*d;
        }
    }
}






void addMiMi(vector<vector<AAF>> &x, vector<vector<AAF>> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<jacdim; j++) {
            x[i][j]+=y[i][j];
        }
    }
}

void addViVi(vector<interval> &x, vector<interval> &y) {
    for (int i=0 ; i<sysdim ; i++) {
            x[i]+=y[i];
    }
}


void addViVi(vector<AAF> &x, vector<double> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        x[i]+=y[i];
    }
}

void addViVi(vector<AAF> &x, vector<interval> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        x[i]+=y[i];
    }
}

void intersectViVi(vector<interval> &x, vector<interval> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        x[i]=intersect(x[i],y[i]);
    }
}

void intersectViVi(vector<interval> &x, vector<AAF> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        x[i]=intersect(x[i],y[i].convert_int());
    }
}

void emptyVi(vector<interval> &x) {
     for (int i=0 ; i<sysdim ; i++)
         x[i] = empty();
}

void setVi(vector<interval> &x,vector<AAF> &y) {
    for (int i=0 ; i<sysdim ; i++) 
        x[i]= y[i].convert_int();
}

// with constraints
void setVi(vector<interval> &x,vector<AAF> &y, vector<AAF> &constraints) {
    for (int i=0 ; i<sysdim ; i++)
        x[i]= y[i].convert_int(constraints,sysdim);
}





void setId(vector<vector<AAF>> &J) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int k=0 ; k<sysdim ; k++)
            J[i][k] = 0;
        J[i][i] = 1.0;
    }
}

