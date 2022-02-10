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
    for (int i=0; i<M1.size() ; i++) {
        for (int j=0; j<M1[i].size() ; j++) {
            if (! subseteq(M1[i][j].convert_int(),M2[i][j].convert_int()))
                return false;
        }
    }
    return true;
}


bool subseteq(vector<AAF> &x, vector<AAF> &y) {
    for (int i=0; i<x.size() ; i++) {
        if (! subseteq(x[i].convert_int(),y[i].convert_int()))
            return false;
    }
    return true;
}

bool subseteq(vector<interval> &x, vector<interval> &y) {
    for (int i=0; i<x.size() ; i++) {
        if (! subseteq(x[i],y[i]))
            return false;
    }
    return true;
}

// beware that arguments are modified

vector<AAF> hull(vector<AAF> &x, vector<AAF> &y) {
    vector<AAF> res(sysdim);
    for (int i=0; i<res.size() ; i++) {
     //   x[i].sumup(tol_noise);
      //  y[i].sumup(tol_noise);
        res[i] = hull(x[i],y[i]);
    }
    return res;
}


void hull(vector<interval> &res, vector<interval> &x, vector<interval> &y)
{
    for (int i=0; i<res.size() ; i++)
        res[i] = hull(x[i],y[i]);
}

// beware that arguments are modified

void hull(vector<AAF> &res, vector<AAF> &x, vector<AAF> &y) {
    //cout << "hull" << endl;
    for (int i=0; i<res.size() ; i++) {
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
    for (int i=0 ; i<y.size() ; i++) {
        y[i] = 0.0;
        for (int j=0 ; j<A[i].size(); j++)
            y[i]=y[i]+A[i][j]*x[j];
    }
}

// produit d'une matrice (n,n) par une (n,p) -> res (n,p)
/*void multMiMi(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<jacdim; j++) {
            z[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
    }
}*/

// produit d'une matrice (n,n) par une (n,p) -> res (n,p)
void multMiMi(vector<vector<AAF>> &z, const vector<vector<AAF>> x, const vector<vector<AAF>> y) {
    for (int i=0 ; i<z.size() ; i++) {
        for (int j=0 ; j<z[i].size(); j++) {
            z[i][j]=0;
            for (int k=0 ; k<x[i].size() ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
    }
}

// produit d'une matrice (n,n) par une (n,p) -> res (n,p)
void multMiMi(vector<vector<interval>> &z, const vector<vector<double>> x, const vector<vector<interval>> y) {
    for (int i=0 ; i<z.size() ; i++) {
        for (int j=0 ; j<z[i].size(); j++) {
            z[i][j]=0;
            for (int k=0 ; k<x[i].size() ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
    }
}

// produit d'une matrice (n,n) par une (n,p) -> res (n,p)
void multMiMi(vector<vector<AAF>> &z, const  vector<vector<double>> x, const vector<vector<AAF>> y) {
    for (int i=0 ; i<z.size() ; i++) {
        for (int j=0 ; j<z[i].size(); j++) {
            z[i][j]=0;
            for (int k=0 ; k<x[i].size() ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
    }
}

// produit d'une matrice (n,n) par une (n,p) -> res (n,p)
void multMiMi(vector<vector<double>> &z, const vector<vector<double>> x, const vector<vector<double>> y) {
    for (int i=0 ; i<z.size() ; i++) {
        for (int j=0 ; j<z[i].size(); j++) {
            z[i][j]=0;
            for (int k=0 ; k<x[i].size() ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
    }
}


ostream& operator<<(ostream& os, const vector<vector<interval>> &M)
{
    for (int i=0 ; i<M.size() ; i++) {
        for (int j=0 ; j<M[i].size(); j++)
            os << "M["<<i<<"]["<<j<<"]="<<M[i][j]<<" ";
        os << endl;
    }
    
    return os;
}

ostream& operator<<(ostream& os, const vector<interval> &z)
{
    for (int i=0 ; i<z.size() ; i++)
            os << "z["<<i<<"]="<<z[i]<<" ";
        os << endl;
    
    return os;
}



// product (jacdim,jacdim) times (jacdim,jacdim)
// but only the first sysdim lines of x and z  are relevant
void multJacfzJaczz0(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &y) {
    // sysdim is on purpose (only the first sysdim lines of x and z  are relevant)
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<sysdim; j++) {
            z[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
        for (int j=sysdim ; j<jacdim; j++) {
            z[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                z[i][j] += x[i][k]*y[k][j];
            z[i][j] += x[i][j];
        }
    }
}


// product (jacdim,jacdim) times (jacdim,jacdim)
// but only the first sysdim lines of x and z  are relevant
// incliudes the case when u can be a linear function of z (aux2 contains partial u / partial z0)
void multJacfzuJaczz0Jacuz0(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &aux2, vector<vector<AAF>> &y,double offset) {
    // sysdim is on purpose (only the first sysdim lines of x and z  are relevant)
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<sysdim; j++) {
            z[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                z[i][j] += x[i][k]*y[k][j];
        }
        for (int j=sysdim ; j<jacdim; j++) {
            z[i][j]=0;
            for (int k=0 ; k<sysdim ; k++)
                z[i][j] += x[i][k]*y[k][j];
            z[i][j] += x[i][j];
        }
        for (int j=0 ; j<sysdim ; j++)
            for (int k=0 ; k<inputsdim ; k++)
                z[i][j] += x[i][sysdim+index_param_inv[k]+floor(offset*nb_inputs[k])]*aux2[k][j];
    }
}


void scaleM(vector<vector<AAF>> &x, double d) {
    for (int i=0 ; i<x.size() ; i++) {
        for (int j=0 ; j<x[i].size(); j++) {
            x[i][j]=x[i][j]*d;
        }
    }
}

void scaleJacfz(vector<vector<AAF>> &x, double d) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<x[i].size(); j++) {
            x[i][j]=x[i][j]*d;
        }
    }
}




void addMiMi(vector<vector<AAF>> &x, vector<vector<AAF>> &y) {
    for (int i=0 ; i<x.size() ; i++) {
        for (int j=0 ; j<x[i].size(); j++) {
            x[i][j]+=y[i][j];
        }
    }
}

void addJacfzJacfz(vector<vector<AAF>> &x, vector<vector<AAF>> &y) {
    for (int i=0 ; i<sysdim ; i++) {
        for (int j=0 ; j<x[i].size(); j++) {
            x[i][j]+=y[i][j];
        }
    }
}

void scalarproduct(vector<interval> &z, vector<interval> &x, vector<interval> &y) {
    for (int i=0 ; i<x.size() ; i++) {
        z[i]=x[i]*y[i];
    }
}


void addViVi(vector<interval> &x, vector<interval> &y) {
    for (int i=0 ; i<x.size() ; i++) {
            x[i]+=y[i];
    }
}


void addViVi(vector<AAF> &x, vector<double> &y) {
    for (int i=0 ; i<x.size() ; i++) {
        x[i]+=y[i];
    }
}

void addViVi(vector<AAF> &x, vector<interval> &y) {
    for (int i=0 ; i<x.size() ; i++) {
        x[i]+=y[i];
    }
}

void intersectViVi(vector<interval> &x, vector<interval> &y) {
    for (int i=0 ; i<x.size() ; i++) {
        x[i]=intersect(x[i],y[i]);
    }
}

void intersectViVi(vector<interval> &x, vector<AAF> &y) {
    for (int i=0 ; i<x.size() ; i++) {
        x[i]=intersect(x[i],y[i].convert_int());
    }
}

void emptyVi(vector<interval> &x) {
     for (int i=0 ; i<x.size() ; i++)
         x[i] = empty();
}

void setVi(vector<interval> &x,vector<AAF> &y) {
    for (int i=0 ; i<x.size() ; i++)
        x[i]= y[i].convert_int();
}

// with constraints
void setVi(vector<interval> &x,vector<AAF> &y, vector<AAF> &constraints) {
    for (int i=0 ; i<x.size() ; i++)
        x[i]= y[i].convert_int(constraints,sysdim);
}





void setId(vector<vector<AAF>> &J) {
    for (int i=0 ; i<J.size() ; i++) {
        for (int k=0 ; k<J[i].size() ; k++)
            J[i][k] = 0;
        J[i][i] = 1.0;
    }
}

 // resulting quadrilatere A * inner is an inner approximation of the range of f
void compute_print_skewbox(interval &temp_inner_x, interval &temp_inner_y, vector<vector<double>> &A, int varx, int vary, const char *approxtype)
{
    vector<double> output_skewedbox(8);
    
    output_skewedbox[0] = inf(temp_inner_x)*A[varx][varx] + inf(temp_inner_y)*A[varx][vary];
    output_skewedbox[1] = inf(temp_inner_x)*A[vary][varx] + inf(temp_inner_y)*A[vary][vary];
    output_skewedbox[2] = inf(temp_inner_x)*A[varx][varx] + sup(temp_inner_y)*A[varx][vary];
    output_skewedbox[3] = inf(temp_inner_x)*A[vary][varx] + sup(temp_inner_y)*A[vary][vary];
    output_skewedbox[4] = sup(temp_inner_x)*A[varx][varx] + sup(temp_inner_y)*A[varx][vary];
    output_skewedbox[5] = sup(temp_inner_x)*A[vary][varx] + sup(temp_inner_y)*A[vary][vary];
    output_skewedbox[6] = sup(temp_inner_x)*A[varx][varx] + inf(temp_inner_y)*A[varx][vary];
    output_skewedbox[7] = sup(temp_inner_x)*A[vary][varx] + inf(temp_inner_y)*A[vary][vary];
    
    out_approx << YAML::Key << approxtype;
    out_approx << YAML::Value << output_skewedbox;
}

vector<vector<double>> compute_skewbox_3d(vector<interval> &temp_inner, vector<vector<double>> &A, int varx, int vary, int varz)
{
    vector<vector<double>> output_skewedbox(8);
    for (int i=0; i<8; i++)
        output_skewedbox[i] = vector<double>(3);
    
    output_skewedbox[0][0] = inf(temp_inner[varx])*A[varx][varx] + inf(temp_inner[vary])*A[varx][vary] + inf(temp_inner[varz])*A[varx][varz];
    output_skewedbox[0][1] = inf(temp_inner[varx])*A[vary][varx] + inf(temp_inner[vary])*A[vary][vary] + inf(temp_inner[varz])*A[vary][varz];
    output_skewedbox[0][2] = inf(temp_inner[varx])*A[varz][varx] + inf(temp_inner[vary])*A[varz][vary] + inf(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[1][0] = inf(temp_inner[varx])*A[varx][varx] + sup(temp_inner[vary])*A[varx][vary] + inf(temp_inner[varz])*A[varx][varz];
    output_skewedbox[1][1] = inf(temp_inner[varx])*A[vary][varx] + sup(temp_inner[vary])*A[vary][vary] + inf(temp_inner[varz])*A[vary][varz];
    output_skewedbox[1][2] = inf(temp_inner[varx])*A[varz][varx] + sup(temp_inner[vary])*A[varz][vary] + inf(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[2][0] = sup(temp_inner[varx])*A[varx][varx] + sup(temp_inner[vary])*A[varx][vary] + inf(temp_inner[varz])*A[varx][varz];
    output_skewedbox[2][1] = sup(temp_inner[varx])*A[vary][varx] + sup(temp_inner[vary])*A[vary][vary] + inf(temp_inner[varz])*A[vary][varz];
    output_skewedbox[2][2] = sup(temp_inner[varx])*A[varz][varx] + sup(temp_inner[vary])*A[varz][vary] + inf(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[3][0] = sup(temp_inner[varx])*A[varx][varx] + inf(temp_inner[vary])*A[varx][vary] + inf(temp_inner[varz])*A[varx][varz];
    output_skewedbox[3][1] = sup(temp_inner[varx])*A[vary][varx] + inf(temp_inner[vary])*A[vary][vary] + inf(temp_inner[varz])*A[vary][varz];
    output_skewedbox[3][2] = sup(temp_inner[varx])*A[varz][varx] + inf(temp_inner[vary])*A[varz][vary] + inf(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[4][0] = inf(temp_inner[varx])*A[varx][varx] + inf(temp_inner[vary])*A[varx][vary] + sup(temp_inner[varz])*A[varx][varz];
    output_skewedbox[4][1] = inf(temp_inner[varx])*A[vary][varx] + inf(temp_inner[vary])*A[vary][vary] + sup(temp_inner[varz])*A[vary][varz];
    output_skewedbox[4][2] = inf(temp_inner[varx])*A[varz][varx] + inf(temp_inner[vary])*A[varz][vary] + sup(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[5][0] = inf(temp_inner[varx])*A[varx][varx] + sup(temp_inner[vary])*A[varx][vary] + sup(temp_inner[varz])*A[varx][varz];
    output_skewedbox[5][1] = inf(temp_inner[varx])*A[vary][varx] + sup(temp_inner[vary])*A[vary][vary] + sup(temp_inner[varz])*A[vary][varz];
    output_skewedbox[5][2] = inf(temp_inner[varx])*A[varz][varx] + sup(temp_inner[vary])*A[varz][vary] + sup(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[6][0] = sup(temp_inner[varx])*A[varx][varx] + sup(temp_inner[vary])*A[varx][vary] + sup(temp_inner[varz])*A[varx][varz];
    output_skewedbox[6][1] = sup(temp_inner[varx])*A[vary][varx] + sup(temp_inner[vary])*A[vary][vary] + sup(temp_inner[varz])*A[vary][varz];
    output_skewedbox[6][2] = sup(temp_inner[varx])*A[varz][varx] + sup(temp_inner[vary])*A[varz][vary] + sup(temp_inner[varz])*A[varz][varz];
    
    output_skewedbox[7][0] = sup(temp_inner[varx])*A[varx][varx] + inf(temp_inner[vary])*A[varx][vary] + sup(temp_inner[varz])*A[varx][varz];
    output_skewedbox[7][1] = sup(temp_inner[varx])*A[vary][varx] + inf(temp_inner[vary])*A[vary][vary] + sup(temp_inner[varz])*A[vary][varz];
    output_skewedbox[7][2] = sup(temp_inner[varx])*A[varz][varx] + inf(temp_inner[vary])*A[varz][vary] + sup(temp_inner[varz])*A[varz][varz];
    
    return output_skewedbox;
}


// builds conditionner for skewbox computation
// A is center of Jaux on components i and k (otherwise diagonal), C is inverse of A
void build_2dpreconditionner(vector<vector<double>> &A, vector<vector<double>> &C, vector<vector<interval>> Jaux, int i, int k)
{
    double determinant;
    
    for (int p=0 ; p<sysdim; p++) {
        for (int q=0 ; q<sysdim; q++) {
            A[p][q] = 0.0;
            C[p][q] = 0.0;
        }
        A[p][p] = 1.0;
        C[p][p] = 1.0;
    }
    
    A[i][i] = mid(Jaux[i][i]);
    A[k][k] = mid(Jaux[k][k]);
    A[i][k] = mid(Jaux[i][k]);
    A[k][i] = mid(Jaux[k][i]);
    
    // C is inverse of A
    determinant = 1.0/(A[i][i]*A[k][k]-A[i][k]*A[k][i]);
    C[i][i] = determinant*A[k][k];
    C[i][k] = - determinant*A[i][k];
    C[k][i] = - determinant*A[k][i];
    C[k][k] = determinant*A[i][i];
}
