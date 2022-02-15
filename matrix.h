/* ============================================================================
 File   : matrix.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file defines rather standard vector/matrix operations over intervals and affine forms (AAF)
 ============================================================================ */

#ifndef MATRIX_H
#define MATRIX_H

#include "filib_interval.h"
#include "ode_def.h"
#include "fadbad_aa.h"

using namespace std;

bool subseteq (vector<vector<AAF>> &M1, vector<vector<AAF>> &M2);
bool subseteq(vector<interval> &x, vector<interval> &y);
bool subseteq(vector<AAF> &x, vector<AAF> &y);
vector<AAF> hull(vector<AAF> &x, vector<AAF> &y);
void hull(vector<interval> &res, vector<interval> &x, vector<interval> &y);
void hull(vector<AAF> &res, vector<AAF> &x, vector<AAF> &y);

void multMiMi(vector<vector<AAF>> &z, const vector<vector<AAF>> x, const vector<vector<AAF>> y);
void multMiMi(vector<vector<interval>> &z, const vector<vector<double>> x, const vector<vector<interval>> y);
void multMiMi(vector<vector<AAF>> &z, const vector<vector<double>> x, const vector<vector<AAF>> y);
void multMiMi(vector<vector<double>> &z, const vector<vector<double>> x, const vector<vector<double>> y);
void multJacfzJaczz0(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &y);
void multJacfzuJaczz0Jacuz0(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &aux2, vector<vector<AAF>> &y,double offset);

void multMiVi(vector<interval> &y, vector<vector<interval>> &A, vector<interval> &x) ;

ostream& operator<<(ostream& os, const vector<vector<interval>> &M);
ostream& operator<<(ostream& os, const vector<interval> &z);

void scaleM(vector<vector<AAF>> &x, double d);
void scaleJacfz(vector<vector<AAF>> &x, double d);

void addMiMi(vector<vector<AAF>> &x, vector<vector<AAF>> &y);
void addJacfzJacfz(vector<vector<AAF>> &x, vector<vector<AAF>> &y);

void scalarproduct(vector<interval> &z, vector<interval> &x, vector<interval> &y);

void addViVi(vector<interval> &x, vector<interval> &y);
void addViVi(vector<AAF> &x, vector<double> &y);
void addViVi(vector<AAF> &x, vector<interval> &y);
void intersectViVi(vector<interval> &x, vector<interval> &y);
void intersectViVi(vector<interval> &x, vector<AAF> &y);

void emptyVi(vector<interval> &x);
void setVi(vector<interval> &x,vector<AAF> &y);
void setVi(vector<interval> &x,vector<AAF> &y, vector<AAF> &constraints);

void setId(vector<vector<AAF>> &J);

void compute_print_skewbox(interval &temp_inner_x, interval &temp_inner_y, vector<vector<double>> &A, int varx, int vary, const char *approxtype);
vector<vector<double>> compute_skewbox_3d(vector<interval> &temp_inner, vector<vector<double>> &A, int varx, int vary, int varz);

// builds conditionner for skewbox computation
// A is center of Jaux on components i and k (otherwise diagonal), C is inverse of A
void build_2dpreconditionner(vector<vector<double>> &A, vector<vector<double>> &C, vector<vector<interval>> Jaux, int i, int k);
void build_3dpreconditionner(vector<vector<double>> &A, vector<vector<double>> &C, vector<vector<interval>> J, int i, int k, int l);

#endif



