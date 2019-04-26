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

void multMiMi(vector<vector<AAF>> &z, vector<vector<AAF>> &x, vector<vector<AAF>> &y);

void multMiVi(vector<interval> &y, vector<vector<interval>> &A, vector<interval> &x) ;

void scaleM(vector<vector<AAF>> &x, double d);

void addMiMi(vector<vector<AAF>> &x, vector<vector<AAF>> &y);
void addViVi(vector<interval> &x, vector<interval> &y);
void addViVi(vector<AAF> &x, vector<double> &y);
void addViVi(vector<AAF> &x, vector<interval> &y);
void intersectViVi(vector<interval> &x, vector<interval> &y);
void intersectViVi(vector<interval> &x, vector<AAF> &y);

void emptyVi(vector<interval> &x);
void setVi(vector<interval> &x,vector<AAF> &y);
void setVi(vector<interval> &x,vector<AAF> &y, vector<AAF> &constraints);

void setId(vector<vector<AAF>> &J);

#endif



