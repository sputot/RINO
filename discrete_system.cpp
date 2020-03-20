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


vector<vector<vector<interval>>> constr_eps;  // constraints on noise symbols when partitioning 2D region in 4
vector<vector<vector<interval>>> eps_loc;     // consequence on [x]-x0  when partitioning 2D region in 4

vector<vector<vector<vector<interval>>>> constr_eps_discr;  // same but with additional discretization in each direction
vector<vector<vector<vector<interval>>>> eps_loc_discr;     // consequence on [x]-x0  when partitioning 2D region in 4
vector<vector<vector<vector<double>>>> extremity_eps_loc_discr;

int nb_discr, nb_discr1, nb_discr2;

void range_discrete_system(void) {
    
    if (syschoice == 1) {
        jacdim = 1;  // number of input components
        sysdim = 1; // number of output components
    }
    else if (syschoice == 2) {
        jacdim = 2;
        sysdim = 1;
    }
    else if (syschoice == 3) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 4) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 5) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 6) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 7) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 8) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 9) {
        jacdim = 2;
        sysdim = 2;
    }
    else if (syschoice == 10) {
        jacdim = 3;
        sysdim = 3;
    }
    else if (syschoice == 11) {
        jacdim = 3;
        sysdim = 3;
    }
    initial_values = vector<AAF>(jacdim);
        
    if (syschoice == 1) {
        initial_values[0] = interval(2,3);
    }
    else if (syschoice == 2) {
        initial_values[0] = interval(2,3);
        initial_values[1] = interval(2,3);
    }
    else if (syschoice == 3) { // example 3.5 Goldztejn
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
    }
    else if (syschoice == 4) { // example 3 CDC
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
    }
    else if (syschoice == 5) { // example 5.1 Goldstzjen
        initial_values[0] = interval(0.99,1.01);
        initial_values[1] = interval(0.99,1.01);
    }
    else if (syschoice == 6) { //
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
    }
    else if (syschoice == 7) { //
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
    }
    else if (syschoice == 8) { //
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
    }
    else if (syschoice == 9) { //
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
    }
    else if (syschoice == 10) {
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
        initial_values[2] = interval(0.9,1.1);
    }
    else if (syschoice == 11) {
        initial_values[0] = interval(0.9,1.1);
        initial_values[1] = interval(0.9,1.1);
        initial_values[2] = interval(0.9,1.1);
    }
    else if (syschoice == 12) {
        initial_values[0] = interval(-1.0,1.0);
        initial_values[1] = interval(-1.0,1.0);
        initial_values[2] = interval(-1.0,1.0);
    }
    
    nb_discr = 10;
    nb_discr1 = nb_discr;
    nb_discr2 = nb_discr;
    
    
    vector<F<AAF>> x(jacdim);
    vector<vector<interval>> Jacf(sysdim);
    vector<vector<AAF>> JacAff(sysdim);
    
    for (int i=0; i < sysdim ; i++) {
        Jacf[i] = vector<interval>(jacdim);
        JacAff[i] = vector<AAF>(jacdim);
    }
    vector<F<AAF>> z(sysdim), z1(sysdim), z2(sysdim);
    vector<interval> x0(jacdim);    // center
    
    system("rm -r output");
    system("mkdir output");
    
    DiscreteFunc f;
    
    
    
    // estimation of the range of f
    vector<double> input(jacdim);
    vector<double> output(sysdim);
    vector<double> max_output(sysdim);
    vector<double> min_output(sysdim);
    int discr = 100;
    for (int i=0; i < jacdim ; i++)
        input[i] = initial_values[i].convert_int().mid();
    for (int i=0; i < sysdim ; i++) {
        max_output = f(input);
        min_output = f(input);
    }
    ofstream outFile_xi;
    stringstream file_name;
    file_name.str("");
    file_name << "output/xi.out";
    outFile_xi.open(file_name.str().c_str());
    for (int i1=0; i1 <= discr ; i1++) {
        input[0] = initial_values[0].convert_int().inf() + (2.0*i1*initial_values[0].convert_int().rad())/discr;
        for (int i2=0; i2 <= discr ; i2++) {
            input[1] = initial_values[1].convert_int().inf() + (2.0*i2*initial_values[1].convert_int().rad())/discr;
            output = f(input);
            for (int i=0; i < sysdim ; i++)
                outFile_xi << output[i] <<  "\t";
            outFile_xi << endl;
            for (int i=0; i < sysdim ; i++) {
                if (output[i] < min_output[i])
                    min_output[i] = output[i];
                else if (output[i] > max_output[i])
                    max_output[i] = output[i];
            }
        }
    }
    outFile_xi.close();
    cout << "Estimated range of f:"<<endl;
    for (int i=0; i < sysdim ; i++)
        cout << "z["<<i << "]=[" << min_output[i] << ", " << max_output[i] <<"]  ";
    cout << endl;
    
    eps = vector<interval>(jacdim);
    
    
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]); //+(eps[i].sup()-mid(eps[i]))/2.0;
        eps[i] = eps[i] - x0[i];
    }
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    
    z = f(x);
    
  
    vector<ofstream> outFile_outer_direct(sysdim);
    
    for (int i=0; i < jacdim ; i++) {
        cout << "outer range of f(x) direct evaluation:" << z[i].x().convert_int() << endl;
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_direct.out";
        outFile_outer_direct[i].open(file_name.str().c_str());
        outFile_outer_direct[i] << inf(z[i].x().convert_int()) << "\t" << sup(z[i].x().convert_int()) << endl;
        outFile_outer_direct[i].close();
    }
    
    // begin for subdivisions
    constr_eps = vector<vector<vector<interval>>>(jacdim);
    eps_loc = vector<vector<vector<interval>>>(jacdim);
    for (int i=0; i < jacdim ; i++) {
        constr_eps[i] = vector<vector<interval>>(jacdim);
        eps_loc[i] = vector<vector<interval>>(jacdim);
        for (int j=0; j < jacdim ; j++) {
            constr_eps[i][j] = vector<interval>(jacdim);
            eps_loc[i][j] = vector<interval>(jacdim);
            for (int k=0; k < jacdim ; k++) {
                constr_eps[i][j][k] = interval(-1,1);
                eps_loc[i][j][k] = eps[k];
            }
        }
    }
    
    constr_eps[0][0][0] = interval(-1,0);
    constr_eps[0][0][1] = interval(-1,0);
    constr_eps[0][1][0] = interval(-1,0);
    constr_eps[0][1][1] = interval(0,1);
    constr_eps[1][0][0] = interval(0,1);
    constr_eps[1][0][1] = interval(-1,0);
    constr_eps[1][1][0] = interval(0,1);
    constr_eps[1][1][1] = interval(0,1);
    
    for (int i=0; i < jacdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            for (int k=0; k < jacdim ; k++)
                eps_loc[i][j][k] = eps[k].mid() + eps[k].rad() * constr_eps[i][j][k];
        }
    }
    // end subdiv
    
    // additional discretization of each component
    constr_eps_discr = vector<vector<vector<vector<interval>>>>(jacdim);
    eps_loc_discr = vector<vector<vector<vector<interval>>>>(jacdim);
    extremity_eps_loc_discr = vector<vector<vector<vector<double>>>>(jacdim);
    for (int i=0; i < jacdim ; i++) {
        constr_eps_discr[i] = vector<vector<vector<interval>>>(jacdim);
        eps_loc_discr[i] = vector<vector<vector<interval>>>(jacdim);
        extremity_eps_loc_discr[i] = vector<vector<vector<double>>>(jacdim);
        for (int j=0; j < jacdim ; j++) {
            constr_eps_discr[i][j] = vector<vector<interval>>(nb_discr);
            eps_loc_discr[i][j] = vector<vector<interval>>(nb_discr);
            extremity_eps_loc_discr[i][j] = vector<vector<double>>(nb_discr);
            for (int m=0; m < nb_discr ; m++) {
                constr_eps_discr[i][j][m] = vector<interval>(jacdim);
                eps_loc_discr[i][j][m] = vector<interval>(jacdim);
                extremity_eps_loc_discr[i][j][m] = vector<double>(jacdim);
                for (int k=0; k < jacdim ; k++) {
                    constr_eps_discr[i][j][m][k] = interval(-1,1);
                    eps_loc_discr[i][j][m][k] = eps[k];
                   // completer pour extremity_eps_loc_discr ?
                }
            }
        }
    }
    
    for (int m=0; m < nb_discr ; m++) {
        constr_eps_discr[1][1][m][0] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        constr_eps_discr[1][1][m][1] = constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][0][m][0] = -constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][0][m][1] = -constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][1][m][0] = -constr_eps_discr[1][1][m][0];
        constr_eps_discr[0][1][m][1] = constr_eps_discr[1][1][m][0];
        constr_eps_discr[1][0][m][0] = constr_eps_discr[1][1][m][0];
        constr_eps_discr[1][0][m][1] = -constr_eps_discr[1][1][m][0];
        extremity_eps_loc_discr[1][1][m][0] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[1][1][m][1] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][0][m][0] = -(m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][0][m][1] = -(m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][1][m][0] = -(m+1.0)/nb_discr;
        extremity_eps_loc_discr[0][1][m][1] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[1][0][m][0] = (m+1.0)/nb_discr;
        extremity_eps_loc_discr[1][0][m][1] = -(m+1.0)/nb_discr;
    }
    
    
    for (int i=0; i < jacdim ; i++) {
        for (int j=0; j < jacdim ; j++) {
            for (int m=0; m < nb_discr ; m++) {
                for (int k=0; k < jacdim ; k++) {
                    eps_loc_discr[i][j][m][k] = eps[k].mid() + eps[k].rad() * constr_eps_discr[i][j][m][k];
                    extremity_eps_loc_discr[i][j][m][k] = eps[k].mid() + eps[k].rad() * extremity_eps_loc_discr[i][j][m][k];
              //      cout << "constr_eps_discr[i][j][m][k]=" << constr_eps_discr[i][j][m][k] << endl;
               //     cout << "eps_loc_discr[i][j][m][k]=" << eps_loc_discr[i][j][m][k] << endl;
                }
            }
        }
    }
        
    
    
    // end discretization
    
    
   
    cout << endl;
    
    for (int i=0; i < sysdim ; i++)
        for (int j=0; j < jacdim ; j++) {
            Jacf[i][j] = z[i].d(j).convert_int();
            JacAff[i][j] = z[i].d(j);
        }
    cout << "Jacf evaluated on all [x]" << endl;
    cout << Jacf;
    
    
   
    vector<interval> z0 = f(x0);
 
    evaluate_projections(z0, Jacf);
    evaluate_projections_subdiv(z0, JacAff);
  evaluate_projections_subdiv_discretize(z0, JacAff);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf,0,1);
        joint_ranges_subdiv(z0,JacAff,0,1);
        joint_ranges_subdiv_discretize(z0,JacAff,0,1);
        preconditioned_joint_ranges(z0,Jacf,0,1);
        preconditioned_joint_ranges_subdiv(z0,JacAff,0,1);
        preconditioned_joint_ranges_subdiv_discretize(z0,JacAff,0,1);
    }
  
    
    
    /*********************************/
/*
    cout << endl;
   // cout << "Improved mean-value:" << endl;
    vector<vector<interval>> Jacf_refined(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jacf_refined[i] = vector<interval>(jacdim);
    
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    x[jacdim-1] = mid(initial_values[jacdim-1].convert_int());
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    z1 = f(x);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim-1 ; j++) { // first column evaluated on ([x1],x2_0)
            Jacf_refined[i][j] = z1[i].d(j).convert_int();
            JacAff[i][j] = z1[i].d(j);
        }
        for (int j=jacdim-1; j < jacdim ; j++) { // rest evaluated on box
            Jacf_refined[i][j] = z[i].d(j).convert_int();
            JacAff[i][j] = z[i].d(j);
        }
    }
    cout << "Jacf: 1st column evaluated on ([x1],x2_0), 2nd column on [x]" << endl;
    cout << Jacf_refined;

    evaluate_projections(z0, Jacf_refined);
    evaluate_projections_subdiv(z0, JacAff);
    evaluate_projections_subdiv_discretize(z0, JacAff);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        joint_ranges_subdiv(z0,JacAff,0,1);
        joint_ranges_subdiv_discretize(z0,JacAff,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges_subdiv(z0,JacAff,0,1);
        preconditioned_joint_ranges_subdiv_discretize(z0,JacAff,0,1);
    }
 */
 
    
    /*
    cout << endl;
   // cout << "Improved mean-value:" << endl;
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    x[0] = mid(initial_values[0].convert_int());
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    z2 = f(x);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim-1 ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z2[i].d(j).convert_int();
    }
    cout << "Jacf: 1st column evaluated on [x], 2nd column on (x0_0,[x1])" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
  */

    /*
    cout << endl;
 //   cout << "Improved mean-value:" << endl;
    
    for (int i=0; i < sysdim-1 ; i++) {
        for (int j=0; j < jacdim-1 ; j++)
            Jacf_refined[i][j] = z1[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
    }
    for (int i=sysdim-1; i < sysdim ; i++) {
        for (int j=0; j < jacdim ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z2[i].d(j).convert_int();
    }
    cout << "Jacf: f1: 1st column evaluated on ([x0],x1_0), 2nd column on [x], f2: 1st column evaluated on [x], 2nd column on (x0_0,[x1])" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
  */
   
}

void evaluate_projections(vector<interval> &z0,  vector<vector<interval>> &Jacf) {
    interval inner_impro;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = evaluate_outerrange_x(z0, Jacf, i);
        z_inner[i] = evaluate_innerrange_x(z0, Jacf, true, aux, i);
    }
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value: ";
    cout << z_inner;
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}

// index1 is the first index that we subdivide, n1 the number of associated subdivisions
// evaluate projections on the 4 quadrant and join them
void evaluate_projections_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = evaluate_outerrange_x_subdiv(z0, JacAff, i);
        z_inner[i] = z0[i];
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++) {
                temp_inner = evaluate_innerrange_x_subdiv(z0, JacAff, true, aux, i, i1, i2); // maximal inner approx on quadrant (i1,i2)
                // on maximal inner approx we can do interval hull
                z_inner[i] = interval_hull(z_inner[i],temp_inner);
            }
        }
    }
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value with 2 subdiv each input:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value with 2 subdiv each input: ";
    cout << z_inner;
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}

// index1 is the first index that we subdivide, n1 the number of associated subdivisions
// evaluate projections on the 4 quadrant and join them
void evaluate_projections_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff) {
    interval inner_impro, temp_inner, temp_outer;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    vector<int> aux;
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = evaluate_outerrange_x_subdiv_discretize(z0, JacAff, i);
        z_inner[i] = z0[i];
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++) {
                temp_inner = evaluate_innerrange_x_subdiv_discretize(z0, JacAff, true, aux, i, i1, i2); // maximal inner approx on quadrant (i1,i2)
                // on maximal inner approx we can do interval hull
                z_inner[i] = interval_hull(z_inner[i],temp_inner);
            }
        }
    }
    
    
    vector<ofstream> outFile_outer_mean_value(sysdim);
    vector<ofstream> outFile_inner(sysdim);
    stringstream file_name;
    
    cout << "Outer range of f(x) by mean-value with 2 subdiv each input + discretization:               ";
    cout << z_outer;
    cout << "Projection of inner range of f(x) by mean-value with 2 subdiv each input + discretization: ";
    cout << z_inner;
    
    
    for (int i=0; i < sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_mean_value.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str(),fstream::app);
        outFile_outer_mean_value[i] << inf(z_outer[i]) << "\t" << sup(z_outer[i]) << endl;
        outFile_outer_mean_value[i].close();
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_outer_mean_value[i].open(file_name.str().c_str());
        outFile_outer_mean_value[i] << inf(z_inner[i]) << "\t" << sup(z_inner[i]) << endl;
        outFile_outer_mean_value[i].close();
    }
}









interval evaluate_outerrange_x(vector<interval> &z0,  vector<vector<interval>> &Jacf, int i) {
    interval z_outer = z0[i];
    
    for (int j=0; j < jacdim ; j++) {
        z_outer += Jacf[i][j]*eps[j];
    }
    
    return z_outer;
}

interval evaluate_outerrange_x_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp, z_outer = z0[i];
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
            z_temp = z0[i];
            for (int j=0; j < jacdim ; j++) {
                z_temp += JacAff[i][j].convert_int(constr_eps[i1][i2],jacdim) * eps_loc[i1][i2][j];
            }
            
            z_outer = interval_hull(z_outer,z_temp);
        }
        
    }
    
    return z_outer;
}

// subdivise n times in the 2 inputs
// I think trouble is I progress in crabe (synchronously on the 2 components)
/*
interval evaluate_outerrange_x_subdiv_discretize_old(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, z_outer = z0[i];
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
            z_temp_partial = z0[i];
            // evaluated at the extremity on the n-1 first subdivisions
            for (int m=0; m < nb_discr-1 ; m++) {
                z_temp = z_temp_partial;
                for (int j=0; j < jacdim ; j++) {
                    if (m>0)
                    {
                        z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
                        z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (extremity_eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
                    }
                    else
                    {
                        z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * eps_loc_discr[i1][i2][m][j];
                        z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * extremity_eps_loc_discr[i1][i2][m][j];
                    }
                }
                // hull on all subdivisions
                z_outer = interval_hull(z_outer,z_temp);
            }
            // evaluated on interval on the last subdivision
            z_temp = z_temp_partial;
            for (int j=0; j < jacdim ; j++) {
                z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][nb_discr-1],jacdim) * (eps_loc_discr[i1][i2][nb_discr-1][j]-extremity_eps_loc_discr[i1][i2][nb_discr-2][j]);
            }
            
            z_outer = interval_hull(z_outer,z_temp);
        }
    }
    return z_outer;
}
 */

// the vector of constraints on eps
vector<interval> constraint_eps(int i1, int i2, int m1, int m2) {
    vector<interval> res(jacdim);
    
    for (int j=0; j < jacdim ; j++)
        res[j] = interval(-1,1);
    
    if (constr_eps[i1][i2][0].sup() > 0)
        res[0] = interval(((double)m1)/nb_discr1,(m1+1.0)/nb_discr1);
    else
        res[0] = -interval(((double)m1)/nb_discr1,(m1+1.0)/nb_discr1);
    
    if (constr_eps[i1][i2][1].sup() > 0)
        res[1] = interval(((double)m2)/nb_discr2,(m2+1.0)/nb_discr2);
    else
        res[1] = -interval(((double)m2)/nb_discr2,(m2+1.0)/nb_discr2);
    
    return res;
}

//interval concretize_eps(int i1, int i2, int m1, int m2, int k) {
//    return eps[k].mid() + eps[k].rad() * constr_eps_discr[i][j][m][k];
//}

double extremity_eps_max(int i1, int i2, int m1, int m2, int k) {
    if (constr_eps[i1][i2][k].sup() > 0) {
        if (k == 0)
            return eps[k].mid() + eps[k].rad() * (m1+1.0)/nb_discr1;
        else
            return eps[k].mid() + eps[k].rad() * (m2+1.0)/nb_discr2;
    }
    else {
        if (k == 0)
            return eps[k].mid() - eps[k].rad() * (m1+1.0)/nb_discr1;
        else
            return eps[k].mid() - eps[k].rad() * (m2+1.0)/nb_discr2;
    }
}

double extremity_eps_min(int i1, int i2, int m1, int m2, int k) {
    if (constr_eps[i1][i2][k].sup() > 0) {
        if (k == 0)
            return eps[k].mid() + eps[k].rad() * (double(m1))/nb_discr1;
        else
            return eps[k].mid() + eps[k].rad() * (double(m2))/nb_discr2;
    }
    else {
        if (k == 0)
            return eps[k].mid() - eps[k].rad() * (double(m1))/nb_discr1;
        else
            return eps[k].mid() - eps[k].rad() * (double(m2))/nb_discr2;
    }
}

/*
// subdivise n times in the 2 inputs
// I think trouble is I progress in crabe (synchronously on the 2 components)
interval evaluate_outerrange_x_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, loc_eps, extremity_loc_eps_min, extremity_loc_eps_max, z_outer = z0[i];
    vector<interval> c_eps;
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
         //   z_temp_partial = z0[i];
            z_temp = z0[i]; // z_temp_partial;
            for (int m1=0; m1 < nb_discr1 ; m1++) {
                for (int m2=0; m2 < nb_discr2 ; m2++) {
                    c_eps = constraint_eps(i1,i2,m1,m2);
                    cout << "m1=" << m1 << " m2=" << m2 << " c_eps=" << c_eps;
                    for (int j=0; j < jacdim ; j++) {
                        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
                       // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                        extremity_loc_eps_min = extremity_eps_min(i1,i2,m1,m2,j);
                        z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
                        cout << "j=" << j << " z_temp=" << z_temp << " JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    }
                    // hull on all subdivisions
                    z_outer = interval_hull(z_outer,z_temp);
                    cout << "z_outer=" << z_outer << endl;
                }
            }
        }
    }
    return z_outer;
}
 */

// subdivise n times in the 2 inputs
interval evaluate_outerrange_x_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, loc_eps, extremity_loc_eps_min, extremity_loc_eps_max, z_outer = z0[i];
    vector<interval> c_eps;
    
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++) {
            //   z_temp_partial = z0[i];
         //   if (i == 0)
          //      cout << "i1=" << i1 << " i2=" << i2 << endl;
            z_temp = z0[i]; // z_temp_partial;
            int m1 = 0; int m2 = 0;
            
            c_eps = constraint_eps(i1,i2,m1,m2);
            for (int j=0; j < jacdim ; j++) {
                loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
                // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                extremity_loc_eps_min = extremity_eps_min(i1,i2,m1,m2,j);
                z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
           /*     if (i == 0) {
                    cout << "j=" << j << " m1=" << m1 << " m2=" << m2 <<endl;
                    cout << "loc_eps=" << loc_eps << " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    cout << "JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " z_temp=" << z_temp << endl;
                }*/
            }
            // hull on all subdivisions
            //    z_outer = interval_hull(z_outer,z_temp);
            
            // first translate in x
            for (int m1=1; m1 < nb_discr1 ; m1++) {
                int j = 0;
                c_eps = constraint_eps(i1,i2,m1,m2);
                //for (int j=0; j < jacdim ; j++) {
                loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
                // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                extremity_loc_eps_min = extremity_eps_min(i1,i2,m1,m2,j);
                z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
         /*       if (i == 0) {
                    cout << "j=" << j << " m1=" << m1 << endl;
                    cout << "loc_eps=" << loc_eps<< " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    cout << "JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " z_temp=" << z_temp << endl;
                } */
            }
            
            // then translate in y but for whole [x]: Jacobian should be evaluated on all m1=0..nb_discr
            for (int m2=1; m2 < nb_discr2 ; m2++) {
                int j = 1;
                c_eps = constraint_eps(i1,i2,m1,m2);
                c_eps[0] = constr_eps[i1][i2][0];
                //for (int j=0; j < jacdim ; j++) {
                loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
                // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
                extremity_loc_eps_min = extremity_eps_min(i1,i2,m1,m2,j);
                z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
           /*     if (i == 0) {
                    cout << "j=" << j << " m2=" << m2 << endl;
                    cout << "loc_eps=" << loc_eps<< " loc_eps - extremity_loc_eps_min" << loc_eps - extremity_loc_eps_min << endl;
                    cout << "JacAff[i][j].convert_int(c_eps,jacdim)=" << JacAff[i][j].convert_int(c_eps,jacdim) << " z_temp=" << z_temp << endl;
                } */
            }
            
            // hull on all 4 subdivisions
            z_outer = interval_hull(z_outer,z_temp);
        //    cout << "z_outer=" << z_outer << endl;
        }
        
    }
    return z_outer;
}


// subdivise n times in the 2 inputs with simultaneous progress in the 2 dimensions around the center
interval evaluate_outerrange_x_subdiv_discretize_simultaneous(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int i) {
    interval z_temp_partial, z_temp, loc_eps, extremity_loc_eps_min, extremity_loc_eps_max, z_outer = z0[i];
    vector<interval> c_eps;
    
    z_temp = z0[i]; // z_temp_partial;
    
    int m = 0;
    
    int m1 = 0; int m2 = 0;
    c_eps = constraint_eps(m1,m2);
    for (int j=0; j < jacdim ; j++) {
        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
        extremity_loc_eps_min = extremity_eps_min(m1,m2,j);
        z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    // hull on all subdivisions
    //    z_outer = interval_hull(z_outer,z_temp);
    
    // first translate in x
    for (int m1=1; m1 < nb_discr1 ; m1++) {
        int j = 0;
        c_eps = constraint_eps(m1,m2);
        //for (int j=0; j < jacdim ; j++) {
        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
        // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
        extremity_loc_eps_min = extremity_eps_min(m1,m2,j);
        z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    // then translate in y but for whole [x]: Jacobian should be evaluated on all m1=0..nb_discr
    for (int m2=1; m2 < nb_discr2 ; m2++) {
        int j = 1;
        c_eps = constraint_eps(i1,i2,m1,m2);
        c_eps[0] = constr_eps[i1][i2][0];
        //for (int j=0; j < jacdim ; j++) {
        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
        // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
        extremity_loc_eps_min = extremity_eps_min(i1,i2,m1,m2,j);
        z_temp += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
        
    }
    
    // hull on all 4 subdivisions
    z_outer = interval_hull(z_outer,z_temp);
    
    
    
    return z_outer;
}




interval evaluate_innerrange_x(vector<interval> &z0,  vector<vector<interval>> &Jacf, bool maximal, vector<int> &exist_quantified, int i) {
    interval inner_impro, inner_pro;
    interval z_inner;
    
        inner_impro = 0;
        inner_pro = z0[i];
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],eps[j]);
            else
                inner_pro += Jacf[i][j]*eps[j];
        }
        z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}


// evaluation on one subdivision (defined by index1 and index2)
// defined per subdivision contrarily to over-approx because we cannot do the hull for inner approx
interval evaluate_innerrange_x_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2) {
    interval inner_impro, inner_pro;
    interval z_inner;
    
    inner_impro = 0;
    inner_pro = z0[i];
    for (int j=0; j < jacdim ; j++) {
        if (maximal || (exist_quantified[j] == i)) {// output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(constr_eps[index1][index2],jacdim),eps_loc[index1][index2][j]);
        }
        else {
            inner_pro += JacAff[i][j].convert_int(constr_eps[index1][index2],jacdim) * eps_loc[index1][index2][j];
        }
    }
    //   cout << endl;
    z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}

// evaluation on one subdivision (defined by index1 and index2)
// defined per subdivision contrarily to over-approx because we cannot do the hull for inner approx
/* interval evaluate_innerrange_x_subdiv_discretize_old(vector<interval> &z0,  vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2) {
    interval inner_impro, inner_pro;
    interval z_inner;
    
    inner_impro = 0;
    inner_pro = z0[i];
    for (int m=0 ; m<nb_discr; m++)
    {
        for (int j=0; j < jacdim ; j++) {
            if (maximal || (exist_quantified[j] == i)) {// output component i the one where j is existentially quantified {
                if (m>0)
                    inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim),(eps_loc_discr[index1][index2][m][j]-extremity_eps_loc_discr[index1][index2][m-1][j]));
                else
                    inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim),eps_loc_discr[index1][index2][m][j]);
            }
            else {
                if (m>0)
                    inner_pro += JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim) * (eps_loc_discr[index1][index2][m][j]-extremity_eps_loc_discr[index1][index2][m-1][j]);
                else
                    inner_pro += JacAff[i][j].convert_int(constr_eps_discr[index1][index2][m],jacdim) * eps_loc_discr[index1][index2][m][j];
            }
        }
    }
    //   cout << endl;
    z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}
*/
    
/*    z_temp_partial = z0[i];
    // evaluated at the extrimity on the n-1 first subdivisions
    for (int m=0; m < nb_discr-1 ; m++) {
        z_temp = z_temp_partial;
        for (int j=0; j < jacdim ; j++) {
            if (m>0)
            {
                z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
                z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * (extremity_eps_loc_discr[i1][i2][m][j]-extremity_eps_loc_discr[i1][i2][m-1][j]);
            }
            else
            {
                z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * eps_loc_discr[i1][i2][m][j];
                z_temp_partial += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][m],jacdim) * extremity_eps_loc_discr[i1][i2][m][j];
            }
        }
        // hull on all subdivisions
        z_outer = interval_hull(z_outer,z_temp);
        // shoudl be doing hull here I guess
    }
    // evaluated on interval on the last subdivision
    z_temp = z_temp_partial;
    for (int j=0; j < jacdim ; j++) {
        z_temp += JacAff[i][j].convert_int(constr_eps_discr[i1][i2][nb_discr-1],jacdim) * (eps_loc_discr[i1][i2][nb_discr-1][j]-extremity_eps_loc_discr[i1][i2][nb_discr-2][j]);
    }
    
    z_outer = interval_hull(z_outer,z_temp);*/




// evaluation on one subdivision (defined by index1 and index2)
// defined per subdivision contrarily to over-approx because we cannot do the hull for inner approx
interval evaluate_innerrange_x_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff, bool maximal, vector<int> &exist_quantified, int i, int index1, int index2) {
    interval inner_impro, inner_pro, loc_eps;
    interval z_inner;
    vector<interval> c_eps(jacdim);
    double extremity_loc_eps_min;
    
    inner_impro = 0;
    inner_pro = z0[i];
    
    // mean-value on first block
    int m1 = 0, m2 = 0;
    c_eps = constraint_eps(index1,index2,m1,m2);
    for (int j=0; j < jacdim ; j++) {
        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
        extremity_loc_eps_min = extremity_eps_min(index1,index2,m1,m2,j); // always 0 ?
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(c_eps,jacdim),loc_eps - extremity_loc_eps_min);
        else
            inner_pro += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    // first translate in first component x
    for (int m1=1; m1 < nb_discr1 ; m1++) {
        int j = 0;
        c_eps = constraint_eps(index1,index2,m1,m2);
        //for (int j=0; j < jacdim ; j++) {
        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
        // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
        extremity_loc_eps_min = extremity_eps_min(index1,index2,m1,m2,j);
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(c_eps,jacdim),loc_eps - extremity_loc_eps_min);
        else
            inner_pro += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    // then translate in y but for whole [x]: Jacobian should be evaluated on all m1=0..nb_discr
    for (int m2=1; m2 < nb_discr2 ; m2++) {
        int j = 1;
        c_eps = constraint_eps(index1,index2,m1,m2);
        c_eps[0] = constr_eps[index1][index2][0];
        //for (int j=0; j < jacdim ; j++) {
        loc_eps = eps[j].mid() + eps[j].rad() * c_eps[j];
        // extremity_loc_eps_max = extremity_eps_max(i1,i2,m1,m2,j);
        extremity_loc_eps_min = extremity_eps_min(index1,index2,m1,m2,j);
        if (maximal || (exist_quantified[j] == i))  // output component i the one where j is existentially quantified {
            inner_impro += Kaucher_multeps(JacAff[i][j].convert_int(c_eps,jacdim),loc_eps - extremity_loc_eps_min);
        else
            inner_pro += JacAff[i][j].convert_int(c_eps,jacdim) * (loc_eps - extremity_loc_eps_min);
    }
    
    //   cout << endl;
    z_inner = Kaucher_add_pro_impro(inner_pro,inner_impro);
    
    return z_inner;
}


    
    

void print_pi(vector<int> &exist_quantified) {
    cout << "(pi : ";
    for (int j = 0 ; j < jacdim ; j++)
        cout << j+1 << " -> " << exist_quantified[j]+1 << ", " ;
    cout << ") :  ";
}


// z0 = f(x0)
void joint_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "joint_inner.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    interval inner_x = evaluate_innerrange_x(z0,Jacf,false,exist_quantified,varx);
    interval inner_y = evaluate_innerrange_x(z0,Jacf,false,exist_quantified,vary);
    cout << "inner box ";
    print_pi(exist_quantified);
    cout << inner_x << "  " << inner_y << endl;
    outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    inner_x = evaluate_innerrange_x(z0,Jacf,false,exist_quantified,varx);
    inner_y = evaluate_innerrange_x(z0,Jacf,false,exist_quantified,vary);
    cout << "inner box ";
    print_pi(exist_quantified);
    cout << inner_x << "  " << inner_y << endl;
    outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    outFile_jointinner.close();
}




// z0 = f(x0)

void joint_ranges_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    interval temp_x, tempy, inner_x, inner_y;
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "joint_inner.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    
    cout << "inner box "; print_pi(exist_quantified); cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    // TODO : en deriver une boite gloabel ?
   //  inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,varx);
   //  inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,vary);
    
   // cout << inner_x << "  " << inner_y << endl;
   // outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    cout << "inner box ";
    print_pi(exist_quantified);
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    outFile_jointinner.close();
}


void joint_ranges_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<bool> is_existential(jacdim);
    interval temp_x, tempy, inner_x, inner_y;
    
    ofstream outFile_jointinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "joint_inner.out";
    outFile_jointinner.open(file_name.str().c_str(),fstream::app);
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[varx] = vary;
    exist_quantified[vary] = varx;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    
    cout << "inner box "; print_pi(exist_quantified); cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv_discretize(z0,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv_discretize(z0,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    // TODO : en deriver une boite gloabel ?
    //  inner_x = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,varx);
    //  inner_y = evaluate_innerrange_x_subdiv(z0,JacAff,exist_quantified,vary);
    
    // cout << inner_x << "  " << inner_y << endl;
    // outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    cout << "inner box ";
    print_pi(exist_quantified);
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            inner_x = evaluate_innerrange_x_subdiv_discretize(z0,JacAff,false,exist_quantified,varx,i1,i2);
            inner_y = evaluate_innerrange_x_subdiv_discretize(z0,JacAff,false,exist_quantified,vary,i1,i2);
            cout << inner_x << "  " << inner_y << "       ";
            outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
        }
    }
    cout << endl;
    outFile_jointinner.close();
}



// f0 = f(x0),
// tests by Hansen Sengupta if a candidate box (fix,fixy) is in the joint range for f_varx, f_fary
/*void Hansen_Sengupta(vector<interval> &f0, vector<vector<interval>> &Jacf, interval &fix, interval &fi2, int varx, int vary)
{
    interval xtilde = x0 + gamma(Jacf)
}*/



vector<vector<double>> print_skewbox(interval &temp_inner_x, interval &temp_inner_y, vector<vector<double>> &A, int varx, int vary, ofstream &outFile) {
    // resulting quadrilatere A * inner is an inner approximation of the range of f
    vector<vector<double>> output_skewedbox(4);
    for (int i=0; i<4; i++) {
        output_skewedbox[i] = vector<double>(2);
    }
    
    output_skewedbox[0][0] = inf(temp_inner_x)*A[varx][varx] + inf(temp_inner_y)*A[varx][vary];
    output_skewedbox[0][1] = inf(temp_inner_x)*A[vary][varx] + inf(temp_inner_y)*A[vary][vary];
    output_skewedbox[1][0] = inf(temp_inner_x)*A[varx][varx] + sup(temp_inner_y)*A[varx][vary];
    output_skewedbox[1][1] = inf(temp_inner_x)*A[vary][varx] + sup(temp_inner_y)*A[vary][vary];
    output_skewedbox[2][0] = sup(temp_inner_x)*A[varx][varx] + sup(temp_inner_y)*A[varx][vary];
    output_skewedbox[2][1] = sup(temp_inner_x)*A[vary][varx] + sup(temp_inner_y)*A[vary][vary];
    output_skewedbox[3][0] = sup(temp_inner_x)*A[varx][varx] + inf(temp_inner_y)*A[varx][vary];
    output_skewedbox[3][1] = sup(temp_inner_x)*A[vary][varx] + inf(temp_inner_y)*A[vary][vary];
    
  //  cout << "inner skewed box: (pi : 1 -> " << exist_quantified[0]+1 << ", 2 -> " << exist_quantified[1]+1 << ")" << endl;
    for (int i=0; i<4; i++)
        cout << "(" << output_skewedbox[i][0] <<", " << output_skewedbox[i][1] << ") ";
    cout << endl;
    
    for (int i=0; i<3; i++)
        outFile << output_skewedbox[i][0] << "\t" << output_skewedbox[i][1] << "\t" ;
    outFile << output_skewedbox[3][0] << "\t" << output_skewedbox[3][1] << endl ;
    
    return output_skewedbox;
}




void preconditioned_joint_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim), C(sysdim);
    for (int i=0 ; i<sysdim; i++) {
        A[i] = vector<double> (sysdim);
        C[i] = vector<double> (sysdim);
    }
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(Jacf[varx][varx]);
    A[vary][vary] = mid(Jacf[vary][vary]);
    A[varx][vary] = mid(Jacf[varx][vary]);
    A[vary][varx] = mid(Jacf[vary][varx]);
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // CJacf = C * Jacf
    vector<vector<interval>> CJacf(sysdim);
    for (int i=0; i < sysdim ; i++)
        CJacf[i] = vector<interval>(jacdim);
   
    multMiMi(CJacf,C,Jacf);
    cout << "CJacf=" << endl;
    cout << CJacf;
  
    
    // outer range
    interval temp_outer_x = evaluate_outerrange_x(f0,CJacf,varx);
    interval temp_outer_y = evaluate_outerrange_x(f0,CJacf,vary);
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, outFile_skewedouter);
    
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    interval temp_inner_x, temp_inner_y;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
     temp_inner_x = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,varx);
     temp_inner_y = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,vary);
    
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    temp_inner_x = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,varx);
    temp_inner_y = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,vary);
    
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        temp_inner_x = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,varx);
        temp_inner_y = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,vary);
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        temp_inner_x = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,varx);
        temp_inner_y = evaluate_innerrange_x(f0,CJacf,false,exist_quantified,vary);
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}



void preconditioned_joint_ranges_subdiv(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim), C(sysdim);
    for (int i=0 ; i<sysdim; i++) {
        A[i] = vector<double> (sysdim);
        C[i] = vector<double> (sysdim);
    }
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(JacAff[varx][varx].convert_int());
    A[vary][vary] = mid(JacAff[vary][vary].convert_int());
    A[varx][vary] = mid(JacAff[varx][vary].convert_int());
    A[vary][varx] = mid(JacAff[vary][varx].convert_int());
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // CJacf = C * Jacf
    vector<vector<AAF>> CJacAff(sysdim);
    for (int i=0; i < sysdim ; i++)
        CJacAff[i] = vector<AAF>(jacdim);
    
    multMiMi(CJacAff,C,JacAff);
  //  cout << "CJacf=" << endl;
  //  cout << CJacf;
    
    
    // outer range
    interval temp_outer_x = evaluate_outerrange_x_subdiv(f0,CJacAff,varx);
    interval temp_outer_y = evaluate_outerrange_x_subdiv(f0,CJacAff,vary);
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, outFile_skewedouter);
    
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    interval temp_inner_x, temp_inner_y;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            temp_inner_x = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,vary,i1,i2);
      //      cout << "f0= " << f0;
      //      cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
        }
    }
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
    //        cout << "f0= " << f0 ;
            temp_inner_x = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,vary,i1,i2);
   //         cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
        }
    }
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,vary,i1,i2);
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
            }
        }
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv(f0,CJacAff,false,exist_quantified,vary,i1,i2);
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
            }
        }
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}



void preconditioned_joint_ranges_subdiv_discretize(vector<interval> &z0,  vector<vector<AAF>> &JacAff, int varx, int vary) {
    
    vector<vector<double>> output_skewedbox;
    
    ofstream outFile_skewedinner;
    stringstream file_name;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_inner.out";
    outFile_skewedinner.open(file_name.str().c_str(),fstream::app);
    
    ofstream outFile_skewedouter;
    file_name.str("");
    file_name << "output/x" << varx+1 << "x" << vary+1 << "skewed_outer.out";
    outFile_skewedouter.open(file_name.str().c_str(),fstream::app);
    
    // precobditionning: A is center of Jacobian, C its inverse
    vector<vector<double>> A(sysdim), C(sysdim);
    for (int i=0 ; i<sysdim; i++) {
        A[i] = vector<double> (sysdim);
        C[i] = vector<double> (sysdim);
    }
    
    for (int i=0 ; i<sysdim; i++) {
        A[i][i] = 1.0;
        C[i][i] = 1.0;
    }
    // 2D preconditioner - just on components varx and vary
    A[varx][varx] = mid(JacAff[varx][varx].convert_int());
    A[vary][vary] = mid(JacAff[vary][vary].convert_int());
    A[varx][vary] = mid(JacAff[varx][vary].convert_int());
    A[vary][varx] = mid(JacAff[vary][varx].convert_int());
    
    // C is inverse of A
    double determinant = 1.0/(A[varx][varx]*A[vary][vary]-A[varx][vary]*A[vary][varx]);
    C[varx][varx] = determinant*A[vary][vary];
    C[varx][vary] = - determinant*A[varx][vary];
    C[vary][varx] = - determinant*A[vary][varx];
    C[vary][vary] = determinant*A[varx][varx];
    
    // f0 = C * z0
    vector<interval> f0(sysdim);
    for (int i=0 ; i<sysdim; i++)
        f0[i] = z0[i];
    f0[varx] = C[varx][varx]*z0[varx] + C[varx][vary]*z0[vary];
    f0[vary] = C[vary][vary]*z0[vary] + C[vary][varx]*z0[varx];
    
    // CJacf = C * Jacf
    vector<vector<AAF>> CJacAff(sysdim);
    for (int i=0; i < sysdim ; i++)
        CJacAff[i] = vector<AAF>(jacdim);
    
    multMiMi(CJacAff,C,JacAff);
    //  cout << "CJacf=" << endl;
    //  cout << CJacf;
    
    
    // outer range
    interval temp_outer_x = evaluate_outerrange_x_subdiv_discretize(f0,CJacAff,varx);
    interval temp_outer_y = evaluate_outerrange_x_subdiv_discretize(f0,CJacAff,vary);
    cout << "outer skewed box: ";
    output_skewedbox = print_skewbox(temp_outer_x, temp_outer_y, A,  varx,  vary, outFile_skewedouter);
    
    
    // for each input, index of the output in which it is existentially quantified
    vector<int> exist_quantified(jacdim);
    interval temp_inner_x, temp_inner_y;
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = vary;
    exist_quantified[1] = varx;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,vary,i1,i2);
            //      cout << "f0= " << f0;
            //      cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
        }
    }
    
    // TO BE generalized for more then 2 variables/components
    exist_quantified[0] = varx;
    exist_quantified[1] = vary;
    for (int i=2; i < jacdim ; i++)
    {
        if (i % 2 == 0)
            exist_quantified[i] = varx;
        else
            exist_quantified[i] = vary;
    }
    cout << "inner skewed box: ";
    print_pi(exist_quantified);
    cout << endl;
    for (int i1=0 ; i1<=1 ; i1++) {
        for (int i2=0 ; i2<=1 ; i2++)  {
            //        cout << "f0= " << f0 ;
            temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,varx,i1,i2);
            temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,vary,i1,i2);
            //         cout << "temp_inner_x = " << temp_inner_x << " temp_inner_y = " << temp_inner_y << endl;
            output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
        }
    }
    
    if (jacdim >= 3)
    {
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = vary;
        exist_quantified[1] = varx;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,vary,i1,i2);
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
            }
        }
        
        // TO BE generalized for more then 2 variables/components
        exist_quantified[0] = varx;
        exist_quantified[1] = vary;
        for (int i=2; i < jacdim ; i++)
        {
            if (i % 2 == 0)
                exist_quantified[i] = vary;
            else
                exist_quantified[i] = varx;
        }
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        
        for (int i1=0 ; i1<=1 ; i1++) {
            for (int i2=0 ; i2<=1 ; i2++)  {
                temp_inner_x = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,varx,i1,i2);
                temp_inner_y = evaluate_innerrange_x_subdiv_discretize(f0,CJacAff,false,exist_quantified,vary,i1,i2);
                
                output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
            }
        }
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}
