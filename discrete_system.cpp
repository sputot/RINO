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
    
    vector<F<AAF>> x(jacdim);
    vector<vector<interval>> Jacf(sysdim);
    for (int i=0; i < sysdim ; i++)
        Jacf[i] = vector<interval>(jacdim);
    vector<F<AAF>> z(sysdim), z1(sysdim), z2(sysdim);
    vector<interval> x0(jacdim);    // center
    
    DiscreteFunc f;
    
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
    
    system("rm -r output");
    system("mkdir output");
    vector<ofstream> outFile_outer_direct(sysdim);
    stringstream file_name;
    
    for (int i=0; i < jacdim ; i++) {
        cout << "outer range of f(x) direct evaluation:" << z[i].x().convert_int() << endl;
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_direct.out";
        outFile_outer_direct[i].open(file_name.str().c_str());
        outFile_outer_direct[i] << inf(z[i].x().convert_int()) << "\t" << sup(z[i].x().convert_int()) << endl;
        outFile_outer_direct[i].close();
    }
    
  
   
    cout << endl;
    
    for (int i=0; i < sysdim ; i++)
        for (int j=0; j < jacdim ; j++)
            Jacf[i][j] = z[i].d(j).convert_int();
    cout << "Jacf evaluated on all [x]" << endl;
    cout << Jacf;
   
    vector<interval> z0 = f(x0);
    
    evaluate_projections(z0, Jacf);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf,0,1);
        preconditioned_joint_ranges(z0,Jacf,0,1);
    }
    
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
    x[0] = mid(initial_values[0].convert_int());
    
    for (int i=0; i < jacdim ; i++) {
        x[i].diff(i,jacdim);    // differentiate to x
    }
    z1 = f(x);
    
    for (int i=0; i < sysdim ; i++) {
        for (int j=0; j < jacdim-1 ; j++)
            Jacf_refined[i][j] = z1[i].d(j).convert_int();
        for (int j=jacdim-1; j < jacdim ; j++)
            Jacf_refined[i][j] = z[i].d(j).convert_int();
    }
    cout << "Jacf: 1st column evaluated on x0, 2nd column on [x]" << endl;
    cout << Jacf_refined;

    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
    
    cout << endl;
   // cout << "Improved mean-value:" << endl;
    for (int i=0; i < jacdim ; i++) {
        x[i] = initial_values[i];
        eps[i] = initial_values[i].convert_int();
        x0[i] = mid(eps[i]);
        eps[i] = eps[i] - x0[i];
    }
    x[1] = mid(initial_values[1].convert_int());
    
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
    cout << "Jacf: 1st column evaluated on [x], 2nd column on all x0" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
    
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
    cout << "Jacf: f1: 1st column evaluated on x0, 2nd column on [x], f2: 1st column evaluated on [x], 2nd column on x0" << endl;
    cout << Jacf_refined;
    
    evaluate_projections(z0, Jacf_refined);
    if (sysdim >= 2) {
        joint_ranges(z0,Jacf_refined,0,1);
        preconditioned_joint_ranges(z0,Jacf_refined,0,1);
    }
    
    
}

void evaluate_projections(vector<interval> &z0,  vector<vector<interval>> &Jacf) {
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


void evaluate_ranges(vector<interval> &z0,  vector<vector<interval>> &Jacf, vector<bool> &is_existential) {
    interval inner_impro, inner_pro, outer_impro;
    vector<interval> z_inner(sysdim);
    vector<interval> z_outer(sysdim);
    
    for (int i=0; i < sysdim ; i++) {
        z_outer[i] = z0[i];
        inner_impro = 0;
        inner_pro = z0[i];
        for (int j=0; j < jacdim ; j++) {
            if (is_existential[j]) {
                z_outer[i] += Jacf[i][j]*eps[j];
                inner_impro += Kaucher_multeps(Jacf[i][j],eps[j]);
            }
            else {
                outer_impro += Kaucher_multeps(Jacf[i][j],eps[j]);
                inner_pro += Jacf[i][j]*eps[j];
            }
        }
        z_inner[i] = Kaucher_add_pro_impro(inner_pro,inner_impro);
        z_outer[i] = Kaucher_add_pro_impro_resultpro(z_outer[i],outer_impro);
    }
    
    string str = "";
    for (int i=0; i < jacdim ; i++)
        if (!is_existential[i])
            str = str+to_string(i+1);
        
    for (int i=0; i < sysdim ; i++) {
        cout << "outer range of f(x,"<< str <<") by mean-value:" << z_outer[i] << endl;
        cout << "inner range of f(x,"<< str <<") by mean-value:" << z_inner[i] << endl;
    }
}

interval evaluate_outerrange_x(vector<interval> &z0,  vector<vector<interval>> &Jacf, int i) {
    interval z_outer = z0[i];
    
    for (int j=0; j < jacdim ; j++) {
        z_outer += Jacf[i][j]*eps[j];
    }
    
    return z_outer;
}



interval evaluate_innerrange_x(vector<interval> &z0,  vector<vector<interval>> &Jacf, vector<int> &exist_quantified, int i) {
    interval inner_impro, inner_pro;
    interval z_inner;
    
        inner_impro = 0;
        inner_pro = z0[i];
        for (int j=0; j < jacdim ; j++) {
            if (exist_quantified[j] == i) // output component i the one where j is existentially quantified
                inner_impro += Kaucher_multeps(Jacf[i][j],eps[j]);
            else
                inner_pro += Jacf[i][j]*eps[j];
        }
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
    interval inner_x = evaluate_innerrange_x(z0,Jacf,exist_quantified,varx);
    interval inner_y = evaluate_innerrange_x(z0,Jacf,exist_quantified,vary);
    cout << "inner box ";
    print_pi(exist_quantified);
    cout << inner_x << "  " << inner_y << endl;
    outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
    exist_quantified[varx] = varx;
    exist_quantified[vary] = vary;
    if (jacdim >=3)
        exist_quantified[2] = 0;
    inner_x = evaluate_innerrange_x(z0,Jacf,exist_quantified,varx);
    inner_y = evaluate_innerrange_x(z0,Jacf,exist_quantified,vary);
    cout << "inner box ";
    print_pi(exist_quantified);
    cout << inner_x << "  " << inner_y << endl;
    outFile_jointinner << inf(inner_x) << "\t" << sup(inner_x) << "\t" << inf(inner_y) << "\t" << sup(inner_y) <<  endl;
    
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
     temp_inner_x = evaluate_innerrange_x(f0,CJacf,exist_quantified,varx);
     temp_inner_y = evaluate_innerrange_x(f0,CJacf,exist_quantified,vary);
    
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
    temp_inner_x = evaluate_innerrange_x(f0,CJacf,exist_quantified,varx);
    temp_inner_y = evaluate_innerrange_x(f0,CJacf,exist_quantified,vary);
    
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
        temp_inner_x = evaluate_innerrange_x(f0,CJacf,exist_quantified,varx);
        temp_inner_y = evaluate_innerrange_x(f0,CJacf,exist_quantified,vary);
        
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
        temp_inner_x = evaluate_innerrange_x(f0,CJacf,exist_quantified,varx);
        temp_inner_y = evaluate_innerrange_x(f0,CJacf,exist_quantified,vary);
        
        cout << "inner skewed box: ";
        print_pi(exist_quantified);
        output_skewedbox = print_skewbox(temp_inner_x, temp_inner_y, A,  varx,  vary, outFile_skewedinner);
    }
    
    
    
    outFile_skewedinner.close();
    outFile_skewedouter.close();
}
