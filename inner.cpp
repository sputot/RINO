/* ============================================================================
 File   : inner.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This file defines operations specific to Kaucher arithmetic for computing inner ranges
 ============================================================================ */

#include "filib_interval.h"
#include "fadbad_aa.h"
#include "matrix.h"
#include "inner.h"
#include <iostream>
#include <ostream>
#include <fstream>
using namespace std;



// Kaucher multiplication Jacobian * (z-z0) where we consider the improper dual(eps=(z-z0))
// eps in (dual Z) => inf(y)=sup(eps) and sup(y)=inf(eps) in Kaucher multiplication formula
interval Kaucher_multeps(interval Jac, interval eps)
{
    interval res;
    if ((inf(Jac) <= 0) && (sup(Jac) >= 0))
        res = 0.0;
    else if (sup(Jac) < 0) // (-P) . dual Z = [sup(x)sup(y),sup(x)inf(y)] => we take the dual [sup(x)inf(y),sup(x)sup(y)] because we need to return a proper interval
        res = interval(sup(Jac)*sup(eps),sup(Jac)*inf(eps));
    else  // (P) . dual Z = [inf(x)inf(y),inf(x)sup(y)] => we take the dual  [inf(x)sup(y),inf(x)inf(y)] because we need to return a proper interval
        res = interval(inf(Jac)*inf(eps),inf(Jac)*sup(eps));
    return res;
}

// Kaucher addition where we consider pro as proper and impro as improper and we want an improper result
//  (though they all are represented as proper intervals)
interval Kaucher_add_pro_impro(interval pro, interval impro)
{
    interval res;
    if  (inf(impro)+sup(pro) <= sup(impro)+inf(pro)) // result is improper
        res = interval(inf(impro)+sup(pro),sup(impro)+inf(pro));
    else
        res = empty();
    return res;
}

// Kaucher addition where we consider pro as proper and impro as improper and we want an proper result
//  (though they all are represented as proper intervals)
interval Kaucher_add_pro_impro_resultpro(interval pro, interval impro)
{
    interval res;
    if  (inf(impro)+sup(pro) >= sup(impro)+inf(pro)) // result is proper
        res = interval(sup(impro)+inf(pro),inf(impro)+sup(pro));
    else
        res = empty();
    return res;
}



// computer inner and outer-approx by mean-value theorem
// Xinner_robust : when the 'uncontrolled' last parameters of parameter set are considered as not controllable => proper part that goes agains the improper part of the inner-approx on the other parameters
// Application of the Generalized Mean-Value Theorem
void InnerOuter(vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xinner_minimal, vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<interval> &Xouter_minimal, vector<AAF> &x0p1, vector<vector<AAF>> &Jtau, vector<interval> &eps, double tnp1)
{
    vector<vector<interval>> Jaux(sysdim, vector<interval>(jacdim)); // will copy only the sysdim relevant lines of Jaux
    vector<interval> ix0(sysdim);
    interval initialcondition_impro, initialcondition_pro, uncontrolled_impro, uncontrolled_pro, controlled_impro, controlled_pro, aux_impro;
    
    for (int i=0 ; i<sysdim ; i++) {
        ix0[i] = x0p1[i].convert_int();
        for (int k=0 ; k<jacdim ; k++)
            Jaux[i][k] = Jtau[i][k].convert_int();
    }
    
    // inner approx: special addition corresponding to addition of proper and improper intervals: the result is an inner range for the solution of ODE system
    for (int i=0 ; i<sysdim; i++) {
        
        initialcondition_impro = 0;
        initialcondition_pro = 0;
        uncontrolled_impro = 0;
        uncontrolled_pro = 0;
        controlled_impro = 0;
        controlled_pro = 0;
   //     innerjoint_pro = 0;
   //     innerjoint_impro = 0;
        
        
        Xinner_robust[i] = 0;
        Xinner_minimal[i] = 0;
        Xouter_robust[i] = 0;
        
        for (int j=0 ; j<jacdim ; j++)
        { // adding the controlled (exists) part
            if (j<sysdim) // initial conditions //(is_initialcondition[j])
            {
                initialcondition_impro = initialcondition_impro + Kaucher_multeps(Jaux[i][j],eps[j]);
                initialcondition_pro = initialcondition_pro + Jaux[i][j]*eps[j];
         /*       if (i == j) {
                    innerjoint_impro = innerjoint_impro + Kaucher_multeps(Jaux[i][j],eps[j]);
                }
                else {
                    innerjoint_pro = innerjoint_pro + Jaux[i][j]*eps[j];
                } */
            }
            else if (!is_uncontrolled[j-sysdim])
            {
                controlled_impro = controlled_impro + Kaucher_multeps(Jaux[i][j],eps[j]);
                controlled_pro = controlled_pro + Jaux[i][j]*eps[j];
            }
            else // uncontrolled
            {
                uncontrolled_impro = uncontrolled_impro + Kaucher_multeps(Jaux[i][j],eps[j]) ;
                uncontrolled_pro = uncontrolled_pro + Jaux[i][j]*eps[j] ;
            }
        }
        
        
        aux_impro = initialcondition_impro + controlled_impro;
    //    innerjoint_impro = innerjoint_impro + controlled_impro;
        
        
        Xinner[i] = Kaucher_add_pro_impro(ix0[i], aux_impro + uncontrolled_impro);
        
   //     innerjoint_pro = ix0[i] + innerjoint_pro;
    //    Xinner_joint[i] = Kaucher_add_pro_impro(innerjoint_pro, innerjoint_impro + uncontrolled_impro);
        // added here so that we can store the different joint inner-approximations
//        outFile_inner_joint[i] << tnp1 << "\t" << inf(Xinner_joint[i]) << "\t" << sup(Xinner_joint[i]) << endl;
        
        if (uncontrolled > 0)
        {
            Xinner_robust[i]  = ix0[i] + uncontrolled_pro;
            Xinner_robust[i] = Kaucher_add_pro_impro(Xinner_robust[i], aux_impro);
            Xouter_robust[i] = ix0[i] + initialcondition_pro + controlled_pro;
            Xouter_robust[i] = Kaucher_add_pro_impro_resultpro(Xouter_robust[i],uncontrolled_impro);
        }
        if (controlled > 0 || uncontrolled > 0)
        {
            Xinner_minimal[i] = ix0[i] + uncontrolled_pro + controlled_pro;
            Xinner_minimal[i] = Kaucher_add_pro_impro(Xinner_minimal[i], initialcondition_impro);
            Xouter_minimal[i] = ix0[i] + initialcondition_pro;
            Xouter_minimal[i] =Kaucher_add_pro_impro_resultpro(Xouter_minimal[i],controlled_impro + uncontrolled_impro);
        }
        
    }
    
    // outer-approx computation : Xouter = x0p1 + Jtau*eps;
    multMiVi(Xouter,Jaux,eps);
    addViVi(Xouter,ix0);
    
    
    // computing joint inner-approximations
    // for any pair (i,k) of system outputs which we wish to jointly inner-approximate
    //    component j for j < sysdim/2 is existential in z_i and universal in z_k
    //    except i is existential in z_i and k is existential in z_k
    //    same for inputs
    interval range_i, range_k, range_i_impro, range_k_impro;
    for (int i=0 ; i<sysdim; i++) {
        for (int k=i+1 ; k<sysdim; k++) {
            range_i = ix0[i];   // center and uncontrolled (forall) part
            range_i_impro = 0;      // controlled (existential) part
            range_k = ix0[k];
            range_k_impro = 0;
            for (int j=0 ; j<sysdim/2 ; j++) {
                range_k = range_k + Jaux[k][j]*eps[j];
                range_i_impro = range_i_impro + Kaucher_multeps(Jaux[i][j],eps[j]);
            }
            for (int j=sysdim/2 ; j<sysdim ; j++) {
                range_i = range_i + Jaux[i][j]*eps[j];
                range_k_impro = range_k_impro + Kaucher_multeps(Jaux[k][j],eps[j]);
            }
            for (int j=0 ; j<fullinputsdim/2 ; j++) {
                if (is_uncontrolled[j])
                    range_i = range_i + Jaux[i][j+sysdim]*eps[j+sysdim];
                else
                    range_i_impro = range_i_impro + Kaucher_multeps(Jaux[i][j+sysdim],eps[j+sysdim]);
                range_k = range_k + Jaux[k][j+sysdim]*eps[j+sysdim];
            }
            for (int j=fullinputsdim/2 ; j<fullinputsdim ; j++) {
                range_i = range_i + Jaux[i][j+sysdim]*eps[j+sysdim];
                if (is_uncontrolled[j])
                    range_k = range_k + Jaux[k][j+sysdim]*eps[j+sysdim];
                else
                    range_k_impro = range_k_impro + Kaucher_multeps(Jaux[k][j+sysdim],eps[j+sysdim]);
            }
           
            range_i = Kaucher_add_pro_impro(range_i,range_i_impro);
            range_k = Kaucher_add_pro_impro(range_k,range_k_impro);
            
            outFile_joint_inner[i][k] << tnp1 << "\t" << inf(range_i) << "\t" << sup(range_i) << "\t" << inf(range_k) << "\t" << sup(range_k) << endl;
        }
    }
    
    // note that we could enumerate all possible mappings and superpose all corresponding boxes in the same file ?
    
    
}

