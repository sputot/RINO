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
//#include "discrete_system.h"
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

void compute_print_jointinnerranges(interval &range_i, interval &range_k, vector<vector<interval>> &Jaux, vector<interval> &eps, int i, int k, const bool skewed, vector<vector<double>> &A)
{
    interval range_i_impro, range_k_impro, robust_i, robust_k, robust_i_impro, robust_k_impro, max_i, max_k, max_i_impro, max_k_impro;
    vector<double> temp(4);
    
    range_i_impro = 0;      // controlled (existential) part
    range_k_impro = 0;
    for (int j=0 ; j<sysdim/2 ; j++) {
        range_k = range_k + Jaux[k][j]*eps[j];
        range_i_impro = range_i_impro + Kaucher_multeps(Jaux[i][j],eps[j]);
    }
    for (int j=sysdim/2 ; j<sysdim ; j++) {
        range_i = range_i + Jaux[i][j]*eps[j];
        range_k_impro = range_k_impro + Kaucher_multeps(Jaux[k][j],eps[j]);
    }
    
    robust_i = 0; max_i = 0;
    robust_i_impro = 0; max_i_impro = 0;
    robust_k = 0; max_k = 0;
    robust_k_impro = 0; max_k_impro = 0;
    for (int j=0 ; j<fullinputsdim/2 ; j++) {
        if (is_uncontrolled[index_param[j]])
            robust_i += Jaux[i][j+sysdim]*eps[j+sysdim];
        else
            robust_i_impro += Kaucher_multeps(Jaux[i][j+sysdim],eps[j+sysdim]);
        max_i_impro += Kaucher_multeps(Jaux[i][j+sysdim],eps[j+sysdim]);
        
        robust_k = robust_k + Jaux[k][j+sysdim]*eps[j+sysdim];
        max_k += Jaux[k][j+sysdim]*eps[j+sysdim];
    }
    for (int j=fullinputsdim/2 ; j<fullinputsdim ; j++) {
        robust_i += Jaux[i][j+sysdim]*eps[j+sysdim];
        max_i += Jaux[i][j+sysdim]*eps[j+sysdim];
        
        if (is_uncontrolled[index_param[j]])
            robust_k += Jaux[k][j+sysdim]*eps[j+sysdim];
        else
            robust_k_impro += Kaucher_multeps(Jaux[k][j+sysdim],eps[j+sysdim]);
        max_k_impro += Kaucher_multeps(Jaux[k][j+sysdim],eps[j+sysdim]);
    }
    
    robust_i = Kaucher_add_pro_impro(robust_i + range_i, range_i_impro + robust_i_impro);
    robust_k = Kaucher_add_pro_impro(robust_k + range_k,range_k_impro + robust_k_impro);
    
    max_i = Kaucher_add_pro_impro(max_i + range_i, range_i_impro + max_i_impro);
    max_k = Kaucher_add_pro_impro(max_k + range_k, range_k_impro + max_k_impro);
    
    if (skewed) {
        compute_print_skewbox(max_i,max_k,A,i,k,"maxskew");
        
        if (uncontrolled > 0)
            compute_print_skewbox(robust_i,robust_k,A,i,k,"robskew");
    }
    else
    {
        temp[0] = inf(max_i); temp[1] = sup(max_i); temp[2] = inf(max_k); temp[3] = sup(max_k);
        out_approx << YAML::Key << "maxbox";
        out_approx << YAML::Value << temp;
    
        if (uncontrolled > 0) {
            temp[0] = inf(robust_i); temp[1] = sup(robust_i); temp[2] = inf(robust_k); temp[3] = sup(robust_k);
            out_approx << YAML::Key << "robbox";
            out_approx << YAML::Value << temp;
        }
    }
}




// computer inner and outer-approx by mean-value theorem
// Xinner_robust : when the 'uncontrolled' last parameters of parameter set are considered as not controllable => proper part that goes agains the improper part of the inner-approx on the other parameters
// Application of the Generalized Mean-Value Theorem
void InnerOuter(vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<AAF> &x0p1, vector<vector<AAF>> &Jtau, vector<interval> &eps, vector<AAF> &xp1, double tnp1)
{
    vector<vector<interval>> Jloc(sysdim, vector<interval>(jacdim)); // will copy only the sysdim relevant lines of Jaux
    vector<interval> ix0(sysdim);
    interval initialcondition_impro, initialcondition_pro, uncontrolled_impro, uncontrolled_pro, controlled_impro, controlled_pro, aux_impro;
    
    for (int i=0 ; i<sysdim ; i++) {
        ix0[i] = x0p1[i].convert_int();
        for (int k=0 ; k<jacdim ; k++)
            Jloc[i][k] = Jtau[i][k].convert_int();
    }
    
    // inner approx: special addition corresponding to addition of proper and improper intervals: the result is an inner range for the solution of ODE system
    for (int i=0 ; i<sysdim; i++) {
        
        initialcondition_impro = 0;
        initialcondition_pro = 0;
        uncontrolled_impro = 0;
        uncontrolled_pro = 0;
        controlled_impro = 0;
        controlled_pro = 0;
        
        Xinner_robust[i] = 0;
        Xouter_robust[i] = 0;
        
        for (int j=0 ; j<jacdim ; j++)
        { // adding the controlled (exists) part
            if (j<sysdim) // initial conditions //(is_initialcondition[j])
            {
                initialcondition_impro = initialcondition_impro + Kaucher_multeps(Jloc[i][j],eps[j]);
                initialcondition_pro = initialcondition_pro + Jloc[i][j]*eps[j];
            }
            else if (!is_uncontrolled[index_param[j-sysdim]])
            {
                controlled_impro = controlled_impro + Kaucher_multeps(Jloc[i][j],eps[j]);
                controlled_pro = controlled_pro + Jloc[i][j]*eps[j];
            }
            else // uncontrolled part
            {
                uncontrolled_impro = uncontrolled_impro + Kaucher_multeps(Jloc[i][j],eps[j]) ;
                uncontrolled_pro = uncontrolled_pro + Jloc[i][j]*eps[j] ;
            }
        }
        
        aux_impro = initialcondition_impro + controlled_impro;
        Xinner[i] = Kaucher_add_pro_impro(ix0[i] , aux_impro + uncontrolled_impro);

        if (uncontrolled > 0)
        {
            Xinner_robust[i]  = ix0[i] + uncontrolled_pro;
            Xinner_robust[i] = Kaucher_add_pro_impro(Xinner_robust[i], aux_impro);
            Xouter_robust[i] = ix0[i] + initialcondition_pro + controlled_pro;
            Xouter_robust[i] = Kaucher_add_pro_impro_resultpro(Xouter_robust[i],uncontrolled_impro);
        }
    }
    
    // outer-approx computation : Xouter = x0p1 + Jtau*eps;
    multMiVi(Xouter,Jloc,eps);
    addViVi(Xouter,ix0);
    
    
    // computing joint inner-approximations
    // for any pair (i,k) of system outputs which we wish to jointly inner-approximate
    //    component j for j < sysdim/2 is existential in z_i and universal in z_k
    //    same for inputs
    // note that we could enumerate all possible mappings and superpose all corresponding boxes in the same file ?
    if (print_debug)
    {
        vector<interval> range(3);
        interval range_i, range_k;
        interval robust_i, max_i, robust_i_impro, max_i_impro, range_i_impro;
        interval robust_k, max_k, robust_k_impro, max_k_impro, range_k_impro;
        
        // precobditionning: A is center of Jacobian, C its inverse
        vector<vector<double>> A(sysdim), C(sysdim);
        for (int i=0 ; i<sysdim; i++) {
            A[i] = vector<double> (sysdim);
            C[i] = vector<double> (sysdim);
        }
        
        vector<vector<interval>> CJac(sysdim, vector<interval>(jacdim));
        vector<vector<double>> output_skewbox;
        vector<double> temp(4);
        vector<double> temp2(8);
        vector<double> temp3d(6);
        
        out_approx << YAML::Key << "inner2d";
        out_approx << YAML::Value << YAML::BeginSeq;
    
        for (int i=0 ; i<sysdim; i++) {
            for (int k=i+1 ; k<sysdim; k++) {
                
                out_approx << YAML::BeginMap;
                
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << i;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << k;
                
                // computing and printing inner boxes (max and robust)
                range[0] = ix0[i];   // center and uncontrolled (forall) part
                range[1] = ix0[k];
                compute_print_jointinnerranges(range[0], range[1], Jloc, eps, i, k, false,A);
           
                // computing and printing skew inner boxes  (max and robust)
                build_2dpreconditionner(A,C,Jloc,i,k); // A is center of Jaux on components i and k (otherwise diagonal), C is inverse of A
                multMiMi(CJac,C,Jloc); // CJac = C*Jaux

                range[0] = C[i][i]*ix0[i] + C[i][k]*ix0[k];   //
                range[1] = C[k][k]*ix0[k] + C[k][i]*ix0[i];
                compute_print_jointinnerranges(range[0], range[1], CJac, eps, i, k, true, A);

                 
                out_approx << YAML::EndMap;
                
            }
        }
        
        out_approx << YAML::EndSeq;
       
        
        out_approx << YAML::Key << "outer2d";
        out_approx << YAML::Value << YAML::BeginSeq;
    
        for (int i=0 ; i<sysdim; i++) {
            for (int k=i+1 ; k<sysdim; k++) {
                out_approx << YAML::BeginMap;
                
                out_approx << YAML::Key << "x1";
                out_approx << YAML::Value << i;
                out_approx << YAML::Key << "x2";
                out_approx << YAML::Value << k;
                
                // computing and printing skew inner boxes  (max and robust)
                build_2dpreconditionner(A,C,Jloc,i,k); // A is center of Jaux on components i and k (otherwise diagonal), C is inverse of A
        
                multMiMi(CJac,C,Jloc); // CJac = C*Jaux
                
                //range[0] = C[i][i]*ix0[i] + C[i][k]*ix0[k];   //
                //range[1] = C[k][k]*ix0[k] + C[k][i]*ix0[i];
                range[0] = (C[i][i]*x0p1[i] + C[i][k]*x0p1[k]).convert_int();
                range[1] = (C[k][k]*x0p1[k] + C[k][i]*x0p1[i]).convert_int();
                
                vector<int> indices(2);
                indices[0] = i;
                indices[1] = k;
                vector<interval> max(2);
                
                for (int l=0 ; l<=1; l++)
                {
                    max[l] = range[l];
                    for (int j=0 ; j<jacdim ; j++)
                        max[l] += Jloc[indices[l]][j]*eps[j];
                }
                
                max[0] = intersect(max[0],(C[i][i]*xp1[i] + C[i][k]*xp1[k]).convert_int());
                max[1] = intersect(max[1],(C[k][k]*xp1[k] + C[k][i]*xp1[i]).convert_int());
                compute_print_skewbox(max[0],max[1],A,i,k,"maxskewJac");
                
                // skewed  direct solution
                // attention faux de faire l'intersection de chaque composante, il faudrait faire l'intersection des objets...
                max[0] = (C[i][i]*xp1[i] + C[i][k]*xp1[k]).convert_int();
                max[1] = (C[k][k]*xp1[k] + C[k][i]*xp1[i]).convert_int();
                compute_print_skewbox(max[0],max[1],A,i,k,"maxskew");
                //cout << "max0" << C[i][i]*xp1[i] + C[i][k]*xp1[k];
                //cout << "max1" << C[k][k]*xp1[k] + C[k][i]*xp1[i];
                
                //cout << "xp1[i]=" << xp1[i];
                //cout << "xp1[k]=" << xp1[k];
                //build_2dpreconditionner_frominit(A, C, xp1, i, k);
                //max[0] = (C[i][i]*xp1[i] + C[i][k]*xp1[k]).convert_int();
                //max[1] = (C[k][k]*xp1[k] + C[k][i]*xp1[i]).convert_int();
                //compute_print_skewbox(max[0],max[1],A,i,k,"maxskew");
                
                //compute_print_outerskewbox(range, CJac, eps, i, k, A);
                
                if (uncontrolled > 0)
                {
                    vector<interval> robust(2);
                    vector<interval> initialcondition_pro(2), uncontrolled_impro(2), controlled_pro(2);
                    
                    for (int l=0 ; l<=1; l++)
                    {
                        initialcondition_pro[l] = 0;
                        uncontrolled_impro[l] = 0;
                        controlled_pro[l] = 0;
                    
                        for (int j=0 ; j<jacdim ; j++)
                        {
                            if (j<sysdim) // initial conditions (controlled)
                                initialcondition_pro[l] = initialcondition_pro[l] + Jloc[indices[l]][j]*eps[j];
                            else if (!is_uncontrolled[index_param[j-sysdim]])
                                controlled_pro[l] = controlled_pro[l] + Jloc[indices[l]][j]*eps[j];
                            else // uncontrolled part
                                uncontrolled_impro[l] = uncontrolled_impro[l] + Kaucher_multeps(Jloc[indices[l]][j],eps[j]) ;
                        }
             
                        robust[l] = range[l] + initialcondition_pro[l] + controlled_pro[l];
                        robust[l] = Kaucher_add_pro_impro_resultpro(robust[l],uncontrolled_impro[l]);
                    }
                    compute_print_skewbox(robust[0],robust[1],A,i,k,"robskew");
                }
                
                
                
                out_approx << YAML::EndMap;
            }
        }
        
        out_approx << YAML::EndSeq;
        
        
        
        
        out_approx << YAML::Key << "inner3d";
        out_approx << YAML::Value << YAML::BeginSeq;
        
        // computing joint inner-approximations
        // for any pair (i,k,l) of system outputs which we wish to jointly inner-approximate
        //    component j for j < sysdim/3 is existential in z_i and universal in z_k and z_l
        //    component j for sysdim/3 <= j < 2*sysdim/3 is existential in z_k and universal in z_i and z_l
        //    component j for 2sysdim/3 <= j  is existential in z_l and universal in z_i and z_k
        //    same for inputs
        // note that we could enumerate all possible mappings and superpose all corresponding boxes in the same file ?
       
            interval range_l, range_l_impro;
            interval robust_l, max_l, robust_l_impro, max_l_impro;
            for (int i=0 ; i<sysdim; i++) {
                for (int k=i+1 ; k<sysdim; k++) {
                    for (int l=k+1 ; l<sysdim; l++) {
                        range_i = ix0[i];   // center and uncontrolled (forall) part
                        range_i_impro = 0;      // controlled (existential) part
                        range_k = ix0[k];
                        range_k_impro = 0;
                        range_l = ix0[l];
                        range_l_impro = 0;
                        for (int j=0 ; j<sysdim/3 ; j++) {
                            range_k = range_k + Jloc[k][j]*eps[j];
                            range_l = range_l + Jloc[l][j]*eps[j];
                            range_i_impro = range_i_impro + Kaucher_multeps(Jloc[i][j],eps[j]);
                        }
                        for (int j=sysdim/3 ; j<2*sysdim/3 ; j++) {
                            range_i = range_i + Jloc[i][j]*eps[j];
                            range_l = range_l + Jloc[l][j]*eps[j];
                            range_k_impro = range_k_impro + Kaucher_multeps(Jloc[k][j],eps[j]);
                        }
                        for (int j=2*sysdim/3 ; j<sysdim ; j++) {
                            range_i = range_i + Jloc[i][j]*eps[j];
                            range_k = range_k + Jloc[k][j]*eps[j];
                            range_l_impro = range_l_impro + Kaucher_multeps(Jloc[l][j],eps[j]);
                        }
                        
                        robust_i = 0; max_i = 0;
                        robust_i_impro = 0; max_i_impro = 0;
                        robust_k = 0; max_k = 0;
                        robust_k_impro = 0; max_k_impro = 0;
                        
                        for (int j=0 ; j<fullinputsdim/3 ; j++) {
                            if (is_uncontrolled[index_param[j]])
                                robust_i += Jloc[i][j+sysdim]*eps[j+sysdim];
                            else
                                robust_i_impro += Kaucher_multeps(Jloc[i][j+sysdim],eps[j+sysdim]);
                            max_i_impro += Kaucher_multeps(Jloc[i][j+sysdim],eps[j+sysdim]);
                            
                            robust_k += Jloc[k][j+sysdim]*eps[j+sysdim];
                            max_k += Jloc[k][j+sysdim]*eps[j+sysdim];
                            
                            robust_l += Jloc[l][j+sysdim]*eps[j+sysdim];
                            max_l += Jloc[l][j+sysdim]*eps[j+sysdim];
                        }
                        
                        for (int j=fullinputsdim/3 ; j<2*fullinputsdim/3 ; j++) {
                            robust_i += Jloc[i][j+sysdim]*eps[j+sysdim];
                            max_i += Jloc[i][j+sysdim]*eps[j+sysdim];
                            
                            if (is_uncontrolled[index_param[j]])
                                robust_k += Jloc[k][j+sysdim]*eps[j+sysdim];
                            else
                                robust_k_impro += Kaucher_multeps(Jloc[k][j+sysdim],eps[j+sysdim]);
                            max_k_impro += Kaucher_multeps(Jloc[k][j+sysdim],eps[j+sysdim]);
                            
                            robust_l += Jloc[l][j+sysdim]*eps[j+sysdim];
                            max_l += Jloc[l][j+sysdim]*eps[j+sysdim];
                        }
                        for (int j=2*fullinputsdim/2 ; j<fullinputsdim ; j++) {
                            robust_i += Jloc[i][j+sysdim]*eps[j+sysdim];
                            max_i += Jloc[i][j+sysdim]*eps[j+sysdim];
                            
                            robust_k += Jloc[k][j+sysdim]*eps[j+sysdim];
                            max_k += Jloc[k][j+sysdim]*eps[j+sysdim];
                            
                            if (is_uncontrolled[index_param[j]])
                                robust_l += Jloc[l][j+sysdim]*eps[j+sysdim];
                            else
                                robust_l_impro += Kaucher_multeps(Jloc[l][j+sysdim],eps[j+sysdim]);
                            max_l_impro += Kaucher_multeps(Jloc[l][j+sysdim],eps[j+sysdim]);
                        }
                        
                        robust_i = Kaucher_add_pro_impro(robust_i + range_i, range_i_impro + robust_i_impro);
                        robust_k = Kaucher_add_pro_impro(robust_k + range_k, range_k_impro + robust_k_impro);
                        robust_l = Kaucher_add_pro_impro(robust_l + range_l, range_l_impro + robust_l_impro);
                        
                        max_i = Kaucher_add_pro_impro(max_i + range_i, range_i_impro + max_i_impro);
                        max_k = Kaucher_add_pro_impro(max_k + range_k, range_k_impro + max_k_impro);
                        max_l = Kaucher_add_pro_impro(max_l + range_l, range_l_impro + max_l_impro);
                        
                        // range_i = Kaucher_add_pro_impro(range_i,range_i_impro);
                        // range_k = Kaucher_add_pro_impro(range_k,range_k_impro);
                        // range_l = Kaucher_add_pro_impro(range_l,range_l_impro);
                        
                        out_approx << YAML::BeginMap;
                        
                        out_approx << YAML::Key << "x1";
                        out_approx << YAML::Value << i;
                        out_approx << YAML::Key << "x2";
                        out_approx << YAML::Value << k;
                        out_approx << YAML::Key << "x3";
                        out_approx << YAML::Value << l;
                        
                        out_approx << YAML::Key << "maxbox";
                        
                        temp3d[0] = inf(max_i); temp3d[1] = sup(max_i); temp3d[2] = inf(max_k); temp3d[3] = sup(max_k); temp3d[4] = inf(max_l); temp3d[5] = sup(max_l);
                        out_approx << YAML::Value << temp3d;
                        
                        
                        if (uncontrolled > 0) {
                            temp3d[0] = inf(robust_i); temp3d[1] = sup(robust_i); temp3d[2] = inf(robust_k); temp3d[3] = sup(robust_k); temp3d[4] = inf(robust_l); temp3d[5] = sup(robust_l);
                            out_approx << YAML::Key << "robbox";
                            out_approx << YAML::Value << temp3d;
                        }
                      
                        
                        out_approx << YAML::EndMap;
                    }
                }
            }
        
        out_approx << YAML::EndSeq;
        
    }
    
    
    
}



// specializing JacAff given that we are between [x]^m and [x]^{m+1} in the subdivision
void constraint_eps(vector<vector<interval>> &Jac_m, vector<vector<AAF>> &JacAff, int m, int nb_discr)
{
    vector<interval> c_eps(jacdim);
    
    int index1 = 0, index2 = 1;
    if (component_to_subdiv > -1)
        index1 = component_to_subdiv;
    if (component_to_subdiv2 > -1)
        index2 = component_to_subdiv2;
    
    for (int j=0; j < jacdim ; j++)
        c_eps[j] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr); // interval(-1,1);
    
    if (m == 0)
    {
        /*    c_eps[index1] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        if (jacdim > 1)
            c_eps[index2] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr); */
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++) {
         //              cout << "JacAff[i][j]=" << JacAff[i][j] << endl;
                Jac_m[i][j] = JacAff[i][j].convert_int(c_eps,jacdim);
            }
    }
    else
    {
        
        c_eps[index1] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        if (jacdim > 1)
            c_eps[index2] = interval(-(m+1.0)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = JacAff[i][j].convert_int(c_eps,jacdim);
        
        c_eps[index1] = -interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
        for (int i=0 ; i<sysdim ; i++)
            for (int j=0 ; j<jacdim ; j++)
                Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
        
        if (jacdim > 1) {
            c_eps[index1] = interval(-((double)m)/nb_discr,((double)m)/nb_discr);
            if (jacdim > 1)
                c_eps[index2] = interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
            
            if (jacdim > 1)
                c_eps[index2] = -interval(((double)m)/nb_discr,(m+1.0)/nb_discr);
            for (int i=0 ; i<sysdim ; i++)
                for (int j=0 ; j<jacdim ; j++)
                    Jac_m[i][j] = interval_hull(Jac_m[i][j],JacAff[i][j].convert_int(c_eps,jacdim));
        }
    }
}





// same as InnerOuter but with quadrature
void InnerOuter_discretize(vector<interval> &Xinner, vector<interval> &Xinner_robust, vector<interval> &Xouter, vector<interval> &Xouter_robust, vector<AAF> &x0p1, vector<vector<AAF>> &Jtau, vector<interval> &eps, double tnp1)
{
    vector<vector<interval>> Jaux(sysdim, vector<interval>(jacdim)); // will copy only the sysdim relevant lines of Jaux
    vector<interval> ix0(sysdim);
    vector<interval> initialcondition_impro(sysdim), initialcondition_pro(sysdim), uncontrolled_impro(sysdim), uncontrolled_pro(sysdim), controlled_impro(sysdim), controlled_pro(sysdim), aux_impro(sysdim);
    vector<interval> loc_eps(jacdim);
    
    for (int i=0 ; i<sysdim ; i++) {
        ix0[i] = x0p1[i].convert_int();
        for (int k=0 ; k<jacdim ; k++)
            Jaux[i][k] = Jtau[i][k].convert_int();
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        initialcondition_impro[i] = 0;
        initialcondition_pro[i] = 0;
        uncontrolled_impro[i] = 0;
        uncontrolled_pro[i] = 0;
        controlled_impro[i] = 0;
        controlled_pro[i] = 0;
        
        Xinner_robust[i] = 0;
        Xouter_robust[i] = 0;
        Xouter[i] = ix0[i];
    }
    
//    int nb_discr = 1;
    int nb_discr = 10;
    
    for (int j=0 ; j<jacdim; j++)
        loc_eps[j] = eps[j]/((double)nb_discr);
    
    
    for (int m=0; m<nb_discr; m++)
    {
        constraint_eps(Jaux,Jtau,m,nb_discr);
        
        // inner approx: special addition corresponding to addition of proper and improper intervals: the result is an inner range for the solution of ODE system
        for (int i=0 ; i<sysdim; i++) {
            
            for (int j=0 ; j<jacdim ; j++)
            { // adding the controlled (exists) part
                if (j<sysdim) // initial conditions //(is_initialcondition[j])
                {
                    initialcondition_impro[i] = initialcondition_impro[i] + Kaucher_multeps(Jaux[i][j],loc_eps[j]);
                    initialcondition_pro[i] = initialcondition_pro[i] + Jaux[i][j]*loc_eps[j];
                }
                else if (!is_uncontrolled[index_param[j-sysdim]])
                {
                    controlled_impro[i] = controlled_impro[i] + Kaucher_multeps(Jaux[i][j],loc_eps[j]);
                    controlled_pro[i] = controlled_pro[i] + Jaux[i][j]*loc_eps[j];
                }
                else // uncontrolled
                {
                    uncontrolled_impro[i] = uncontrolled_impro[i] + Kaucher_multeps(Jaux[i][j],loc_eps[j]) ;
                    uncontrolled_pro[i] = uncontrolled_pro[i] + Jaux[i][j]*loc_eps[j] ;
                }
                
                Xouter[i] += Jaux[i][j]*loc_eps[j] ;
            }
           
        }
    }
    
    for (int i=0 ; i<sysdim ; i++) {
        aux_impro[i] = initialcondition_impro[i] + controlled_impro[i];
        Xinner[i] = Kaucher_add_pro_impro(ix0[i] , aux_impro[i] + uncontrolled_impro[i]);
        
        if (uncontrolled > 0)
        {
            Xinner_robust[i]  = ix0[i] + uncontrolled_pro[i];
            Xinner_robust[i] = Kaucher_add_pro_impro(Xinner_robust[i], aux_impro[i]);
            Xouter_robust[i] = ix0[i] + initialcondition_pro[i] + controlled_pro[i];
            Xouter_robust[i] = Kaucher_add_pro_impro_resultpro(Xouter_robust[i],uncontrolled_impro[i]);
        }
        
    }
    
    // outer-approx computation : Xouter = x0p1 + Jtau*eps;
  //  multMiVi(Xouter,Jaux,eps);
  //  addViVi(Xouter,ix0);
}
