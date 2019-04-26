/*
 * aa_aafinnerflow.c -- 
 Sylvie Putot - nov 2016 */


#include "aa.h"

#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <cstdio>
#include <assert.h>

using namespace std;


// Ajout SP
AAF::AAF(const interval iv):
cvalue((inf(iv)+sup(iv))/2),
length(1),
size(1),
#ifdef FAST_RAD
radius(0.0),
#endif
deviations(new double[1]),
indexes(new unsigned[1])
{
    deviations[0] = (sup(iv)-inf(iv))/2;
    
    indexes[0] = inclast();
    
#ifdef FAST_RAD
    radius = fabs(deviations[0]);
#endif
    
#ifdef CLEANUP
    allAAF.push_back(this);
#endif
    
}


// afficher la forme affien sur une ligne, ou n est le numero du dernier symbol eta_n
void AAF::aafprintline(stringstream &TM_stream, int n) const
{
    TM_stream << " ( " << cvalue ;
    for (unsigned i = 0; i < length-1 ; i++)
        TM_stream << " + " << deviations[i] << "*e" << indexes[i];
    TM_stream << " + " << deviations[length-1] << "*n" << n;
    TM_stream << " )";
}

void AAF::aafprintline() const
{
    cout << cvalue ;
    for (unsigned i = 0; i < length ; i++)
        cout << " + " << deviations[i] << "*e" << indexes[i];
    cout << endl;
}

// reduce affine form knowing that last noise symbol is within I instead of [-1,1]
void AAF::reduce_aaf(interval I) {
    interval temp = deviations[length-1]*I;
    cvalue = cvalue + filib::mid(temp);
    deviations[length-1] = filib::rad(temp);
}


interval AAF::convert_int() const
{
    // lower bound == central value of the AAF - the total deviation
    // upper bound == central value of the AAF + the total deviation
    
  interval Temp = interval(cvalue-rad(), cvalue+rad()); //filib::interval<double,filib::native_switched,filib::i_mode_extended>(cvalue-rad(), cvalue+rad());
    return Temp;
}

// suppose we have some constraints on the noise symbol eps_i up to sysdim
interval AAF::convert_int(vector<interval> &constr_eps, int sysdim) const
{
    interval Temp = interval(cvalue,cvalue);
    interval temp2 = interval(-1,1);
    unsigned dimlength = sysdim;
    if (length <sysdim)
        dimlength = length;
    for (unsigned i = 0; i < dimlength; i++)
        Temp = Temp + deviations[i]*constr_eps[i];
    for (unsigned i = sysdim; i < length; i++)
    Temp = Temp + deviations[i]*temp2;
    
    return Temp;
}


// suppose we have some constraints on the noise symbols eps_i, given by some affine forms being <= 0
interval AAF::convert_int(vector<AAF> &constraints, int sysdim) const
{
 //   cout << "this=" << (*this) << endl;
    interval res = convert_int();
  //  cout << "res before=" << res << endl;
    AAF aux;
    
    for (int j=0 ; j<constraints.size(); j++)
    {
        cout << "constraints[j] non empty" << endl; // << constraints[j] << endl;
        unsigned l1 = length;
        unsigned l2 = constraints[j].length;
        
        if (l1 == 0 || l2 == 0)
        { // no common symbols
            continue;
        }
        
        unsigned * id1 = indexes;
        unsigned * id2 = constraints[j].indexes;
        
        double * va1 = deviations;
        double * va2 = constraints[j].deviations;
        
        unsigned * pu1 = id1;
        unsigned * pu2 = id2;
        
        AAF Temp(0.0); // create our resulting AAF and initialize will center
        Temp.indexes = new unsigned [l1+l2+1];
        unsigned * idtemp=Temp.indexes;
        
        // Fill the indexes array
        unsigned * tot = std::set_union(id1, id1 + l1, id2, id2 + l2, idtemp);
        unsigned ltemptot = tot - idtemp;
        
        Temp.deviations = new double [ltemptot + 1];
        double * vatempg = Temp.deviations;
        
        Temp.length = ltemptot + 1;
        Temp.size = Temp.length;
        
        // interval center = hull(P1.cvalue,P2.cvalue);
        // Fill the deviations array
        double deviation = 0.0; // rad(center);
        
        
        for (unsigned int i = 0; i < ltemptot; i++) // iterate on all indexes
        {
            unsigned a = pu1 - id1;
            unsigned b = pu2 - id2;
            
            if (id1[a] == idtemp[i] && id2[b] == idtemp[i])
            {
                // substitute only input noise symbols (?)
                if (idtemp[i] > sysdim)
                    break;
                
              //  cout << "id1[a]=" << id1[a] << " id2[b]=" << id2[b] << endl;
             //   cout << "va1[a]=" << va1[a] << " va2[b]=" << va2[b] << endl;
                
                if (va1[a]*va2[b] < 0.0)
                {
                    aux = (*this) + fabs(va1[a]/va2[b]) * constraints[j]; // such that concretization of this <= concretization of aux
             //      cout << "aux1=" << aux.convert_int() << endl;
                  //  cout << "aux1=" << aux << endl;
                    
                    res = interval(fmax(inf(res),inf(aux.convert_int())),sup(res));
                   // res = interval(inf(res),fmin(sup(res),sup(aux.convert_int())));
                  //  cout << "res1=" << res << endl;
                }
                else
                {
                    aux = (*this) - fabs(va1[a]/va2[b]) * constraints[j]; // such that concretization of this >= concretization of aux
              //      cout << "aux2=" << aux.convert_int() << endl;
               //     cout << "aux2=" << aux << endl;
                    
                    res = interval(inf(res),fmin(sup(res),sup(aux.convert_int())));
                   // res = interval(fmax(inf(res),inf(aux.convert_int())),sup(res));
               //     cout << "res2=" << res << endl;

                }
                
                pu1++;
                pu2++;
            }
            else if (id1[a] == idtemp[i])
            {
                pu1++;
            }
            else if (id2[b] == idtemp[i])
            {
                pu2++;
            }
            
        }
        
        
    }
 //   cout << "res after=" << res << endl;
    return res;
}



// concret (sum alpha_i eps_i) - sans le centre
interval AAF::rad(vector<interval> &constr_eps, int sysdim) const {
    interval temp = interval(-1,1);
    interval concret = 0;
    unsigned dimlength = sysdim;
    if (length <sysdim)
        dimlength = length;
    for (unsigned i = 0; i < dimlength; i++) {
        concret = concret + deviations[i]*constr_eps[i];
    }
    for (unsigned i = sysdim; i < length; i++)
        concret = concret + deviations[i]*temp;
        return concret; //max(sup(concret),-inf(concret));
    
}

/************************************************************
 * INTERPOLATION INFO
 *
 * a is the lower bound of the argument AAF
 * b is the upper bound of the argument AAF
 *
 * alpha is the slope of the line r(x) that interpolates
 * i)  CHEBYSHEV: (a, f(a)) and (b, f(b))
 * ii) MINRANGE:  MIN(f'(a), f'(b))
 *
 * dzeta is the central y-intercept of the parallelogram
 * including (a,f(a) .. b,f(b))
 *
 * delta is the new deviation value being half a hight of
 * the parallelogram
 ************************************************************/

/************************************************************
 multiplication, given some interval constraints on the noise symbols
 ************************************************************/
AAF mult_eps (const AAF & P1, const AAF & P2, vector<interval> &constr_eps, int sysdim) 
{
    unsigned l1 = P1.length;
    unsigned l2 = P2.length;
    
    if (l1+l2 == 0)
    {
        AAF Temp(P1.cvalue*P2.cvalue);
        return Temp;
    }

   if (l1 == 0)
    {
        // if *this is double
        AAF Temp(P2);
        Temp *= P1.cvalue;
        return Temp;
    }
    if (l2 == 0)
    {
        // if P is double
        AAF Temp(P1);
        Temp *= P2.cvalue;
        return Temp;
    }
    
    unsigned * id1 = P1.indexes;
    unsigned * id2 = P2.indexes;
    
    double * va1 = P1.deviations;
    double * va2 = P2.deviations;
    
    unsigned * pu1 = id1;
    unsigned * pu2 = id2;
    
    AAF Temp(P1.cvalue*P2.cvalue);  // Create our resulting AAF
    
    Temp.indexes = new unsigned [l1+l2+1];
    
    unsigned * idtemp=Temp.indexes;
    
    // Fill the indexes array
    
    unsigned * fin = std::set_union(id1, id1 + l1, id2, id2 + l2, idtemp);
    unsigned ltemp = fin - idtemp;
    
    Temp.deviations = new double [ltemp + 1];
    double * vatempg = Temp.deviations;
    
    Temp.length = ltemp + 1;
    Temp.size = Temp.length;
    
    double commonTermCenter = 0.0;
    double commonTermDeviation = 0.0;
    
    // Fill the deviations array
    
    for (unsigned i = 0; i < ltemp; i++)
    {
        unsigned a = pu1 - id1;
        unsigned b = pu2 - id2;
        
        if (a == l1 || id1[a] != idtemp[i])
        {
            vatempg[i] = P1.cvalue*va2[b];  // P1.cvalue*va2[b]+(P2.cvalue)*0
            pu2++;
            continue;
        }
        
        if (b == l2 || id2[b] != idtemp[i])
        {
            vatempg[i] = (P2.cvalue)*va1[a];  // P1.cvalue*0+(P2.cvalue)*va1[a]
            pu1++;
            continue;
        }
        
        vatempg[i] = P1.cvalue*va2[b] + (P2.cvalue)*va1[a];
        commonTermCenter += va2[b]*va1[a];
        commonTermDeviation += fabs(va2[b]*va1[a]);
        pu1++;
        pu2++;
    }
    
    // Compute the error
    // in a new deviation symbol
   // double delta = P1.rad(constr_eps,sysdim)*P2.rad(constr_eps,sysdim);
    interval delta = P1.rad(constr_eps,sysdim)*P2.rad(constr_eps,sysdim); // concret (sum alpha_i eps_i) - sans le centre
    
    Temp.indexes[ltemp] = Temp.inclast();
    
    assert (AAF::approximationType != SECANT);
    
    // consider deviations occuring in both expressions  
   // Temp.cvalue += 0.5*commonTermCenter;
    Temp.cvalue = Temp.cvalue + mid(delta);
    Temp.deviations[ltemp] = rad(delta); // - 0.5*commonTermDeviation;
        return Temp;
}


/* "arg-min" join of affine forms */
AAF hull  (const AAF & P1, const AAF & P2) 
{
    unsigned l1 = P1.length;
    unsigned l2 = P2.length;
    
    if (l1 == 0 || l2 == 0)
    { // no relation left: just join interval concretizations
        AAF Temp(hull(P1.convert_int(),P2.convert_int()));
        return Temp;
    }

    unsigned * id1 = P1.indexes;
    unsigned * id2 = P2.indexes;
    
    double * va1 = P1.deviations;
    double * va2 = P2.deviations;
    
    unsigned * pu1 = id1;
    unsigned * pu2 = id2;

    interval center;
    center = hull(P1.convert_int(),P2.convert_int());
        
    AAF Temp(mid(center)); // create our resulting AAF and initialize will center
    
    Temp.indexes = new unsigned [l1+l2+1];
    unsigned * idtemp=Temp.indexes;
    
    // Fill the indexes array
    unsigned * tot = std::set_union(id1, id1 + l1, id2, id2 + l2, idtemp);
    unsigned ltemptot = tot - idtemp;
    
    Temp.deviations = new double [ltemptot + 1];
    double * vatempg = Temp.deviations;
    
    Temp.length = ltemptot + 1;
    Temp.size = Temp.length;
    
    // interval center = hull(P1.cvalue,P2.cvalue);
    // Fill the deviations array
    double deviation = 0.0; // rad(center);
    double min, max;
    
   // cout << "P1= " << P1.convert_int();
   // cout << "P1= " << P1;
   // cout << "P2= " << P2.convert_int();
   // cout << "P2= " << P2;
   // cout << "deviation= " << deviation << " ltemptot = " << ltemptot << endl;
    
    for (unsigned int i = 0; i < ltemptot; i++) // iterate on all indexes
    {
        unsigned a = pu1 - id1;
        unsigned b = pu2 - id2;
        
   //     cout << "a=" << a << " b=" << b << " i=" << i << endl;
   //     cout << "id1[a]=" << id1[a] << " id2[b]=" << id2[b] << " idtemp[i]" << idtemp[i] << endl;

        if (id1[a] == idtemp[i] && id2[b] == idtemp[i])
        {
            if (va1[a]*va2[b] < 0.0)
            { // no common relation kept
                vatempg[i] = 0.0;
               // deviation += fmax(fabs(va1[a]),fabs(va2[b]));
            }
            else if (va1[a] > 0.0) {
                min = fmin(va1[a],va2[b]);
                max = fmax(va1[a],va2[b]);
                vatempg[i] = min;
                deviation += fabs(min);
               // deviation += max-min;
            }
            else // both negative coeff
            {
                min = fmin(va1[a],va2[b]);
                max = fmax(va1[a],va2[b]);
                vatempg[i] = max;
                deviation += fabs(max);
               // deviation += max-min;
            }
           
            pu1++;
            pu2++;
        }
        else if (id1[a] == idtemp[i])
        {
           // deviation += fabs(va1[a]);
            vatempg[i] = 0.0;
            pu1++;
        }
        else if (id2[b] == idtemp[i])
        {
          //  deviation += fabs(va2[b]);
            pu2++;
            vatempg[i] = 0.0;
        }
        else
            vatempg[i] = 0.0;
        
//        cout << "va1[a]=" << va1[a] << " va2[b]=" << va2[b] << " vatempg[i]=" << vatempg[i] << " deviation=" << deviation << endl ;
    }
    
    
    
    Temp.indexes[ltemptot] = Temp.inclast();
    Temp.deviations[ltemptot] = rad(center)-deviation;
    
    Temp.compact();
    
  //   cout << "Temp=" <<  Temp.convert_int();
  //  cout << "Temp=" << Temp;
    
    return Temp;
}
