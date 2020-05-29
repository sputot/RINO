/*
 * aa_aaftrigo.c -- Trigonometric operations (all non-affine)
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 * Copyright (C) 2009 LUH (Leibniz Universitaet Hannover)
 * Modifications 2018 Franck Djeumou (bug corrections)
 *
 * This file is part of aaflib.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with libaa; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include "aa.h"
#include "aa_aaf.h"
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <iostream>

#ifndef PI
#define PI (4*atan(1.0))
#endif

// number of points for the linear regression approximation
#define NPTS 8 


/************************************************************
 * Method:        sin
 * Author & Date: ??? - ???
 * Description:   
 *   Sine function
 *   sine isn't montonic and the second derivative change its 
 *   sign inside the defined interval so we use the 
 *   approximation by least squares
 *
 *   Input  : double : corresponding center value
 *   Output : -
 ************************************************************/
AAF sin(const AAF & P)
{

  double a, b;
  double fa, fb;
  double alpha, dzeta, delta;

  if ((P.length == 0) || (P.length == 1 && P.deviations[0] == 0.))
  {
    AAF Temp(sin(P.cvalue));
    return Temp;
  }

  double r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

    // Sylvie: j'ai ajouté l'assert car le truc ci-dessous buggait pour un petit intervalle juste sous 2PI.
    // a verifier
// il y a toujours un souci mais j'ai commenté car cette assertion est trop violente, a revoir a l'occasion,
    // assert ( a > -PI && b < PI);
    
    // translation ci-dessous a reverifier a l'occasion
  if (a > PI || a < -PI){
    double temp = floor(a/(2*PI));
    a -= temp * 2 * PI;
    a -= (a > PI) ? 2*PI : 0.0;
    b = a + 2*r;
  }

  //assert ( a > -PI/2.0 && b < PI/2.0);
  /*if (!(a > -PI/2.0 && b < PI/2.0))
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
          std::string("Error <AAF::sin>: a > -PI/2.0 && b < PI/2.0").c_str(), 
      a, b));*/


  fa = sin(a);
  fb = sin(b);

  if (b < 0 || (b < PI && a >0)){ // chebyshev apprix in this case
    if (r > 1.0e-6){
      alpha = (fb - fa) / (b-a);
      double sol = ((a > 0) - (a < 0)) * acos(alpha);
      double fsol = sqrt(1 - alpha*alpha); // fast computation of sin(acos)
      dzeta = (fa + fsol - alpha * (a + sol)) / 2;
      delta = abs(fsol - fa - alpha * (sol - a)) / 2;
    } else {
      alpha = cos(a);
      dzeta =  fa - alpha*a;
      delta =  0.0;
    }
  } else { // min range optimization since the derivative cant be zero because of the radius constraints
    if (a <= 0){
      alpha = 1;
      delta = (-fb + fa - alpha * (a-b)) / 2;
    } else {
      alpha = -1;
      delta = (fb - fa + alpha * (a-b)) / 2;
    }
    dzeta = (fa + fb - alpha * (a + b)) / 2;
  }

/*  // Min range approximation on this interval
  double theta = a;
  alpha = cos(theta);

  dzeta = 0.5 * (fa + fb - alpha * (a + b));
  delta = 0.5 * (fb - fa - alpha * (b - a));*/

  /*AAInterval i = P.convert();

  double w = i.width();

  //i = mintrigo(i); // no more needed

  a = i.getlo();
  b = i.gethi();


  // y' = alpha*x+dzeta , the regression line
  // approximate y = sin(x)

  double x[NPTS];
  double y[NPTS];
  double r[NPTS]; // residues, r[i] = y[i]-y'[i]

  double xm = 0;
  double ym = 0;

  // the trivial case, the interval is larger than 2*PI
  if (w >= 2*PI ) 
  {
    // y' = 0 , delta = 1 cause -1 <= sin(x) <= +1
    alpha = 0.0;
    dzeta = 0.0;
    delta = 1.0;
  }
  else // case of the least squares
  {
    x[0]=a;
    y[0]=sin(a);
    x[NPTS-1]=b;
    y[NPTS-1]=sin(b);
    
    double pas=w/(NPTS-1);
    
    for (unsigned j=1; j< NPTS-1; j++)
    {
      x[j]=x[j-1]+pas;
      y[j]=sin(x[j]);
    }
    
    
    // Calculation of xm and ym , averages of x and y
    
    for (unsigned j=0; j<NPTS; j++)
    {
      xm=xm+x[j];
      ym=ym+y[j];
    }
    
    xm=xm/NPTS;
    ym=ym/NPTS;
    
    // Calculation of alpha and dzeta
    
    double temp1;
    double temp2=0;
    alpha = 0;
    
    for (unsigned j=0; j<NPTS; j++)
    {
      temp1=x[j]-xm;
      alpha+=y[j]*temp1;
      temp2+=temp1*temp1;
    }    
    
    alpha=alpha/temp2;  // final alpha
    dzeta=ym-alpha*xm; // final dzeta
    
    
    // Calculation of the residues
    // We use the absolute value of the residues!
    
    for (unsigned j=0; j<NPTS; j++)
    {
      r[j]=fabs(y[j]-(dzeta+alpha*x[j]));
    }
    
    
    // The error delta is the maximum
    // of the residues (in absolute values)
    
    double *ptr;
    
    ptr = std::max_element(r, r+NPTS);
    
    delta = *ptr;
    
  }*/
    
  // z0 = alpha*x0 + dzeta
  
  AAF Temp(alpha*(P.cvalue)+dzeta);
  
  Temp.length=(P.length)+1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];
  
  // zi = alpha*xi
  
  for (unsigned j=0; j < P.length; j++)
  {
    Temp.indexes[j]=P.indexes[j];
    Temp.deviations[j]=alpha*(P.deviations[j]);
  }
  
  // zk = delta
  
  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.deviations[P.length]=delta;
  
  return Temp;
}


/************************************************************
 * Method:        cos
 * Author & Date: ??? - ???
 * Description:   
 *   Cosine function
 *   we use the identity cos(x)=sin(x+PI/2)
 *
 *   Input  : double : corresponding center value
 *   Output : -
 ************************************************************/
AAF cos(const AAF & P)
{

  return sin( P + PI/2.0);

  /*double a, b;
  double fa, fb;
  double alpha, dzeta, delta;

  if ((P.length == 0) || (P.length == 1 && P.deviations[0] == 0.))
  {
    AAF Temp(cos(P.cvalue));
    return Temp;
  }

  double r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  //assert ( a > -PI/2.0 && b < PI/2.0);
  if (!(a > -PI/2.0 && b < PI/2.0))
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
          std::string("Error <AAF::cos>: a > -PI/2.0 && b < PI/2.0").c_str(), 
      a, b));

  fa = cos(a);
  fb = cos(b);

  // chebyshev approximation on this interval
  if  ( r > 1.0e-6) 
    alpha = (fb - fa) / (b - a);
  else
    alpha = - sin(a);

  double sol = asin(-alpha);
  // std::cout << " a = " << a << " ; sol1 = " << sol << " ; b = " << b << std::endl;
  if ( ! (a <= sol && sol <= b) )
    sol = (a + b) / 2.0;

  double fsol = cos(sol);

  dzeta = 0.5 * (fa + fsol - alpha * (a + sol));
  delta = 0.5 * fabs(fsol - fa - alpha * (sol - a));

  AAF Temp(alpha*(P.cvalue)+dzeta);
  
  Temp.length=(P.length)+1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];
  
  // zi = alpha*xi
  
  for (unsigned j=0; j < P.length; j++)
  {
    Temp.indexes[j]=P.indexes[j];
    Temp.deviations[j]=alpha*(P.deviations[j]);
  }
  
  // zk = delta
  
  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.deviations[P.length]=delta;
  
  return Temp;*/
  // AAF Temp = P;
  // return sin(Temp+PI/2);
}


/************************************************************
 * Method:        tan
 * Author & Date: ??? - ???
 * Description:   
 *   Tangent function
 *   we use the identity tan(x)=sin(x)/cos(x)
 *   Due to the nature of the tan fct remember that
 *   we can have infinite value with small intervals
 *
 *   Input  : double : corresponding center value
 *   Output : -
 ************************************************************/
AAF tan(const AAF & P)
{
  double a, b;
  double fa, fb;
  double alpha, dzeta, delta;

  if ((P.length == 0) || (P.length == 1 && P.deviations[0] == 0.))
  {
    AAF Temp(tan(P.cvalue));
    return Temp;
  }

  double r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  // In case the angles are not between -PI/2 and PI/2
  if (a>= PI/2 || a <= -PI/2){
    double temp = floor(a/PI);
    a -= temp * PI;
    a -= (a > PI/2) ? PI : 0.0;
    b = a + 2*r;
  }

  if ((a <= PI/2 && b >= PI/2) || (a <= -PI/2 && b>= -PI/2)){
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
          std::string("Error <AAF::tan>: (a <= PI/2 && b >= PI/2) || (a <= -PI/2 && b>= -PI/2) ").c_str(), 
      a, b));
  }

  //assert ( a > -PI/2.0 && b < PI/2.0);
  /*if (!(a > -PI/2.0 && b < PI/2.0))
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
          std::string("Error <AAF::tan>: a > -PI/2.0 && b < PI/2.0").c_str(), 
      a, b));*/

  fa = tan(a);
  fb = tan(b);

  double sol;

  if ( r > 1.0e-6)
    alpha = (fb - fa)/(b-a);
  else
    alpha = 1.0 + fa * fb;

  double aux = sqrt(alpha - 1.0);
  double eps = atan(aux);

  if (eps >= a && eps <= b && -eps >= a && -eps <= b){
    dzeta = 0.0;
    delta = abs(aux - alpha * eps);
  } else {
    if (-eps >= a && -eps <= b){
      eps = -eps;
      aux = -aux;
    }
    dzeta = (fa + aux - alpha*(a+eps)) / 2;
    delta = abs(aux - fa - alpha * (eps -a )) / 2;
  }

  /*double sol1 = atan(sqrt(alpha - 1.0));
  double sol2 = atan(-sqrt(alpha - 1.0));
  // std::cout << " a = " << a << " ; sol1 = " << sol1 << " ; b = " << b << std::endl;

  if (a <= sol1 && sol1 <= b)
    sol = sol1;
  else if (a <= sol2 && sol2 <= b)
    sol = sol2;
  else
    sol = (a + b) / 2.0;

  double fsol = tan(sol);

  if (a * b < 0){
    dzeta = 0.0;
    delta = fabs(fsol - alpha * sol);
  } else {
    dzeta = 0.5 * (fa + fsol - alpha * (a + sol));
    delta = 0.5 * fabs(fsol - fa - alpha * (sol - a));
  }*/

  AAF Temp(alpha*(P.cvalue)+dzeta);
  
  Temp.length=(P.length)+1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];
  
  // zi = alpha*xi
  
  for (unsigned j=0; j < P.length; j++)
  {
    Temp.indexes[j]=P.indexes[j];
    Temp.deviations[j]=alpha*(P.deviations[j]);
  }
  
  // zk = delta
  
  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.deviations[P.length]=delta;
  
  return Temp;

  //i = mintrigo(i); // no more needed

  // a = i.getlo();
  // b = i.gethi();

  /*a = a - ((round(a/(2.0*PI))) * 2.0 * PI );
  b = b - ((round(b/(2.0*PI))) * 2.0 * PI );

  if (fabs(a - PI) < 1e-15)
    a = a > 0 ? -a : a;

  if (fabs(b-PI) < 1e-15)
    b = b > 0 ? -b : b;

  double w = (b - a)/2.0;

  if ( w < 0 || w >= PI || (w > 0 && ((a <= PI/2.0 && PI/2.0 <= b) || (a <= -PI/2.0 && -PI/2.0 <= b))) ){
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
          std::string("Error <AAF::tan>: computing tangente function with kPi/2 in the set").c_str(), 
      a, b));
  }

  fa = tan(a);
  fb = tan(b);

  if ( w > 1e-10)
    alpha = (fb - fa)/(b-a);
  else
    alpha = 1 + tan(a)*tan(a);

  double sol1 = acos(1.0/sqrt(alpha));
  double sol2 = - sol1;
  double sol3 = acos(-1.0/sqrt(alpha));
  double sol4 = - sol3;

  double sol;

  if (a <= sol1 && sol1 <= b)
    sol = sol1;
  else if (a <= sol2 && sol2 <= b)
    sol = sol2;
  else if (a <= sol3 && sol3 <= b)
    sol = sol3;
  else if (a <= sol4 && sol4 <= b)
    sol = sol4;
  else 
    assert(false);

  double fsol = tan(sol);

  if ( a * b < 0){ // 0 is in the interval a , b
    dzeta = 0.0;
    delta = fabs(fsol - alpha * sol);
  } else { // The interval doesn't contain 0 so function is concave or convex
    dzeta = 0.5 * (fa + fsol - alpha * (a + sol));
    delta = 0.5 * fabs(fsol - fa - alpha * (sol - a));
  }

  AAF Temp(alpha*(P.cvalue)+dzeta);
  
  Temp.length=(P.length)+1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];
  
  // zi = alpha*xi
  
  for (unsigned j=0; j < P.length; j++)
  {
    Temp.indexes[j]=P.indexes[j];
    Temp.deviations[j]=alpha*(P.deviations[j]);
  }
  
  // zk = delta
  
  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.deviations[P.length]=delta;
  
  return Temp;*/

  // return sin(P)/cos(P);
}


/************************************************************
 * Method:        tan
 * Author & Date: ??? - ???
 * Description:   
 *   Cotangent function
 *   we use the identity cotan(x)=cos(x)/sin(x)
 *   Due to the nature of the cotan fct remember that
 *   we can have infinite value with small intervals
 *
 *   Input  : double : corresponding center value
 *   Output : -
 ************************************************************/
AAF cotan(const AAF & P)
{
  return cos(P)/sin(P);
}
