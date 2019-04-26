/*
 * aa_aafapprox.c -- Standart non-affine operations
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 * Copyright (C) 2009 LUH (Leibniz Universitaet Hannover)
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
#include <cmath>
#include <algorithm>
#include <cstdio>

//#define AAF_DEBUG
#define EPSILON 1E-20
#define AAF_MINRAD 1E-10

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
 * Operator:      *
 * Author & Date: ??? - ???
 * Description:   
 *   computes the affine product (*this * P)
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF AAF::operator * (const AAF & P) const
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    AAF Temp(cvalue*P.cvalue);
    return Temp;
  }
  if (l1 == 0)
  {
    // if *this is double
    AAF Temp(P);
    Temp *= cvalue;
    return Temp;
  }
  if (l2 == 0)
  {
    // if P is double
    AAF Temp(*this);
    Temp *= P.cvalue;
    return Temp;
  }

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  double * va1 = deviations;
  double * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  AAF Temp(cvalue*P.cvalue);  // Create our resulting AAF

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
      vatempg[i] = cvalue*va2[b];  // cvalue*va2[b]+(P.cvalue)*0
      pu2++;
      continue;
    }
    
    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = (P.cvalue)*va1[a];  // cvalue*0+(P.cvalue)*va1[a]
      pu1++;
      continue;
    }
    
    vatempg[i] = cvalue*va2[b] + (P.cvalue)*va1[a];
    commonTermCenter += va2[b]*va1[a];
    commonTermDeviation += fabs(va2[b]*va1[a]);
    pu1++;
    pu2++;
  }

  // Compute the error
  // in a new deviation symbol
  double delta = rad()*P.rad();

  Temp.indexes[ltemp] = inclast();
  if (AAF::approximationType == SECANT)
  {
    commonTermCenter *= 0.5;
    commonTermDeviation *= 0.5;

    delta -= commonTermDeviation;
    // new deviation equals zero
    Temp.deviations[ltemp] = 0.;
    // distribute delta on current deviations
    double r = 0.0;
    for (unsigned i = 0; i < ltemp; i++)
      r += fabs(Temp.deviations[i]);
    double fac;
    if (fabs(r) < EPSILON)
      fac = 1.;
    else
      fac = 1.0 + delta/r;
    for (unsigned i = 0; i < ltemp; i++)
      Temp.deviations[i] *= fac;
#ifdef FAST_RAD
    Temp.radius = r + delta;
#endif
    Temp.cvalue += commonTermCenter;
    return Temp;
  }
  else
  {
    Temp.deviations[ltemp] = delta;
  }
  
  // consider deviations occuring in both expressions  
  commonTermCenter *= 0.5;
  commonTermDeviation *= 0.5;
  Temp.cvalue += commonTermCenter;
  Temp.deviations[ltemp] -= commonTermDeviation;

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < Temp.length; i++)
    Temp.radius += fabs(Temp.deviations[i]);
#endif
  
  return Temp;
}


/************************************************************
 * Operator:      *=
 * Author & Date: Darius Grabowski - 2006
 * Description:   
 *   computes the affine product (*this *= P)
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF &       : AAF result
 ************************************************************/
AAF & AAF::operator *=(const AAF & P)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1 == 0)
  {
    // if *this is double
    double temp = cvalue;
    *this = P;
    *this *= temp;
    return *this;
  }
  if (l2 == 0)
  {
    // if P is double
    *this *= P.cvalue;
    return *this;
  }

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  double * va1 = deviations;
  double * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  unsigned * tempIndexes = new unsigned [l1+l2+1];
  
  unsigned * idtemp = tempIndexes;

  // Fill the indexes array
  unsigned * fin = std::set_union(id1, id1 + l1, id2, id2 + l2, idtemp);
  unsigned ltemp = fin - idtemp;

  double * tempDeviations = new double [ltemp + 1];
  double * vatempg = tempDeviations;
 
  double commonTermCenter = 0.0;
  double commonTermDeviation = 0.0;

  // Fill the deviations array

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1 - id1;
    unsigned b = pu2 - id2;
    
    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = cvalue*va2[b];  // cvalue*va2[b]+(P.cvalue)*0
      pu2++;
      continue;
    }
    
    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = (P.cvalue)*va1[a];  // cvalue*0+(P.cvalue)*va1[a]
      pu1++;
      continue;
    }
    
    vatempg[i] = cvalue*va2[b] + (P.cvalue)*va1[a];
    commonTermCenter += va2[b]*va1[a];
    commonTermDeviation += fabs(va2[b]*va1[a]);
    pu1++;
    pu2++;
  }

  // Compute the error in a new deviation symbol
  double delta = rad()*P.rad();

  // multiply cvalue
  cvalue *= P.cvalue;

  tempIndexes[ltemp] = inclast();
  if (AAF::approximationType == SECANT)
  {
    commonTermCenter *= 0.5;
    commonTermDeviation *= 0.5;

    delta -= commonTermDeviation;
    // new deviation equals zero
    tempDeviations[ltemp] = 0.;
    // distribute delta on current deviations
    double r = 0.0;
    for (unsigned i = 0; i < ltemp; i++)
      r += fabs(tempDeviations[i]);
    double fac;
    if (fabs(r) < EPSILON)
      fac = 1.;
    else
      fac = 1.0 + delta/r;
    for (unsigned i = 0; i < ltemp; i++)
      tempDeviations[i] *= fac;
#ifdef FAST_RAD
    radius = r + delta;
#endif
    cvalue += commonTermCenter;
  }
  else
  {
    tempDeviations[ltemp] = delta;
    // consider deviations occuring in both expressions  
    commonTermCenter *= 0.5;
    commonTermDeviation *= 0.5;
    cvalue += commonTermCenter;
    tempDeviations[ltemp] -= commonTermDeviation;
  }
  
  length = ltemp + 1;
  size = length;

  delete [] indexes;
  delete [] deviations;

  indexes = tempIndexes;
  deviations = tempDeviations;

#ifdef FAST_RAD
  radius = 0.0;
  for (unsigned i = 0; i < length; i++)
    radius += fabs(deviations[i]);
#endif
  
  return *this;
}


/************************************************************
 * Operator:      /
 * Author & Date: ??? - ???
 * Description:   
 *   computes the division, we use the identity x/y = x * (1/y)
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF AAF::operator / (const AAF & P) const
{
  if (this == &P)
    return AAF(1.0);
  else
    return (*this)*inv(P);
}


/************************************************************
 * Operator:      /
 * Author & Date: Darius Grabowski - 10.06.2005
 * Description:   
 *   computes the division <double> / <AAF>
 *
 *   Input  : double    : divisor
 *            AAF       : dividend
 *   Output : AAF       : AAF result
 ************************************************************/
AAF operator / (double cst, const AAF P)
{  
  return inv(P)*cst;
}

/************************************************************
 * Operator:      ^
 * Author & Date: Darius Grabowski - 10.06.2005
 * Description:   
 *   computes the power function
 *
 *   Input  : const int : exponent
 *   Output : AAF       : AAF result
 ************************************************************/
AAF AAF::operator ^ (const int exponent) const
{
  double a, b;
  double fa, fb;
  double r;
  double alpha = 0.0;
  double dzeta = 0.0;
  double delta = 0.0;

  if (length == 0)
  {
    AAF Temp(pow((double)cvalue, double(exponent)));
    return Temp;
  }

  if (exponent == 0)
  {
    return AAF(1.0);
  }
  
  if (exponent == 1)
  {
    return *this;
  }

  if (exponent == -1)
  {
    return inv(*this);
  }
  r = rad();

  a = cvalue - r;
  b = cvalue + r;

  fa = pow(a, (double)exponent);
  fb = pow(b, (double)exponent);

  if ( (b*a < EPSILON) && (exponent < 0) )
  {
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
			std::string("Error <AAF::^>: computing zero inversion").c_str(), 
			a, b));
  }

  if (AAF::approximationType == CHEBYSHEV)
  {
    // use Chebyshev approximation
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      // avoid division by zero
      alpha = exponent*fa/(a+EPSILON);
    }
    
    // we have two points having the slope alpha
    double x_1 = -pow(fabs((double)alpha/exponent), (double)1.0/(exponent-1.0));

    double x_2 = -x_1;

    double fx_1, fx_2;
    
    // check valid points
    if (x_1 > a)
    {
      // x_1 is valid
      fx_1 = pow(x_1, (double)exponent);
    }
    else
    {
      x_1 = a;
      fx_1 = fa;
    }
    
    if (x_2 < b)
    {
      // x_2 is valid
      fx_2 = pow(x_2, (double)exponent);
    }
    else
    {
      x_2 = b;
      fx_2 = fb;
    }
    
    double y_2 = fx_2 - alpha*x_2;
    double y_1 = fx_1 - alpha*x_1;
    
    delta = 0.5*(y_1 - y_2);
    
    dzeta = 0.5*(y_1 + y_2);
#ifdef AAF_DEBUG
    fprintf(stdout, "Function: x^%d\n", exponent);
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_1: %f, f(x_1): %f\n", x_1, fx_1);
    fprintf(stdout, "x_2: %f, f(x_2): %f\n", x_2, fx_2);
    fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    // use MinRange approximation

    // special case: 0.0 in [a,b] : alpha = f'(0.0)
    // exp > 0 and exp even and [a,b] > 0 : MINRANGE: alpha = f'(a), y0_b > y0_a : e1 = 1 e2 = 1 e3 = 1
    // exp > 0 and exp odd  and [a,b] > 0 : MINRANGE: alpha = f'(a), y0_b > y0_a : e1 = 1 e2 = 0 e3 = 1
    // exp > 0 and exp even and [a,b] < 0 : MINRANGE: alpha = f'(b), y0_a > y0_b : e1 = 1 e2 = 1 e3 = 0
    // exp > 0 and exp odd  and [a,b] < 0 : MINRANGE: alpha = f'(b), y0_b > y0_a : e1 = 1 e2 = 0 e3 = 0
    // exp < 0 and exp even and [a,b] > 0 : MINRANGE: alpha = f'(b), y0_a > y0_b : e1 = 0 e2 = 1 e3 = 1
    // exp < 0 and exp odd  and [a,b] > 0 : MINRANGE: alpha = f'(b), y0_a > y0_b : e1 = 0 e2 = 0 e3 = 1
    // exp < 0 and exp even and [a,b] < 0 : MINRANGE: alpha = f'(a), y0_b > y0_a : e1 = 0 e2 = 1 e3 = 0
    // exp < 0 and exp odd  and [a,b] < 0 : MINRANGE: alpha = f'(a), y0_a > y0_b : e1 = 0 e2 = 0 e3 = 0
    
    if (b*a < 0.0)
    {
      alpha = 0.0;
      if (exponent % 2 == 0)
      {
	delta = 0.5*((fa>fb)?fa:fb);
	dzeta = delta;
      }
      else
      {
	delta = 0.5*(fb - fa);
	dzeta = 0.5*(fb + fa);
      }
    }
    else
    {
      if ( ( (exponent > 0) && (a > 0.0) ) || ( (exponent < 0) && (a < 0.0) ) )
      {
	alpha = exponent*fa/(a+EPSILON);
      }
      else
      {
	alpha = exponent*fb/(b+EPSILON);
      }
      
      double y_a = fa-alpha*a;
      double y_b = fb-alpha*b;

      delta = 0.5*fabs(y_a - y_b);
      
      dzeta = 0.5*(y_a + y_b);

#ifdef AAF_DEBUG
      fprintf(stdout, "Function: x^%d\n", exponent);
      fprintf(stdout, "Interval: [%f, %f]\n", a, b);
      fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
      fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
    }
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      // avoid division by zero
      alpha = exponent*fa/(a+EPSILON);
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(cvalue) + dzeta);

  Temp.length = length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi

  for (unsigned i = 0; i < length; i++)
  {
    Temp.indexes[i] = indexes[i];
    Temp.deviations[i] = alpha*(deviations[i]);
  }
  
  
  // Compute the error in a new deviation symbol
  // zk = delta
  Temp.indexes[length] = Temp.inclast();
  Temp.deviations[length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * radius + fabs(delta);
#endif

#ifdef AAF_DEBUG
			fprintf(stdout, "Ergebnis:\n");
			fprintf(stdout, "Mittelwert: %f!\n", Temp.cvalue);

			for(unsigned i=0; i<Temp.length;i++){
				fprintf(stdout, "Abweichung %d: Index: %u Wert: %f!\n", i, Temp.indexes[i], Temp.deviations[i]);
			}
#endif


  return Temp;
}

/************************************************************
 * Operator:      ^
 * Author & Date: Anna Krause - 01.03.2010
 * Description: Operator for the affine exponentiation
 *   computes the power function
 *
 *   Input  : AAF		: exp
 *   Output : AAF       : AAF result
 ************************************************************/
AAF AAF::operator ^ (const AAF exp) const
{
	AAF base=*this;
	AAF temp=aaf_pow(*this,exp);
	return temp;
}



/************************************************************
 * Method:      aaf_pow
 * Author & Date: Anna Krause - 01.03.2010
 * Description:   
 *   computes the power function
 *
 *   Input  : AAF	 	: base
 *   		  AAF		: exponent
 *   Output : AAF       : AAF result
 ************************************************************/
AAF aaf_pow(const AAF & base, const AAF & exp) {
// Properties of the input aafs
	double centre_b=base.getcenter();
	double centre_e=exp.getcenter();
	unsigned length_b=base.getlength();
	unsigned length_e=exp.getlength();
	double min_b=base.getMin();
	
	const double pi=3.14159265;
	
// Variables
	double fx, fy, fc;
	double dmin, dmax, x1, x2, x3, y1, y2, y3, phi, phi0, phi_max, d1, d2, d3, dyc;
	unsigned * merge;
	unsigned * help;
	double * dev_array_b;
	double * dev_array_e;
	double * eps_vec;
	unsigned dev_b, dev_e;
	unsigned i=0, j=0, last_eps;
	double alpha, gamma;
	
// Properties of the resulting aaf
	unsigned length_r;	
	double * devs_r;
	unsigned * ind_r;
		
// Area of validity
	if( min_b<=0 ){
		throw(AAF_Exception(AAF_NEGBASE_EXCEPTION,
			std::string("Error <AAF::aaf_pow>: Exponentiation of negative numbers is only possible for integer exponents and currently not defined in this scope.").c_str(),
			min_b));
	}
	
// Special cases
	if(length_b==0){
		if(length_e==0){
			// b and e are real numbers
			AAF result(pow(centre_b,centre_e));
			std::cout << result << std::endl;
			return result;
		}
		else{// b is a real number, e is an affine form
			AAF result=aaf_pow(centre_b, exp);
			return result;
		}
	}
	else {
		if(length_e==0){ //b is an affine form, e is a real number
			AAF result=aaf_pow(base,centre_e);
			return result;
		}
	}
	
	unsigned * ind_e=exp.indexes;
	double * devs_e=exp.deviations;
	unsigned * ind_b=base.indexes;
	double * devs_b=base.deviations;

// Two dimensional Taylor series at the centre value
	fc=pow(centre_b, centre_e);
	fx=centre_e*pow(centre_b,(centre_e-1.0));
	fy=pow(centre_b,centre_e)*log(centre_b);
	
// Evaluate the edge of the polygon defined by the base and the exponent
	dmin=0.0;
	dmax=0.0;
	x1=centre_b;
	y1=centre_e;
	x2=0.0;
	y2=0.0;
	x3=0.0;
	y3=0.0;
	last_eps=NULL;
	
// Merge vector of indexes and construct array for the deviations and the epsilon
	//TODO: Pr�fen: mehrere unabh�ngige Abweichungssymbole, die nur in einer AAF vorkommen? - sollte egal sein
	help=new unsigned[length_b+length_e+1];
	merge=std::set_union(ind_b,ind_b+length_b,ind_e,ind_e+length_e,help);
	length_r=merge-help+1;
	
	dev_b=0;
	dev_e=0;
	dev_array_b=new double[length_r-1];
	dev_array_e=new double[length_r-1];
	eps_vec=new double[length_r-1];
	for(i=0; i<length_r-1; i++){
		fprintf(stdout, "i: %d, help: %d, desb: %f, devse. %f\n", i, help[i], devs_b[dev_b], devs_e[dev_e] );
		if(dev_b<length_b && ind_b[dev_b]==help[i]){
			dev_array_b[i]=devs_b[dev_b];
			dev_b++;
		}
		else{
			dev_array_b[i]=0.0;
		}
		
		if(dev_e<length_e && ind_e[dev_e]==help[i]){
			dev_array_e[i]=devs_e[dev_e];
			dev_e++;
		}
		else{
			dev_array_e[i]=0.0;
		}
		
		// Write epsilon vector
		if(dev_array_b[i]>0)
			eps_vec[i]=-1.0;
		else if(dev_array_b[i]<0)
			eps_vec[i]=1.0;
		else{
			if(dev_array_e[i]>=0)
				eps_vec[i]=-1.0;
			else
				eps_vec[i]=1.0;
		}
	}
	
// Find minimum and maximum distance between taylor series and exponentiation
	for(i=0; i<(length_r-1); i++){
		x1+=eps_vec[i]*dev_array_b[i];
		y1+=eps_vec[i]*dev_array_e[i];
	}
	
	for(i=0; i<2*(length_r-1); i++){
		// Find the outmost segment
		phi0=atan2((centre_e-y1),(centre_b-x1));
		phi_max=0.0;
		
		for(j=0;j<(length_r-1); j++){ 
			if(dev_array_b[j]==0 && dev_array_e[j]!=0){ //x=const
				x2=x1;
				x3=x1;
				y2=y1+(-1)*eps_vec[j]*dev_array_e[j];
				y3=y1+(-2)*eps_vec[j]*dev_array_e[j];
				
				phi=atan2((y3-y1),(x3-x1))-phi0;
				if(phi>phi_max && phi>0 && phi<pi){
					phi_max=phi;
					last_eps=j;
				}
				
			}
			else if(dev_array_b[j]!=0 && dev_array_e[j]==0){ //y=const
				y2=y1;
				y3=y1;
				x2=x1+(-1)*eps_vec[j]*dev_array_b[j];
				x3=x1+(-2)*eps_vec[j]*dev_array_b[j];
				
				phi=atan2((y3-y1),(x3-x1))-phi0;
				if(phi>phi_max && phi>0 && phi<pi){
					phi_max=phi;
					last_eps=j;
				}
			}
			else{	
				x2=x1+(-1)*eps_vec[j]*dev_array_b[j];
				x3=x1+(-2)*eps_vec[j]*dev_array_b[j];
				y2=y1+(-1)*eps_vec[j]*dev_array_e[j];
				y3=y1+(-2)*eps_vec[j]*dev_array_e[j];
				
				phi=atan2((y3-y1),(x3-x1))-phi0;
				if(phi>phi_max && phi>0 && phi<pi){
					phi_max=phi;
					last_eps=j;
				}
			}
		}
		
		x2=x1+(-1)*eps_vec[last_eps]*dev_array_b[last_eps];
		x3=x1+(-2)*eps_vec[last_eps]*dev_array_b[last_eps];
		y2=y1+(-1)*eps_vec[last_eps]*dev_array_e[last_eps];
		y3=y1+(-2)*eps_vec[last_eps]*dev_array_e[last_eps];
		
		d1=fc+fx*(x1-centre_b)+fy*(y1-centre_e)-pow(x1,y1);
		d2=fc+fx*(x2-centre_b)+fy*(y2-centre_e)-pow(x2,y2);
		d3=fc+fx*(x3-centre_b)+fy*(y3-centre_e)-pow(x3,y3);
		
		if((d1<d2) && (d2<d3)){
			if(i==0) {
				dmin=d1;
				dmax=d3;
			}
			if(d1<dmin) dmin=d1;
			if(d3>dmax) dmax=d3;
		}
		else if((d1>d2) && (d2>d3)){
			if(i==0) {
				dmin=d3;
				dmax=d1;
			}
			if(d3<dmin) dmin=d3;
			if(d1>dmax) dmax=d1;
		}
		else{
			fprintf(stdout, "\t Warning: The distance function is not monotone on the straight line from (%f,%f) to (%f,%f) \n", x1, y1, x3, y3);

			if(i==0){ //initial values for dmin and dmay
				if(d1>d3){
					dmin=d3;
					dmax=d1;
				}
				else if(d1<d3){
					dmin=d1;
					dmax=d3;
				}
				else {
					dmin=d1;
					dmax=d1;
				}
			}

			if(dev_array_b[last_eps]==0){
				if(fy/log(x1)>0){
					dyc=log(fy/log(x1))/log(x1);
					d2=fc+fx*(x1-centre_b)+fy*(dyc-centre_e)-pow(x1,dyc);
				}
			}
			else if(dev_array_e[last_eps]==0){
				if(fx/y1>0){
					dyc=pow(fx/y1,1/(y1-1));
					d2=fc+fx*(dyc-centre_b)+fy*(y1-centre_e)-pow(dyc,y1);
				}
			}
			else{
				//Todo: eine beliebige Gerade
				if((d2<d1 && d2<dmin)||(d2>d1 && d2>dmax)) fprintf(stdout, "Warning: Information may have been lost!\n");
			}
			// write minimum and maximum distance if applicable
			if(d2<=d1){
				if(d2<dmin) dmin=d2;
				if(d1<d3 && d3>dmax) dmax=d3;
				else if(d1>d3 && d1>dmax) dmax=d1;
			}
			else{
				if (d2>dmax) dmax=d2;
				if(d1<d3 && d1<dmin) dmin=d1;
				else if(d1>d3 && d3<dmin) dmin=d3;
			}
		}
		
		x1=x3;
		y1=y3;
		eps_vec[last_eps]=(-1)*eps_vec[last_eps];
		
		for(j=0;j<length_r-1;j++){
		}
	}
	
// Calculate the shift of the centre and the new deviation
	alpha=(dmax+dmin)/2;
	gamma=(dmax-dmin)/2;
	
	fprintf(stdout, "alpha: %f, gamma: %f\n", alpha, gamma);
	
	AAF result(pow(centre_b,centre_e)+alpha);
	
	result.length=length_r;
	result.size=length_r;
	result.indexes=new unsigned[length_r];
	result.deviations=new double[length_r];
	
	devs_r=result.deviations;
	ind_r=result.indexes;
	
	for(i=0; i<(length_r-1); i++){
		ind_r[i]=help[i];
		devs_r[i]=fx*dev_array_b[i]+fy*dev_array_e[i];
		fprintf(stdout, "index: %d, devb: %f, deve: %f\n", help[i], dev_array_b[i], dev_array_e[i]);
	}
	
	result.indexes[length_r-1]=result.inclast();
	result.deviations[length_r-1]=gamma;
	
#ifdef FAST_RAD
	result.radius = 0.0;
	for (unsigned i = 0; i < length_r; i++)
		result.radius += fabs(result.deviations[i]);
#endif
	
	//std::cout << result << std::endl;
	return result;
}

/************************************************************
 * Method:      aaf_pow
 * Author & Date: Anna Krause - 01.03.2010
 * Description:
 *   computes the power function
 *
 *   Input  : AAF	 	: base
 *   		  double	: exponent
 *   Output : AAF       : AAF result
 ************************************************************/
AAF aaf_pow(const AAF & base, const double & exp) {
// AAF properties
	double centre_b=base.getcenter();
	unsigned length_b=base.getlength();
	double min_b=base.getMin();
	double max_b=base.getMax();
	unsigned * ind_b=base.indexes;
	double * devs_b=base.deviations;

// Variables
	double fmin, fmax, alpha, beta, gamma, x2;

	fmin=pow(min_b,exp);
	fmax=pow(max_b,exp);
	beta=0.0;
	alpha=0.0;
	gamma=0.0;

	if(exp>1.0){
		//fprintf(stdout, "Exponent > 1!\n");
		if (AAF::approximationType == MINRANGE){
			beta=exp*pow(min_b,exp-1);
			alpha=0.5*(-beta*(min_b+max_b)+fmin+fmax);
			gamma=0.5*(beta*(min_b-max_b)-fmin+fmax);
		}
		else { //CHEBYSHEV
			beta=(fmax-fmin)/(max_b-min_b);
			x2=pow((beta/exp),1/(exp-1));
			alpha=0.5*(-beta*(min_b+x2)+fmin+pow(x2,exp));
			gamma=0.5*(beta*(x2-min_b)+fmin-pow(x2,exp));
		}
	}
	else if(exp==1.0){
		//fprintf(stdout, "Exponent = 1!\n");
		beta=1.0;
		alpha=0.0;
		gamma=0.0;
	}
	else if(exp<1.0 && exp>0.0){
		//fprintf(stdout, "Exponent < 1 und > 0!\n");
		if (AAF::approximationType == MINRANGE){
			beta=exp*pow(max_b,exp-1);
			alpha=0.5*(-beta*(min_b+max_b)+fmin+fmax);
			gamma=0.5*(beta*(min_b-max_b)-fmin+fmax);
		}
		else{
			beta=(fmax-fmin)/(max_b-min_b);
			x2=pow((beta/exp),1/(exp-1));
			alpha=0.5*(-beta*(min_b+x2)+fmin+pow(x2,exp));
			gamma=0.5*(beta*(-x2+min_b)-fmin+pow(x2,exp));
		}
	}
	else if(exp==0.0){
		//fprintf(stdout, "Exponent = 0!\n");
		beta=0.0;
		alpha=1.0;
		gamma=0.0;
	}
	else if(exp<0.0){
		//fprintf(stdout, "Exponent < 0!\n");
		if (AAF::approximationType == MINRANGE){
			beta=exp*pow(min_b,exp-1);
			alpha=0.5*(-beta*(min_b+max_b)+fmin+fmax);
			gamma=0.5*(beta*(min_b-max_b)-fmin+fmax);
		}
		else { //CHEBYSHEV
			beta=(fmax-fmin)/(max_b-min_b);
			x2=pow((beta/exp),1/(exp-1));
			alpha=0.5*(-beta*(min_b+x2)+fmin+pow(x2,exp));
			gamma=0.5*(beta*(x2-min_b)+fmin-pow(x2,exp));
		}
	}

//build result
	AAF result(beta*centre_b+alpha);
	if (gamma == 0.0){
		result.length = length_b;
		result.size = result.length;
		result.deviations = new double [result.size];
		result.indexes = new unsigned [result.size];

		for (unsigned i = 0; i < length_b; i++){
			result.indexes[i] = ind_b[i];
			result.deviations[i] = beta*(devs_b[i]);
		}
	}
	else{
		result.length = length_b + 1;
		result.size = result.length;
		result.deviations = new double [result.size];
		result.indexes = new unsigned [result.size];

		for (unsigned i = 0; i < length_b+1; i++){
			result.indexes[i] = ind_b[i];
			result.deviations[i] = beta*(devs_b[i]);
		}

		result.indexes[length_b] = result.inclast();
		result.deviations[length_b] = gamma;
	}

	#ifdef FAST_RAD
		result.radius = fabs(beta) * 0.5* (max_b-min_b) + fabs(gamma);
	#endif

	std::cout << result << std::endl;
	return result;
}

/************************************************************
 * Method:      aaf_pow
 * Author & Date: Anna Krause - 01.03.2010
 * Description:
 *   computes the power function
 *
 *   Input  : double 	: base
 *   		  AAF		: exponent
 *   Output : AAF       : AAF result
 ************************************************************/
AAF aaf_pow(const double & base, const AAF & exp) {
// AAF properties
	double centre_e=exp.getcenter();
	unsigned length_e=exp.getlength();
	double min_e=exp.getMin();
	double max_e=exp.getMax();
	unsigned * ind_e=exp.indexes;
	double * devs_e=exp.deviations;

// Variables
	double fmin, fmax, alpha, beta, gamma, x2;

	//fprintf(stdout, "Basis ist reell!\n");

	fmin=pow(base,min_e);
	fmax=pow(base,max_e);
	beta=0.0;
	gamma=0.0;
	alpha=0.0;

	if(base==1.0){
		fprintf(stdout, "\t und gleich 1!\n");
		beta=0.0;
		alpha=1.0;
		gamma=0.0;
	}
	else {
		if (AAF::approximationType == MINRANGE){
			beta=fmin*log(base);
			alpha=-beta*(min_e+max_e)+fmin+fmax;
			gamma=beta*(min_e-max_e)-fmin+fmax;
		}
		else{ //CHEBYCHEV
			beta=(fmax-fmin)/(max_e-min_e);
			x2=log(beta/log(base))/log(base);
			alpha=0.5*(-beta*(min_e+x2)+fmin+pow(base,x2));
			gamma=0.5*(beta*(x2-min_e)+fmin-pow(base,x2));
		}
	}

//build result
	AAF result(beta*centre_e+alpha);
	result.length = length_e + 1;
	result.size = result.length;
	result.deviations = new double [result.size];
	result.indexes = new unsigned [result.size];

	for (unsigned i = 0; i < length_e+1; i++){
		result.indexes[i] = ind_e[i];
		result.deviations[i] = beta*(devs_e[i]);
	}

	result.indexes[length_e] = result.inclast();
	result.deviations[length_e] = gamma;

	#ifdef FAST_RAD
		result.radius = fabs(beta) * 0.5* (max_e-min_e) + fabs(gamma);
	#endif

	std::cout << result << std::endl;
	return result;
}

/************************************************************
 * Method:        sqrt
 * Author & Date: ??? - ???
 * Description:   
 *   This function computes the square root
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF sqrt(const AAF & P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;

  if (P.length == 0)
  {
    AAF Temp(sqrt(P.cvalue));
    return Temp;
  }

  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  // check lower bound
  if (a < EPSILON)
  {
    throw(AAF_Exception(AAF_NEGROOT_EXCEPTION, 
			std::string("Error <AAF::sqrt>: computing negative root").c_str(), 
			a, b));
  }

  double dNeg = 0.0;

  if (a < 0.)
  {
    dNeg = -a/2.;
    a = 0.;
    fa = 0.;
  }
  else
  {
    fa = sqrt(a);
  }
  fb = sqrt(b);


  if (AAF::approximationType == CHEBYSHEV)
  {
    // temporary vars
    double t, u;

    alpha = 1.0/(fa + fb);
    
    t = 0.25*fa*fb*alpha;
    
    u = 0.125*(a + b)*alpha;
    
    // dzeta calculation:
    dzeta = u + 3.0*t;
    
    // Calculation of the error
    delta = u - t;
  }
  else if (AAF::approximationType == MINRANGE)
  {       
    alpha = 1.0/(2.0*fb);

    delta = 0.5*alpha*(a-2*fb*fa+b);
    
    dzeta = 0.5*fb-delta;
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1.0/(2.0*fb);
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue + dNeg) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  alpha -= dNeg/P.rad();
  // zi = alpha*xi
  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }
  
  // Compute the error --> delta in a new deviation symbol 
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}


/************************************************************
 * Method:        isqrt
 * Author & Date: Darius Grabowski - 03.02.2006
 * Description:   
 *   This function computes the inverted square root
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF isqrt(const AAF & P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;

  if (P.length == 0)
  {
    AAF Temp(1.0/sqrt(P.cvalue));
    return Temp;
  }

  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  // check bounds
  if ( (a < EPSILON) || (b < EPSILON) )
  {
    throw(AAF_Exception(AAF_NEGROOT_EXCEPTION, 
			std::string("Error <AAF::isqrt>: computing negative root").c_str(), 
			a, b));
  }

  fa = 1.0/sqrt(a);
  fb = 1.0/sqrt(b);

  if (AAF::approximationType == CHEBYSHEV)
  {
    double x_;  // temporary var

    // compute the slope
    if (r > AAF_MINRAD)
    {
      alpha = (fb - fa)/(b - a);
    }
    else
    {
      alpha = -0.5*fb*fb*fb;
    }
    
    // compute the point having the slope alpha
    x_ = pow((double)0.5*(a/fb + b/fa), (double)2/3);
    
    delta = 0.5*(fa + alpha*(x_ - a) - 1.0/sqrt(x_));
   
    // dzeta calculation:
    dzeta = fa - alpha*a - delta;
  }
  else if (AAF::approximationType == MINRANGE)
  {
    alpha = -0.5*fb*fb*fb;

    delta = 0.5*(fa - fb + alpha*(b - a));
    
    dzeta = fa - alpha*a - delta;
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = -0.5*fb*fb*fb;
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi

  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }
  
  // Compute the error --> delta in a new deviation symbol 

  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}


/************************************************************
 * Method:        inv
 * Author & Date: ??? - ???
 * Description:   
 *   This function computes the inverse 
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF inv(const AAF & P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;

  if (P.length == 0)
  {
    AAF Temp(1.0/(P.cvalue));
    return Temp;
  }

  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  if (a*b < EPSILON)
  {
    throw(AAF_Exception(AAF_DIVZERO_EXCEPTION, 
    			std::string("Error <AAF::inv>: computing zero inversion").c_str(), 
			a, b));
  }

  fa = 1/a;
  fb = 1/b;

  if (AAF::approximationType == CHEBYSHEV)
  {
    alpha = -fa*fb;    
    double u = sqrt(a*b);

    if (a > 0)
    {
      delta = +0.5*(fa+fb-2.0/u);
      dzeta = fa+fb-delta;
    }
    else
    {
      delta = -0.5*(fa+fb+2.0/u);
      dzeta = fa+fb+delta;
    }
    
#ifdef AAF_DEBUG
    fprintf(stdout, "Function: inv(x)\n");
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_1: %f, f(x_1): %f\n", u, 1.0/u);
    fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    double y_a, y_b;
    // Derivative of 1/x is -1/x*x
    if (a > 0.0)
    {
      alpha = -fb/b;
      // y_a = fa - alpha*a;
      // y_b = fb - alpha*b = 2.0*fb;
      y_a = fa - alpha*a;
      y_b = 2.0*fb;     
    }
    else
    {
      alpha = -fa/a;
      // y_a = fa - alpha*a = 2.0*fa;
      // y_b = fb - alpha*b;
      y_a = 2.0*fa;
      y_b = fb - alpha*b;
    }

    delta = 0.5*(y_a - y_b);
    dzeta = 0.5*(y_a + y_b);
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = -fa*fb;
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }
  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi

  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }  
  
  // Compute the error in a new deviation symbol
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}


/************************************************************
 * Method:        exp
 * Author & Date: ??? - ???
 * Description:   
 *   This function computes the exponential fct.
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF exp(const AAF &P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;

  bool chebyshevFailed = false;

  if (P.length == 0)
  {
    AAF Temp(exp(P.cvalue));
    return Temp;
  }

  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  // if upper bound < -50 --> return value is smaller than 1E-20
  if (b < -50.0)
  {
    AAF Temp(exp(P.cvalue));
    return Temp;
  }

  // Derivative of exp(x) is exp(x)
  fa = exp(a); 
  fb = exp(b);  
  
  if (AAF::approximationType == CHEBYSHEV)
  {
    // use Chebyshev approximation
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = fa;
    }
    
    double u = log(alpha);
    
    delta = 0.5*(fa+(u-a-1.0)*alpha);
    
    dzeta = fa-a*alpha-delta;
    
    if ( (fa - 2.0*delta) < 0.0)
      chebyshevFailed = true;

#ifdef AAF_DEBUG
    fprintf(stdout, "Function: exp(x)\n");
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_1: %f, f(x_1): %f\n", u, exp(u));
    fprintf(stdout, "alpha: %g, delta: %g, dzeta: %g\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    // use MinRange approximation
    alpha = fa;
  
    dzeta = 0.5*(fa+fb-alpha*(a+b));

    delta = 0.5*(fb-fa-alpha*(b-a));
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r/fa > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = fa;
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi

  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }
    
  // Compute the error in a new deviation symbol 
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}


/************************************************************
 * Method:        log
 * Author & Date: Darius Grabowski - 10.10.2005
 * Description:   
 *   computes the natural logarithm function
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF log(const AAF &P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;

  if (P.length == 0)
  {
    AAF Temp(log(P.cvalue));
    return Temp;
  }
  
  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  if ( (a < EPSILON) || (b < EPSILON) )
  {
    throw(AAF_Exception(AAF_NEGLOG_EXCEPTION, 
			std::string("Error <log(AAF)>: computing negative log").c_str(), 
			a, b));
  }

  // Derivative of log(x) is 1.0/x  
  
  fa = log(a);
  fb = log(b);

  if (AAF::approximationType == CHEBYSHEV)
  {
    // use Chebyshev approximation
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1.0/a;
    }
   
    double u = 1.0/alpha;
    double fu = log(u);

    double y_1 = fu - alpha*u;
    double y_2 = fb - alpha*b;

    delta = 0.5*(y_1 - y_2);
    
    // dzeta = 0.5*(y_1 + y_2); 
    dzeta = y_1 - delta;
    
    
#ifdef AAF_DEBUG
    fprintf(stdout, "Function: log(x)\n");
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_: %f, f(x_): %f\n", u, fu);
    fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    // use MinRange approximation    
    alpha = 1.0/b;
  
    delta = 0.5*(fb-fa-alpha*(b-a));
    
    dzeta = fb-1.0-delta;    
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1./a;
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi
  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }
  
  // Compute the error in a new deviation symbol  
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;
  
#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}


/************************************************************
 * Method:        logexp
 * Author & Date: Darius Grabowski - 10.10.2005
 * Description:   
 *   computes the chained function log(1.0 + exp(x))
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF logexp(const AAF &P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;

  if (P.length == 0)
  {
    AAF Temp(log(1.0 + exp(P.cvalue)));
    return Temp;
  }
  
  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  if (b < 0)
  {
    fa = log(1.0 + exp(a));
    fb = log(1.0 + exp(b));
  }
  else
  {
    fa = a + log(1.0+exp(-a));
    fb = b + log(1.0+exp(-b));
  }

  // Derivative of log(1+exp(x)) is exp(x)/(1+exp(x))

  if (AAF::approximationType == CHEBYSHEV)
  {
    // use Chebyshev approximation
    if (r/fa > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = exp(a)/(1.0 + exp(a));
    }
   
    if (alpha < 1E-15)
      alpha = 1E-15;
    
    double u, fu;

    if (alpha > 1.0-1E-15)
    {
      u = 0.5*(a+b);
      fu = log(1.0 + exp(u));
    } 
    else
    {
      u = log(alpha/(1.0 - alpha));
      fu = -log(1.0 - alpha);
    }

    delta = 0.5*(fa - fu - alpha*(a - u));
    
    dzeta = fa - alpha*a - delta;

#ifdef AAF_DEBUG
    fprintf(stdout, "Function: log(1.0+exp(x))\n");
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_: %f, f(x_): %f\n", u, fu);
    fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    // use MinRange approximation
    alpha = exp(a)/(1.0 + exp(a));

    double y_a = fa - alpha*a;
    double y_b = fb - alpha*b;
  
    delta = 0.5*(y_b - y_a);
    
    dzeta = 0.5*(y_b + y_a);    
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r/fa > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = exp(a)/(1.0 + exp(a));
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi
  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }
  
  // Compute the error in a new deviation symbol  
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif
  
  return Temp;
}

/************************************************************
 * Method:        atan
 * Author & Date: Darius Grabowski - 10.06.2005
 * Description:   
 *   Computes the arcus tangens function
 *   The atan function is strictly monotonic increasing but it 
 *   has two possible turning points which has to be considered
 *   for the Chebyshev approxiamtion
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF atan(const AAF &P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;  

  if (P.length == 0)
  {
    AAF Temp(atan(P.cvalue));
    return Temp;
  }

  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  // Derivative of atan(x) is 1.0/(1.0+x^2)
  fa = atan(a);
  fb = atan(b);

  if (AAF::approximationType == CHEBYSHEV)
  {
    // use Chebyshev approximation
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1.0/(1.0+a*b);
    }
    
    // we have two points having the slope alpha
    double x_1 = -sqrt((1-alpha)/alpha);
    double x_2 = -x_1;

    double fx_1, fx_2;

    // check valid points
    if (x_1 > a)
    {
      // x_1 is valid
      fx_1 = atan(x_1);
    }
    else
    {
      x_1 = a;
      fx_1 = fa;
    }

    if (x_2 < b)
    {
      // x_2 is valid
      fx_2 = atan(x_2);
    }
    else
    {
      x_2 = b;
      fx_2 = fb;
    }

    double y_2 = fx_2 - alpha*x_2;
    double y_1 = fx_1 - alpha*x_1;
    
    delta = 0.5*(y_2 - y_1);
    
    dzeta = 0.5*(y_2 + y_1);

#ifdef AAF_DEBUG
    fprintf(stdout, "Function: atan(x)\n");
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_1: %f, f(x_1): %f\n", x_1, fx_1);
    fprintf(stdout, "x_2: %f, f(x_2): %f\n", x_2, fx_2);
    fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    // use MinRange approximation    
    if (fabs(b) > fabs(a))
    {
      alpha = 1.0/(1.0+b*b);
    }
    else
    {
      alpha = 1.0/(1.0+a*a);
    }
  
    double y0_a = fa-alpha*a;
    double y0_b = fb-alpha*b;
    
    dzeta = 0.5*(y0_b + y0_a);

    delta = 0.5*(y0_b - y0_a);
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1.0/(1.0+a*a);
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta
  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi
  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }  
  
  // Compute the error in a new deviation symbol
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}

/************************************************************
 * Method:        tanh
 * Author & Date: Darius Grabowski - 10.10.2005
 * Description:   
 *   computes the tangens hyperbolicus function
 *   Definition: tanh(x) := (exp(x)-exp(-x))/(exp(x)+exp(-x))
 *   The tanh function is strictly monotonic increasing but it 
 *   has two possible turning points which has to be considered
 *   for the Chebyshev approximation
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF tanh(const AAF &P)
{
  double a, b;
  double fa, fb;
  double r;
  double alpha, dzeta, delta;  

  if (P.length == 0)
  {
    AAF Temp(tanh(P.cvalue));
    return Temp;
  }

  r = P.rad();

  a = P.cvalue - r;
  b = P.cvalue + r;

  // Derivative of tanh(x) is 1.0-tanh(x)^2

  fa = tanh(a);
  fb = tanh(b);

  if (AAF::approximationType == CHEBYSHEV)
  {
    // use Chebyshev approximation
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1.0-fa*fb;
    }
    
    // we have two points having the slope alpha
    double dSqrtAlpha = sqrt(1-alpha);
    double x_1 = -atanh(dSqrtAlpha);
    double x_2 = -x_1;

    double fx_1, fx_2;

    // check valid points
    if (x_1 > a)
    {
      // x_1 is valid
      fx_1 = -dSqrtAlpha;
    }
    else
    {
      x_1 = a;
      fx_1 = fa;
    }

    if (x_2 < b)
    {
      // x_2 is valid
      fx_2 = dSqrtAlpha;
    }
    else
    {
      x_2 = b;
      fx_2 = fb;
    }

    double y_2 = fx_2 - alpha*x_2;
    double y_1 = fx_1 - alpha*x_1;
    
    delta = 0.5*(y_2 - y_1);
    
    dzeta = 0.5*(y_2 + y_1);

#ifdef AAF_DEBUG
    fprintf(stdout, "Function: tanh(x)\n");
    fprintf(stdout, "Interval: [%f, %f]\n", a, b);
    fprintf(stdout, "f(Interval): [%f, %f]\n", fa, fb);
    fprintf(stdout, "x_1: %f, f(x_1): %f\n", x_1, fx_1);
    fprintf(stdout, "x_2: %f, f(x_2): %f\n", x_2, fx_2);
    fprintf(stdout, "alpha: %f, delta: %f, dzeta: %f\n", alpha, delta, dzeta);
#endif
  }
  else if (AAF::approximationType == MINRANGE)
  {
    // use MinRange approximation    
    if (fabs(fb) > fabs(fa))
    {
      alpha = 1.0-fb*fb;
    }
    else
    {
      alpha = 1.0-fa*fa;
    }
  
    double y0_a = fa-alpha*a;
    double y0_b = fb-alpha*b;
    
    dzeta = 0.5*(y0_b + y0_a);

    delta = 0.5*(y0_b - y0_a);
  }
  else // (AAF::approximationType == SECANT)
  {
    if (r > AAF_MINRAD)
    {
      alpha = (fb-fa)/(b-a);
    }
    else
    {
      alpha = 1.0-fa*fb;
    }
    dzeta = fa - alpha*a;
    delta = 0.;
  }

  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue) + dzeta);

  Temp.length = P.length + 1;
  Temp.size = Temp.length;
  Temp.deviations = new double [Temp.size];
  Temp.indexes = new unsigned [Temp.size];

  // zi = alpha*xi

  for (unsigned i = 0; i < P.length; i++)
  {
    Temp.indexes[i] = P.indexes[i];
    Temp.deviations[i] = alpha*(P.deviations[i]);
  }
  
  
  // Compute the error in a new deviation symbol  
  // zk = delta
  Temp.indexes[P.length] = Temp.inclast();
  Temp.deviations[P.length] = delta;

#ifdef FAST_RAD
  Temp.radius = fabs(alpha) * P.radius + fabs(delta);
#endif

  return Temp;
}


/************************************************************
 * Method:        heaviside
 * Author & Date: Darius Grabowski - 10.12.2005
 * Description:   
 *   computes a simplified heaviside function
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
AAF heaviside(const AAF &P)
{
  AAF Temp(0.0);

  if (P.getcenter() > 0.0)
    Temp += 1.0;

  return Temp;
}

