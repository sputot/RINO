/*
 * aa_aaf.h -- Affine Arithmetic class
 * Copyright (C) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
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


#ifndef AA_AAF_H
#define AA_AAF_H

#include "aa_interval.h"
#include "aa_exceptions.h"

//#include "vnode.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

typedef enum {MINRANGE, CHEBYSHEV, SECANT} tApproximationType;

//#define FAST_RAD

#ifdef USE_QPF
class QPF;
#endif

// Affine Arithmetic Form
class AAF
{

 private:

  // central value
  double cvalue;
  // length of indexes
  unsigned length;
  // array size of indexes and deviations
  unsigned size;

#ifdef FAST_RAD
  double radius;
#endif

  // At creation we don't store null deviations

  // values of parial deviations
  double * deviations;
  // indexes of partial deviations
  unsigned * indexes;

  // current approximation type: <CHEBYSHEV> (default), <MINRANGE> or <SECANT>
  static tApproximationType approximationType;
  // highest deviation symbol in use
  static unsigned last;

#ifdef CLEANUP
  // pointers to all affine variables
  static list <AAF*> allAAF;
#endif

 public:
  
  // constructors
  AAF(double v0 = 0.0);
  AAF(double, const double *, const unsigned *, unsigned);
  AAF(const AAF &);
  AAF(const AAInterval);
  // ajout SP
  AAF(const interval);

  // destructor
  ~AAF();

  unsigned int get_class_size(void);
  unsigned int get_class_size_index(void);
  unsigned int get_class_size_deviations(void);

  double operator[](unsigned) const;

  bool operator < (const AAF &) const;
  bool operator <= (const AAF &) const;
  bool operator > (const AAF &) const;
  bool operator >= (const AAF &) const;
  bool operator == (const AAF &) const;

  AAF & operator = (const AAF &);
  AAF & operator = (const double);
  AAF operator + (const AAF &) const;
  AAF operator - (const AAF &) const;
  AAF operator * (const AAF &) const;
  AAF operator / (const AAF &) const;
  AAF operator ^ (const int) const;

  AAF operator - () const;
  AAF & operator += (double);
  AAF & operator -= (double);
  AAF & operator *= (double);
  AAF & operator /= (double);
  AAF & operator += (const AAF &);
  AAF & operator -= (const AAF &);
  AAF & operator *= (const AAF &);
  AAF operator * (double);
  AAF operator / (double) const;
  AAF operator ^ (const AAF) const;

  friend std::ostream & operator << (std::ostream &, const AAF &);
 
  unsigned getlength() const;
  double getcenter() const;

  AAInterval convert() const;
    
    
    // concretization
    interval convert_int() const;
    interval convert_int(vector<interval> &constr_eps, int sysdim) const;
    // concretization, given constraints on the noise symbols up to sysdim
    
    // suppose we have some constraints on the noise symbols eps_i, given by some affine forms being <= 0
    interval convert_int(vector<AAF> &constraints, int sysdim) const;
    
    void reduce_aaf(interval I);
    interval rad(vector<interval> &constr_eps, int sysdim) const;

    
  double rad() const;
  double getMax() const;
  double getMin() const;
  double getAbsMax() const;
  double getAbsMin() const;
  unsigned getFirstIndex(void) const;
  unsigned getLastIndex(void) const;
  
    void compact(void);
    void sumup(double);
unsigned sumup(unsigned);
    
#ifdef USE_AAF_EXTENSIONS
  // special methods declared in aa_aafspecial.c
#ifdef USE_QPF
  AAF & operator = (const QPF &);
#endif
  unsigned getData(unsigned*, double*) const;
  void get(double *, unsigned, unsigned npp = 1);
  void getPD(double *, unsigned, unsigned, unsigned npp = 1);
  void getAt(double *, unsigned, unsigned, unsigned npp = 1);
  void update(double *, unsigned, unsigned npp = 1);
  
  void sumupall(void);
  void resize(void) {resize(getLastIndex()-1);};
  void resize(unsigned int);
  void resize(unsigned int, unsigned int);
  void ResizeNewSymbol(unsigned int ind);
  void addPD(double);
  double at(unsigned) const;
  void set(unsigned, double);
  void setLast(double);
  void submul(const AAF &, double);
#endif

  // static methods
  static tApproximationType getApproximationType(void);
  static void setApproximationType(tApproximationType);

#ifdef USE_AAF_EXTENSIONS
  // this method must be treaten carefully !!!
  static void setDefault(const unsigned val = 0);
#endif
  static unsigned getDefault(void);
#ifdef CLEANUP
  static void cleanup(double);
#endif

  // friend methods
  friend AAF sqrt(const AAF &);
  friend AAF sqr(const AAF &);
  friend AAF isqrt(const AAF &);
  friend AAF exp(const AAF &);
  friend AAF log(const AAF &);
  friend AAF logexp(const AAF &);
  friend AAF atan(const AAF &);
  friend AAF tanh(const AAF &);
  friend AAF inv(const AAF &);
  friend AAF sin(const AAF &);
  friend AAF cos(const AAF &);
  friend AAF tan(const AAF &);
  friend AAF heaviside(const AAF &);
  friend AAF aaf_pow( const AAF&, const AAF&);
  friend AAF aaf_pow( const AAF&, const double&);
  friend AAF aaf_pow( const double&, const AAF&);
#ifdef USE_AAF_EXTENSIONS
  friend AAF mos(AAF &, AAF &);
  friend AAF arg(AAF &, AAF &);
  friend AAF mag(const AAF &, const AAF &);
  friend AAF magdb(AAF &, AAF &);
#endif
    
    // Sylvie 
  // multiplication, given some constraints on noise symbols
   friend AAF mult_eps (const AAF & P1, const AAF & P2, vector<interval> &constr_eps, int sysdim);
   // hull of 2 affine forms ("argmin" formula) 
   friend AAF hull  (const AAF & P1, const AAF & P2);
   // end Sylvie

private:
  static unsigned inclast();
public:
  void aafprint() const;
    void aafprintline(stringstream &TM_stream, int n) const;
    void aafprintline() const;
};

// binary operators
AAF operator * (double, const AAF);
AAF operator / (double, const AAF);
AAF operator + (double, const AAF);
AAF operator - (double, const AAF);



// unary functions
AAF sqrt(const AAF &);
AAF isqrt(const AAF &);
AAF cotan(const AAF &);
AAF exp(const AAF &);
AAF log(const AAF &);
AAF logexp(const AAF &);
AAF atan(const AAF &);
AAF tanh(const AAF &);
AAF heaviside(const AAF &);

#ifdef USE_AAF_EXTENSIONS
// binary functions
AAF mos(AAF &, AAF &);
AAF arg(AAF &, AAF &);
AAF mag(const AAF &, const AAF &);
AAF magdb(AAF &, AAF &);
#endif

// AAF inline methods

/************************************************************
 * Method:        getcenter
 * Author & Date: ??? - ???
 * Description:   
 *   Returns the central value of the AAF 
 *
 *   Input  : -
 *   Output : double : center value
 ************************************************************/
inline double AAF::getcenter() const
{
  return cvalue;
}


/************************************************************
 * Method:        power
 * Author & Date: Darius Grabowski - 10.05.2005
 * Description:   
 *   Computes the power fct.
 *
 *   Input  : const AAF & : AAF argument
 *            int         : exponent
 *   Output : AAF         : AAF result
 ************************************************************/
inline AAF power(const AAF & P, int exp)
{
  return (P^exp);
}


/************************************************************
 * Method:        sqr
 * Author & Date: Darius Grabowski - 10.05.2005
 * Description:   
 *   Computes the square fct.
 *
 *   Input  : const AAF & : AAF argument
 *   Output : AAF         : AAF result
 ************************************************************/
inline AAF sqr(const AAF & P)
{
  return (P^2);
}


/************************************************************
 * Method:        pow
 * Author & Date: Darius Grabowski - 10.05.2005
 * Description:   
 *   Computes the power fct.
 *
 *   Input  : const AAF & : AAF argument
 *            int         : exponent
 *   Output : AAF         : AAF result
 ************************************************************/
inline AAF pow(const AAF & P, int exp)
{
  return (P^exp);
}


/************************************************************
 * Method:        setDefault
 * Author & Date: ??? - ???
 * Description:   
 *   sets the highest symbol, this function must be treaten 
 *   carefully !!!
 *
 *   Input  : unsigned : highest symbol index
 *   Output : -
 ************************************************************/
#ifdef USE_AAF_EXTENSIONS
inline void AAF::setDefault(const unsigned val)
{
  last = val;
}
#endif

/************************************************************
 * Method:        getDefault
 * Author & Date: Darius Grabowski - 14.06.2005
 * Description:   
 *   returns the highest symbol
 *
 *   Input  : -
 *   Output : unsigned : highest symbol in use
 ************************************************************/
inline unsigned AAF::getDefault(void) 
{
  return last;
}

/************************************************************
 * Method:        inclast
 * Author & Date: ??? - ???
 * Description:   
 *   increases the highest symbol
 *
 *   Input  : -
 *   Output : unsigned : highest symbol to use
 ************************************************************/
inline unsigned AAF::inclast()
{
  return ++last;
}

/************************************************************
 * Method:        getlength
 * Author & Date: ??? - ???
 * Description:   
 *   returns the number of partial deviations
 *
 *   Input  : -
 *   Output : unsigned : number of partial deviations
 ************************************************************/
inline unsigned AAF::getlength() const
{
  return length;
}

/************************************************************
 * Method:        getFirstIndex
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   returns the index of the first deviation
 *
 *   Input  : -
 *   Output : unsigned : first index
 ************************************************************/
inline unsigned AAF::getFirstIndex(void) const
{
  if (length > 0)
  {
    return indexes[0];
  }
  return 0;
}

/************************************************************
 * Method:        getLastIndex
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   returns the index of the last deviation
 *
 *   Input  : -
 *   Output : unsigned : last index
 ************************************************************/
inline unsigned AAF::getLastIndex(void) const
{
  if (length > 0)
  {
    return indexes[length-1];
  }
  return 0;
}

/************************************************************
 * Method:        getApproximationType
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   returns the approximation type
 *
 *   Input  : -
 *   Output : tApproximationType : approximation type
 ************************************************************/
inline tApproximationType AAF::getApproximationType(void)
{
  return approximationType;
}

/************************************************************
 * Method:        setApproximationType
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   sets the approximation type
 *
 *   Input  : tApproximationType : approximation type
 *   Output : -
 ************************************************************/
inline void AAF::setApproximationType(tApproximationType t)
{
  approximationType = t;
}

inline unsigned int AAF::get_class_size(void)
{
  return (sizeof(this) + length * sizeof(indexes) + size * sizeof(deviations));
}

inline unsigned int AAF::get_class_size_index(void)
{
  return (length * sizeof(indexes));
}

inline unsigned int AAF::get_class_size_deviations(void)
{
  return (size * sizeof(deviations));
}


#endif  // AA_AAF_H
