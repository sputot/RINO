/*
 * aa_interval.c -- Implementation of the interval class
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


#include "aa_interval.h"
#include "aa_rounding.h"
#include <cstdio>
#include <iostream>
#include <cmath>

#define PI (4*atan(1.0))


// AAInterval(), modlo(), modhi(), modlohi() are not
// inline functions as they arent used by the lib
// but they are useful for applications


/************************************************************
 * Method:        AAInterval
 * Author & Date: ??? - ???
 * Description:   
 *   Create a null interval
 *   Useful if we only declare a variable
 *
 *   Input  : -
 *   Output : -
 ************************************************************/
AAInterval:: AAInterval():
  lo(0), 
  hi(0)
{
}


/************************************************************
 * Method:        Operator = 
 * Author & Date: ??? - ???
 * Description:   
 *   Affectation operator
 *
 *   Input  : const AAInterval : input Interval
 *   Output : -
 ************************************************************/
AAInterval & AAInterval::operator = (const AAInterval & I)
{
  if (&I != this)
  {
    lo = I.lo;
    hi = I.hi;
  }

  return *this;
}


/************************************************************
 * Method:        modlo
 * Author & Date: ??? - ???
 * Description:   
 *   Modify the lower bound of an Interval
 *
 *   Input  : double : value to modify lo
 *   Output : -
 ************************************************************/
void AAInterval::modlo(const double low)
{
  lo=low;
}


/************************************************************
 * Method:        modhi
 * Author & Date: ??? - ???
 * Description:   
 *   Modify the higher bound of an Interval
 *
 *   Input  : double : value to modify hi
 *   Output : -
 ************************************************************/
void AAInterval::modhi(const double high)
{
  hi=high;
}


/************************************************************
 * Method:        modboth
 * Author & Date: ??? - ???
 * Description:   
 *   Modify the two bounds of an Interval
 *
 *   Input  : double : value to modify lo
 *            double : value to modify hi
 *   Output : -
 ************************************************************/
void AAInterval::modlohi(const double low, const double high)
{
  lo=low;
  hi=high;
}


/************************************************************
 * Method:        mid
 * Author & Date: ??? - ???
 * Description:   
 *   Calculate the midpoint of an Interval
 *
 *   Input  : double : midpoint
 *   Output : -
 ************************************************************/
double AAInterval::mid() const
{

  double t0,t1;

  aa_rnd_t mode = aa_fegetround();

  aa_fesetround(AA_DOWNWARD);
  t0 = lo*0.5;

  aa_fesetround(AA_UPWARD);
  t1 = hi*0.5;

  aa_fesetround(mode);

  return t0+t1;
}


/************************************************************
 * Method:        radius
 * Author & Date: ??? - ???
 * Description:   
 *   Calculate the radius of an Interval
 *
 *   Input  : double : radius
 *   Output : -
 ************************************************************/
double AAInterval::radius() const
{
  double m = mid();

  double t0,t1;

  aa_rnd_t mode = aa_fegetround();

  aa_fesetround(AA_DOWNWARD);
  t0 = m-lo;

  aa_fesetround(AA_UPWARD);
  t1 = hi-m;

  aa_fesetround(mode);

  return (t0 >= t1 ? t0 : t1);
}


/************************************************************
 * Method:        operator >> 
 * Author & Date: ??? - ???
 * Description:   
 *   Istream input of an Interval
 *   Accepts the forms x or [x] or [x,y]
 *
 *   Input  : std::istream & : input stream
 *   Output : AAInterval &   : interval
 ************************************************************/
std::istream & operator >> (std::istream & s, AAInterval & I)
{
  double lo = 0, hi = 0;
  char c = 0;

  s >> c;
  if (c == '[')
    {
      s >> lo >> c;
      if (c == ',') s >> hi >> c;
      //if (c != ']') s.clear(ios_base::badbit);
    }
  else
  {
    //s.putback(c);
    s >> lo;
  }

  if (s) I = AAInterval(lo, hi);

  return s;

}


/************************************************************
 * Method:        intvprint
 * Author & Date: ??? - ???
 * Description:   
 *   Print an Interval to stdout
 *
 *   Input  : -
 *   Output : -
 ************************************************************/
void AAInterval::intvprint() const
{
  printf("[%f,%f]\n", lo, hi);
}

// ajout conversion Sylvie
interval AAInterval::convert_int() const
{
  interval res = interval(lo,hi); // filib::interval<double,filib::native_switched,filib::i_mode_extended>(lo,hi);
    return res;
}




/************************************************************
 * Method:        operator <<
 * Author & Date: ??? - ???
 * Description:   
 *   Ostream output of an Interval
 *
 *   Input  : std::ostream & : output stream
 *            AAInterval &   : interval
 *   Output : -
 ************************************************************/
std::ostream & operator << (std::ostream & s, const AAInterval &I)
{
  // s.setf(0, ios_base::floatfield);
  // cause we don't want to display in scientific format

  s << "[" << I.getlo() << "," << I.gethi() << "]\n";
  return s;
}


/************************************************************
 * Method:        mintrigo
 * Author & Date: ??? - ???
 * Description:   
 *   Caculate the minimal 2PI periodic interval
 *   of an interval, e.g. : [5*PI, 6*PI] -> [PI, 2PI]
 *   In fact it is the closer interval to 0
 *
 *   Input  : const AAInterval & : interval
 *   Output : AAInterval         : interval
 ************************************************************/
AAInterval mintrigo(const AAInterval & I)
{

  // This function is no more needed
  // for our new algorithm of the sine of an AAF

  double a = I.getlo();
  double b = I.gethi();
  double t1, t2 ;

  t1 = floor(a/(2*PI));
  t2 = b-a;

  a=a-(t1*2*PI);
  b=a+t2;

  return AAInterval(a,b);
}
