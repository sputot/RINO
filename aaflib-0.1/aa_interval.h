/*
 * aa_interval.h -- A simple interval class
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


#ifndef AA_INTERVAL_H
#define AA_INTERVAL_H

//#include "vnode.h"
#include "filib_interval.h"

#include <iostream>

// A class for Interval representation
// the class is used by our AAF class

class AAInterval
{
 private:

  double lo, hi;

 public:

  AAInterval();
  AAInterval(double, double);
  AAInterval(interval);
  AAInterval & operator = (const AAInterval &);

  friend std::istream & operator >> (std::istream &, AAInterval &);

  double getlo() const;
  double gethi() const;
  void modlo(const double);
  void modhi(const double);
  void modlohi(const double, const double);
  double mid() const;
  double radius() const;
  double width() const;
  void intvprint() const;
  // ajout SP
  interval convert_int() const;

};

std::ostream & operator << (std::ostream &, const AAInterval &);
AAInterval mintrigo(const AAInterval &);


// Interval inline functions

/************************************************************
 * Method:        AAInterval
 * Author & Date: ??? - ???
 * Description:   
 *   Create an Interval object
 *
 *   Input  : double : lower bound
 *            double : upper bound
 *   Output : -
 ************************************************************/
inline AAInterval:: AAInterval(double l, double h):
     lo(l), hi(h)
{
}

// Ajout SP
inline AAInterval:: AAInterval(interval I):
		  lo(inf(I)), hi(inf(I))
{
}

/************************************************************
 * Method:        getlo
 * Author & Date: ??? - ???
 * Description:   
 *   Get the lower bound of an Interval
 *
 *   Input  : 
 *   Output : double : lower bound
 ************************************************************/
inline double AAInterval::getlo() const
{
  return lo;
}


/************************************************************
 * Method:        getho
 * Author & Date: ??? - ???
 * Description:   
 *   Get the higher bound of an Interval
 *
 *   Input  : 
 *   Output : double : higher bound
 ************************************************************/
inline double AAInterval::gethi() const
{
  return hi;
}


/************************************************************
 * Method:        width
 * Author & Date: ??? - ???
 * Description:   
 *   Calculate the width of an Interval
 *
 *   Input  : 
 *   Output : double : width
 ************************************************************/
inline double AAInterval::width() const
{
  return (hi-lo);
}

#endif  // AA_INTERVAL_H
