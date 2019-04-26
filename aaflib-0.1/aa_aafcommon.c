/*
 * aa_aafcommon.c -- Common functions used to manipulate AAF
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
#include <iostream>
#include <ostream>
#include <fstream>
#include <cstdio>
#include <assert.h>

using namespace std;


// static variables
// defualt starting index for deviations
unsigned AAF::last = 0;
// default approximation type
tApproximationType AAF::approximationType = CHEBYSHEV; // MINRANGE; //CHEBYSHEV;

#ifdef CLEANUP
list <AAF*> AAF::allAAF;
#endif

/************************************************************
 * Method:        AAF
 * Author & Date: ??? - ???
 * Description:   
 *   Default constructor creating an AAF without deviations
 *
 *   Input  : double : corresponding center value
 *   Output : -
 ************************************************************/
AAF::AAF(double v0):
  cvalue(v0), 
  length(0),
  size(0),
#ifdef FAST_RAD
  radius(0.0),
#endif 
  deviations(NULL),
  indexes(NULL)
{
#ifdef CLEANUP
  allAAF.push_back(this);
#endif
}


/************************************************************
 * Method:        AAF
 * Author & Date: ??? - ?AAF operator + (double, const AAF);??
 * Description:   
 *   Default constructor creating an AAF with deviations
 *
 *   Input  : double     : corresponding center value
 *            double *   : array of index values
 *            unsigned * : array of indices
 *            unsigned   : # of indices
 *   Output : -
 ************************************************************/
AAF::AAF(double v0, const double * t1, const unsigned * t2, unsigned T):
  cvalue(v0), 
  length(T),
  size(T),
#ifdef FAST_RAD
  radius(0.0),
#endif
  deviations(new double[T]),
  indexes(new unsigned[T])
{
  for (unsigned i = 0; i < length; i++)
  {
    deviations[i] = t1[i];
#ifdef FAST_RAD
    radius += fabs(deviations[i]);
#endif
    indexes[i] = t2[i];
  }

  if (indexes[length-1] > last) 
    last = indexes[length-1];
#ifdef CLEANUP
  allAAF.push_back(this);
#endif
}


/************************************************************
 * Method:        AAF
 * Author & Date: ??? - ???
 * Description:   
 *   Copy constructor
 *
 *   Input  : AAF : AAF to be copied
 *   Output : -
 ************************************************************/
AAF::AAF(const AAF &P):
  cvalue(P.cvalue), 
  length(P.length),
  size(P.size),
#ifdef FAST_RAD
  radius(P.radius),
#endif 
  deviations(NULL),
  indexes(NULL)
{
  if (size)
  {
    deviations = new double [size];
    indexes = new unsigned [size];
  }

  for (unsigned i = 0; i < length; i++)
  {
    deviations[i] = P.deviations[i];
    indexes[i] = P.indexes[i];
  }

#ifdef CLEANUP
  allAAF.push_back(this);
#endif
}


/************************************************************
 * Method:        AAF
 * Author & Date: ??? - ???
 * Description:   
 *   Constructor creating AAF from interval data
 *
 *   Input  : AAInterval           : Interval
 *   Output : -
 ************************************************************/
AAF::AAF(const AAInterval iv):
  cvalue((iv.gethi()+iv.getlo())/2), 
  length(1),
  size(1),
#ifdef FAST_RAD
  radius(0.0),
#endif 
  deviations(new double[1]),
  indexes(new unsigned[1])
{  
  deviations[0] = (iv.gethi()-iv.getlo())/2;

  indexes[0] = inclast();

#ifdef FAST_RAD
  radius = fabs(deviations[0]);
#endif

#ifdef CLEANUP
  allAAF.push_back(this);
#endif
}


/************************************************************
 * Method:        ~AAF
 * Author & Date: ??? - ???
 * Description:   
 *   Destructor
 *
 *   Input  : -
 *   Output : -
 ************************************************************/
AAF::~AAF()
{
  if (size)
  {
    delete [] deviations;
    delete [] indexes;
  }
#ifdef CLEANUP
  allAAF.remove(this);
#endif
}


/************************************************************
 * Method:        =
 * Author & Date: ??? - ???
 * Description:   
 *   Affectation operator
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
AAF & AAF::operator = (const double d)
{
  length = 0;
  cvalue = d;

#ifdef FAST_RAD
  radius = 0.0;
#endif

  return *this;
}

/************************************************************
 * Method:        =
 * Author & Date: ??? - ???
 * Description:   
 *   Affectation operator
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
AAF & AAF::operator = (const AAF & P)
{
  unsigned plength = P.getlength();

  if (&P != this)
  {
    if (size < plength)
    {
      if (size)
      {
	delete [] deviations;
	delete [] indexes;
      }
      size = plength;
      if (size)
      {
	deviations = new double [size];
	indexes = new unsigned [size];
      }
    }

    cvalue = P.cvalue;
    length = plength;

    for (unsigned i = 0; i < length; i++)
    {
      deviations[i] = P.deviations[i];
      indexes[i] = P.indexes[i];
    }

#ifdef FAST_RAD
    radius = P.radius;
#endif
  }  

  return *this;
}


/************************************************************
 * Method:        <
 * Author & Date: ??? - ???
 * Description:   
 *   
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
bool AAF::operator < (const AAF & P) const
{
  if ((cvalue+rad()) < (P.cvalue-P.rad()))
    return true; 
  else
    return false;
}

/************************************************************
 * Method:        <=
 * Author & Date: ??? - ???
 * Description:   
 *   
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
bool AAF::operator <= (const AAF & P) const
{
  if ((cvalue+rad()) <= (P.cvalue-P.rad()))
    return true;
  else
    return false;
}

/************************************************************
 * Method:        >
 * Author & Date: ??? - ???
 * Description:   
 *   
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
bool AAF::operator > (const AAF & P) const
{
  if ((cvalue-rad()) > (P.cvalue+P.rad()))
    return true; 
  else
    return false;
}

/************************************************************
 * Method:        >=
 * Author & Date: ??? - ???
 * Description:   
 *   
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
bool AAF::operator >= (const AAF & P) const
{
  if ((cvalue-rad()) >= (P.cvalue+P.rad()))
    return true; 
  else
    return false;
}


/************************************************************
 * Method:        ==
 * Author & Date: ??? - ???
 * Description:   
 *   Compare if two affine expressions are equal
 *
 *   Input  : AAF
 *   Output : -
 ************************************************************/
bool AAF::operator == (const AAF & P) const
{
  bool result = true;
  
  // if even the lengths are different we do not have to check further
  if (length != P.length)
    return false;
  
  // if the central value is not equal - there will be no equivalence
  if (fabs(cvalue) < 1 && fabs(P.cvalue) < 1)
  {
    if (fabs(cvalue - P.cvalue) > 1E-15)
      return false; 
  }
  else
  {
    if ( fabs(cvalue - P.cvalue) / (fabs(cvalue) + fabs(P.cvalue)) > 1E-15)
      return false;
  }
  // compare if all deviations and their respective indexes both are equal
  for (unsigned i = 0; i < length; i++) 
  {	
    if ( fabs(deviations[i]) < 1 &&  fabs(P.deviations[i]) < 1)
    {
      if ( fabs( deviations[i] - P.deviations[i]) > 1E-15 ||
	   (indexes[i] != P.indexes[i]) )
      {
	result = false;
	break;
      }
    }
    else
    {
      if ( fabs(deviations[i] - P.deviations[i]) /  
	   (fabs(deviations[i]) + fabs(P.deviations[i]) )  > 1E-15 ||
	   (indexes[i] != P.indexes[i]) )
      {
	result = false;
	break;
      }
      
    }
  }
  return result;
}


/************************************************************
 * Method:        <<
 * Author & Date: ??? - ???
 * Description:   
 *   Ostream output of an AAF
 *
 *   Input  : ostream
 *            AAF
 *   Output : -
 ************************************************************/
std::ostream & operator << (std::ostream & s, const AAF & P)
{

  // s.setf(0, ios_base::floatfield);
  s << "-------------\n";
  s << "Length = " << P.length << "\n";
  s << "v0     = " << P.cvalue << "\n";
#ifdef FAST_RAD
  s << "Radius = " << P.rad() << "\n";
#endif

  for (unsigned i = 0; i < P.length ; i++)
  {
    s << "e" << P.indexes[i] << " -> " << P.deviations[i] << "\n";
  }
  s << "-------------\n";

  return s;
}


/************************************************************
 * Method:        aafprint
 * Author & Date: ??? - ???
 * Description:   
 *   Print length and deviations of an AAF to stdout
 *
 *   Input  : -
 *   Output : -
 ************************************************************/
void AAF::aafprint() const
{
  std::cout << "-------------\n";
  std::cout << "Length = " << length << "\n";
  std::cout << "v0 = " << cvalue << "\n";

  for (unsigned i = 0; i < length ; i++)
  {
    std::cout << "e" << indexes[i] << " -> " << deviations[i] << "\n";
  }  
  std::cout << "-------------\n";
}

/************************************************************
 * Method:        convert
 * Author & Date: ??? - ???
 * Description:   
 *   Convert an AAF to an Interval representation
 *
 *   Input  : -
 *   Output : AAInterval : Interval representation
 ************************************************************/
AAInterval AAF::convert() const
{
  // lower bound == central value of the AAF - the total deviation
  // upper bound == central value of the AAF + the total deviation

  AAInterval Temp(cvalue-rad(), cvalue+rad());
  return Temp;
}


/************************************************************
 * Method:        rad
 * Author & Date: ??? - ???
 * Description:   
 *   Get the total deviation of an AAF
 *   i.e. the sum of all deviations (their abs value)
 *
 *   Input  : -
 *   Output : double : total deviation
 ************************************************************/
double AAF::rad() const
{
#ifdef FAST_RAD
  return radius;
#else
  double sum = 0.0;

  for (unsigned i = 0; i < length; i++)
  {
    if (deviations[i] >= 0.0)
      sum+=deviations[i];
    else
      sum+=-deviations[i];
  }

  return sum;
#endif
}


/************************************************************
 * Method:        getMax
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   Get the total maximum of an AAF
 *
 *   Input  : -
 *   Output : double : total maximum
 ************************************************************/
double AAF::getMax() const
{
  return getcenter()+rad();
}

/************************************************************
 * Method:        getMin
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   Get the total minimum of an AAF
 *
 *   Input  : -
 *   Output : double : total minimum
 ************************************************************/
double AAF::getMin() const
{
  return getcenter()-rad();
}

/************************************************************
 * Method:        getAbsMax
 * Author & Date: Michael Kaergel - 08/2009
 * Description:   
 *   Get the total absolute maximum of an AAF
 *
 *   Input  : -
 *   Output : double : total maximum
 ************************************************************/
double AAF::getAbsMax() const
{
  if ( fabs(getcenter()-rad()) > fabs(getcenter()+rad()))
    return fabs(getcenter()-rad());
  else
    return fabs(getcenter()+rad());
}

/************************************************************
 * Method:        getAbsMin
 * Author & Date: Michael Kaergel - 08/2009
 * Description:   
 *   Get the total absolute minimum of an AAF
 *
 *   Input  : -
 *   Output : double : total minimum
 ************************************************************/
double AAF::getAbsMin() const
/*{
  if ((getcenter()-rad()) < 0)
    return 0;
  else
    return getcenter()-rad();
    }*/
{
  if ( ( ((getcenter()+rad()) > 0) && ((getcenter()-rad()) > 0) ) || ( ((getcenter()+rad()) < 0) && ((getcenter()-rad()) < 0) ) )
    if ( fabs(getcenter()-rad()) > fabs(getcenter()+rad()))
      return fabs(getcenter()+rad());
    else
      return fabs(getcenter()-rad());
  else
    return 0;
}


/************************************************************
 * Method:        []
 * Author & Date: Darius Grabowski - 05/2005
 * Description:   
 *   Index operator
 *
 *   Input  : unsigned : index
 *   Output : double   : reference to the corresponding index
 ************************************************************/
double AAF::operator[](unsigned ind) const
{  
  if (ind < 0 || ind > length)
  {
    printf("Index out of bounds [0..%d]...\n", length);
    // return cvalue;
    return 0.0;
  }
  if (ind == 0)
    return cvalue;
  return deviations[ind-1];
}


//last waechst und dadurch auch der Aufwand - das muss noch geaendert werden
/************************************************************
 * Method:        cleanup
 * Author & Date: Wilhelm Heupke - ???
 * Description:   
 *   cleanup ("garbage collection")
 *
 *   Input  : double : cleanup level
 *   Output : -
 ************************************************************/
#ifdef CLEANUP
void AAF::cleanup(double level) 
{
#ifdef CLEANUP_DEBUG
  cout << "Call to cleanup():" << endl;
  cout << "highest deviation term index is " << last << "." << endl;
  // the "last" member variable has got to be changed to prevent an uncontrolled increase
#endif
  
  // the maximum amplitudes found
  double * amplitudes;
  amplitudes = new double [last];
  
  // clear the amplitudes
  unsigned i;
  for (i = 0; i < last; i++) 
  {
    amplitudes[i] = 0;
  }
  
  unsigned max = 0; // the highest index that can be found
  unsigned myindex;
  double mydeviation;
  
  // iteration over all AAF variables in use at the moment to 
  // find the significant deviation symbol in _all_ variables
  list <AAF*>::const_iterator iter;
  for (iter = allAAF.begin(); iter != allAAF.end(); iter++)
  {
    // is there a higher index - if so change max    
    if (max < ((*iter)->length))
      max = ((*iter)->length); 

    // iteration across all used indexes in this variable
    for (i = 0; i<((*iter)->length); i++) 
    {
      myindex = (*iter)->indexes[i];
      mydeviation = (*iter)->deviations[i];
	
      if ( fabs(amplitudes[myindex]) < fabs(mydeviation) ) 
      {
	amplitudes[myindex] = fabs(mydeviation);
      }
    }
  }
  
  last = max;  // prevent last from growing uncontrolled
  unsigned number_of_signifcant;
  double * new_deviations;
  unsigned * new_indexes;
  // now reconstruct all affine variables to contain only relevant deviation terms
  for (iter = allAAF.begin(); iter != allAAF.end(); iter++)
  {
#ifdef CLEANUP_DEBUG
    // document how the expression looked like before
    cout << (*iter) << " L=" << (*iter)->length << " ; " << *(*iter) << endl;
#endif
    
    // check the number of significant deviation symbols in _this_ variable
    // to know how large the two new arrays have to be
    number_of_signifcant = 0;
    for (i = 0; i < ((*iter)->length); i++) 
    {
      if (fabs((*iter)->deviations[i]) > level) 
	number_of_signifcant++;
    }
    
#ifdef CLEANUP_DEBUG
    // debugging information for the step above
    cout << "Found " << number_of_signifcant << " significant deviation symbol(s)." << endl;
#endif
    // don't reserve space if there are no significant left, else start to clean
    if (number_of_signifcant == 0) 
    {
      new_deviations = NULL;
      new_indexes = NULL;
    }
    else 
    {
      // reserve space for the new indexes and deviations if needed
      new_deviations = new double [number_of_signifcant];
      new_indexes = new unsigned [number_of_signifcant];
      
      // ... and copy the significant deviation symbols in _this_ variable
      myindex = 0;
      for (i = 0; i < ((*iter)->length); i++) 
      {
	if (fabs((*iter)->deviations[i]) > level) 
	{
	  new_deviations[myindex] = (*iter)->deviations[i];
	  new_indexes[myindex] = (*iter)->indexes[i];
	  myindex++;
	}
      }
    }
    
    // free space of the old indexes and deviations	  
    delete [] (*iter)->deviations;
    delete [] (*iter)->indexes;
    
    // assign the new array to this variable
    (*iter)->deviations = new_deviations;
    (*iter)->indexes = new_indexes;
    
    // correct the length
    (*iter)->length = number_of_signifcant;
    
#ifdef CLEANUP_DEBUG
    cout << "new compacted:" << endl;
    cout << (*iter) << " L=" << (*iter)->length << " ; " << *(*iter) << endl;
    cout << "-----------------------------------------" << endl;
#endif
  }
  
  delete [] amplitudes;
#ifdef CLEANUP_DEBUG
  cout << "Maximum number of _used_ terms is " << max << endl;
#endif
}
#endif
