/*
 * aa_aafarithm.c -- Affine arithmetical operations
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

/************************************************************
 * Operator:      +=
 * Author & Date: ??? - ???
 * Description:   
 *   Addition is an affine operation
 *
 *   Input  : double : double value to be added
 *   Output : AAF
 ************************************************************/
AAF & AAF::operator += (double cst)
{
  cvalue += cst;
  return (*this);
}


/************************************************************
 * Operator:      -=
 * Author & Date: ??? - ???
 * Description:   
 *   Subtraction is an affine operation
 *
 *   Input  : double : double value to be subtracted
 *   Output : AAF
 ************************************************************/
AAF & AAF::operator -= (double cst)
{
  cvalue -= cst;
  return (*this);
}


/************************************************************
 * Operator:      *=
 * Author & Date: ??? - ???
 * Description:   
 *   Scalar multiplication is an affine operation
 *
 *   Input  : double : double value to be multiplied by
 *   Output : AAF
 ************************************************************/
AAF & AAF::operator *= (double cst)
{
  cvalue *= cst;
  for (unsigned int i = 0; i < length; i++)
    deviations[i] *= cst;

#ifdef FAST_RAD
  radius *= fabs(cst);
#endif

  return (*this);
}

/************************************************************
 * Operator:      /=
 * Author & Date: ??? - ???
 * Description:   
 *   Scalar division is an affine operation
 *
 *   Input  : double : double value to be divided by
 *   Output : AAF
 ************************************************************/
AAF & AAF::operator /= (double cst)
{
  cst = 1.0/cst;
  cvalue *= cst;
  for (unsigned int i = 0; i < length; i++)
    deviations[i] *= cst;

#ifdef FAST_RAD
  radius *= fabs(cst);
#endif

  return (*this);
}


/************************************************************
 * Operator:      +
 * Author & Date: ??? - ???
 * Description:   
 *   Affine addition
 *
 *   Input  : const AAF : AAF to be added
 *   Output : AAF
 ************************************************************/
AAF AAF::operator + (const AAF & P) const
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    AAF Temp(cvalue+P.cvalue);
    return (Temp);
  }

  if (l1 == 0)
  {
    AAF Temp(P);
    Temp += cvalue;
    return (Temp);
  }

  if (l2 == 0)
  {
    AAF Temp(*this);
    Temp += P.cvalue;
    return (Temp);
  }

  // Create our resulting AAF
  AAF Temp(cvalue + P.cvalue);

  unsigned * id1=indexes;
  unsigned * id2=P.indexes;

  double * va1=deviations;
  double * va2=P.deviations;

  unsigned * pu1=id1;
  unsigned * pu2=id2;

  if (l1+l2)
    Temp.indexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  if (ltemp)
    Temp.deviations = new double [ltemp];
  double * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting AAF

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i]=va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]+0
      pu1++;
      continue;
    }
    
    vatempg[i]=va1[a] + va2[b];
    pu1++;
    pu2++;
  }

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif
  
  return Temp;
}


/************************************************************
 * Operator:      +=
 * Author & Date: Darius Grabowski - 07/2007
 * Description:   
 *   Affine addition
 *
 *   Input  : const AAF : AAF to be added
 *   Output : AAF &     : *this
 ************************************************************/
AAF & AAF::operator += (const AAF & P)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    cvalue += P.cvalue;
    return (*this);
  }

  if (l1 == 0)
  {
    double c = cvalue;
    *this = P;
    cvalue += c;
    return (*this);
  }

  if (l2 == 0)
  {
    cvalue += P.cvalue;
    return (*this);
  }

  // Create our resulting AAF

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  double * va1 = deviations;
  double * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  unsigned * tempIndexes = NULL;
  double * tempDeviations = NULL;

  if (l1+l2)
    tempIndexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp = tempIndexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  if (ltemp)
    tempDeviations = new double [ltemp];
  double * vatempg = tempDeviations;

  // Fill the deviations array
  // of the resulting AAF

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i]=va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]+0
      pu1++;
      continue;
    }
    
    vatempg[i]=va1[a] + va2[b];
    pu1++;
    pu2++;
  }

  // set new properties
  length = ltemp;
  size = ltemp;
  delete [] deviations;
  delete [] indexes;
  deviations = tempDeviations;
  indexes = tempIndexes;
  cvalue += P.cvalue;

#ifdef FAST_RAD
  radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    radius += fabs(vatempg[i]);
#endif
  
  return *this;
}


/************************************************************
 * Operator:      -
 * Author & Date: ??? - ???
 * Description:   
 *   Affine subtraction
 *
 *   Input  : const AAF : AAF to be subtracted
 *   Output : AAF
 ************************************************************/
AAF AAF::operator - (const AAF & P) const
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    AAF Temp(cvalue-P.cvalue);
    return (Temp);
  }

  if (l1 == 0)
  {
    AAF Temp(P);
    Temp *= -1.0;
    Temp += cvalue;
    return (Temp);
  }

  if (l2 == 0)
  {
    AAF Temp(*this);
    Temp -= P.cvalue;
    return (Temp);
  }

  // Create our resulting AAF
  AAF Temp(cvalue-P.cvalue);

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  double * va1 = deviations;
  double * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  Temp.indexes = new unsigned [l1+l2];
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  Temp.deviations = new double [ltemp];
  double * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting AAF
  
  for (unsigned i=0;i<ltemp;i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = -va2[b];  // 0-va2[b]
      pu2++;
      continue;
    }
    
    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = va1[a];  // va1[a]-0
      pu1++;
      continue;
    }

    vatempg[i]=va1[a]-va2[b];
    pu1++;
    pu2++;
  }

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif
  
  return Temp;
}


/************************************************************
 * Operator:      -=
 * Author & Date: Darius Grabowski - 07/2007
 * Description:   
 *   Affine addition
 *
 *   Input  : const AAF : AAF to be subtracted
 *   Output : AAF &     : *this
 ************************************************************/
AAF & AAF::operator -= (const AAF & P)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    cvalue -= P.cvalue;
    return (*this);
  }

  if (l1 == 0)
  {
    double c = cvalue;
    *this = -P;
    cvalue += c;
    return (*this);
  }

  if (l2 == 0)
  {
    cvalue -= P.cvalue;
    return (*this);
  }

  // Create our resulting AAF

  unsigned * id1 = indexes;
  unsigned * id2 = P.indexes;

  double * va1 = deviations;
  double * va2 = P.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  unsigned * tempIndexes = NULL;
  double * tempDeviations = NULL;

  if (l1+l2)
    tempIndexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp = tempIndexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin-idtemp;

  if (ltemp)
    tempDeviations = new double [ltemp];
  double * vatempg = tempDeviations;

  // Fill the deviations array
  // of the resulting AAF

  for (unsigned i = 0; i < ltemp; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i]=-va2[b];  // -va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]-0
      pu1++;
      continue;
    }
    
    vatempg[i]=va1[a] - va2[b];
    pu1++;
    pu2++;
  }

  // set new properties
  length = ltemp;
  size = ltemp;
  delete [] deviations;
  delete [] indexes;
  deviations = tempDeviations;
  indexes = tempIndexes;
  cvalue -= P.cvalue;

#ifdef FAST_RAD
  radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    radius += fabs(vatempg[i]);
#endif
  
  return *this;
}


/************************************************************
 * Operator:      -
 * Author & Date: ??? - ???
 * Description:   
 *   Unary operator
 *
 *   Input  : const AAF : AAF to be multiplied by -1
 *   Output : AAF
 ************************************************************/
AAF AAF::operator - () const
{
  AAF Temp(*this);

  Temp.cvalue = -(Temp.cvalue);
  for (unsigned i = 0; i < length; i++)
  {
    Temp.deviations[i] = -(Temp.deviations[i]);
  }

  return Temp;
}


/************************************************************
 * Operator:      *
 * Author & Date: ??? - ???
 * Description:   
 *   Mul by a constant (on right), Affine operation
 *
 *   Input  : double : double to be multiplyed by
 *   Output : AAF
 ************************************************************/
AAF AAF::operator * (double cst)
{
  AAF Temp(*this);
  Temp.cvalue = cst*cvalue;

  for (unsigned i = 0; i < length; i++)
  {
    Temp.deviations[i] = cst*(Temp.deviations[i]);
  }

#ifdef FAST_RAD
  Temp.radius = fabs(cst)*radius;
#endif

  return Temp;
}


/************************************************************
 * Operator:      /
 * Author & Date: ??? - ???
 * Description:   
 *   Division by a constant (on right), Affine operation
 *
 *   Input  : double : double to be divided by
 *   Output : AAF
 ************************************************************/
AAF AAF::operator / (double cst) const
{
  AAF Temp(*this);
  Temp.cvalue = cvalue/cst;

  for (unsigned i = 0; i < length; i++)
  {
    Temp.deviations[i] = (Temp.deviations[i])/cst;
  }

#ifdef FAST_RAD
  Temp.radius = radius/fabs(cst);
#endif
  return Temp;
}

// -- Non member AAF functions --


/************************************************************
 * Operator:      *
 * Author & Date: ??? - ???
 * Description:   
 *   Mul by a constant (the left case)
 *
 *   Input  : double : double factor
 *            AAF    : AAF factor
 *   Output : AAF
 ************************************************************/
AAF operator * (double cst, const AAF P)
{
  AAF Temp(P);
  return Temp*cst;
}


/************************************************************
 * Operator:      +
 * Author & Date: ??? - ???
 * Description:   
 *   Add a constant (the left case)
 *
 *   Input  : double : double summand
 *            AAF    : AAF summand
 *   Output : AAF
 ************************************************************/
AAF operator + (double cst, const AAF P)
{
  AAF Temp(P);
  Temp += cst;
  return (Temp);
}


/************************************************************
 * Operator:      +
 * Author & Date: ??? - ???
 * Description:   
 *   Sub a constant (the left case)
 *
 *   Input  : double : double 
 *            AAF    : AAF 
 *   Output : AAF
 ************************************************************/
AAF operator - (double cst, AAF P)
{
  AAF Temp = -P;
  Temp += cst;
  return (Temp);
}
