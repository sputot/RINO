/*
 * aa_aafspecial.c -- Special operations
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

/**
 * Class: AAF
 * This class includes special operations on affine forms
 * used for analog circuit simulation with parameter uncertainties
 *
 * @author Darius Grabowski <darius.grabowski@ims.uni-hannover.de>
 * @date 2009-04-01
 * @version 1
 *
 * @see http://www.ims.uni-hannover.de/
 *
 * (c) 2008-2009
 * Institute of Microelectronic Systems (IMS)
 * Leibniz Universitaet Hannover
 *
 * This program comes without any warranty;
 * without even the implied warranty of
 * merchantability or fitness for a particular purpose.
 * Proprietary and confidential.
 * Distribution only by express authority of the IMS.
 */


#include "aa.h"
#ifdef USE_QPF
#include "qpf.h"
#endif

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cmath>

#define AAF_THRESHOLD 1E-15

#ifdef USE_AAF_EXTENSIONS

#ifdef USE_QPF
/************************************************************
 * Method:        =
 * Author & Date: Darius Grabowski - 23.08.2006
 * Description:
 *   Affectation operator: AAF = QPF
 *
 *   Input  : QPF
 *   Output : this
 ************************************************************/
AAF & AAF::operator = (const QPF & P)
{
  cvalue = P.getCenter();

  if (size)
  {
    delete [] deviations;
    delete [] indexes;
  }

  length = P.getAAFLength();
  size = length;

  if (size)
  {
    indexes = new unsigned[size];
    deviations = new double[size];

    P.getAAF(indexes, deviations);
    // check last index
    if (indexes[length-1] > last)
      last = indexes[length-1];
#ifdef FAST_RAD
    radius = P.rad();
#endif
  }
  else
  {
    deviations = NULL;
    indexes = NULL;
#ifdef FAST_RAD
    radius = 0.0;
#endif
  }
  return *this;
}
#endif

/************************************************************
 * Method:        submul
 * Author & Date: Darius Grabowski - 05/2008
 * Description:
 *   Affine substraction with scaling: this = this - d*P
 *
 *   Input  : const AAF : AAF to be subtracted
 *          : double    : scaling
 *   Output : -
 ************************************************************/
void AAF::submul(const AAF & P, double d)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  if (l1+l2 == 0)
  {
    cvalue -= d*P.cvalue;
    return;
  }

  if (l1 == 0)
  {
    double c = cvalue;
    *this = -P;
    *this *= d;
    cvalue += c;
    return;
  }

  if (l2 == 0)
  {
    cvalue -= d*P.cvalue;
    return;
  }

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
      vatempg[i]=-d*va2[b];  // -va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i]=va1[a];  // va1[a]-0
      pu1++;
      continue;
    }

    vatempg[i]=va1[a] - d*va2[b];
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
  cvalue -= d*P.cvalue;

#ifdef FAST_RAD
  radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    radius += fabs(vatempg[i]);
#endif

  return;
}


/************************************************************
 * Method:        getData
 * Author & Date: Darius Grabowski - 05/2005
 * Description:
 *   Get the total data of the AAF
 *
 *   Input  : unsigned * : pointer to allocated index array
 *            double   * : pointer to allocated deviation array
 *   Output : unsigned * : filled index array
 *            double   * : filled deviation array
 *            unsigned   : array length
 ************************************************************/
unsigned AAF::getData(unsigned * destIndexes, double * destDeviations) const
{
  if (length == 0)
    return length;

  for (unsigned i = 0; i < length; i++)
  {
    destIndexes[i] = indexes[i];
    destDeviations[i] = deviations[i];
  }

  return length;
}


/************************************************************
 * Method:        get
 * Author & Date: Darius Grabowski - 01/2008
 * Description:
 *   get a part of the AAF
 *
 *   Input  : double   * : pointer to allocated deviation array
 *            unsigned   : number of deviations to be stored
 *            unsigned   : increment
 *   Output : double   * : pointer to allocated deviation array
 ************************************************************/
void AAF::get(double * data, unsigned n, unsigned npp)
{
  if (length+1 < n)
    return;

  unsigned index = 0;
#ifdef FAST_RAD
  double dev = 0.;
#endif

  data[index] = cvalue;
  cvalue = 0.;
  index += npp;

  for (unsigned i = 0; i < n-1; i++)
  {
#ifdef FAST_RAD
    dev += fabs(deviations[i]);
#endif
    data[index] = deviations[i];
    deviations[i] = 0.;
    index += npp;
  }

#ifdef FAST_RAD
  radius -= dev;
#endif

  return;
}

/************************************************************
 * Method:        getPD
 * Author & Date: Darius Grabowski - 11/2008
 * Description:
 *   get partial deviations of the AAF
 *
 *   Input  : double   * : pointer to allocated deviation array
 *            unsigned   : first PD to be stored
 *            unsigned   : number of deviations to be stored
 *            unsigned   : increment
 *   Output : double   * : pointer to allocated deviation array
 ************************************************************/
void AAF::getPD(double * data, unsigned first, unsigned n, unsigned npp)
{
  if (length < n+first)
    return;

  unsigned index = 0;

  for (unsigned i = first; i < first+n; i++)
  {
    data[index] = deviations[i];
    index += npp;
  }

  return;
}


/************************************************************
 * Method:        getAt
 * Author & Date: Darius Grabowski - 01/2008
 * Description:
 *   get a part of the AAF
 *
 *   Input  : double   * : pointer to allocated deviation array
 *            unsigned   : first epsilon to be stored
 *            unsigned   : number of deviations to be stored
 *            unsigned   : increment
 *   Output : double   * : pointer to allocated deviation array
 ************************************************************/
void AAF::getAt(double * data, unsigned ind, unsigned n, unsigned npp)
{
  if (length == 0)
    return;

  unsigned iter;
  unsigned index = 0;
  unsigned min = 0;
  unsigned max = length-1;
  unsigned mid;
#ifdef FAST_RAD
  double rad_ = 0.0;
#endif

  // do a binary search for the epsilon given by ind
  while (max > min+1)
  {
    mid = (max+min)/2;
    if (indexes[mid] > ind)
      max = mid;
    else
      min = mid;
  }
  if (indexes[min] >= ind)
    iter = min;
  else
    iter = max;

  for (unsigned i = 0; i < n; i++)
  {
    if (indexes[iter] == ind)
    {
      data[index] = deviations[iter];
#ifdef FAST_RAD
      rad_ += fabs(deviations[iter]);
#endif
      deviations[iter] = 0.0;
      if (iter < length-1)
	iter++;
    }
    else
    {
      data[index] = 0.0;
    }
    index += npp;
    ind++;
  }

#ifdef FAST_RAD
  radius -= rad_;
#endif
  return;
}


/************************************************************
 * Method:        update
 * Author & Date: Darius Grabowski - 01/2008
 * Description:
 *   update a part of the AAF
 *
 *   Input  : double   * : pointer to allocated deviation array
 *            unsigned   : number of deviations to be stored
 *            unsigned   : increment
 *   Output : -
 ************************************************************/
void AAF::update(double * data, unsigned n, unsigned npp)
{
  if (length+1 < n)
    return;

#ifdef FAST_RAD
  double radOld = 0.0;
  double radNew = 0.0;
#endif
  unsigned index = 0;

  cvalue -= data[index];
  index += npp;

  for (unsigned i = 0; i < n-1; i++)
  {
#ifdef FAST_RAD
    radOld += fabs(deviations[i]);
#endif
    deviations[i] -= data[index];
#ifdef FAST_RAD
    radNew += fabs(deviations[i]);
#endif
    index += npp;
  }

#ifdef FAST_RAD
  radius += radNew - radOld;
#endif

  return;
}


/************************************************************
 * Method:        at
 * Author & Date: ??? - ???
 * Description:
 *   Get the deviation symbol at "index"
 *
 *   Input  : unsigned : index
 *   Output : double   : deviation symbol at "index"
 ************************************************************/
#ifdef LINEAR_SEARCH
double AAF::at(unsigned index) const
{
  for (unsigned i = 0; i < length; i++)
  {
    if (indexes[i] == index)
      return (deviations[i]);
  }

  return 0.0;
}
#else
double AAF::at(unsigned index) const
{
  if (length == 0)
    return 0.0;
  unsigned min = 0;
  unsigned max = length-1;
  unsigned mid;

  while (max > min+1)
  {
    mid = (max+min)/2;
    if (indexes[mid] > index)
      max = mid;
    else
      min = mid;
  }
  if (indexes[min] == index)
    return deviations[min];
  if (indexes[max] == index)
    return deviations[max];

  return 0.0;
}
#endif


/************************************************************
 * Method:        set
 * Author & Date: Darius Grabowski - 05/2005
 * Description:
 *   Set the deviation symbol at "index"
 *
 *   Input  : unsigned : index
 *            double   : deviation symbol at "index"
 *   Output : -
 ************************************************************/
#ifdef LINEAR_SEARCH
void AAF::set(unsigned index, double value)
{
  for (unsigned i = 0; i < length; i++)
  {
    if (indexes[i] == index)
    {
#ifdef FAST_RAD
      radius -= fabs(deviations[i]) - fabs(value);
#endif
      deviations[i] = value;
      break;
    }
  }
}
#else
void AAF::set(unsigned index, double value)
{
  unsigned min = 0;
  unsigned max = length-1;
  unsigned mid;

  while (max > min+1)
  {
    mid = (max+min)/2;
    if (indexes[mid] > index)
      max = mid;
    else
      min = mid;
  }
  if (indexes[min] == index)
  {
#ifdef FAST_RAD
    radius -= fabs(deviations[min]) - fabs(value);
#endif
    deviations[min] = value;
    return;
  }
  if (indexes[max] == index)
  {
#ifdef FAST_RAD
    radius -= fabs(deviations[max]) - fabs(value);
#endif
    deviations[max] = value;
    return;
  }
  std::cout << "Warning <AAF::set>: Could not find index " << index << "\n";
}
#endif


/************************************************************
 * Method:        setLast
 * Author & Date: Darius Grabowski - 08/2006
 * Description:
 *   Set the deviation symbol at last index
 *
 *   Input  : double   : deviation symbol at "index"
 *   Output : -
 ************************************************************/
void AAF::setLast(double value)
{
  if (length)
  {
#ifdef FAST_RAD
    radius -= fabs(deviations[length-1]) - fabs(value);
#endif
    deviations[length-1] = value;
  }
}


/************************************************************
 * Method:        sumup
 * Author & Date: Darius Grabowski - 07/2005
 * Description:
 *   Sum up all partial deviations with index higher than ind
 *
 *   Input  : unsigned int : sumup index
 *   Output : new highest index
 ************************************************************/
unsigned AAF::sumup(unsigned ind)
{
  if (ind == 0)
    return 0;

  if (length == 0)
    return length;

  if (indexes[length-1] <= ind)
    return indexes[length-1];

  double sumUpDeviation = 0.0;
  unsigned nIndex;

  // sum up deviations until index falls below ind
  for (nIndex = length-1; nIndex >= 0; nIndex--)
  {
    if (indexes[nIndex] <= ind)
    {
      nIndex++;
      break;
    }
    sumUpDeviation += fabs(deviations[nIndex]);
  }
  if (nIndex == length-1)
    return indexes[length-1];

  length = nIndex+1;

  // replace last deviation
  deviations[nIndex] = sumUpDeviation;
  // replace last index
  indexes[nIndex] = inclast();

  return indexes[nIndex];
}

/************************************************************
 * Method:        sumup
 * Author & Date: Darius Grabowski - 08/2005
 * Description:
 *   sums up "small" deviations depending on "level", which is
 *   scaled by the total deviation (radius)
 *
 *   Input  : double : sumup level
 *   Output : -
 ************************************************************/
void AAF::sumup(double level)
{
  if (!length)
    return;

  int nSumUpTerms = 0;
  int nIndex = 0;
  double sumUpDeviation = 0.0;
  double dThreshold = fabs(level*rad());

  if (dThreshold < AAF_THRESHOLD)
    dThreshold = AAF_THRESHOLD;

  for (unsigned int i = 0; i < length; i++)
  {
    if (fabs(deviations[i]) >= dThreshold)
    {
      // deviation
      deviations[nIndex] = deviations[i];
      indexes[nIndex] = indexes[i];
      nIndex++;
    }
    else
    {
      sumUpDeviation += fabs(deviations[i]);
      nSumUpTerms++;
    }
  }

  if (nSumUpTerms)
  {
    length = nIndex+1;
    indexes[length-1] = inclast();
    deviations[length-1] = sumUpDeviation;
  }

  return;
}

/************************************************************
 * Method:        sumupall
 * Author & Date: Michael Kaergel - 09/2009
 * Description:
 *   Sum up all partial deviations
 *
 *   Input  :
 *   Output :
 ************************************************************/
void AAF::sumupall(void)
{
  if (length == 0)
    return;

  double sumUpDeviation = 0.0;

  // sum up deviations until index falls below ind
  for (unsigned int nIndex = 0; nIndex < length; nIndex++)
  {
    sumUpDeviation += fabs(deviations[nIndex]);
  }

  length = 1;

  // replace last deviation
  deviations[length - 1] = sumUpDeviation;
  // replace last index
  indexes[length - 1] = inclast();

  return;
}

/************************************************************
 * Method:        compact
 * Author & Date: Darius Grabowski - 05/2008
 * Description:
 *   removes zero deviations
 *
 *   Input  : -
 *   Output : -
 ************************************************************/
void AAF::compact(void)
{
  if (!length)
    return;

  int newLength = 0;

  for (unsigned int i = 0; i < length; i++)
  {
    if (deviations[i] != 0.)
    {
      deviations[newLength] = deviations[i];
      indexes[newLength] = indexes[i];
      newLength++;
    }
  }

  length = newLength;
}


/************************************************************
 * Method:        resize
 * Author & Date: Darius Grabowski - 04/2007
 * Description:
 *   Resizes the AAF by distributing all PDs with index > ind
 *   on PDs with index <= ind
 *
 *   Input  : unsigned int : marginal index
 *   Output : -
 ************************************************************/
void AAF::resize(unsigned int ind)
{
  double dDeviation = 0.0;
  double dScale;
  int i;

  if (!length)
    return;

  if ((ind == 0) or (ind > indexes[length-1]))
    return;

  // get total deviation of indexes <= ind
  for (i = (int)length-1; i >= 0; i--)
  {
    if (indexes[i] > ind)
      dDeviation += fabs(deviations[i]);
    else
      break;
  }

  // set new length
  length = unsigned (i+1);

  // set scaling factor --> radius remains the same
  dScale = rad()/(rad() - dDeviation + AAF_THRESHOLD);

  // update deviations
  for (i = 0; i < (int)length; i++)
  {
    deviations[i] *= dScale;
  }
}

/************************************************************
 * Method:        ResizeNewSymbol
 * Author & Date: Michael Kaergel - 09/2009
 * Description:
 *   Resizes the AAF by putting the Sum of all PDs with
 *   index > ind on new PD
 *   This is needed for FFT, else the result is wrong, because
 *   of wrong correlation
 * ! Should be used with extreme care, dirty implementation !!
 * ! there is no new Array for deviations and index created !!
 * ! so memory is not freed                                 !!
 *
 *   Input  : unsigned int : marginal index
 *   Output : -
 ************************************************************/
void AAF::ResizeNewSymbol(unsigned int ind)
{
  //std::cout << "Starting Resize of AAF" << std::endl;
  double dDeviation = 0.0;
  unsigned int newLength = 0;

  if (!length)
    return;

  if ((ind == 0) or (ind > indexes[length-1]))
    return;

  // get total deviation of indexes <= ind
  for (int i = (int)length - 1; i >= 0; i--)
  {
    if (indexes[i] > ind)
      dDeviation += fabs(deviations[i]);
    else
      {
	newLength = (unsigned int)(i);
	break;
      }
  }
  //std::cout << "Done Calculating Sum of Deviations to remove" << std::endl;

  // set new length, +1 because of array starting at 0, +1 because of new PD
  //length += 2;

  //indexes[length - 1] = inclast();
  //deviations[length - 1] = dDeviation;
  //increase length by one for new Symbol
  if (dDeviation != 0) newLength += 1;

  //std::cout << "New Length is: "<< newLength << ", Old Length was: "<< length <<std::endl;
  // create new arrays
  double * tempDeviations = new double[newLength];
  unsigned * tempIndexes = new unsigned[newLength];

  // copy deviations and indexes
  for (unsigned int i = 0; i < newLength - 1; i++)
  {
    tempDeviations[i] = deviations[i];
    tempIndexes[i] = indexes[i];
  }

  // add new deviation and index
  if (dDeviation != 0)
    {
      tempIndexes[newLength - 1] = inclast();
      tempDeviations[newLength - 1] = dDeviation;
    }

  if (length > 1)
    {
      delete [] deviations;
      delete [] indexes;
    }
  else if (length == 1)
    {
      delete deviations;
      delete indexes;
    }

  // update AAF data
  deviations = tempDeviations;
  indexes = tempIndexes;
  length = newLength;
  //std::cout << "Done Calculating new reduced AAF" << std::endl;
}


/************************************************************
 * Method:        resize
 * Author & Date: Darius Grabowski - 04/2007
 * Description:
 *   Resizes the AAF by distributing all PDs with index > ind
 *   on PDs with index <= ind
 *
 *   Input  : unsigned int : marginal index
 *   Output : -
 ************************************************************/
void AAF::resize(unsigned int ind, unsigned int n)
{
  double dDeviation = 0.0;
  double dScale;
  unsigned int i;

  if(!length)
    return;

  if (ind+n > indexes[length-1])
    return;

  // get total deviation of indexes <= ind
  for (i = ind; i < ind+n; i++)
  {
    dDeviation += fabs(deviations[i]);
  }

  // set scaling factor --> radius remains the same
  dScale = rad()/(rad() - dDeviation + AAF_THRESHOLD);

  // update deviations with index smaller than ind
  for (i = 0; i < ind; i++)
  {
    deviations[i] *= dScale;
  }

  // udate deviations with index larger than ind+n
  for (i = ind+n; i < length; i++)
  {
    deviations[i-n] = dScale*deviations[i];
    indexes[i-n] = indexes[i];
  }

  // set new length
  length -= n;
}


/************************************************************
 * Method:        addPD
 * Author & Date: Darius Grabowski - 04/2007
 * Description:
 *   Adds a new partial deviation to the AAF
 *
 *   Input  : double : new partial deviation
 *   Output : -
 ************************************************************/
void AAF::addPD(double pd)
{
#ifdef FAST_RAD
  radius += fabs(pd);
#endif

  if (size > length)
  {
    // array size > array length: just add new deviation
    deviations[length]  = pd;
    indexes[length] = inclast();
    length++;
    return;
  }

  // resize the arrays
  size = length+1;

  // create new arrays
  double * tempDeviations = new double[size];
  unsigned * tempIndexes = new unsigned[size];

  // copy deviations and indexes
  for (unsigned int i = 0; i < length; i++)
  {
    tempDeviations[i] = deviations[i];
    tempIndexes[i] = indexes[i];
  }

  // add new deviation and index
  tempIndexes[length] = inclast();
  tempDeviations[length] = pd;

  if (length > 1)
  {
    delete [] deviations;
    delete [] indexes;
  }

  // update AAF data
  deviations = tempDeviations;
  indexes = tempIndexes;
  length++;
}


/************************************************************
 * Method:        ID
 * Author & Date: Darius Grabowski - 10.04.2007
 * Description:
 *   computes simplified drain current
 *
 *   Input  : double : VGS
 *            double : VDS
 *   Output : double : ID
 ************************************************************/
double ID(double VGS, double VDS)
{
  if (VDS > VGS)
    return (0.5*VGS*VGS);
  else
    return ((VGS - 0.5*VDS)*VDS);
}


/************************************************************
 * Method:        mos
 * Author & Date: Darius Grabowski - 10.04.2007
 * Description:
 *   computes simplified drain current
 *
 *   Input  : AAF : VGS
 *            AAF : VDS
 *   Output : AAF : ID
 ************************************************************/
AAF mos(AAF & VGS, AAF & VDS)
{
  AAF dV = VDS - VGS;

#ifdef MOSDEBUG
  fprintf(stdout, "dv: %f | [%f..%f]\n", dV.getcenter(), dV.getMin(), dV.getMax());
  // dV.aafprint();
#endif

  if (dV.getMin() > 0.0)
  {
#ifdef MOSDEBUG
    fprintf(stdout, "Saturation\n");
#endif
    // saturation
    AAF temp = 0.5*(VGS^2);
    return (temp);
  }
  if (dV.getMax() < 0.0)
  {
#ifdef MOSDEBUG
    fprintf(stdout, "Active region\n");
#endif
    // active region
    AAF temp = (VGS - 0.5*VDS)*VDS;
    return (temp);
  }

  //#ifdef MOSDEBUG
  fprintf(stdout, "Saturation <--> Active region\n");
  //#endif

  double VDSmin = VDS.getMin();
  double VDSmax = VDS.getMax();
  double VGSmin = VGS.getMin();
  double VGSmax = VGS.getMax();

  double I_VGSmin_VDSmin = ID(VGSmin, VDSmin);
  double I_VGSmax_VDSmin = ID(VGSmax, VDSmin);
  double I_VGSmin_VDSmax = ID(VGSmin, VDSmax);
  double I_VGSmax_VDSmax = ID(VGSmax, VDSmax);

  double dIDdVGS = 0.5*(I_VGSmax_VDSmin - I_VGSmin_VDSmin + I_VGSmax_VDSmax - I_VGSmin_VDSmax)/(VGSmax - VGSmin);
  double dIDdVDS = 0.5*(I_VGSmin_VDSmax - I_VGSmin_VDSmin + I_VGSmax_VDSmax - I_VGSmax_VDSmin)/(VDSmax - VDSmin);

  //  double dMin = -0.5*dIDdVGS*dIDdVGS - dIDdVDS*VDSmax;
  //  double dMax = 0.5*(VGSmax*VGSmax + dIDdVDS*dIDdVDS) - (dIDdVGS + dIDdVDS)*VGSmax;

  //  double dOffs = -0.5*(dMin + dMax);
  //  double dDelta = 0.5*(dMax - dMin);

#ifdef MOSDEBUG
  fprintf(stdout, "dIDdVGS: %f, dIDdVDS: %f\n", dIDdVGS, dIDdVDS);
  // fprintf(stdout, "dMin: %f, dMax: %f\n", dMin, dMax);
  fprintf(stdout, "offset: %f, delta: %f\n", dOffs, dDelta);
#endif

  // compute id_min and id_max at central value
  double IDmin = ID(dIDdVGS, VDSmax);
  double IDmax = ID(VGSmax, -dIDdVDS + VGSmax);

  IDmin += dIDdVGS*(VGS.getcenter() - dIDdVGS) + dIDdVDS*(VDS.getcenter() - VDSmax);
  IDmax += dIDdVGS*(VGS.getcenter() - VGSmax) + dIDdVDS*(VDS.getcenter() + dIDdVDS - VGSmax);

  double dOffs = 0.5*(IDmin + IDmax);
  double dDelta = 0.5*(IDmax - IDmin);

#ifdef MOSDEBUG
  fprintf(stdout, "Corrected:\n");
  fprintf(stdout, "IDmin: %f, IDmax: %f\n", IDmin, IDmax);
  fprintf(stdout, "Offset: %f, Delta: %f\n", 0.5*(IDmin + IDmax), 0.5*(IDmax - IDmin));
#endif

  AAF Temp1(dOffs);

  Temp1.length = VGS.length + 1;
  Temp1.size = Temp1.length;
  Temp1.deviations = new double [Temp1.size];
  Temp1.indexes = new unsigned [Temp1.size];

  for (unsigned i = 0; i < VGS.length; i++)
  {
    Temp1.indexes[i] = VGS.indexes[i];
    Temp1.deviations[i] = dIDdVGS*(VGS.deviations[i]);
  }

  // Compute the error in a new deviation symbol
  // zk = delta
  Temp1.indexes[Temp1.length-1] = Temp1.inclast();
  Temp1.deviations[Temp1.length-1] = dDelta;

  AAF Temp2(0.0);

  Temp2.length = VDS.length;
  Temp2.size = Temp2.length;
  Temp2.deviations = new double [Temp2.size];
  Temp2.indexes = new unsigned [Temp2.size];

  for (unsigned i = 0; i < VDS.length; i++)
  {
    Temp2.indexes[i] = VDS.indexes[i];
    Temp2.deviations[i] = dIDdVDS*(VDS.deviations[i]);
  }

  AAF Temp = Temp1 + Temp2;

  AAF test = 0.5*(VGS^2);

#ifdef MOSDEBUG
  Temp.aafprint();
  test.aafprint();

  VGS.aafprint();
  VDS.aafprint();
  Temp.aafprint();
#endif

  return (Temp);
}


/************************************************************
 * Method:        arg
 * Author & Date: Darius Grabowski - 06/2007
 * Description:
 *   computes the arg of (re|im)
 *
 *   Input  : const AAF & : AAF real part
 *            const AAF & : AAF imaginary part
 *   Output : AAF         : AAF result
 ************************************************************/
AAF arg(AAF & re, AAF & im)
{
  double dRe = re.getcenter();
  double dIm = im.getcenter();

  double dMag_2 = dRe*dRe + dIm*dIm;
  double dxx = 57.29578*dIm / dMag_2;
  double dyy = -57.29578*dRe / dMag_2;

  double dPhase = 57.29578*atan(dIm/dRe);
  if (dRe < 0.)
  {
    if (dIm < 0.)
      dPhase -= 180.0;
    else
      dPhase += 180.0;
  }

  unsigned l1 = re.length;
  unsigned l2 = im.length;

  if (l1+l2 == 0)
  {
    AAF Temp(dPhase);
    return (Temp);
  }

  if (l1 == 0)
  {
    // todo
  }

  if (l2 == 0)
  {
    // todo
  }

  // Create our resulting AAF
  AAF Temp(dPhase);

  unsigned * id1 = re.indexes;
  unsigned * id2 = im.indexes;

  double * va1 = re.deviations;
  double * va2 = im.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  if (l1+l2)
    Temp.indexes = new unsigned [l1+l2 + 1]; // the indexes of the result
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin - idtemp + 1;

  if (ltemp)
    Temp.deviations = new double [ltemp];
  double * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting AAF

  for (unsigned i = 0; i < ltemp-1; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = dyy*va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = dxx*va1[a];  // va1[a]+0
      pu1++;
      continue;
    }

    vatempg[i] = dxx*va1[a] + dyy*va2[b];
    pu1++;
    pu2++;
  }

  // todo
  Temp.indexes[ltemp-1] = Temp.inclast();
  Temp.deviations[ltemp-1] = 0.0;

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif

  //  std::cout << "Radius: " << Temp.radius << std::endl;

  return Temp;
}


/************************************************************
 * Method:        mag
 * Author & Date: Darius Grabowski - 06/2007
 * Description:
 *   computes the magnitude of (re|im)
 *   This is the first order taylor approximation of the mag
 *
 *   Input  : const AAF & : AAF real part
 *            const AAF & : AAF imaginary part
 *   Output : AAF         : AAF result
 ************************************************************/
AAF mag(const AAF & re, const AAF & im)
{
  double dRe = re.getcenter();
  double dIm = im.getcenter();

  double dMag = sqrt(dRe*dRe + dIm*dIm);
  double dxx = dRe / dMag;
  double dyy = dIm / dMag;

  unsigned l1 = re.length;
  unsigned l2 = im.length;

  if (l1+l2 == 0)
  {
    AAF Temp(dMag);
    return (Temp);
  }

  if (l1 == 0)
  {
    // todo
  }

  if (l2 == 0)
  {
    // todo
  }

  // Create the resulting AAF
  AAF Temp(dMag);

  unsigned * id1 = re.indexes;
  unsigned * id2 = im.indexes;

  double * va1 = re.deviations;
  double * va2 = im.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  if (l1+l2)
    Temp.indexes = new unsigned [l1+l2 + 1]; // the indexes of the result
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin - idtemp + 1;

  if (ltemp)
    Temp.deviations = new double [ltemp];
  double * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting AAF

  for (unsigned i = 0; i < ltemp-1; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = dyy*va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = dxx*va1[a];  // va1[a]+0
      pu1++;
      continue;
    }

    vatempg[i] = dxx*va1[a] + dyy*va2[b];
    pu1++;
    pu2++;
  }

  // todo
  Temp.indexes[ltemp-1] = Temp.inclast();
  Temp.deviations[ltemp-1] = 0.0;

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif

  //  std::cout << "Radius: " << Temp.radius << std::endl;

  return Temp;
}


/************************************************************
 * Method:        magdb
 * Author & Date: Darius Grabowski - 06/2007
 * Description:
 *   computes the magnitude (dB) of (re|im)
 *   This is the first order taylor approximation of the log mag
 *
 *   Input  : const AAF & : AAF real part
 *            const AAF & : AAF imaginary part
 *   Output : AAF         : AAF result
 ************************************************************/
AAF magdb(AAF & re, AAF & im)
{
  double dRe = re.getcenter();
  double dIm = im.getcenter();

  double dMag = dRe*dRe + dIm*dIm;
  double dxx = 2*4.342945*dRe / dMag;
  double dyy = 2*4.342945*dIm / dMag;

  unsigned l1 = re.length;
  unsigned l2 = im.length;

  if (l1+l2 == 0)
  {
    AAF Temp(2*4.342945*log(dMag));
    return (Temp);
  }

  if (l1 == 0)
  {
    // todo
  }

  if (l2 == 0)
  {
    // todo
  }

  // Create the resulting AAF
  AAF Temp(2*4.342945*log(dMag));

  unsigned * id1 = re.indexes;
  unsigned * id2 = im.indexes;

  double * va1 = re.deviations;
  double * va2 = im.deviations;

  unsigned * pu1 = id1;
  unsigned * pu2 = id2;

  if (l1+l2)
    Temp.indexes = new unsigned [l1+l2 + 1]; // the indexes of the result
  unsigned * idtemp = Temp.indexes;

  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp = fin - idtemp + 1;

  if (ltemp)
    Temp.deviations = new double [ltemp];
  double * vatempg = Temp.deviations;

  Temp.length = ltemp;
  Temp.size = ltemp;

  // Fill the deviations array
  // of the resulting AAF

  for (unsigned i = 0; i < ltemp-1; i++)
  {
    unsigned a = pu1-id1;
    unsigned b = pu2-id2;

    if (a == l1 || id1[a] != idtemp[i])
    {
      vatempg[i] = dyy*va2[b];  // va2[b]+0
      pu2++;
      continue;
    }

    if (b == l2 || id2[b] != idtemp[i])
    {
      vatempg[i] = dxx*va1[a];  // va1[a]+0
      pu1++;
      continue;
    }

    vatempg[i] = dxx*va1[a] + dyy*va2[b];
    pu1++;
    pu2++;
  }

  // todo
  Temp.indexes[ltemp-1] = Temp.inclast();
  Temp.deviations[ltemp-1] = 0.0;

#ifdef FAST_RAD
  Temp.radius = 0.0;
  for (unsigned i = 0; i < ltemp; i++)
    Temp.radius += fabs(vatempg[i]);
#endif

  return Temp;
}

#endif
