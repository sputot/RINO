/*
 * aa_exceptions.h -- Affine Arithmetic Exceptions Library
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
 * Class: AAF_Exception
 * This class includes exceptions for operations on affine forms 
 *
 * @author Darius Grabowski <darius.grabowski@ims.uni-hannover.de>
 * @date 2009-04-01
 * @version 1
 *
 * @see http://www.ims.uni-hannover.de/
 *
 * Last changed by Anna Krause, 2010-03-19
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


#ifndef AAF_EXCEPTIONS_H
#define AAF_EXCEPTIONS_H

#define AAF_ERROR_EXCEPTION     0
#define AAF_DIVZERO_EXCEPTION   1
#define AAF_NEGROOT_EXCEPTION   2
#define AAF_NEGLOG_EXCEPTION    3
#define AAF_NEGBASE_EXCEPTION    4

#include <string>
#include <exception>

class AAF_Exception: public std::exception
{
public:
  std::string message;

  int errorCode;

  double lValue;

  double hValue;

public:

  AAF_Exception(const int, const char *, const double);

  AAF_Exception(const int, const char *, const double, const double);

  AAF_Exception(const int, const char *);

  ~AAF_Exception() throw();

  const char* what() const throw ();

  void report() const throw ();
};

#endif
