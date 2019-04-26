/* 
 * aa_rounding.h -- Platform dependent rounding
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


#ifndef AA_ROUNDING_H
#define AA_ROUNDING_H

#ifdef SOLARIS
#define fesetround(x) x
#define fegetround() 0
#endif

#include <fenv.h>
#define AA_UPWARD FE_UPWARD
#define AA_DOWNWARD FE_DOWNWARD
typedef unsigned int aa_rnd_t;

unsigned int aa_fesetround(aa_rnd_t);
aa_rnd_t aa_fegetround(void);

#endif
