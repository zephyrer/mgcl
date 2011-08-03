/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDefaultVector_HH_
#define _MGDefaultVector_HH_

#include "mg/MGCL.h"
#include "mg/Position.h"

//  MGDefault.h
//  header for class MGDefault

// Difines default values of each class.
MGEXTERN const MGPosition mgNULL_Pos;
MGEXTERN const MGPosition mgORIGIN;
MGEXTERN const MGPosition mgORIGIN_2D;

MGEXTERN const MGVector mgX_UVEC;
MGEXTERN const MGVector mgY_UVEC;
MGEXTERN const MGVector mgZ_UVEC;

MGEXTERN const MGVector mgX_UVEC_2D;
MGEXTERN const MGVector mgY_UVEC_2D;

MGEXTERN const MGVector mgNULL_VEC;
MGEXTERN const MGVector mgZERO_VEC;
MGEXTERN const MGVector mgZERO_VEC_2D;

#endif
