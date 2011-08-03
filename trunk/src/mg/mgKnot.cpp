/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Knot.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implement Class Knot.

//Constructor
MGKnot::MGKnot(	double t,			//Knot value
				int multiplicity)	//Multiplicity
		:m_value(t),m_multiplicity(multiplicity)
{
	assert(multiplicity>0);
}

bool MGKnot::operator== (const MGKnot& kt2)const{
	if(m_multiplicity!=kt2.m_multiplicity) return false;
	return MGREqual(m_value,kt2.m_value);
}
