/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/KnotArray.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implemetation of MGKnotArray Class.

//Friend Function

//Constructor
MGKnotArray::MGKnotArray(const MGKnot& knot)		//From Knot.
:m_ktArray(1,knot){;}

MGKnotArray::MGKnotArray(double knot, int mult)
//From knot and the multiplicity.
:m_ktArray(1,MGKnot(knot,mult))
{assert(mult>0);}

//	MGKnotArray(const MGKnotArray&);	//Copy Constructor.

//Destructor
//	~MGKnotArray();	

//Member Function
MGKnotArray& MGKnotArray::add(const MGKnot& knot)
	//Add to the end of list.
{
	push_back(knot);
	return *this;
}

MGKnotArray& MGKnotArray::add(double t, int mult)
//Add to the end of list.
{
	MGKnot knot(t, mult);
	return add(knot);
}

//Operator overload.
//	MGKnotArray& operator =(MGKnotArray&);//Assignment operator overload.
