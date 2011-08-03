/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Vector.h"
#include "mg/Curve.h"
#include "mg/CParam_list.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCParam_list defines linked list of parameter values of a curve.
// Used to represent list of curve's parameter.

// Constructor
MGCParam_list::MGCParam_list(const MGCurve* c):m_curve(c)
{
	double tol=(MGTolerance::rc_zero())*2.;
	double er=0.;
	if(c) er=(c->param_e() - c->param_s())*tol;
	m_error=er;	//set error.
}

//Copy Constructor.

// Destructor.

// Operator overload.

//Assignment.

// Member Function.

void MGCParam_list::append(double param){
// Adds the param to the end of the list.
	Citerator itr;
	for(itr=begin(); itr!=end(); itr++){
		if(fabs((*itr)-param)<=m_error)
			return;
	}
	push_back(param);
}

void MGCParam_list::append(const MGCParam_list& list)
// Adds the parameter list to the end of the list.
{
	const_Citerator i;
	for(i=list.begin(); i!=list.end(); i++)
		append(*i);
}

double MGCParam_list::removeAt(Citerator i){
//Remove the param and return the param. If i is no valid, 
// behavior is undefined.
	double param=*i;
	erase(i);
	return param;
}

double MGCParam_list::removeFirst(){
//Remove the first param int the list and return the param.
//If i is not valid, behavior is undefined.
	double param=front();
	pop_front();
	return param;
}

double MGCParam_list::removeLast(){
//Remove the last param in the list and return the param.
//If i is not valid, behavior is undefined.
	double param=back();
	pop_back();
	return param;
}
