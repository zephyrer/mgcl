/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/LPoint.h"
#include "topo/LEPoint.h"
#include "topo/LCisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGLPoint Class.
//MGLPoint is to represent Loop's point. This is represented as
//(i, t), where i is the pcell id(i.e. edge number) of in the loop, 
//and t is the parameter value of the curve of the pcell(edge).

///////Constructor////////

//Conversion constructor from MGLEPoint.
MGLPoint::MGLPoint(const MGLEPoint& le)
:m_i(le.edge_num()), m_t(le.param()){;}

///////Operator oveload///////

bool MGLPoint::operator< (const MGLPoint& lp)const{
	if(m_i<lp.m_i) return true;
	if(m_i>lp.m_i) return false;
	return (m_t<lp.m_t);
}
bool MGLPoint::operator> (const MGLPoint& lp)const{
	return lp<*this;
}
bool MGLPoint::operator<= (const MGLPoint& lp)const{
	return !(*this>lp);
}
bool MGLPoint::operator>= (const MGLPoint& lp)const{
	return !(*this<lp);
}
bool MGLPoint::operator== (const MGLPoint& lp)const{
	return (m_i==lp.m_i && m_t==lp.m_t);
}

///////Member function///////

//Debug Function
std::ostream& operator<< (std::ostream& ostrm, const MGLPoint& lp){
	ostrm<<"MGLPoint::m_i="<<lp.m_i<<", m_t="<<lp.m_t;
	return ostrm;
}
