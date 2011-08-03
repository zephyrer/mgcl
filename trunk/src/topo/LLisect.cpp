/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/LLisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGLLisect Class.
//MGLLisect is to represent two loops intersection point of
//a parent face parameter space.
//Holds two MGLEPoint data of intersection points.

///////Constructor////////

MGLLisect::MGLLisect(){;}

MGLLisect::MGLLisect(
	const MGPosition& uv,	//Intersection point data.
	const MGLEPoint& lp1,	//First loop's LPoint data.
	const MGLEPoint& lp2)	//Second loop's LPoint data.
	:m_uv(uv),m_is1(lp1),m_is2(lp2){;}

///////Operator oveload///////

bool MGLLisect::operator< (const MGLLisect& li2)const{
	if(m_is1.edge_num()<li2.m_is1.edge_num()) return true;
	if(m_is1.edge_num()>li2.m_is1.edge_num()) return false;
	return m_is1.param()<li2.m_is1.param();
}

bool MGLLisect::operator== (const MGLLisect& li2)const{
	if(m_is1 != li2.m_is1) return false;
	if(m_is2 != li2.m_is2) return false;
	return m_uv==li2.m_uv;
}

///////Member function///////

//Compute distance square of two isect.
double  MGLLisect::distance_square(const MGLLisect& is2) const{
	MGVector dif=m_uv-is2.m_uv;
	return dif%dif;
}

//Debug Function
std::ostream& operator<< (std::ostream& ostrm, const MGLLisect& lli){
	ostrm<<"MGLLisect::m_uv="<<lli.m_uv<<", m_is1="<<lli.m_is1
		<<", m_is2="<<lli.m_is2;
	return ostrm;
}
