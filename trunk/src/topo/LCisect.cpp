/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/Tolerance.h"
#include "topo/LCisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Define MGLCisect Class.
//MGLCisect is to represent Loop and curve intersection point of
//a parent face parameter space.
//Holds (lp, t, uv), where lp=loop point, t=curve parameter value, and
//uv=face aprameter value.

///////Constructor////////

MGLCisect::MGLCisect():m_t(0){;}

/*MGLCisect::MGLCisect(size_t i, const MGCCisect isect)
:m_i(i), m_isect(isect){;}*/

MGLCisect::MGLCisect(
	const MGLEPoint& lp,		//loop's parameter with edge id.
	double t,					//Curve's parameter value.
	const MGPosition& uv)		//Face's parameter value(u,v) data.
	:m_lp(lp),m_t(t), m_uv(uv){;}

///////Operator oveload///////

bool MGLCisect::operator< (const MGLCisect& lci2)const{
	if(m_lp<lci2.m_lp) return true;
	if(m_lp>lci2.m_lp) return false;

	return m_t<lci2.m_t;
}

bool MGLCisect::operator== (const MGLCisect& lci2)const{
	if(m_lp!=lci2.m_lp) return false;
	if(!MGREqual(m_t, lci2.m_t)) return false;
	return m_uv==lci2.m_uv;
}

///////Member function///////

//Compute distance square of two isect.
double MGLCisect::distance_square(const MGLCisect& is2) const
{
	MGVector dif=m_uv-is2.m_uv;
	return dif%dif;
}

//Debug Function
std::ostream& operator<< (std::ostream& ostrm, const MGLCisect& lcis)
{
	ostrm<<"MGLCisect::m_lp="<<lcis.m_lp
		<<", m_t="<<lcis.m_t<<", m_uv="<<lcis.m_uv;
	return ostrm;
}
