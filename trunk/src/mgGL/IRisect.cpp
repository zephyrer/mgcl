/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include <iostream>
#include "mgGL/IRisect.h"

//mgIRisect is a proprietry class for the class mgImageRect.
//mgIRisect represents one intesection of an mgImageRect perimeter and u=const(or v=const)
//line. When intersection with u=const, m_t is v value, or vise versa.

//Comparison with mgIRisect.
bool mgIRisect::operator> (const mgIRisect& is)const{
	if(m_perimeter!=is.m_perimeter)
		return m_perimeter>is.m_perimeter;
	if(m_perimeter<=1)
		return m_t>is.m_t;
	return m_t<is.m_t;
}

bool mgIRisect::operator== (const mgIRisect& is)const{
	if(m_perimeter!=is.m_perimeter)
		return false;
	return MGAEqual(m_t,is.m_t);
}

std::ostream& operator<< (std::ostream& out, const mgIRisect& isect){
	out<<"perimeter="<<isect.m_perimeter<<", t="<<isect.m_t<<", tri_id="<<isect.m_tri_id;
	out<<std::endl;
	return out;
}
