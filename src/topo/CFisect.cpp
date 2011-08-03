/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "topo/Face.h"
#include "topo/CFisect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCFisect Class.
//MGCFisect is to represent an intersection of a face and a curve.
//(MGCSisect csi, MGFace* f) where csi consists of world point, curve parameter,
//and face(surface) parameter, and f is a face pointer.

///////Constructor////////

//Construct from all the necessary data.
MGCFisect::MGCFisect(
	const MGPosition& point,	//World coordinate point data of the isect.
	const double& t,			//curve parameter value of the isect.
	const MGPosition& uv,		//Face(Surface) parameter value of the isect.
	const MGFace& face)			//face.
	:m_csi(point,t,uv), m_face(&face){;}

///////Operator oveload///////

bool MGCFisect::operator< (const MGCFisect& fp)const{
	return m_csi.param_curve()<fp.m_csi.param_curve();
}

//Ordering functions.
bool MGCFisect::operator< (const MGisect& is)const{
	const MGCFisect* cis=dynamic_cast<const MGCFisect*>(&is);
	if(cis) return operator<(*cis);
	return is>(*this);
}

bool MGCFisect::operator== (const MGCFisect& fp)const{
	if(m_face!=fp.m_face) return false;
	return m_csi==fp.m_csi;
}

bool MGCFisect::operator== (const MGisect& is)const{
	const MGCFisect* cis=dynamic_cast<const MGCFisect*>(&is);
	if(!cis) return false;
	return operator==(*cis);
}

///////Member function///////

// Output virtual function.
std::ostream& MGCFisect::out(std::ostream& ostrm)const{
	ostrm<<"MGCFisect::m_face="<<m_face<<", m_csi="<<m_csi;
	return ostrm;
}
