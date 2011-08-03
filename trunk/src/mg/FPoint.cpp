/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/FSurface.h"
#include "mg/FPoint.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGFPoint Class.
//MGFPoint is to represent a Face's point. The expression is:
//(MGFace* f, MGPosition uv), where f is a face pointer of interest, 
//and uv is the parameter value of the face f.

///////Constructor////////

///////Operator oveload///////

bool MGFPoint::operator< (const MGFPoint& fp)const{
	if(fp.m_face==m_face){
		if(m_uv[0]==fp.m_uv[0]) return (m_uv[1]<fp.m_uv[1]);
		else return (m_uv[0]<fp.m_uv[0]);
	}else{
		return false;
	}
}

//Two Fpoints are equal if they belong to one face and their distance in parameter
//(u1,v2) and (u2,v2) is less than parameter_error() of the face.
bool MGFPoint::operator== (const MGFPoint& fp)const{
	if(!m_face)
		return false;
	if(fp.m_face!=m_face)
		return false;
	if(fabs(m_uv[0]-fp.m_uv[0])>m_face->param_error_u())
		return false;
	if(fabs(m_uv[1]-fp.m_uv[1])>m_face->param_error_v())
		return false;

	return true;
}

///////Member function///////

//Evaluation of the Face at the FPoint.
//When nderi=0, get the positional data at the point.
MGVector MGFPoint::eval(size_t ndu, size_t ndv)const{
	return fsurface().eval(uv(), ndu, ndv);
}

//Debug Function
std::ostream& operator<< (std::ostream& ostrm, const MGFPoint& fp){
	ostrm<<"MGFPoint::m_face="<<fp.m_face<<", m_uv="<<fp.m_uv;
	return ostrm;
}

/////////// GLOBAL FUNCTIONS FOR SHELL ///////////////////

//Evaluate face's (shell's) point from the MGFPoint.
MGVector eval(const MGFPoint& fp,	//Face point.
	  size_t ndu, size_t ndv)		//Order of derivative along u and v direction.
{	return fp.fsurface().eval(fp.uv(), ndu,ndv);}

//Compute normal vector(not unit) at uv.
MGVector normal(const MGFPoint& fp){ return fp.fsurface().normal(fp.uv());}

//Compute unit normal vector at uv.
MGUnit_vector unit_normal(const MGFPoint& fp)
{return fp.fsurface().unit_normal(fp.uv());}
