/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/FPline.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;

#endif
// MGFPline.cpp
// MGFPline の実装ファイル
//MGFPline is to represent an parameter (u,v) line of a face.
//(MGFace* f, MGCurve uvline) where f is a face pointer, and uvline is
//a parameter (u,v) line of the face f.
//
//MGFPline is used to express Shell's intersection lines.
//The behavior of MGFPline is like an auto_ptr. Copy or assignment
//of MGFPline means transfer of the ownership of all the included curves
//to copied or assigned MGFPline and original MGFPline does not have the
//ownership of the curves any more. Users should be aware of it.

//
// Constructor
//

// Copy Constructor;
MGFPline::MGFPline(const MGFPline& fpl)
:m_face(fpl.m_face){
	m_uvline=fpl.m_uvline;
	fpl.m_uvline=0;
}

////////// Destructor //////////
MGFPline::~MGFPline(){delete m_uvline;};

// Operator overload
//
//Assignment
MGFPline& MGFPline::operator= (const MGFPline& fpl){
	m_face=fpl.m_face;
	delete m_uvline;
	m_uvline=fpl.m_uvline;
	fpl.m_uvline=0;
	return *this;
}

bool MGFPline::operator< (const MGFPline& fpl2)const{
	return (*m_face)<(*(fpl2.m_face));
}

bool MGFPline::operator== (const MGFPline& fpl2)const{
	if(m_face!=fpl2.m_face) return false;
	return ((*m_uvline)== *(fpl2.m_uvline));
}

//
// メンバ関数
//

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
void MGFPline::change_range(
	double t0,		//Parameter value for the start of original. 
	double t1	//Parameter value for the end of original. 
){
	if(m_uvline) m_uvline->change_range(t0,t1);
}

//Release each curve pointer from this.
//After the use of release_line(), MGFPline does not have the ownership of
//the each curve.
MGCurve* MGFPline::release_line(){MGCurve* a=m_uvline; m_uvline=0; return a;}

//Reverse the direction of this line.
void MGFPline::reverse_direction(){
	if(m_uvline) m_uvline->negate();
}

//
// デバッグ関数
//
std::ostream& operator<< (std::ostream& out, const MGFPline& fpl){
//	out.setf ( ios::scientific, ios::floatfield );
//	out.precision ( 10 );
	out <<"MGFPline::m_face="<<(fpl.m_face)<<", m_uvline="<<(fpl.m_uvline);
	if(fpl.m_uvline) out<<(*(fpl.m_uvline));
	return out;
}
