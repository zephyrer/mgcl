/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/FFisect.h"

// MGFFisect.h
// Header for MGFFisect
//
//MGFFisect represents one intersection line of a MGFace with MGFace or MGSurface..
//The behavior of MGFFisect is like a auto_ptr. Copy or assignment
//of MGFFisect means transfer of the ownership of all the included curve
//to copied or assigned MGFFisect and original MGFFisect does not have the
//ownership of the curves any more. User should be aware of it.

////////// Constructor //////////

//Construct providing all the raw data.
//Copy version. Copy of the three curves will take place.
MGFFisect::MGFFisect(
	const MGCurve& iline,
	const MGFace* face1,
	const MGCurve& param1,
	const MGFace* face2,
	const MGCurve& param2
):m_curve(iline.clone()),
m_face1line(face1,param1.clone()),m_face2line(face2,param2.clone()){;
}

//Construct providing all the raw data.
MGFFisect::MGFFisect(
	MGCurve* iline,	//Pointer of a newed object.
	MGFPline face1uv,
	MGFPline face2uv
):m_curve(iline),m_face1line(face1uv),m_face2line(face2uv){;
}

// Copy Constructor;
// ffi's ownership of all the three curves will be released.
MGFFisect::MGFFisect(const MGFFisect& ffi)
:m_curve(ffi.m_curve),m_face1line(ffi.m_face1line),m_face2line(ffi.m_face2line){
	MGFFisect* uffi=const_cast<MGFFisect*>(&ffi);
	uffi->m_curve=0;
}

////////// Operator overload //////////

//Assignment
// ssi's ownership of all the three curves will be released.
MGFFisect& MGFFisect::operator= (const MGFFisect& ffi){
	MGFFisect* uffi=const_cast<MGFFisect*>(&ffi);
	m_curve=uffi->m_curve; uffi->m_curve=0;
	m_face1line=uffi->m_face1line;
	m_face2line=uffi->m_face2line;
	return *this;
}

//Comparison operator.
bool MGFFisect::operator< (const MGFFisect& ffi2)const{
	const MGFSurface* f11=m_face1line.face();
	const MGFSurface* f21=ffi2.m_face1line.face();
	if(f11!=f21) return f11<f21;
	const MGFSurface* f12=m_face2line.face();
	const MGFSurface* f22=ffi2.m_face2line.face();
	if(f12!=f22) return f12<f22;
	double t10=m_face1line.uvline().param_s();
	double t20=ffi2.m_face1line.uvline().param_s();
	if(t10!=t20) return t10<t20;
	return m_face1line.uvline().param_s()<ffi2.m_face1line.uvline().param_s();
}
bool MGFFisect::operator== (const MGFFisect& ffi2)const{
	const MGFSurface* f11=m_face1line.face();
	const MGFSurface* f21=ffi2.m_face1line.face();
	if(f11!=f21) false;
	const MGFSurface* f12=m_face2line.face();
	const MGFSurface* f22=ffi2.m_face2line.face();
	if(f12!=f22) return false;
	double t10=m_face1line.uvline().param_s();
	double t20=ffi2.m_face1line.uvline().param_s();
	if(t10!=t20) return false;
	return m_face1line.uvline().param_s()==ffi2.m_face1line.uvline().param_s();
}
bool MGFFisect::operator< (const MGisect& is)const{
	const MGFFisect* cis=dynamic_cast<const MGFFisect*>(&is);
	if(cis) return operator<(*cis);
	return is>(*this);
}
bool MGFFisect::operator== (const MGisect& is)const{
	const MGFFisect* cis=dynamic_cast<const MGFFisect*>(&is);
	if(!cis) return false;
	return operator==(*cis);
}

////////// Memeber Function //////////

// Output virtual function.
std::ostream& MGFFisect::out(std::ostream& ostrm)const{
//	ostrm.setf ( ios::scientific, ios::floatfield );
//	ostrm.precision ( 10 );
	ostrm <<"MGFFisect::m_curve="<<*m_curve
		<<",m_face1line="<<m_face1line
		<<", m_face2line="<<m_face2line;
	return ostrm;
}

//Replace 1st and 2nd order of the parameter line representation.
void MGFFisect::exchange12(){
	MGFPline face1=m_face1line;
	MGFPline face2=m_face2line;
	m_face2line=face1;
	m_face1line=face2;
}
