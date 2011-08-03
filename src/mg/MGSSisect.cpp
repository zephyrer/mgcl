/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/SSisect.h"
#include "mg/CCisect_list.h"
#include "mg/Curve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSSisect.cc
// MGSSisect の実装ファイル

//
// Constructor
//

// Copy Constructor;
MGSSisect::MGSSisect(const MGSSisect& ssi)
:m_rel(ssi.m_rel){
	MGSSisect* ssip=const_cast<MGSSisect*>(&ssi);
	m_iline=ssip->m_iline; m_param1=ssip->m_param1; m_param2=ssip->m_param2;
	ssip->m_iline=ssip->m_param1=ssip->m_param2=0;
}

//Construct providing all the raw data.
//Copy version. Copy of the three curves will take place.
MGSSisect::MGSSisect (
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION rel)
:m_rel(rel){
	m_iline=iline.clone();
	m_param1=param1.clone();
	m_param2=param2.clone();
}

//Destructor
//
MGSSisect::~MGSSisect(){
	delete m_iline; delete m_param1; delete m_param2;
}

// Operator overload
//
//Assignment
MGSSisect& MGSSisect::operator= (const MGSSisect& ssi){
	delete m_iline; delete m_param1; delete m_param2;
	m_rel=ssi.m_rel;
	MGSSisect* ssip=const_cast<MGSSisect*>(&ssi);
	m_iline=ssip->m_iline; m_param1=ssip->m_param1; m_param2=ssip->m_param2;
	ssip->m_iline=ssip->m_param1=ssip->m_param2=0;
	return *this;
}

bool MGSSisect::operator< (const MGSSisect& ssi2)const{
	return m_param1->start_point()<ssi2.m_param1->start_point();
}

bool MGSSisect::operator== (const MGSSisect& ssi2)const{
	if((*m_param1)!= *(ssi2.m_param1)) return false;
	if((*m_param2)!= *(ssi2.m_param2)) return false;
	return true;
}

bool MGSSisect::operator< (const MGisect& is)const{
	const MGSSisect* cis=dynamic_cast<const MGSSisect*>(&is);
	if(cis) return operator<(*cis);
	return is>(*this);
}
bool MGSSisect::operator== (const MGisect& is)const{
	const MGSSisect* cis=dynamic_cast<const MGSSisect*>(&is);
	if(!cis) return false;
	return operator==(*cis);
}

//
// メンバ関数
//

//Test if two ssi's world curve have common parts (in line_zero()).
//Fucntion's return value is 
//		1:have commonpart.
//		0>=:no common part(except a point).
int MGSSisect::has_common(const MGSSisect& ssi2)const{
	std::vector<double> dvec;
	MGCCisect_list isects;
	int retcode=line().common(ssi2.line(),dvec,isects);
	if(retcode==1 || retcode==3) return 1;
	return 0;
}

//negate the direction of the intersection line.
void MGSSisect::negate(){
	m_iline->negate();	//(x,y,z)coordinate representaion of the line.
	m_param1->negate();	//(u,v) representaion of the line of the first
	m_param2->negate();	//(u,v) representaion of the line of the second
}

// Output virtual function.
std::ostream& MGSSisect::out(std::ostream& ostrm)const{
//	ostrm.setf ( ios::scientific, ios::floatfield );
//	ostrm.precision ( 10 );
	if(m_iline){
		ostrm<<"MGSSisect::m_iline="<<(*m_iline);
		if(m_param1) ostrm<<", m_param1="<<(*m_param1);
		if(m_param2) ostrm<<", m_param2="<<(*m_param2);
		ostrm<<", m_rel="   <<(m_rel);
	}else{
		ostrm<<"MGSSisect::m_iline="<<m_iline;
	}
	return ostrm;
}

//Release each curve pointer from this.
//After the use of release_xxxx(), MGSSisect does not have the ownership of
//the each curve.
MGCurve* MGSSisect::release_line(){MGCurve* a=m_iline; m_iline=0; return a;}
MGCurve* MGSSisect::release_param1(){MGCurve* a=m_param1; m_param1=0; return a;}
MGCurve* MGSSisect::release_param2(){MGCurve* a=m_param2; m_param2=0; return a;}

//Replace 1st and 2nd order of the parameter line representation.
MGSSisect& MGSSisect::replace12(){
	MGCurve* save=m_param1;
	m_param1=m_param2;
	m_param2=save;
	return *this;
}

void MGSSisect::set_null(){
	if(!m_iline)
		return;

	delete m_iline; delete m_param1; delete m_param2;
	m_iline=m_param1=m_param2=0;
}
