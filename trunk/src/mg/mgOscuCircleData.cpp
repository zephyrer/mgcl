/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/OscuCircleData.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implement Class MGOscuCircleData.

//Constructor
MGOscuCircleData::MGOscuCircleData(
	size_t index,	//Index
	double radius)	//Rdius
:m_index(index),m_radius(radius){
	assert(radius>0.);
}

bool MGOscuCircleData::operator< (const MGOscuCircleData& ocd2)const{
	if(m_index==ocd2.m_index) return m_radius<ocd2.m_radius;
	return m_index<ocd2.m_index;
}

bool MGOscuCircleData::operator== (const MGOscuCircleData& ocd2)const{
	if(m_index!=ocd2.m_index) return false;
	return MGREqual(m_radius, ocd2.m_radius);
}
