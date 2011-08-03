/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// MGPickObject.cpp : Implements MGPickObject
//
/////////////////////////////////////////////////////////////////////////////

#include "MGCLStdAfx.h"
#include "mg/PickObject.h"
#include "mg/Group.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//constructor.
MGPickObject::MGPickObject(const MGGelPosition& gp2){
	MGGelPosition* gp=const_cast<MGGelPosition*>(&gp2);
	set_gel(gp->gel());
	set_group(gp->group());
}
//Generate a newed clone object.
MGPickObject* MGPickObject::clone()const{
	return new MGPickObject(*this);
}

//Assignment operator.
MGPickObject& MGPickObject::operator=(const MGPickObject& pobj){
	MGGelPosition::operator =(pobj);
	m_parameter=pobj.m_parameter;
	return *this;
}
