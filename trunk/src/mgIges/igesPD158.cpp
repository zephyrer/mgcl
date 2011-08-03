/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD158.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Sphere.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD158.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD158 is the class for Iges parameter data type 158(SPHERE).

// Constructors.

//! Constructs an object of class MGIgesPD158.
MGIgesPD158::MGIgesPD158(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(SPHERE,DEpointer){
}

//! Constructs an object of class MGIgesPD158 from a MGPlane
MGIgesPD158::MGIgesPD158(const MGSphere& sphere)
:MGIgesPD(SPHERE){
	m_radius=sphere.radius();
	const MGPosition& center=sphere.C();
	for(int i=0; i<3; i++)
		m_center[i]=center[i];//Center point.
}

//Read in parameter data from string stream data.
void MGIgesPD158::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_real(pDelimeter,pdstream,m_radius);
	get_real(pDelimeter,pdstream,m_center[0]);
	get_real(pDelimeter,pdstream,m_center[1]);
	get_real(pDelimeter,pdstream,m_center[2]);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD158::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_real(m_radius,gsec,plines);
	put_real(m_center[0],gsec,plines);
	put_real(m_center[1],gsec,plines);
	put_real(m_center[2],gsec,plines);
}

//Convert de to MGObject(a newed object). de must be of type 158(SPHERE).
MGSphere* MGIgesIfstream::convert_sphere158(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD158* pd158=static_cast<const MGIgesPD158*>(pd.get());
	MGPosition cntr(3,pd158->m_center);
	
	return new MGSphere(cntr,pd158->m_radius);
}
