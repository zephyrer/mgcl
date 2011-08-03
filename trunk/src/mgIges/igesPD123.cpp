/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD123.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD123.h"
#include "mg/Vector.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD123 is the class for Iges parameter data type 123(DIRECTION).

// Constructors.

//! Constructs an object of class MGIgesPD123.
MGIgesPD123::MGIgesPD123(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(DIRECTION,DEpointer){
	for(int i=0; i<3; i++)
		m_xyz[i]=0.;
}

//! Constructs an object of class MGIgesPD123.
MGIgesPD123::MGIgesPD123(const MGVector& vec)
:MGIgesPD(DIRECTION){
	for(int i=0; i<3; i++)
		m_xyz[i]=vec[i];
}

//Read in parameter data from string stream data.
void MGIgesPD123::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	MGIges::get_real(pDelimeter,pdstream,m_xyz[0]);
	MGIges::get_real(pDelimeter,pdstream,m_xyz[1]);
	MGIges::get_real(pDelimeter,pdstream,m_xyz[2]);
}

//Convert the direction data to MGVector direction.
void MGIgesPD123::convert_to_vector(MGVector& direction)const{
	direction=MGVector(3,m_xyz);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD123::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	MGIges::put_real(m_xyz[0],gsec,plines);
	MGIges::put_real(m_xyz[1],gsec,plines);
	MGIges::put_real(m_xyz[2],gsec,plines);
}
