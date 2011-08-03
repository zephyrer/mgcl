/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD116.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD116.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD116 is the class for Iges parameter data type 116(POINT).

// Constructors.

//! Constructs an object of class MGIgesPD116.
MGIgesPD116::MGIgesPD116(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(MGIges::POINT,DEpointer), m_display_symbolDE(0){
}

//! Constructs an object of class MGIgesPD116.
MGIgesPD116::MGIgesPD116(const MGPosition& P,int display_symbolDE)
:MGIgesPD(MGIges::POINT), m_display_symbolDE(display_symbolDE){
	for(int i=0; i<3; i++)
		m_coordinates[i]=P[i];
}
MGIgesPD116::MGIgesPD116(const MGPoint& P,int display_symbolDE)
:MGIgesPD(MGIges::POINT), m_display_symbolDE(display_symbolDE){
	for(int i=0; i<3; i++)
		m_coordinates[i]=P[i];
}

//! Constructs an object of class MGIgesPD116.
MGIgesPD116::MGIgesPD116(const double coordinates[3],int display_symbolDE)
:MGIgesPD(MGIges::POINT), m_display_symbolDE(display_symbolDE){
	for(int i=0; i<3; i++)
		m_coordinates[i]=coordinates[i];
}

	//Convert the point data to MGPosition position.
void MGIgesPD116::convert_to_position(MGPosition& position)const{
	position.resize(3);
	position(0)=m_coordinates[0];
	position(1)=m_coordinates[1];
	position(2)=m_coordinates[2];
}

//Read in parameter data from string stream data.
void MGIgesPD116::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_real(pDelimeter,pdstream,m_coordinates[0]);
	get_real(pDelimeter,pdstream,m_coordinates[1]);
	get_real(pDelimeter,pdstream,m_coordinates[2]);
	get_DEpointer(pDelimeter,pdstream,m_display_symbolDE);
}

//PD116 write_out_intostring
//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD116::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_real(m_coordinates[0],gsec,plines);
	put_real(m_coordinates[1],gsec,plines);
	put_real(m_coordinates[2],gsec,plines);
	put_DEpointer(m_display_symbolDE,gsec,plines);
}

//Convert de(type=116: point) to MGPoint.
//Returned is a newed object.
MGPoint* MGIgesIfstream::convert_point(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD116* pd116=static_cast<const MGIgesPD116*>(pd.get());
	return new MGPoint(MGPosition(3,pd116->m_coordinates));
}
