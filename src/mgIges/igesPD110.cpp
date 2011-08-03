/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Implementation for class MGIgesPD110(LINE).
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD110.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD110 is the class for Iges parameter data type 110(LINE).

// Constructors.

//! Constructs an object of class MGIgesPD110.
MGIgesPD110::MGIgesPD110(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(LINE,DEpointer){
}

	//! Constructs an object of class MGIgesPD110.
MGIgesPD110::MGIgesPD110(
	const MGPosition& start,
	const MGPosition& terminate
):MGIgesPD(LINE){
	for(int i=0; i<3; i++){
		m_start[i]=start[i];
		m_terminate[i]=terminate[i];
	}
}

//Read in parameter data from string stream data.
void MGIgesPD110::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_real(pDelimeter,pdstream,m_start[0]);
	get_real(pDelimeter,pdstream,m_start[1]);
	get_real(pDelimeter,pdstream,m_start[2]);

	get_real(pDelimeter,pdstream,m_terminate[0]);
	get_real(pDelimeter,pdstream,m_terminate[1]);
	get_real(pDelimeter,pdstream,m_terminate[2]);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD110::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	char pdel=gsec.paramDelimeter();
	put_real(m_start[0],gsec,plines);
	put_real(m_start[1],gsec,plines);
	put_real(m_start[2],gsec,plines);

	put_real(m_terminate[0],gsec,plines);
	put_real(m_terminate[1],gsec,plines);
	put_real(m_terminate[2],gsec,plines);
}

//Convert de to MGObject(a newed object). de must be of type 110(line).
//Output MGObject is a MGStraight(finite, semifinite, and infinit).
MGStraight* MGIgesIfstream::convert_line(
	const MGIgesDirectoryEntry& de
)const{
	int fnum=de.FormNumber();
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD110* pd110=static_cast<const MGIgesPD110*>(pd.get());
	MGStraight* line;
	MGPosition E(3,pd110->m_terminate), S(3,pd110->m_start);
	if(fnum==0){//When finite line segment.
		line=new MGStraight(E,S);
		line->change_range(0.,1.);
	}else{
		MGSTRAIGHT_TYPE stype;
		if(fnum==1)
			stype=MGSTRAIGHT_HALF_LIMIT;
		else
			stype=MGSTRAIGHT_UNLIMIT;
		MGVector direction(E,S);
		line=new MGStraight(stype,direction,S);
	}
	return line;
}
