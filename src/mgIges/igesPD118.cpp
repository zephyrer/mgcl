/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD118.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD118.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD118 is the class for Iges parameter data type 118(Ruled Surface).

// Constructors.

//! Constructs an object of class MGIgesPD118.
MGIgesPD118::MGIgesPD118(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(RULED_SURFACE,DEpointer), m_1st_Curve_DE(0),m_2nd_Curve_DE(0),
m_direction_flag(0), m_developable_flag(0){
}

//Read in parameter data from string stream data.
void MGIgesPD118::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_1st_Curve_DE);
	get_DEpointer(pDelimeter,pdstream,m_2nd_Curve_DE);
	get_integer(pDelimeter,pdstream,m_direction_flag);
	get_integer(pDelimeter,pdstream,m_developable_flag);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD118::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_1st_Curve_DE,gsec,plines);
	put_DEpointer(m_2nd_Curve_DE,gsec,plines);
	put_integer(m_direction_flag,gsec,plines);
	put_integer(m_developable_flag,gsec,plines);
}

//Convert de(type=118: Ruled surface) to MGSurface.
//Returned is a newed object.
MGSurface* MGIgesIfstream::convert_ruled_surface(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD118* pd118=static_cast<const MGIgesPD118*>(pd.get());
	int de1=pd118->m_1st_Curve_DE, de2=pd118->m_2nd_Curve_DE;
	std::auto_ptr<MGGel> obj1=std::auto_ptr<MGGel>(convert_to_gel(de1));
	std::auto_ptr<MGGel> obj2=std::auto_ptr<MGGel>(convert_to_gel(de2));
	MGPoint* P1=dynamic_cast<MGPoint*>(obj1.get());
	MGPoint* P2=dynamic_cast<MGPoint*>(obj2.get());
	if(P1 && P2)
		return 0;
	
	std::auto_ptr<MGCurve> crv1; std::auto_ptr<MGCurve> crv2;
	if(P1){
		const MGPosition& pos1=P1->position();
		crv1=std::auto_ptr<MGCurve>(new MGStraight(pos1,pos1));
	}else{
		MGCurve* crv1p=dynamic_cast<MGCurve*>(obj1.get());
		if(crv1p)
			crv1=std::auto_ptr<MGCurve>(dynamic_cast<MGCurve*>(obj1.release()));
		else
			return 0;
	}
	if(P2){
		const MGPosition& pos2=P2->position();
		crv1=std::auto_ptr<MGCurve>(new MGStraight(pos2,pos2));
	}else{
		MGCurve* crv2p=dynamic_cast<MGCurve*>(obj2.get());
		if(crv2p)
			crv2=std::auto_ptr<MGCurve>(dynamic_cast<MGCurve*>(obj2.release()));
		else
			return 0;
	}
	if(pd118->m_direction_flag)
		crv2->negate();
	std::auto_ptr<MGSurface> surf=MGCL::create_ruled_surface(*crv1, *crv2);
	surf->change_range(true,0.,1.);
	surf->change_range(false,0.,1.);
	return surf.release();
}