/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD100.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD100.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD100 is the class for Iges parameter data type 100(circular arc).

// Constructors.

//! Constructs an object of class MGIgesPD100.
MGIgesPD100::MGIgesPD100(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(CIRCULAR_ARC,DEpointer), m_Zt(0.){
	m_center[0]=m_center[1]=0.;//(x1, y1) of the center.
	m_start[0]=m_start[1]=0.;//(x2, y2) of the start point.
	m_terminate[0]=m_terminate[1]=0.;//(x3, y3) of the terminate point.
}

//Construct PD100, supplying 2D coordinate data in each array.
MGIgesPD100::MGIgesPD100(
	const double center[2], const double start[2], const double terminate[2],
	double Z	//Z coord
):MGIgesPD(CIRCULAR_ARC), m_Zt(Z){
	m_center[0]=center[0];m_center[1]=center[1];//(x1, y1) of the center.
	m_start[0]=start[0];m_start[1]=start[1];//(x2, y2) of the start point.
	m_terminate[0]=terminate[0];m_terminate[1]=terminate[1];//(x3, y3) of the terminate point.
}

//Read in parameter data from string stream data.
void MGIgesPD100::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_real(pDelimeter,pdstream,m_Zt);
	get_real(pDelimeter,pdstream,m_center[0]);
	get_real(pDelimeter,pdstream,m_center[1]);
	get_real(pDelimeter,pdstream,m_start[0]);
	get_real(pDelimeter,pdstream,m_start[1]);
	get_real(pDelimeter,pdstream,m_terminate[0]);
	get_real(pDelimeter,pdstream,m_terminate[1]);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD100::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_real(m_Zt,gsec,plines);
	put_real(m_center[0],gsec,plines);
	put_real(m_center[1],gsec,plines);
	put_real(m_start[0],gsec,plines);
	put_real(m_start[1],gsec,plines);
	put_real(m_terminate[0],gsec,plines);
	put_real(m_terminate[1],gsec,plines);
}

//Convert de to MGEllipse(a newed object). de must be of type 100(Circular Arc).
MGEllipse* MGIgesIfstream::convert_arc(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD100* pd100=static_cast<const MGIgesPD100*>(pd.get());
	MGPosition start(pd100->m_start[0],pd100->m_start[1],pd100->m_Zt);
	MGPosition terminate(pd100->m_terminate[0],pd100->m_terminate[1],pd100->m_Zt);
	MGPosition center(pd100->m_center[0],pd100->m_center[1],pd100->m_Zt);
	MGVector normal(0.,0.,1.);
	MGVector V2(start,center), V3(terminate,center);
	const MGVector& Xaxis=MGDefault::x_unit_vector();
	double r=V2.len();//radius.
	double t2=Xaxis.angle2pai(V2,normal), t3=V2.angle2pai(V3,normal);
	return new MGEllipse(center, Xaxis*r,MGDefault::y_unit_vector()*r,MGInterval(t2,t2+t3));
}
