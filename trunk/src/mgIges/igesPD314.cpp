/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD314.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgGL/Color.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD314.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD314 is the class for Iges parameter data type 314(Color Definition).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD314.
MGIgesPD314::MGIgesPD314(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(COLOR_DEFINITION,DEpointer){
}

//! Constructs an object of class MGIgesPD314.
MGIgesPD314::MGIgesPD314(const MGColor& color)
:MGIgesPD(COLOR_DEFINITION,0){
	float rgba[4];
	color.get_color(rgba);
	for(size_t i=0; i<3; i++)
		m_rgb[i]=rgba[i]*float(100.);//m_rgb[] are percentage.
	m_color_name="Color";
}

//Read in parameter data from string stream data.
void MGIgesPD314::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_real(pDelimeter,pdstream,m_rgb[0]);
	get_real(pDelimeter,pdstream,m_rgb[1]);
	get_real(pDelimeter,pdstream,m_rgb[2]);
	get_Hollerith_string(pDelimeter,pdstream,m_color_name);//color name.
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD314::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_real(m_rgb[0],gsec,plines);
	put_real(m_rgb[1],gsec,plines);
	put_real(m_rgb[2],gsec,plines);
	put_Hollerith_string(m_color_name,gsec,plines);
}

//Convert de(type=314: Color definition entity) to MGColor.
//Returned is a newed MGColor.
MGColor* MGIgesIfstream::convert_color(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD314* pd314=static_cast<const MGIgesPD314*>(pd.get());
	const float* rgb=pd314->m_rgb;
	float rgbOne[3];
	for(size_t i=0; i<3; i++){
		rgbOne[i]=float(rgb[i]/100.);
		if(rgbOne[i]>1.)
			rgbOne[i]=1.;
	}
	return new MGColor(rgbOne[0],rgbOne[1],rgbOne[2]);
}
