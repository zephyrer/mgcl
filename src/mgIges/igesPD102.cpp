/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD102.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/CompositeCurve.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD102.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD102 is the class for Iges parameter data type 102(composite curve).

// Constructors.

//! Constructs an object of class MGIgesPD102.
MGIgesPD102::MGIgesPD102(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(COMPOSITE_CURVE,DEpointer){
}

//Read in parameter data from string stream data.
void MGIgesPD102::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int num_curves;
	get_integer(pDelimeter,pdstream,num_curves);
	m_curve_DEs.resize(num_curves);
	for(int i=0; i<num_curves; i++){
		get_DEpointer(pDelimeter,pdstream,m_curve_DEs[i]);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD102::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	int num_curves=m_curve_DEs.size();
	put_integer(num_curves,gsec,plines);
	for(int i=0; i<num_curves; i++)
		put_DEpointer(m_curve_DEs[i],gsec,plines);
}

//Convert de to MGCompositeCurve(a newed object). de must be of type 102
//(Composite curve).
MGCompositeCurve* MGIgesIfstream::convert_composite(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD102* pd102=static_cast<const MGIgesPD102*>(pd.get());

	MGCompositeCurve* compcrv=new MGCompositeCurve;
	int n=pd102->m_curve_DEs.size();
	for(int i=0; i<n; i++){
		int dei=pd102->m_curve_DEs[i];
		MGGel* obj=convert_to_gel(dei);
		if(obj){
			MGCurve* crv=static_cast<MGCurve*>(obj);//obj must be MGCurve.
			compcrv->connect(crv);
		}
	}

	return compcrv;
}
