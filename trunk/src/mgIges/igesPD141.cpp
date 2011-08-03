/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD141.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD141.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD141 is the class for Iges parameter data type 141(BOUNDARY entity).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD141.
MGIgesPD141::MGIgesPD141(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(BOUNDARY,DEpointer),m_type(0),m_prefered(2),
m_surface_DE(0){
}

//Read in parameter data from string stream data.
void MGIgesPD141::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_type);
	get_integer(pDelimeter,pdstream,m_prefered);
	get_DEpointer(pDelimeter,pdstream,m_surface_DE);

	int n;
	get_integer(pDelimeter,pdstream,n);//Number of curves.
	for(int i=0; i<n; i++){
		int curve_DE;
		get_DEpointer(pDelimeter,pdstream,curve_DE);
		int sense;
		get_integer(pDelimeter,pdstream,sense);

		MGIges141Edge edgei(curve_DE,sense);
		int npcrv;
		get_integer(pDelimeter,pdstream,npcrv);
		for(int j=0; j<npcrv; j++){
			int pcurve_DE;
			get_DEpointer(pDelimeter,pdstream,pcurve_DE);
			edgei.push_back_pcurve(pcurve_DE);
		}
		m_edges.push_back(edgei);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD141::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_integer(m_type,gsec,plines);
	put_integer(m_prefered,gsec,plines);
	put_DEpointer(m_surface_DE,gsec,plines);

	int n=m_edges.size();
	put_integer(n,gsec,plines);
	for(int i=0; i<n; i++){
		MGIges141Edge edgei=m_edges[i];
		put_DEpointer(edgei.m_pcurves[0],gsec,plines);
		put_integer(edgei.m_sense,gsec,plines);
		int m=edgei.m_pcurves.size();
		for(int j=0; j<m; j++){
			put_DEpointer(edgei.m_pcurves[j],gsec,plines);
		}
	}
}
