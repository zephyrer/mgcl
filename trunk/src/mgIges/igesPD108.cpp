/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD108.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD108.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD108 is the class for Iges parameter data type 108(PLANE).

// Constructors.

//! Constructs an object of class MGIgesPD108.
MGIgesPD108::MGIgesPD108(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(PLANE,DEpointer){
}

//! Constructs an object of class MGIgesPD108 from a MGPlane
MGIgesPD108::MGIgesPD108(const MGPlane& plane)
:MGIgesPD(PLANE),m_boundCurve_DE(0),m_symbol_size(1.0){
	plane.abcd(m_coef);
	const MGPosition& rootP=plane.root_point();
	for(int i=0; i<3; i++)
		m_ref_point[i]=rootP[i];//Reference point on the plane(at which symbol be displayed).
}

//Read in parameter data from string stream data.
void MGIgesPD108::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	for(int i=0; i<4; i++)
		get_real(pDelimeter,pdstream,m_coef[i]);
	get_DEpointer(pDelimeter,pdstream,m_boundCurve_DE);
	get_real(pDelimeter,pdstream,m_ref_point[0]);
	get_real(pDelimeter,pdstream,m_ref_point[1]);
	get_real(pDelimeter,pdstream,m_ref_point[2]);
	get_real(pDelimeter,pdstream,m_symbol_size);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD108::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	for(int i=0; i<4; i++)
		put_real(m_coef[i],gsec,plines);
	put_DEpointer(m_boundCurve_DE,gsec,plines);
	put_real(m_ref_point[0],gsec,plines);
	put_real(m_ref_point[1],gsec,plines);
	put_real(m_ref_point[2],gsec,plines);
	put_real(m_symbol_size,gsec,plines);
}

//Convert de to MGObject(a newed object). de must be of type 108(PLANE).
//Output MGObject is a MGPlane(unbounded infinite plane) or a MGFace(bounded plane).
//MGCL does not treat unbounded infinite plane with a hole. de of FomNumber -1 will be
//converted to a MGPlane(unbounded infinite plane).
MGObject* MGIgesIfstream::convert_plane(
	const MGIgesDirectoryEntry& de
)const{
	int fnum=de.FormNumber();
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD108* pd108=static_cast<const MGIgesPD108*>(pd.get());
	MGPlane* pl=new MGPlane(pd108->m_coef,pd108->m_ref_point);
	MGObject* obj=pl;
	if(fnum==1){//When positive bounded planar portion.
		MGGel* crvobj=convert_to_gel(pd108->m_boundCurve_DE);
		MGCurve* crv=dynamic_cast<MGCurve*>(crvobj);
		if(crv){
			MGFace* f=new MGFace(pl);
			std::auto_ptr<MGLoop> loop=f->build_loop(*crv);
			if(loop->area()<0.)
				loop->negate();
			f->prepend_boundary(loop.release());
			obj=f;
		}
		delete crvobj;
	}
	return obj;
}
