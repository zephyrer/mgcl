/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD144.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD142.h"
#include "mgiges/IgesPD144.h"
#include "mgiges/IgesGsec.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD144 is the class for Iges parameter data type 144(Trimmed Surface).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD144.
MGIgesPD144::MGIgesPD144(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(TRIMMED_SURFACE,DEpointer),m_surface_DE(0),m_outer_boundary_type(0),
m_outer_boudary_DE(0){
}

//! Constructs an object of class MGIgesPD144.
MGIgesPD144::MGIgesPD144(
	int surfaceDE,		//Base surface DE.
	int outerboundaryDE	//if =0, no outer boundary.
):MGIgesPD(TRIMMED_SURFACE),m_surface_DE(surfaceDE),m_outer_boundary_type(0),
m_outer_boudary_DE(outerboundaryDE){
	if(outerboundaryDE)
		m_outer_boundary_type=1;
}

//Read in parameter data from string stream data.
void MGIgesPD144::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_surface_DE);
	get_integer(pDelimeter,pdstream,m_outer_boundary_type);

	int num_inner;
	get_integer(pDelimeter,pdstream,num_inner);
	get_DEpointer(pDelimeter,pdstream,m_outer_boudary_DE);
	for(int i=0; i<num_inner; i++){
		int boundary_DE;
		get_DEpointer(pDelimeter,pdstream,boundary_DE);
		m_inner_boundaries.push_back(boundary_DE);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD144::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_surface_DE,gsec,plines);
	put_integer(m_outer_boundary_type,gsec,plines);

	int n=m_inner_boundaries.size();
	put_integer(n,gsec,plines);
	put_DEpointer(m_outer_boudary_DE,gsec,plines);
	for(int i=0; i<n; i++){
		put_DEpointer(m_inner_boundaries[i],gsec,plines);
	}
}

//Convert de(type=144: trimmed surface) to MGFace.
//Returned is a newed object.
MGFace* MGIgesIfstream::convert_trimmed_surface(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD144* pd144=static_cast<const MGIgesPD144*>(pd.get());
	const MGIgesDirectoryEntry* surfDE=directoryEntry(pd144->m_surface_DE);
	std::auto_ptr<MGGel> surfObj(convert_to_gel(*surfDE));
	MGSurface* surfp=dynamic_cast<MGSurface*>(surfObj.get());
	if(!surfp)
		return 0;

	//std::cout<<(*surfp)<<std::endl;
	std::auto_ptr<MGSurface> surf(static_cast<MGSurface*>(surfObj.release()));
	std::auto_ptr<MGFace> face(new MGFace(surf.release()));

	//1. Make the outer boundary.
	if(pd144->m_outer_boundary_type){
		const MGIgesDirectoryEntry* outBDE=directoryEntry(pd144->m_outer_boudary_DE);
		assert(outBDE->EntityTypeNumber() == CURVE_ON_PARAMETRIC_SURFACE);//Must be curve on a surface.
		const MGIgesPD142* pd142=static_cast<const MGIgesPD142*>(outBDE->paramData().get());
		pd142->trim_face(*this,face);
		if(face->number_of_loops())
			face->loop(0)->make_close();
	}else
		face->make_outer_boundary();

	//2. Make inner boundaries.
	const std::vector<int>& innerBoundaries=pd144->m_inner_boundaries;
	size_t n=innerBoundaries.size();
	for(size_t i=0; i<n; i++){
		const MGIgesDirectoryEntry* boundaryDE=directoryEntry(innerBoundaries[i]);
		if(boundaryDE->EntityTypeNumber() != CURVE_ON_PARAMETRIC_SURFACE)
			continue;

		const std::auto_ptr<MGIgesPD>& pd=boundaryDE->paramData();
		const MGIgesPD142* pd142=static_cast<const MGIgesPD142*>(pd.get());
		pd142->trim_face(*this,face,false);
	}
	return face.release();
}
