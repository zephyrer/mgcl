/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD510.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD510.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesOfstream.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD510 is the class for Iges parameter data type 510(FACE Entity).

// Constructors.

//! Constructs an object of class MGIgesPD510.
MGIgesPD510::MGIgesPD510(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(FACE,DEpointer),m_outer_loop_identified(true){
}

//Read in parameter data from string stream data.
void MGIgesPD510::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_surface_DE);
	int num_loops;
	get_integer(pDelimeter,pdstream,num_loops);
	m_loops.resize(num_loops);
	int outer;
	get_integer(pDelimeter,pdstream,outer);
	m_outer_loop_identified=outer ? true:false;
	for(int i=0; i<num_loops; i++){
		get_DEpointer(pDelimeter,pdstream,m_loops[i]);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD510::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_surface_DE,gsec,plines);
	int num_loops=m_loops.size();
	put_integer(num_loops,gsec,plines);
	int outer=m_outer_loop_identified ? 1:0;
	put_integer(outer,gsec,plines);
	for(int i=0; i<num_loops; i++)
		put_DEpointer(m_loops[i],gsec,plines);
}

//consvert m_surface_DE surface to MGSurface.
//Returned is a newed MGSurface object.
MGSurface* MGIgesPD510::convert_to_surface(
	const MGIgesIfstream& igesIstream
)const{
	MGGel* gel=igesIstream.convert_to_gel(m_surface_DE);
	MGSurface* srf=dynamic_cast<MGSurface*>(gel);
	return srf;
}

//Convert de(type=510: FACE) to MGFace.
//Returned is a newed MGFace object.
MGFace* MGIgesIfstream::convert_face(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD510* pd510=static_cast<const MGIgesPD510*>(pd.get());
	MGSurface* srf=pd510->convert_to_surface(*this);

/////////******
//	bool this_is_it=false;
//	MGSBRep* rsb=dynamic_cast<MGSBRep*>(srf);
//	if(rsb){
//		MGKnotVector& tu=rsb->knot_vector_u();
//		MGKnotVector& tv=rsb->knot_vector_v();
//		if(tu.order()==6 && tu.bdim()==6){
//		if(tv.order()==2 && tv.bdim()==2){
//		if(4.625<tu(0)&&tu(0)<4.626){
//		if(1.661<tv(3)&&tv(3)<1.662){
//		this_is_it=true;
//		std::cout<<(*rsb)<<std::endl;
//		}		
//		}
//		}
//		}
//	}
/////////******
	MGFace* face;
	if(srf){
		face=new MGFace(srf);
		int inner_loop_id;
		if(pd510->m_outer_loop_identified){
			MGLoop* loop0=convert_loop(*directoryEntry(pd510->m_loops[0]),*srf);
			if(loop0){
				if(loop0->area()<0.)
					loop0->negate();
				face->prepend_boundary(loop0);
			}else
				face->make_outer_boundary();
			inner_loop_id=1;
		}else{
			face->make_outer_boundary();
			inner_loop_id=0;
		}
		size_t nloop=pd510->m_loops.size();
		for(size_t i=inner_loop_id; i<nloop; i++){
			MGLoop* loopi=convert_loop(*directoryEntry(pd510->m_loops[i]),*srf);
			if(loopi){
				if(loopi->area()>0.)
					loopi->negate();
				face->append_boundary(loopi);		
			}
		}
	}else
		face=0;
	return face;
}
