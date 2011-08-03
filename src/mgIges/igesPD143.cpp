/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD143.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD141.h"
#include "mgiges/IgesPD143.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD143 is the class for Iges parameter data type 143(Bounded Surface).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD143.
MGIgesPD143::MGIgesPD143(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(BOUNDED_SURFACE,DEpointer),m_type(0),m_surface_DE(0){
}

//Read in parameter data from string stream data.
void MGIgesPD143::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_type);
	get_DEpointer(pDelimeter,pdstream,m_surface_DE);

	int n;
	get_integer(pDelimeter,pdstream,n);//Number of boundaries.
	for(int i=0; i<n; i++){
		int boundary_DE;
		get_DEpointer(pDelimeter,pdstream,boundary_DE);
		m_boundaries.push_back(boundary_DE);
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD143::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_integer(m_type,gsec,plines);
	put_DEpointer(m_surface_DE,gsec,plines);

	int n=m_boundaries.size();
	put_integer(n,gsec,plines);
	for(int i=0; i<n; i++)
		put_DEpointer(m_boundaries[i],gsec,plines);
}

//Convert de(type=143: bounded surface) to MGFace.
//Returned is a newed object.
MGFace* MGIgesIfstream::convert_bounded_surface(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD143* pd143=static_cast<const MGIgesPD143*>(pd.get());
	const MGIgesDirectoryEntry* surfDE=directoryEntry(pd143->m_surface_DE);
	std::auto_ptr<MGGel> surfObj(convert_to_gel(*surfDE));
	MGSurface* surfp=dynamic_cast<MGSurface*>(surfObj.get());
	if(!surfp)
		return 0;

	bool curve_is_not_parameter_space=
		(pd143->m_type==0 || surfDE->EntityTypeNumber() != RATIONAL_BSPLINE_SURFACE);
	std::auto_ptr<MGSurface> surf(static_cast<MGSurface*>(surfObj.release()));
	std::auto_ptr<MGFace> face(new MGFace(surf.release()));
	const std::vector<int>& boundaries=pd143->m_boundaries;
	size_t n=boundaries.size();
	for(size_t i=0; i<n; i++){
		const MGIgesDirectoryEntry* boundaryDE=directoryEntry(boundaries[i]);
		if(boundaryDE->EntityTypeNumber() != 141)
			continue;

		std::auto_ptr<MGFace> faceSave(face->clone());
		const std::auto_ptr<MGIgesPD>& pd=boundaryDE->paramData();
		const MGIgesPD141* pd141=static_cast<const MGIgesPD141*>(pd.get());
		const std::vector<MGIges141Edge>& edges=pd141->m_edges;
		size_t nedges=edges.size(), j=0;
		for(; j<nedges; j++){//Loop for one edge.
			int error=0;
			if(curve_is_not_parameter_space){
			//128==NURBS Surface. We can handle only NURBS surface that has the boundary curve of
			//the parameter space as its boundary.
				std::auto_ptr<MGGel> objp(convert_to_gel(edges[j].m_curve_DE));
				MGCurve* crvp=dynamic_cast<MGCurve*>(objp.get());
					//crvp is a model space curve.
				if(!crvp)
					break;
				if(edges[j].m_sense==2)
					crvp->negate();
				std::auto_ptr<MGLoop> loop=face->build_loop(*crvp);
				if(loop->area()>0.)
					face->prepend_boundary(loop.release());
				else
					face->append_boundary(loop.release());
			}else{
			//The case that the base surface is NURBS Surface and the boudaries are expressed
			//in the parameter space.
				const std::vector<int>& pcurves=edges[j].m_pcurves;
				size_t npcurves=pcurves.size();
				size_t k=0;
				for(; k<npcurves; k++){
					std::auto_ptr<MGGel> objp(convert_to_gel(pcurves[k]));
					MGCurve* crvp=dynamic_cast<MGCurve*>(objp.get());
						//crvp is in the surface parameter space(u,v).
					if(!crvp)
						break;
					error=face->trim(*crvp);
					if(error)
						break;
				}
				if(k!=npcurves)
					break;
			}
			if(error)
				break;
		}//Loop for one edge.
		if(j!=nedges){
			face=faceSave;//Restore the face without i-th boundary.
			continue;
		}
	}
	face->make_outer_boundary();//If outer boundary is not built, will be added.
	return face.release();
}
