/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#include "mg/PickObjectSB.h"
#include "topo/Face.h"
#include "topo/Edge.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGPickObjectSB is a MGPickObject that includes the perimeter information of
//a MGSurface.
//SB stands for surface boundary.
//MGPickObjectSB object is generated when users spedified 2-manifold and boundary
//selection, and the result is the boundary of a MGSurface.

// MGPickObject is a class to locate where a picked object is in a group
// hierarchy. Generally, A group includes other groups, and the included groups
// include other groups. In that way the groups make a group hierachy.
// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
// When MGPickObject represents an MGObject, gel() returns MGObject
// pointer and gel_is_object() returns true.
// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
// and gel_is_object() returns false.

//Assignment operator.
MGPickObjectSB& MGPickObjectSB::operator=(const MGPickObject& pobj){
	MGPickObject::operator=(pobj);
	const MGPickObjectSB* psb=dynamic_cast<const MGPickObjectSB*>(&pobj);
	if(psb)
		m_perimeter=psb->m_perimeter;
	else
		m_perimeter=0;
	return *this;
}

////////////////オペレーション////////////

//Generate a newed clone object.
MGPickObjectSB* MGPickObjectSB::clone()const{
	return new MGPickObjectSB(*this);
}

//Return the face of the edge.
MGSurface* MGPickObjectSB::surface(){
	if(m_perimeter<0) return 0;
	return static_cast<MGSurface*>(object());
}

//Return the face of the edge.
const MGSurface* MGPickObjectSB::surface()const{
	if(m_perimeter<0) return 0;
	return static_cast<const MGSurface*>(object());
}
