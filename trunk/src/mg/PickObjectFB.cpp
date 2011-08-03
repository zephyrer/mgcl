/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#include "mg/PickObjectFB.h"
#include "topo/Face.h"
#include "topo/Edge.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGPickObjectFB is MGPickObject that includes the edge information.
//FB stands for face boundary.
//MGPickObjectFB object is generated when users spedified 2-manifold and boundary
//selection, and the result is the boundary of a face.

// MGPickObject is a class to locate where a picked object is in a group
// hierarchy. Generally, A group includes other groups, and the included groups
// include other groups. In that way the groups make a group hierachy.
// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
// When MGPickObject represents an MGObject, gel() returns MGObject
// pointer and gel_is_object() returns true.
// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
// and gel_is_object() returns false.

//Assignment operator.
MGPickObjectFB& MGPickObjectFB::operator=(const MGPickObject& pobj){
	MGPickObject::operator=(pobj);
	const MGPickObjectFB* pfb=dynamic_cast<const MGPickObjectFB*>(&pobj);
	if(pfb)
		m_edge=pfb->m_edge;
	else
		m_edge=0;
	return *this;
}

////////////////オペレーション////////////

//Generate a newed clone object.
MGPickObjectFB* MGPickObjectFB::clone()const{
	return new MGPickObjectFB(*this);
}

//Return the face of the edge.
MGFace* MGPickObjectFB::face(){
	if(m_edge) return static_cast<MGFace*>(object());
	return 0;
}
