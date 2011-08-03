/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#include "mg/PickObjectCB.h"
#include "mg/Curve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGPickObjectCB is a MGPickObject that includes the boundary information of
//a MGCurve.
//CB stands for curve boundary.
//MGPickObjectCB object is generated when users spedified 1-manifold and boundary
//selection.

// MGPickObject is a class to locate where a picked object is in a group
// hierarchy. Generally, A group includes other groups, and the included groups
// include other groups. In that way the groups make a group hierachy.
// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
// When MGPickObject represents an MGObject, gel() returns MGObject
// pointer and gel_is_object() returns true.
// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
// and gel_is_object() returns false.

//Assignment operator.
MGPickObjectCB& MGPickObjectCB::operator=(const MGPickObject& pobj){
	MGPickObject::operator=(pobj);
	const MGPickObjectCB* pcb=dynamic_cast<const MGPickObjectCB*>(&pobj);
	if(pcb)
		m_start_end=pcb->m_start_end;
	else
		m_start_end=0;
	return *this;
}

////////////////オペレーション////////////

//Generate a newed clone object.
MGPickObjectCB* MGPickObjectCB::clone()const{
	return new MGPickObjectCB(*this);
}

//Return the face of the edge.
const MGCurve* MGPickObjectCB::curve()const{
	if(m_start_end<0) return 0;
	return static_cast<const MGCurve*>(object());
}
