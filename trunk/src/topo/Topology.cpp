/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/Topology.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implement MGTopology Class.
//MGTopology is an abstract class which represents a whole Topology,
//Complex, Cell, and Boundary.

//Constructor

//Void constructor(初期化なしでオブジェクトを作成する。)
MGTopology::MGTopology(){;};

//Virtual Destructor
MGTopology::~MGTopology(){;};

MGTopology& MGTopology::operator=(const MGTopology& gel2){
	MGObject::operator=(gel2);
	return *this;
}

//Member Function

//Compute the intersections of two objects.
MGisects MGTopology::intersection(const MGObject& obj2)const{
	return MGisects();
}
MGisects MGTopology::intersection(const MGCurve& obj2)const{
	return MGisects();
}
MGisects MGTopology::intersection(const MGFSurface& obj2)const{
	return MGisects();
}
MGisects MGTopology::intersection(const MGSurface& obj2)const{
	return MGisects();
}
MGisects MGTopology::intersection(const MGFace& obj2)const{
	return MGisects();
}
MGisects MGTopology::intersection(const MGShell& obj2)const{
	return MGisects();
}
