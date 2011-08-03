/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mg/Point.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/MGStl.h"
#include "mg/Tolerance.h"

#include "topo/Topology.h"
#include "topo/Complex.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/Boundary.h"
#include "topo/Loop.h"
#include "topo/Shell.h"
#include "Tl/TLData.h"
#include "mgGL/Appearance.h"
#include "mg/GelFactory.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGObject
// Implementation of MGObject.

//Void constructor(初期化なしでオブジェクトを作成する。)
MGObject::MGObject():m_appearance(0){
}

//Copy constructor.
MGObject::MGObject(const MGObject& obj2):m_appearance(0){
	if(obj2.m_appearance)
		m_appearance=obj2.m_appearance->clone();
}

//Virtual Destructor
MGObject::~MGObject(){
	delete m_appearance;
}

//Assignment.
//When the leaf object of this and gel2 are not equal, this assignment
//does nothing.
MGObject& MGObject::set_object(const MGObject& gel2){
	MGAttribedGel::operator=(gel2);
	delete m_appearance;
	if(gel2.m_appearance){
		m_appearance=gel2.m_appearance->clone();
	}else{
		m_appearance=0;
	}
	return *this;
}

//Construct a null newed MGObject from the type id TID.
//Object handled by MGIfstream or MGOfstream is only the following objects.
MGObject* MGNullObj(long TID){
	MGGelFactoryRegistry* reg = MGGelFactoryRegistry::get_instance();
	return static_cast<MGObject*>(reg->create_gel(TID));
}

AUTO_GEL_REGISTER(MGPoint, MGPOINT_TID);

AUTO_GEL_REGISTER(MGStraight, MGSTRAIGHT_TID);
AUTO_GEL_REGISTER(MGEllipse, MGELLIPSE_TID);
AUTO_GEL_REGISTER(MGLBRep, MGLBREP_TID);
AUTO_GEL_REGISTER(MGRLBRep, MGRLBREP_TID);
AUTO_GEL_REGISTER(MGSurfCurve, MGSRFCRV_TID);
AUTO_GEL_REGISTER(MGTrimmedCurve, MGTRMCRV_TID);
AUTO_GEL_REGISTER(MGCompositeCurve, MGCOMPCRV_TID);
AUTO_GEL_REGISTER(MGBSumCurve, MGBSUMCRV_TID);

AUTO_GEL_REGISTER(MGPlane, MGPLANE_TID);
AUTO_GEL_REGISTER(MGSphere, MGSPHERE_TID);
AUTO_GEL_REGISTER(MGSBRep, MGSBREP_TID);
AUTO_GEL_REGISTER(MGRSBRep, MGRSBREP_TID);
AUTO_GEL_REGISTER(MGCylinder, MGCYLINDER_TID);
AUTO_GEL_REGISTER(MGBSumSurf, MGBSUMSURF_TID);

AUTO_GEL_REGISTER(MGPVertex, MGPVERTEX_TID);
AUTO_GEL_REGISTER(MGBVertex, MGBVERTEX_TID);
AUTO_GEL_REGISTER(MGEdge, MGEDGE_TID);
AUTO_GEL_REGISTER(MGFace, MGFACE_TID);
AUTO_GEL_REGISTER(MGLoop, MGLOOP_TID);
AUTO_GEL_REGISTER(MGShell, MGSHELL_TID);
AUTO_GEL_REGISTER(MGComplex, MGCOMPLEX_TID);
AUTO_GEL_REGISTER(MGStl, MGSTL_TID);

//make this group has appearance and get the MGAppearance pointer.
MGAppearance* MGObject::ensure_appearance(){
	if(!m_appearance)
		m_appearance=new MGAppearance();
	return m_appearance;
}

//Test if this and 2nd object has common area about their box(),
//taking error into account.
bool MGObject::has_common(const MGObject& obj2) const{
	double err=MGTolerance::wc_zero(); err*=.5;
	MGBox bx=box(); bx.expand(err);
	MGBox bx2=obj2.box(); bx2.expand(err);
	bx&=bx2;
	return !bx.empty();
}

//Remove the MGAppearance of this MGAttribedGel.
void MGObject::remove_appearance(){
	if(!m_appearance)
		return;
	delete m_appearance;
	m_appearance=0;
}

//Read all member data.
void MGObject::ReadMembers(MGIfstream& buf){
	MGAttribedGel::ReadMembers(buf);
	const char* iver=buf.version();
	m_appearance=static_cast<MGAppearance*>(buf.ReadPointer());
}

//Write all member data
void MGObject::WriteMembers(MGOfstream& buf)const{
	MGAttribedGel::WriteMembers(buf);
	buf.WritePointer(m_appearance);
}

