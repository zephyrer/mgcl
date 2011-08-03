/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

//! @file
//!	@brief  Implementaion for class MGIgesOfstream.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"

#if defined(MGCL_NO_MFC)
#include "mg/Point.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/CompositeCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/Sphere.h"
#include "mg/RSBRep.h"
#include "mg/SurfCurve.h"
#include "mg/AttribedGel.h"
#include "mg/Group.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGL/Color.h"
int MGPosition::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGVector::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGTransf::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGPoint::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGStraight::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGEllipse::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGLBRep::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGRLBRep::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGCompositeCurve::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGTrimmedCurve::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGSurfCurve::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGBSumCurve::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGPlane::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGCylinder::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGSBRep::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGRSBRep::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGSphere::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGFace::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGGroup::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGColor::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
int MGShell::out_to_IGES(MGIgesOfstream& igesfile,int SubordinateEntitySwitch)const{return 0;}
#endif // defined MGCL_NO_MFC
