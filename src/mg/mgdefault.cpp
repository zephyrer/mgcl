/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Default.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//  Implementation of class MGDefault
//
// Defines default values of each class.
//
MGEXTERN const MGPosition mgNULL_Pos(0);
MGEXTERN const MGPosition mgORIGIN(0.,0.,0.);
MGEXTERN const MGPosition mgORIGIN_2D(0.,0.);

MGEXTERN const MGVector mgX_UVEC(1.,0.,0.);
MGEXTERN const MGVector mgY_UVEC(0.,1.,0.);
MGEXTERN const MGVector mgZ_UVEC(0.,0.,1.);

MGEXTERN const MGVector mgX_UVEC_2D(1.,0.);
MGEXTERN const MGVector mgY_UVEC_2D(0.,1.);

MGEXTERN const MGInterval mgEMP_INTERV(MGINTERVAL_EMPTY);
MGEXTERN const MGInterval mgZERO_TO_DBLPAI(0.,mgDBLPAI);

MGEXTERN const MGBox mgNULL_BOX(0);
MGEXTERN const MGBox mgEMP_BOX = MGBox(
	MGInterval(),
	MGInterval(),
	MGInterval()
);
MGEXTERN const MGBox mgEMP_BOX_2D = MGBox(
	MGInterval(),
	MGInterval()
);

MGEXTERN const MGMatrix mgNULL_MATR(0);
MGEXTERN const MGMatrix mgUNIT_MATR(3,1.0);
MGEXTERN const MGMatrix mgUNIT_MATR_2D(2,1.0);

MGEXTERN const MGTransf mgNULL_TRANSF(0);
MGEXTERN const MGTransf mgID_TRANSF(MGMatrix(3,1.0), MGVector(0.,0.,0.));
MGEXTERN const MGTransf mgID_TRANSF_2D(MGMatrix(2,1.), MGVector(0.,0.));

MGEXTERN const MGVector mgNULL_VEC(0);
MGEXTERN const MGVector mgZERO_VEC(0.0, 0.0, 0.0);
MGEXTERN const MGVector mgZERO_VEC_2D(0.0, 0.0);

MGEXTERN const MGKnotVector mgNULL_KNOT_VECTOR(0);
MGEXTERN const MGKnotVector mgINFINITE_KNOT_VECTOR(2,2,-mgInfiniteVal, 2.*mgInfiniteVal);

//////////// Abstraction of Gell kind////////////
MGEXTERN const MGAbstractGel mgAll_Gell(MGALL_GELL,MGALL_TID);

MGEXTERN const MGAbstractGel mgAll_Group(MGTOP_KIND,MGGROUP_TID);
MGEXTERN const MGAbstractGel mgAll_Attrib(MGTOP_KIND,MGATTRIB_TID);

////////////////////////// MGObject ///////////////////
MGEXTERN const MGAbstractGel mgAll_Object(MGTOP_KIND,MGOBJECT_TID);

MGEXTERN const MGAbstractGel mgAll_0Manifold(MGMANIFOLD,MG0MANIFOLD);
MGEXTERN const MGAbstractGel mgAll_1Manifold(MGMANIFOLD,MG1MANIFOLD);
MGEXTERN const MGAbstractGel mgAll_2Manifold(MGMANIFOLD,MG2MANIFOLD);

MGEXTERN const MGAbstractGel mgAll_Geo(MGGEO_TOPO,MGGEOMETRY_TID);
MGEXTERN const MGAbstractGel mgAll_Topo(MGGEO_TOPO,MGTOPOLOGY_TID);
MGEXTERN const MGAbstractGel mgAll_STL(MGGEO_TOPO,MGSTL_TID);

MGEXTERN const MGAbstractGel mgAll_Point(MGGEO_KIND,MGPOINT_TID0);
MGEXTERN const MGAbstractGel mgAll_Curve(MGGEO_KIND,MGCURVE_TID);
MGEXTERN const MGAbstractGel mgAll_Surface(MGGEO_KIND,MGSURFACE_TID);
MGEXTERN const MGAbstractGel mgAll_FSurface(MGFSURFACE_KIND,MGFSURFACE_TID);

MGEXTERN const MGAbstractGel mgAll_Straight(MGLEAF_KIND,MGSTRAIGHT_TID);
MGEXTERN const MGAbstractGel mgAll_Ellipse(MGLEAF_KIND,MGELLIPSE_TID);
MGEXTERN const MGAbstractGel mgAll_LBRep(MGLEAF_KIND,MGLBREP_TID);
MGEXTERN const MGAbstractGel mgAll_RLBRep(MGLEAF_KIND,MGRLBREP_TID);
MGEXTERN const MGAbstractGel mgAll_SurfCurve(MGLEAF_KIND,MGSRFCRV_TID);
MGEXTERN const MGAbstractGel mgAll_TrimmedCurve(MGLEAF_KIND,MGTRMCRV_TID);
MGEXTERN const MGAbstractGel mgAll_CompositeCurve(MGLEAF_KIND,MGCOMPCRV_TID);
MGEXTERN const MGAbstractGel mgAll_Plane(MGLEAF_KIND,MGPLANE_TID);
MGEXTERN const MGAbstractGel mgAll_SPhere(MGLEAF_KIND,MGSPHERE_TID);
MGEXTERN const MGAbstractGel mgAll_SBRep(MGLEAF_KIND,MGSBREP_TID);
MGEXTERN const MGAbstractGel mgAll_RSBRep(MGLEAF_KIND,MGRSBREP_TID);
MGEXTERN const MGAbstractGel mgAll_Cylinder(MGLEAF_KIND,MGCYLINDER_TID);
MGEXTERN const MGAbstractGel mgAll_PVertex(MGLEAF_KIND,MGPVERTEX_TID);
MGEXTERN const MGAbstractGel mgAll_BVertex(MGLEAF_KIND,MGBVERTEX_TID);
MGEXTERN const MGAbstractGel mgAll_Edge(MGLEAF_KIND,MGEDGE_TID);
MGEXTERN const MGAbstractGel mgAll_Face(MGLEAF_KIND,MGFACE_TID);
MGEXTERN const MGAbstractGel mgAll_Loop(MGLEAF_KIND,MGLOOP_TID);
MGEXTERN const MGAbstractGel mgAll_Shell(MGLEAF_KIND,MGSHELL_TID);
