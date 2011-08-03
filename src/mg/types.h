/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGEL_TID_HH_
#define _MGGEL_TID_HH_

#include <utility>

/** @addtogroup GelRelated
 *  @{
 */

///Tyoe id of subclasses of MGGel.
enum MGGEL_TID{
///The mm in 0xmm000000L of the following TID's are subclass id of MGGel.
///When a new subclass is necessary to add, this mm will have a new number.
	MGALL_TID =			0x00000000L,	///<all of specified kind of MGGEL_KIND.
	MGOBJECT_TID =		0x00000000L,
	MGGROUP_TID =		0x01000000L,
	MGATTRIB_TID =		0x02000000L,

	MG0MANIFOLD=		0x00000000L,
	MG1MANIFOLD=		0x00000100L,
	MG2MANIFOLD=		0x00000200L,
	MG3MANIFOLD=		0x00000300L,

///*************Define MGOBJECT_TID***************
///	Following nn is the manifold dimension, xx is the name id.
///						0x0010nnxxL
	MGGEOMETRY_TID =	0x00100000L,
	MGPOINT_TID0 =		0x00100000L,
	MGPOINT_TID =		0x00100001L,

	MGCURVE_TID =		0x00100100L,
	MGSTRAIGHT_TID=		0x00100101L,
	MGELLIPSE_TID =		0x00100102L,
	MGLBREP_TID =		0x00100103L,
	MGRLBREP_TID =		0x00100104L,
	MGSRFCRV_TID =		0x00100105L,
	MGTRMCRV_TID =		0x00100106L,
	MGCOMPCRV_TID =		0x00100107L,
	MGBSUMCRV_TID =		0x00100108L,

	MGSURFACE_TID =		0x00100200L,
	MGPLANE_TID =		0x00100201L,
	MGSPHERE_TID =		0x00100203L,
	MGSBREP_TID =		0x00100205L,
	MGRSBREP_TID =		0x00100206L,
	MGCYLINDER_TID =	0x00100207L,
	MGBSUMSURF_TID =	0x00100208L,

/////////////////////////////////////////////

///	Following nn is the manifold dimension, xx is the name id,
///	and m is 0 for cells, and m is 1 for complexes.
///						0x002mnnxxL
	MGTOPOLOGY_TID =	0x00200000L,///< is a topology
	MGCELL_TID=			0x00200000L,///< is a cell.
	MGCOMPLEX_TID=		0x00210000L,///< is a complex.

///	MGCELLBASE_TID =	0x0020nnxxL,
	MGPVERTEX_TID =		0x00200001L,
	MGBVERTEX_TID =		0x00200002L,
	MGEDGE_TID =		0x00200101L,
	MGFACE_TID =		0x00200201L,
	MGSOLID_TID =		0x00200301L,

///	MGBOUNDARYND_TID =	0x0021nnxxL,
	MGLOOP_TID =		0x00210101L,
	MGSHELL_TID =		0x00210201L,

///////////////////////////////////
///id for MGSurface or MGFace
	MGFSURFACE_TID =	0x00000200L,

///////////////////////////////////
///id for MGStl
	MGSTL_TID =	0x00300200L,

/////////////////////////////////////////////
///  MGGROUP id.
	MGAPPEARANCE_TID =	0x01000002L,	///<MGAppearance(attributes).
	MGLLOBJECTS_TID =	0x01000003L,	///<Locally lighted objects.

/////////////////////////////////////////////
///  MGAttrib id.
	MGCONTEXT_TID =		0x02010010L,

	MGGLATTRIBUTE_TID=	0x02010000L,
	MGLIGHTS_TID =		0x02010020L,
	MGLIGHT_TID =		0x02010030L,
	MGDIRECTIONAL_LIGHT_TID =
						0x02010031L,
	MGPOINT_LIGHT_TID = 0x02010032L,
	MGSPOT_LIGHT_TID =	0x02010033L,

	MGFOG_TID =			0x02010040L,

	MGMATERIAL_TID =	0x02010050L,
	MGALPHA_FUNC_TID =	0x02010060L,
	MGBLEND_FUNC_TID =	0x02010070L,
	MGCOLOR_TID =		0x02010080L,
	MGCOLOR_MASK_TID =	0x02010090L,
	MGDEPTH_FUNC_TID =	0x020100A0L,
	MGDEPTH_MASK_TID =	0x020100B0L,
	MGLIGHT_ENABLE_TID=	0x020100C0L,
	MGLINE_STIPPLE_TID=	0x020100D0L,
	MGLINE_WIDTH_TID =	0x020100E0L,
	MGPOLYGON_MODE_TID=	0x020100F0L,
	MGRENDER_ATTR_TID=	0x02010100L,
	MGSHADE_MODEL_TID=	0x02010110L,
	MGTRANSP_MODE_TID=	0x02010120L,
	MGTEXTURE_TID=		0x02010200L
};

///MGGEL_KIND_TID is used to specify what kind of group is used to identify gels.
enum MGGEL_KIND{
	MGALL_GELL =	0x00000000L,///<all of the gels
	MGTOP_KIND=		0xff000000L,///<subkind is MGOBJECT_TID, MGGROUP_TID, MGATTRIB_TID
	MGMANIFOLD=		0xff00ff00L,///<subkind is MG0MANIFOLD, MG1MANIFOLD, MG2MANIFOLD, MG3MANIFOLD
	MGFSURFACE_KIND=0xff0fff00L,///<subkind is MGFSURFACE_TID(MGFace or MGSurface, MG2MANIFOLD that is not
		///<MGShell)
	MGGEO_TOPO=		0xfff00000L,///<MGGeometry, or MGTopology.
		///<subkind is MGGEOMETRY_TID, MGTOPOLOGY_TID
	MGGEO_KIND=		0xffffff00L,///<subkind is MGPOINT_TID0, MGCURVE_TID, MGSURFACE_TID
	MGLEAF_KIND=	0xffffffffL	///<subkind is all of the leaf class MGGEL_TID, that is:
		///<MGPOINT_TID,MGSTRAIGHT_TID,MGELLIPSE_TID,MGLBREP_TID,MGRLBREP_TID,
		///<MGSRFCRV_TID,MGTRMCRV_TID,MGCOMPCRV_TID,MGPLANE_TID,MGSPHERE_TID,
		///<MGSBREP_TID,MGRSBREP_TID,MGCYLINDER_TID,MGPVERTEX_TID,
		///<MGBVERTEX_TID,MGEDGE_TID,MGFACE_TID,MGLOOP_TID,MGSHELL_TID
};

///MGAbstractGel is a class to specify what kind of abstract gel group.
///Let MGAbstractGel agel(gel_kind, gel_tid), then
///gel_kind is either MGALL_GELL, MGTOP_KIND, MGMANIFOLD, MGGEO_TOPO,
///MGGEO_KIND, or MGLEAF_KIND. And gel_tid is the specific type of
///the gel_kind. For the value of gel_tid, see MGGEL_KIND above.
///Possible combinations are defined in MGDefault.h(as mgAll_xxxx).
///See the definition.
typedef std::pair<MGGEL_KIND,MGGEL_TID> MGAbstractGel;

/** @} */ // end of GelRelated group
#endif
