/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDefault_HH_
#define _MGDefault_HH_
/** @addtogroup BASE
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/types.h"
#include "mg/DefaultVector.h"
#include "mg/Box.h"
#include "mg/Transf.h"
#include "mg/KnotVector.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"

//  MGDefault.h
//  Header for class MGDefault

MGEXTERN const MGInterval mgEMP_INTERV;
MGEXTERN const MGInterval mgZERO_TO_DBLPAI;

MGEXTERN const MGBox mgNULL_BOX;
MGEXTERN const MGBox mgEMP_BOX;
MGEXTERN const MGBox mgEMP_BOX_2D;

MGEXTERN const MGMatrix mgNULL_MATR;
MGEXTERN const MGMatrix mgUNIT_MATR;
MGEXTERN const MGMatrix mgUNIT_MATR_2D;

MGEXTERN const MGTransf mgNULL_TRANSF;
MGEXTERN const MGTransf mgID_TRANSF;
MGEXTERN const MGTransf mgID_TRANSF_2D;

MGEXTERN const MGKnotVector mgNULL_KNOT_VECTOR;
MGEXTERN const MGKnotVector mgINFINITE_KNOT_VECTOR;

MGEXTERN const MGAbstractGel mgAll_Gell;

MGEXTERN const MGAbstractGel mgAll_Object;
MGEXTERN const MGAbstractGel mgAll_Group;
MGEXTERN const MGAbstractGel mgAll_Attrib;

MGEXTERN const MGAbstractGel mgAll_0Manifold;
MGEXTERN const MGAbstractGel mgAll_1Manifold;
MGEXTERN const MGAbstractGel mgAll_2Manifold;

MGEXTERN const MGAbstractGel mgAll_Geo;
MGEXTERN const MGAbstractGel mgAll_Topo;
MGEXTERN const MGAbstractGel mgAll_STL;

MGEXTERN const MGAbstractGel mgAll_Point;
MGEXTERN const MGAbstractGel mgAll_Curve;
MGEXTERN const MGAbstractGel mgAll_Surface;
MGEXTERN const MGAbstractGel mgAll_FSurface;

MGEXTERN const MGAbstractGel mgAll_Straight;
MGEXTERN const MGAbstractGel mgAll_Ellipse;
MGEXTERN const MGAbstractGel mgAll_LBRep;
MGEXTERN const MGAbstractGel mgAll_RLBRep;
MGEXTERN const MGAbstractGel mgAll_SurfCurve;
MGEXTERN const MGAbstractGel mgAll_TrimmedCurve;
MGEXTERN const MGAbstractGel mgAll_CompositeCurve;
MGEXTERN const MGAbstractGel mgAll_Plane;
MGEXTERN const MGAbstractGel mgAll_SPhere;
MGEXTERN const MGAbstractGel mgAll_SBRep;
MGEXTERN const MGAbstractGel mgAll_RSBRep;
MGEXTERN const MGAbstractGel mgAll_Cylinder;
MGEXTERN const MGAbstractGel mgAll_PVertex;
MGEXTERN const MGAbstractGel mgAll_BVertex;
MGEXTERN const MGAbstractGel mgAll_Edge;
MGEXTERN const MGAbstractGel mgAll_Face;
MGEXTERN const MGAbstractGel mgAll_Loop;
MGEXTERN const MGAbstractGel mgAll_Shell;

/// Defines default values of each class.
class MGDefault{

public:

///String stream function.
MGDECL friend std::ostream& operator<<(std::ostream&, const MGDefault&);

///void constructor
///MGDefault();

///Return null box.
static const MGBox& null_box(){return mgNULL_BOX;};

///Return empty 3D box.  m_empty_box ‚ð•Ô‹pB
static const MGBox& empty_box(){return mgEMP_BOX;};

///Return empty 2D box.  m_empty_box_2D ‚ð•Ô‹pB
static const MGBox& empty_box_2D(){return mgEMP_BOX_2D;};

///Return empty interval.  m_empty_interval ‚ð•Ô‹pB
static const MGInterval& empty_interval(){return mgEMP_INTERV;};

///Return null transformation.
static const MGTransf& null_transf(){return mgNULL_TRANSF;};

///Return 4 by 4 unit transformation.  m_identity_transf ‚ð•Ô‹pB
static const MGTransf& identity_transf(){return mgID_TRANSF;};

///Return 3 by 3 unit transformation.  m_identity_transf_2D ‚ð•Ô‹pB
static const MGTransf& identity_transf_2D(){return mgID_TRANSF_2D;};

///Return 3D origin point(0,0,0).  m_origin ‚ð•Ô‹pB
static const MGPosition& origin(){return mgORIGIN;};

///Return 2D origin point(0,0).  m_origin_2D ‚ð•Ô‹pB
static const MGPosition& origin_2D(){return mgORIGIN_2D;};

///Return null matrix.
static const MGMatrix& null_matrix(){return mgNULL_MATR;};

///Return 3 by 3 unit matrix.  m_unit_matrix ‚ð•Ô‹pB
static const MGMatrix& unit_matrix(){return mgUNIT_MATR;};

///Return 2 by 2 unit matrix.  m_unit_matrix_2D ‚ð•Ô‹pB
static const MGMatrix& unit_matrix_2D(){return mgUNIT_MATR_2D;};

///Return x axis unit vector.  m_x_unit_vector_2D ‚ð•Ô‹pB
static const MGVector& x_unit_vector_2D(){return mgX_UVEC_2D;};

///Return y axis unit vector.  m_y_unit_vector_2D ‚ð•Ô‹pB
static const MGVector& y_unit_vector_2D(){return mgY_UVEC_2D;};

///Return x axis unit vector.  m_x_unit_vector ‚ð•Ô‹pB
static const MGVector& x_unit_vector(){return mgX_UVEC;};

///Return y axis unit vector.  m_y_unit_vector ‚ð•Ô‹pB
static const MGVector& y_unit_vector(){return mgY_UVEC;};

///Return z axis unit vector.  m_z_unit_vector ‚ð•Ô‹pB
static const MGVector& z_unit_vector(){return mgZ_UVEC;};

///Return null position.
static const MGPosition& null_position(){return mgNULL_Pos;};

///Return null vector.
static const MGVector& null_vector(){return mgNULL_VEC;};

///Return 3D zero vector.  m_zero_vector ‚ð•Ô‹pB
static const MGVector& zero_vector(){return mgZERO_VEC;};

///Return 2D zero vector.  m_zero_vector_2D ‚ð•Ô‹pB
static const MGVector& zero_vector_2D(){return mgZERO_VEC_2D;};

///Return null knot vector.  m_null_knot_vector ‚ð•Ô‹pB
static const MGKnotVector& null_knot_vector(){return mgNULL_KNOT_VECTOR;};

///Return infinite knot vector.  m_infinite_knot_vector ‚ð•Ô‹pB
static const MGKnotVector& infinite_knot_vector(){return mgINFINITE_KNOT_VECTOR;};

};

/** @} */ // end of BASE group
#endif
