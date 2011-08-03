/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/// MGPickObjectFB.h : MGPickObjectFB クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _MGPickObjectFB_HH_
#define _MGPickObjectFB_HH_

#include "mg/PickObject.h"

class MGEdge;
class MGFace;

/** @addtogroup MGObjectRelated
 *  @{
 */

///MGPickObjectFB is MGPickObject that includes the edge information.
///FB stands for face boundary.
///MGPickObjectFB object is generated when users spedified 2-manifold and boundary
///selection, and the result is the boundary of a face.
/// MGPickObject is a class to locate where a picked object is in a group
/// hierarchy. Generally, A group includes other groups, and the included groups
/// include other groups. In that way the groups make a group hierachy.
/// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
/// When MGPickObject represents an MGObject, gel() returns MGObject
/// pointer and gel_is_object() returns true.
/// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
/// and gel_is_object() returns false.
class MGCLASS MGPickObjectFB:public MGPickObject{

public:

///////////////Constructor//////////////

MGPickObjectFB():MGPickObject(),m_edge(0){;};

///Conversion constructor from MGGelPosition and MGEdge.
MGPickObjectFB(
	MGGelPosition& gelp,
	const MGEdge* edge
):MGPickObject(gelp),m_edge(edge){;};

///Conversion constructor from MGPickObject and start/end.
MGPickObjectFB(
	MGPickObject& pobj,
	const MGEdge* edge
):MGPickObject(pobj),m_edge(edge){;};

///Copy constructor.
///MGPickObjectFB(const MGPickObjectFB& pobj2);

////////////Destructor////////////////
virtual ~MGPickObjectFB(){;};

///Assignment operator.
MGPickObjectFB& operator=(const MGPickObject& pobj);

//////////////////オペレーション//////////////

///Generate a newed clone object.
virtual MGPickObjectFB* clone()const;

///Return the edge pointer.
const MGEdge* edge()const{return m_edge;};

///Return the face of the edge.
MGFace* face();

///Highlightthe object using the display list of this object.
void hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Set the object pointer.
void set_edge(const MGEdge* edge){m_edge=edge;};

private:

	const MGEdge* m_edge;	///MGEdge pointer of the face picked.
};

/** @} */ // end of MGObjectRelated group
#endif
