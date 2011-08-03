/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/// MGPickObjectSB.h : MGPickObjectSB クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _MGPickObjectSB_HH_
#define _MGPickObjectSB_HH_

#include "mg/PickObject.h"

class MGSurface;

/** @addtogroup MGObjectRelated
 *  @{
 */

///MGPickObjectSB is a MGPickObject that includes the perimeter information of
///a MGSurface.
///SB stands for surface boundary.
///MGPickObjectSB object is generated when users spedified 2-manifold and boundary
///selection, and the result is the boundary of a MGSurface.
/// MGPickObject is a class to locate where a picked object is in a group
/// hierarchy. Generally, A group includes other groups, and the included groups
/// include other groups. In that way the groups make a group hierachy.
/// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
/// When MGPickObject represents an MGObject, gel() returns MGObject
/// pointer and gel_is_object() returns true.
/// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
/// and gel_is_object() returns false.
class MGCLASS MGPickObjectSB:public MGPickObject{

public:

///////////////Constructor//////////////

MGPickObjectSB():MGPickObject(),m_perimeter(-1){;};

///Conversion constructor from MGGelPosition and perimeter.
MGPickObjectSB(
	MGGelPosition& gelp,
	int perimeter
):MGPickObject(gelp),m_perimeter(perimeter){;};

///Conversion constructor from MGPickObject and start/end.
MGPickObjectSB(
	const MGPickObject& pobj,
	int perimeter
):MGPickObject(pobj),m_perimeter(perimeter){;};

///Copy constructor.
///MGPickObjectSB(const MGPickObjectSB& pobj2);

////////////Destructor////////////////
virtual ~MGPickObjectSB(){;};

///Assignment operator.
MGPickObjectSB& operator=(const MGPickObject& pobj);

//////////////////オペレーション//////////////

///Generate a newed clone object.
virtual MGPickObjectSB* clone()const;

///Highlightthe object using the display list of this object.
void hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Return the edge pointer.
int perimeter()const{return m_perimeter;};

///Return the face of the edge.
MGSurface* surface();
const MGSurface* surface()const;

///Set the object pointer.
void set_perimeter(int perimeter){m_perimeter=perimeter;};

private:

	int m_perimeter;	///perimeter number of the MGSurface picked.
};

/** @} */ // end of MGObjectRelated group

#endif
