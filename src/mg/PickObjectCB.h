/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/// MGPickObjectCB.h : MGPickObjectCB クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _MGPickObjectCB_HH_
#define _MGPickObjectCB_HH_

#include "mg/PickObject.h"

class MGCurve;

/** @addtogroup MGObjectRelated
 *  @{
 */

///MGPickObjectCB is a MGPickObject that includes the boundary information of
///a MGCurve.
///CB stands for curve boundary.
///MGPickObjectCB object is generated when users spedified 1-manifold and boundary
///selection.
/// MGPickObject is a class to locate where a picked object is in a group
/// hierarchy. Generally, A group includes other groups, and the included groups
/// include other groups. In that way the groups make a group hierachy.
/// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
/// When MGPickObject represents an MGObject, gel() returns MGObject
/// pointer and gel_is_object() returns true.
/// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
/// and gel_is_object() returns false.
class MGCLASS MGPickObjectCB:public MGPickObject{

public:

///////////////Constructor//////////////

MGPickObjectCB():MGPickObject(),m_start_end(-1){;};

///Conversion constructor from MGGelPosition and MGEdge.
MGPickObjectCB(
	MGGelPosition& gelp,
	int start_end
):MGPickObject(gelp),m_start_end(start_end){;};

///Conversion constructor from MGPickObject and start/end.
MGPickObjectCB(
	MGPickObject& pobj,
	int start_end
):MGPickObject(pobj),m_start_end(start_end){;};

//Copy constructor.
//MGPickObjectCB(const MGPickObjectCB& pobj2);

////////////Destructor////////////////
virtual ~MGPickObjectCB(){;};

///Assignment operator.
MGPickObjectCB& operator=(const MGPickObject& pobj);

//////////////////オペレーション//////////////

///Generate a newed clone object.
virtual MGPickObjectCB* clone()const;

///Highlightthe object using the display list of this object.
void hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Return the edge pointer.
int start_end(){return m_start_end;};

///Return the curve of the target.
const MGCurve* curve()const;

///Set the object pointer.
void set_start_end(int start_end){m_start_end=start_end;};

private:

	int m_start_end; ///<0:start, 1:end point.
					///<-1:undefined.
};

/** @} */ // end of MGObjectRelated group
#endif
