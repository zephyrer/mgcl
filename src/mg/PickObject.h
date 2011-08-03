/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/// MGPickObject.h : MGPickObject クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _MGPickObject_HH_
#define _MGPickObject_HH_

#include "mg/MGCL.h"
#include "mg/GelPosition.h"
#include "mg/Position.h"
#include <vector>

class MGGel;
class MGGroup;
class MGObject;
class MGColor;

/** @addtogroup MGObjectRelated
 *  @{
 */

/// MGPickObject is a class to locate where a picked object is in a group
/// hierarchy. Generally, A group includes other groups, and the included groups
/// include other groups. In that way the groups make a group hierachy.
/// MGPickObject represents this hierarcy, an MGObject or hierarchied MGGroup's.
/// When MGPickObject represents an MGObject, gel() returns MGObject
/// pointer and gel_is_object() returns true.
/// When MGPickObject represents an MGGroup, gel() returns MGGroup pointer,
/// and gel_is_object() returns false.
class MGCLASS MGPickObject:public MGGelPosition{

public:

///////////////Constructor//////////////

MGPickObject():MGGelPosition(){;};

///constructor.
MGPickObject(const MGGelPosition& gp2);

///constructor.
MGPickObject(const MGPickObject& obj2):MGGelPosition(obj2){
	m_parameter=obj2.m_parameter;
}

////////////Destructor////////////////
virtual ~MGPickObject(){;};

////////////////Operator overload/////////////////

///Assignment operator.
virtual MGPickObject& operator=(const MGPickObject& pobj);

bool operator<(const MGPickObject& po2)const;
bool operator>(const MGPickObject& po2)const{return po2<(*this);};
bool operator<=(const MGPickObject& po2)const{return !((*this)>po2);};
bool operator>=(const MGPickObject& po2)const{return !(po2>(*this));};

//////////////////オペレーション//////////////

///Generate a newed clone object.
virtual MGPickObject* clone()const;

///Highlightthe object using the display list of this object.
virtual void hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Make the display list of this object as a highlighted one.
virtual void make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Get the parameter value of the object at the picked position.
MGPosition& parameter(){return m_parameter;};
const MGPosition& parameter()const{return m_parameter;};

///Set the object parameter value.
void set_parameter(const MGPosition& param){m_parameter=param;};

private:
	MGPosition m_parameter;///<m_object's parameter value at the picked position.
		///<m_parameter.sdim()=m_object->manifold_dimension().
};

/** @} */ // end of MGObjectRelated group
#endif
