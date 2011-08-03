/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGObject_HH_
#define _MGObject_HH_

#include "mg/MGCL.h"
#include "mg/Pvector.h"
#include "mg/Plist.h"
#include "mg/Position.h"
#include "mg/isects.h"
#include "mg/AttribedGel.h"

//
//Define MGObject Class.
class MGIfstream;
class MGOfstream;
class MGBox;
class MGVector;
class MGMatrix;
class MGTransf;
class MGObject;
class MGGeometry;
class MGTopology;
class MGPoint;
class MGStraight;
class MGCurve;
class MGFSurface;
class MGSurface;
class MGFace;
class MGShell;

/** @defgroup MGObjectRelated Object Related class
 *  MGObject is top abstract class for MGPoint, MGCurve, and MGSurface.
 *  @{
 */

///MGObject is an abstract class which represents a whole geometry and 
///a topology.
class MGCLASS MGObject:public MGAttribedGel{

public:

////////////Constructor////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGObject();

///Copy constructor.
MGObject(const MGObject& obj2);

///Virtual Destructor
virtual ~MGObject();

///Assignment.
///When the leaf object of this and gel2 are not equal, this assignment
///does nothing.
virtual MGObject& operator=(const MGObject& obj2){return set_object(obj2);};

///Object transformation.
virtual MGObject& operator+=(const MGVector& v)=0;
virtual MGObject& operator-=(const MGVector& v)=0;
virtual MGObject& operator*=(double scale)=0;
virtual MGObject& operator*=(const MGMatrix& mat)=0;
virtual MGObject& operator*=(const MGTransf& tr)=0;

////////////Member Function////////////

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

///Get the MGAppearance pointer of this object. If not defined, null will be
///returned.
///See ensure_appearance().
MGAppearance* appearance(){return m_appearance;};
const MGAppearance* appearance()const{return m_appearance;};

///Get the box of the object.
virtual const MGBox& box() const=0;

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGObject* clone() const=0;

///Delete a display list of this gel.
virtual void delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	///<mgSysGL pointer will be apppended.
)const;

///Draw the object in wire mode, in the world coordinates.
///The object is converted to curve(s) and is drawn.
virtual void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const=0;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void draw3DVertex(
)const=0;

///Shade the object in world coordinates.
virtual void shade(
	double span_length	///Line segment span length.
)const{drawWire(span_length);};

///make this group has appearance and get the MGAppearance pointer.
///See appearance().
MGAppearance* ensure_appearance();

///Make 2 types of display list of this gel(wire and shading).
///Return is the display list name.
virtual size_t make_display_list(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Make a display list without color of this gel.
///Return is the display list name.
virtual size_t make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Test if this and 2nd object has common area about their box(),
///taking error into account.
bool has_common(const MGObject& obj2) const;

///Test if this gel includes an object.
const MGObject* includes_object()const{return this;};
///Test if this gel includes an object.
MGObject* includes_object(){return this;};

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
virtual MGisects intersection(const MGObject& obj2)const=0;
virtual MGisects intersection(const MGPoint& obj2)const{return MGisects();};
virtual MGisects intersection(const MGCurve& obj2)const=0;
virtual MGisects intersection(const MGFSurface& obj2)const=0;
virtual MGisects intersection(const MGSurface& obj2)const=0;
virtual MGisects intersection(const MGFace& obj2)const=0;
virtual MGisects intersection(const MGShell& obj2)const=0;

///Get manifold dimension.
virtual unsigned manifold_dimension() const=0;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
virtual MGPosition pick_closest(const MGStraight& sl)const{
	return MGPosition();
}

///Remove the MGAppearance of this MGAttribedGel.
void remove_appearance();

///Return MGObject pointer if this MGGel is an MGObject, else return null.
virtual MGObject* object(){return this;};
virtual const MGObject* object()const{return this;};

///Get the MGFSurface pointer if this is MGSurface or MGFace.
virtual const MGFSurface* fsurface()const{return (const MGFSurface*)0;};
virtual MGFSurface* fsurface(){return (MGFSurface*)0;};

///Transform the gel by the argument.
virtual void transform(const MGVector& v){(*this)+=v;};///translation
virtual void transform(double scale){(*this)*=scale;};///scaling.
virtual void transform(const MGMatrix& mat){(*this)*=mat;};///matrix transformation.
virtual void transform(const MGTransf& tr){(*this)*=tr;};///general transformation.

protected:

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

///Assignment.
///When the leaf object of this and gel2 are not equal, this assignment
///does nothing.
MGObject& set_object(const MGObject& gel2);

private:

MGAppearance* m_appearance;///<MGAppearance pointer. This is a newed object.

friend class MGIfstream;
friend class MGOfstream;

};

///Construct a null newed MGObject from the type id TID.
MGDECL MGObject* MGNullObj(long TID);

///Draw an object in its parameter space.
///This is valid only for Surface, Face, Loop, Edge.
MGDECL void MGDraw_in_parameter_space(
	const MGObject& obj,///<The object to draw.
	double span_length	///<Line segment span length.
);

/** @} */ // end of MGObjectRelated group
#endif
