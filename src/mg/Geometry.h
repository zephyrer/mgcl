/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGeometry_HH_
#define _MGGeometry_HH_

#include "mg/Object.h"
#include "mg/Box.h"
#include "mg/Unit_vector.h"

//Define MGGeometry Class.
class MGIfstream;
class MGOfstream;
class MGInterval;
class MGBox;
class MGVector;
class MGPosition;
class MGPosition_list;
class MGMatrix;
class MGTransf;
class MGCurve;
class MGSurface;
class MGPoint;

/** @defgroup GEO Geometry (sub) classes
 *  MGGeometry is top abstract class for MGPoint, MGCurve, and MGSurface.
 *  @{
 */

///MGGeometry is an abstract class which represents a whole geometry.
///Geometry represents point(MGPoint), curve(MGCurve), and surface
///(MGSurface).
class MGCLASS MGGeometry: public MGObject{

public:

////////////Constructor////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGGeometry();

///Copy constructor.
MGGeometry(const MGGeometry& geo2);

///Virtual Destructor
virtual ~MGGeometry();

////////////Operator overload(演算子多重定義)////////////

///Assignment.
///When the leaf object of this and obj2 are not equal, this assignment
///does nothing.
virtual MGGeometry& operator=(const MGGeometry& gel2){return set_geometry(gel2);};

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

///Return MGGeometry pointer if this MGGel is an MGGeometry, else return null.
virtual MGGeometry* geometry(){return this;};
virtual const MGGeometry* geometry()const{return this;};

////////////Member Function////////////

///Return minimum box that includes whole of the geometry.
const MGBox& box() const;

///Obtain ceter coordinate of the geometry.
virtual MGPosition center() const=0;

///Obtain ceter parameter value of the geometry.
virtual MGPosition center_param() const=0;

///Changing this object's space dimension.
virtual MGGeometry& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
)=0;

///Construct new geometry object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGGeometry* clone()const=0;

///Construct new geometry object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
virtual MGGeometry* copy_change_dimension(
	size_t sdim,		///< new space dimension
	size_t start1=0, 	///< Destination order of new object.
	size_t start2=0		///< Source order of this object.
)const=0;

///Compute direction unit vector of the geometry.
virtual MGUnit_vector direction(const MGPosition& param) const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void draw3DVertex()const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
virtual MGVector evaluate(
	const MGPosition& t,	///< Parameter value,
				///<t's space dimension is geometry's manifold dimension.
	const size_t* nderiv=0	///<Order of derivative of i-th parameter in nderiv[i],
				///<When nderiv=null, nderiv[i]=0 is assumed for all i.
) const=0;

///Test if input parameter value t is inside parameter range of the line.
///t's space dimension is geometry's manifold dimension.
virtual bool in_range(const MGPosition& t)const{return false;}

///Test if this is null.
bool is_null()const{return sdim()==0;};

///Negate direction of this geometry.
virtual void negate(){return;};

///Transform the coordinates of boundary of this geometry so that
///new coordinate of boundary is the same coordinate as the new one of
///this geometry after negate() of this geometry is done.
///That is, boundary coordinates are of parameter world of this geometry.
virtual void negate_transform(MGGeometry& boundary)const=0;

///Error allowed in the parameter space of the geometry.
double parameter_error() const;

///Return parameter range of the geometry(パラメータ範囲を返す)
///Returned box's space dimension is geometry's manifold dimension.
virtual MGBox parameter_range() const=0;

///Return space dimension
virtual size_t sdim() const=0;

virtual std::string whoami()const{return "Geometry";};

protected:

mutable MGBox* m_box;

///Compute box of the geometry.
///Returned is a newed MGBox pointer. MGGeometry takes the ownership.
virtual MGBox* compute_box() const=0;

/// ファイル入力関数
/// 戻り値はオブジェクトのPID(識別ID)
virtual void ReadMembers(MGIfstream& buf);

///Assignment.
MGGeometry& set_geometry(const MGGeometry& geo2);

///Mark this as updated.
virtual void update_mark(){if(m_box){delete m_box;m_box=0;}};
	
/// ファイル出力関数
/// 戻り値はオブジェクトのPID(識別ID)
virtual void WriteMembers(MGOfstream& buf)const;

private:

};

/** @} */ // end of GEO group
#endif
