/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPoint_HH_
#define _MGPoint_HH_

#include "mg/Geometry.h"
#include "mg/Position.h"
#include "mg/isects.h"

class MGIfstream;
class MGOfstream;
class MGInterval;
class MGBox;
class MGVector;
class MGPosition_list;
class MGMatrix;
class MGTransf;

/** @addtogroup GEO
 *  @{
 */

///MGPoint represents one dimensional manifold, a point in a space.
///See MGPosition class, whic is recommended to use usually.
class MGCLASS MGPoint:public MGGeometry{

public:

//////////// Constructor ////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGPoint();

///Construct a point from a position.
MGPoint(const MGPosition& P);

///Construct a point by changing the space dimension or ordering
///the space dimension element.
MGPoint(
	size_t sdim,		///<new space dimension.
	const MGPoint& P		///<original point.
	, size_t start1=0		///<start position coordinate of new point.
	, size_t start2=0		///<start position coordinate of the original.
);

//////////// Destructor ////////////
~MGPoint();

//////////// Operator overload(演算子多重定義) ////////////

///Assignment.
///When the leaf object of this and obj2 are not equal, this assignment
///does nothing.
MGPoint& operator=(const MGGel& gel2);
MGPoint& operator=(const MGPoint& gel2);

///Return i-th element of the position.
double operator[] (size_t i) const{return m_point.ref(i);}
double operator() (size_t i) const{return m_point.ref(i);}

///Access to i-th element.
double& operator()(size_t i){return m_point(i);};

///Object transformation.
MGPoint& operator+=(const MGVector& v);
MGPoint& operator-=(const MGVector& v);
MGPoint& operator*=(double scale);
MGPoint& operator*=(const MGMatrix& mat);
MGPoint& operator*=(const MGTransf& tr);

////////////Logical operator overload/////////

///comparison
bool operator==(const MGPoint& point)const;
bool operator<(const MGPoint& gel2)const;
bool operator==(const MGGel& gel2)const{return gel2==(*this);};
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGPoint& gel2)const{return !(gel2==(*this));};
bool operator<(const MGGel& gel2)const{return gel2>(*this);};

////////// Member Function ////////////

///Return minimum box that includes the geometry of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
/*MGBox box(
	const MGBox& bx	/// Parameter Range of the geometry.
					///bx's space dimension is manifold dimension
					///of the geometry.
) const;
*/

///Obtain ceter coordinate of the geometry.
MGPosition center() const{return m_point;};

///Obtain ceter parameter value of the geometry.
MGPosition center_param() const{ return MGPosition();};

///Changing this object's space dimension.
MGPoint& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Construct new geometry object by copying to newed area.
///User must delete this copied object by "delete".
MGPoint* clone() const;

///Construct new geometry object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGPoint* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Compute direction unit vector of the geometry.
///For the point, this is undefined.
MGUnit_vector direction(const MGPosition& param) const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void draw3DVertex()const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector evaluate(
	const MGPosition& t,	///< Parameter value,
				///<t's space dimension is geometry's manifold dimension.
	const size_t* nderiv=0	///<Order of derivative of i-th parameter
				///<in nderiv[i],
				///<When nderiv=null, nderiv[i]=0 is assumed for all i.
)const{return m_point;}

/// Return This object's typeID
long identify_type()const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(const MGPosition& t)const;

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
MGisects intersection(const MGObject& obj2)const{return MGisects();};
MGisects intersection(const MGCurve& obj2)const{return MGisects();};
MGisects intersection(const MGFSurface& obj2)const{return MGisects();};
MGisects intersection(const MGSurface& obj2)const{return MGisects();};
MGisects intersection(const MGFace& obj2)const{return MGisects();};
MGisects intersection(const MGShell& obj2)const{return MGisects();};

///Return manifold dimension, i.e. 0:point, 1:curve, 2:surface.
size_t manifold_dimension() const{return 0;};

///Negate direction of this geometry.
void negate(){;};

///Transform the coordinates of boundary of this geometry so that
///new coordinate of boundary is the same coordinate as the new one of
///this geometry after negate() of this geometry is done.
///That is, boundary coordinates are parameter world of this geometry.
void negate_transform(MGGeometry& boundary)const{;};

///Test if given point is on the geometry or not. If yes, return parameter
///value of the geometry. Even if not, return nearest point's parameter.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the geometry,
/// and 0 if the point is not on the geometry.
bool on(
	const MGPosition& P,///<Point(指定点)
	MGPosition&	param	///<Parameter of the geometry(パラメータ)
) const;

///Compute parameter value of given point.
/// 自身の上の指定点を表すパラメータ値を返す。
/// If input point is not on the geometry, return the nearest point on the
/// geometry.
MGPosition parameter(
	const MGPosition& P	///<Point(指定点)
) const;

///Return parameter range of the geometry(パラメータ範囲を返す)
MGBox parameter_range() const;

///Return point pointer if this MGGel is an MGPoint, else return null.
MGPoint* point(){return this;};
const MGPoint* point()const{return this;};

const MGPosition& position() const{return m_point;}

///Round t into geometry's parameter range.
/// 入力パラメータをパラメータ範囲でまるめて返却する。
MGPosition range(const MGPosition& t) const;

///Return i-th element of the point.
double ref(size_t i) const{return m_point.ref(i);};

///Return space dimension
size_t sdim() const;

///IGES output function.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Output function.
std::ostream& out(std::ostream&) const;

std::string whoami()const{return "Point";};

protected:

///メンバデータを読み出す関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void WriteMembers(MGOfstream& buf) const;

private:
	MGPosition m_point;

///Return minimum box that includes whole of the geometry.
///曲線部分を囲むボックスを返す。
MGBox* compute_box() const;

};

/** @} */ // end of GEO group
#endif
