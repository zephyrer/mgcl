/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Unit_vector.h"
#include "mg/Point.h"
#include "mg/Box.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implement MGPoint Class.
//MGPoint is an abstract class which represents zero dimensional manifold.
//

//Constructor

//Void constructor(初期化なしでオブジェクトを作成する。)
MGPoint::MGPoint(){;}

//Construct a point from a position.
MGPoint::MGPoint(const MGPosition& P)
:m_point(P){;}

//Construct a point by changing the space dimension or ordering
//the space dimension element.
MGPoint::MGPoint(
	size_t sdim,		//new space dimension.
	const MGPoint& P	//original point.
	, size_t start1		//start position coordinate of new point.
	, size_t start2)	//start position coordinate of the original.
:m_point(sdim,P.m_point,start1,start2){;}

//Virtual Destructor
MGPoint::~MGPoint(){;}

//Operator overload(演算子多重定義)
//Logical operator overload(論理演算子多重定義)
bool MGPoint::operator<(const MGPoint& gel2)const{
	return m_point.len()<gel2.m_point.len();
}
bool MGPoint::operator==(const MGPoint& geo)const{
	return m_point==(geo.m_point);
}

//Member Function

//Return minimum box that includes whole of the geometry.
//曲線部分を囲むボックスを返す。
MGBox* MGPoint::compute_box() const{return new MGBox(m_point, m_point);}

//Return minimum box that includes the geometry of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
/*MGBox MGPoint::box(
	const MGBox& bx	// Parameter Range of the geometry.
					//bx's space dimension is manifold dimension
					//of the geometry.
	) const{return MGGeometry::box();}
*/

//Changing this object's space dimension.
MGPoint& MGPoint::change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2 		// Source order of this object.
){
	m_point=MGPosition(sdim,m_point,start1,start2);
	update_mark();
	return *this;
}

//Construct new geometry object by copying to newed area.
//User must delete this copied object by "delete".
MGPoint* MGPoint::clone() const{return new MGPoint(*this);}

//Construct new geometry object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGPoint* MGPoint::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1,	 	// Destination order of new point.
	size_t start2		// Source order of this point.
)const{
	return new MGPoint(sdim, *this, start1, start2);
}

//Compute direction unit vector of the geometry.
//For the point, this is undefined.
MGUnit_vector MGPoint::direction(
	const MGPosition& param
)const{
	return MGDefault::z_unit_vector();
}

//Test if input parameter value is inside parameter range of the line.
bool MGPoint::in_range(const MGPosition& t) const{return 1;}

//Test if given point is on the geometry or not. If yes, return parameter
//value of the geometry. Even if not, return nearest point's parameter.
// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the geometry,
// and 0 if the point is not on the geometry.
bool MGPoint::on(
	const MGPosition& P,//Point(指定点)
	MGPosition&	param	//Parameter of the geometry(パラメータ)
)const{
	param=MGPosition();
	return MGAZero((P-m_point).len());
}

//Compute parameter value of given point.
// 自身の上の指定点を表すパラメータ値を返す。
// If input point is not on the geometry, return the nearest point on the
// geometry.
MGPosition MGPoint::parameter(
	const MGPosition& P	//Point(指定点)
)const{
	return MGPosition();
}

//Return parameter range of the geometry(パラメータ範囲を返す)
MGBox MGPoint::parameter_range()const{
	return MGBox();
}

//Round t into geometry's parameter range.
// 入力パラメータをパラメータ範囲でまるめて返却する。
MGPosition MGPoint::range(const MGPosition& t)const{
	return MGPosition();
}

//Return space dimension
size_t MGPoint::sdim() const{return m_point.sdim();}

//Assignment.
//When the leaf object of this and obj2 are not equal, this assignment
//does nothing.
MGPoint& MGPoint::operator=(const MGPoint& gel2){
	if(this==&gel2)
		return *this;

	MGGeometry::operator=(gel2);
	m_point=gel2.m_point;
	return *this;
}
MGPoint& MGPoint::operator=(const MGGel& gel2){
	const MGPoint* gel2_is_this=dynamic_cast<const MGPoint*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Update the geometry by translation.
// 与ベクトルだけ曲線を平行移動して自身とする。
MGPoint& MGPoint::operator+=(const MGVector& v){
	m_point+=v;
	return *this;
}
MGPoint& MGPoint::operator-=(const MGVector& v){
	m_point-=v;
	return *this;
}

//Update the geometry by multiplying scale.
// 与えられたスケールを曲線にかける。
MGPoint& MGPoint::operator*=(double scale){
	m_point*=scale;
	return *this;
}

//Update the geometry by transformation of matrix.
// 与えられた変換で直線の変換を行い自身の直線とする。
MGPoint& MGPoint::operator*=(const MGMatrix& mat){
	m_point*=mat;
	return *this;
}

//Update the geometry by transformation of transf.
// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGPoint& MGPoint::operator*=(const MGTransf& tr){
	m_point*=tr;
	return *this;
}

