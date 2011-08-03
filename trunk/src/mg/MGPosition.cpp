/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Transf.h"
#include "mg/Position.h"
#include "mg/Point.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGPosition.cpp
// Implemetation of class MGPosition
//
// Represent a positional data. MGPosition is the same class as MGposition,
// only is an alias name.

// Constructor
// Conversion Constructor from a point
MGPosition::MGPosition(const MGPoint& point)
:m_element(point.position()){;}

//Copy constructor
//	MGPosition ( const MGPosition& ); We use default copy constructor.

// Destructor
//	~MGPosition(); We use default destructor.

//Member Function

//Let this be the center of the rotation, then compute the angle rotated 
//around the normal from start to end.
//Although normal is assumed to be parallel to N2=V1*v2, normal may not perpendicular
//to v1 and v2, in which case, the projected normal to N2 is used to measure the angle.
//Here v1=start-*this, and v2=end-*this.
//this->angle(start,end,normal)+this->angle(end,start,normal)=2*pai always holds.
double MGPosition::angle(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& normal
)const{
	MGVector V1(start-*this), V2(end-*this);
	return V1.angle2pai(V2,normal);
}

//Clear all the elements by the value init.
MGPosition& MGPosition::clear(double init){
	m_element.clear(init);
	return *this;
}

//Return the distance of this and P2.
double MGPosition::distance(const MGPosition& P2)const{
	return (m_element-P2.m_element).len();
}

// Generate a Position by interpolating two Position. Input scalar is
// a ratio and when zero, output position is a copy of the own vector.
MGPosition MGPosition::interpolate(double t2, const MGPosition& vec2) const{
	return MGPosition(m_element.interpolate(t2,vec2.m_element));
}

//Store vec2 data into *this.
//Store length is minimum of len() and vec2.len().
//Storing will be done rap-around. That is, if id i or j reached to
//each sdim(), the id will be changed to 0.
void MGPosition::store_at(
	size_t i,				//Displacement of *this.
	const MGVector& vec2,	//Vector 2.
	size_t j)				//Displacement of vec2.
{	m_element.store_at(i,vec2,j);}

//Store vec2 data into *this.
//Storing will be done rap-around. That is, if id i or j reached to
//each sdim(), the id will be changed to 0.
void MGPosition::store_at(
	size_t i,				//Displacement of *this.
	const MGVector& vec2,	//Vector 2.
	size_t j,				//Displacement of vec2.
	size_t len)				//Length to store 
{	m_element.store_at(i,vec2,j,len);}

//Operator Oveload

//Assignment
//	MGposition& operator =(const MGposition &); We use default assignment;
	
//Update position data by array of double.
MGPosition& MGPosition::operator=(const double* a){
	m_element = a;
	return *this;
}

// 自身のPositionに与えられたVectorを加算して自身のPositionとする 
MGPosition& MGPosition::operator+= (const MGVector& vec)
{	m_element+=vec;	return *this;}
MGPosition& MGPosition::operator+= (const MGPosition& pos)
{	m_element+=pos.m_element;	return *this;}

// 単項マイナス。自身のPositionを反転し、Positionを生成
MGPosition MGPosition::operator- () const{
	return MGPosition(-m_element);
}

// 自身のPositionから与えられたVectorを減算し自身のPositionとする
MGPosition& MGPosition::operator-= (const MGVector& vec){
	m_element -= vec;
	return *this;
}

// Scalarの乗算を行い自身のPositionとする
MGPosition& MGPosition::operator*= (double a){
	m_element*= a;
	return *this;
}
	
// MatrixによるPositionの変換を行いPositionを生成
MGPosition operator*(const MGPosition& p1,const MGMatrix& mat){
	MGPosition temp(p1);
	return temp*= mat;
}

// MatrixによるPositionの変換を行い自身のPositionとする
MGPosition& MGPosition::operator*= (const MGMatrix& mat){
	m_element *= mat;
	return *this;
}

// PositionのTransformを行いVectorを生成
MGPosition operator*(const MGPosition& p1,const MGTransf& tr){
	return p1.m_element*tr.affine()+tr.translation();
}

// PositionのTransformを行いPositionを生成して，
// 自身のPositionとする
MGPosition& MGPosition::operator*= (const MGTransf& tran){	
	m_element*= tran.affine();
	m_element+= tran.translation();
	return *this;
}

// Scalar除算を行い自身のPositionとする
MGPosition& MGPosition::operator/= (double a){
	m_element/= a;
	return *this;
}

// 与えられたPositionの成分の値を比較し、同じであれば TRUE を返却
bool operator==(const MGPosition& p1,const MGPosition& p2){
	MGVector dif=p2-p1;
	return dif%dif<=MGTolerance::wc_zero_sqr();
}
bool operator==(const MGVector& p1,const MGPosition& p2){
	MGVector dif=p2.m_element-p1;
	return dif%dif<=MGTolerance::wc_zero_sqr();
}
bool operator==(const MGPosition& p1,const MGVector& p2){
	MGVector dif=p2-p1.m_element;
	return dif%dif<=MGTolerance::wc_zero_sqr();
}

//Friend Function

//Test if P1, P2, and P3 are on a single straight line.
//Function's return value is true if the three points are on a straight,
//false if not.
bool is_collinear(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
){return P1.is_collinear(P2,P3);}

