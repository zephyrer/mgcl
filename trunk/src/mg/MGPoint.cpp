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

//Void constructor(�������Ȃ��ŃI�u�W�F�N�g���쐬����B)
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

//Operator overload(���Z�q���d��`)
//Logical operator overload(�_�����Z�q���d��`)
bool MGPoint::operator<(const MGPoint& gel2)const{
	return m_point.len()<gel2.m_point.len();
}
bool MGPoint::operator==(const MGPoint& geo)const{
	return m_point==(geo.m_point);
}

//Member Function

//Return minimum box that includes whole of the geometry.
//�Ȑ��������͂ރ{�b�N�X��Ԃ��B
MGBox* MGPoint::compute_box() const{return new MGBox(m_point, m_point);}

//Return minimum box that includes the geometry of parameter interval.
// ���͂̃p�����[�^�͈͂̋Ȑ��������͂ރ{�b�N�X��Ԃ��B
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
// �w��_�����g��ɂ��邩�𒲂ׂ�B�Ȑ���ɂ���΁C���̃p�����[�^�[�l���C
// �Ȃ��Ă��ŋߖT�_�̃p�����[�^�l��Ԃ��B
// Function's return value is >0 if the point is on the geometry,
// and 0 if the point is not on the geometry.
bool MGPoint::on(
	const MGPosition& P,//Point(�w��_)
	MGPosition&	param	//Parameter of the geometry(�p�����[�^)
)const{
	param=MGPosition();
	return MGAZero((P-m_point).len());
}

//Compute parameter value of given point.
// ���g�̏�̎w��_��\���p�����[�^�l��Ԃ��B
// If input point is not on the geometry, return the nearest point on the
// geometry.
MGPosition MGPoint::parameter(
	const MGPosition& P	//Point(�w��_)
)const{
	return MGPosition();
}

//Return parameter range of the geometry(�p�����[�^�͈͂�Ԃ�)
MGBox MGPoint::parameter_range()const{
	return MGBox();
}

//Round t into geometry's parameter range.
// ���̓p�����[�^���p�����[�^�͈͂ł܂�߂ĕԋp����B
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
// �^�x�N�g�������Ȑ��𕽍s�ړ����Ď��g�Ƃ���B
MGPoint& MGPoint::operator+=(const MGVector& v){
	m_point+=v;
	return *this;
}
MGPoint& MGPoint::operator-=(const MGVector& v){
	m_point-=v;
	return *this;
}

//Update the geometry by multiplying scale.
// �^����ꂽ�X�P�[�����Ȑ��ɂ�����B
MGPoint& MGPoint::operator*=(double scale){
	m_point*=scale;
	return *this;
}

//Update the geometry by transformation of matrix.
// �^����ꂽ�ϊ��Œ����̕ϊ����s�����g�̒����Ƃ���B
MGPoint& MGPoint::operator*=(const MGMatrix& mat){
	m_point*=mat;
	return *this;
}

//Update the geometry by transformation of transf.
// �^����ꂽ�ϊ��ŋȐ��̃g�����X�t�H�[�����s�����g�Ƃ���B
MGPoint& MGPoint::operator*=(const MGTransf& tr){
	m_point*=tr;
	return *this;
}

