/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Vector.h"
#include "mg/Unit_vector.h"
#include "mg/Position.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Object.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGTransf.cc
//Implementation of MGTransf.
//

//
// Constructor

//  void constructor
MGTransf::MGTransf(size_t dim)
	:m_affine(dim), m_translation(dim){;}

//  �S�ẴR���|�[�l���g���w�肵��Transf�𐶐��B
MGTransf::MGTransf(const MGMatrix& mat,
                   const MGVector& vec)
	:m_affine(mat),m_translation(vec)
{
	size_t dim1=mat.sdim(); size_t dim2=vec.sdim();
	if(dim1>dim2) m_translation=MGVector(dim1,vec);
	else if(dim1<dim2) m_affine=MGMatrix(dim2,mat);
}

//  �e�������œ����� Scaling �̂��߂̂�Transf�����B
MGTransf::MGTransf(double scal,
				   const MGVector& vec)
	:m_affine(vec.sdim(),scal), m_translation(vec){;}

//Construct Transf by copying old Transf, changing space dimension,
//and ordering of coordinates.
MGTransf::MGTransf(size_t dim, const MGTransf& transf,
				   size_t start1, size_t start2)
	:m_affine(dim,transf.affine(),start1,start2)
	,m_translation(dim,transf.translation(),start1,start2){;}

//  �e������( x, y )�ňقȂ� Scaling ��2D Transf�𐶐��B
MGTransf::MGTransf(double scalex, double scaley)
	:m_affine(2,scalex), m_translation(0.0, 0.0)
{
	m_affine(1,1)=scaley;
}

// Construct 2D space Transf to transform data for 'origin' to be origin
// and for 'unit' to be x-coordimate.
MGTransf::MGTransf(const MGUnit_vector& unit, //unit vector to be x-coordinate
                   const MGPosition& origin)  //origin to be origin.
{
//	assert(unit.sdim()==2);
    m_affine.set_x_axis(unit);
	m_translation=-origin*m_affine;
}

//  �e������( x,y,z )�ňقȂ� Scaling �̂���Transf�𐶐��B
MGTransf::MGTransf(double scalex, double scaley, double scalez)
	:m_affine(3,scalex)
	,m_translation(0.0, 0.0, 0.0)
{
	m_affine(1,1)=scaley;
	m_affine(2,2)=scalez;
}

//  �^����ꂽ�_�����_�ɂ���悤�ɕ��s�ړ����A�^����ꂽ�Q�̒P��Vector
//  ���e�X�AX���AY���ɂ���悤��]������3D Transf�𐶐�����B�����A
//  �Q�ڂ̒P��Vector���P�ڂ̒P��Vector�ƒ������Ȃ��ꍇ�́A�Q��
//  Vector�̂����ɗ�Vector���܂ޕ��ʓ��Œ�������悤�ϊ�����Vector��
//  �g�p����B
MGTransf::MGTransf ( 
    const MGUnit_vector& uvecx,
    const MGUnit_vector& uvecy,
    const MGPosition& origin)
{
//	assert(uvecx.sdim()<=3 && uvecy.sdim()<=3);
    m_affine.set_xy_axis(uvecx,uvecy);
	m_translation=-origin*m_affine;
}

//Transf to transform the line segment (P0, P1) to the line segment(Q0, Q1)
//P0 is transformed to Q0 and P1 is transformed to Q1.
//Space dimension of each points can be any number more than 1.
MGTransf::MGTransf(
	const MGPosition& P0, const MGPosition& P1,
    const MGPosition& Q0, const MGPosition& Q1){
	m_affine.set_rotate(MGVector(P1,P0), MGVector(Q1,Q0));
	m_translation=-P0*m_affine+Q0;
}

//
//  Member Function
//  
//  Update

//Convert this transf matrix to OpenGL Matrix.
void MGTransf::convert_to_glMatrix(
	double glMat[16]	//OpenGL Matrix will be output.
)const{
	m_affine.convert_to_glMatrix(glMat);
	glMat[12]=m_translation[0];
	glMat[13]=m_translation[1];
	glMat[14]=m_translation[2];
}

// Construct a transf to do the transform of matrix around input point 
// instead of origin of matrix, and replace own transf with it.
MGTransf& MGTransf::set_matrix(
		const MGMatrix& mat,		//Matrix
		const MGPosition& point)	//Point 
{
    m_affine=mat;
	m_translation=-point*mat+point;
	return *this;
}
   
//Set this as a null transf.
void MGTransf::set_null(){
	m_affine.set_null();
	m_translation.set_null();
}

//  �w��_(origin)��ʂ�A�^����ꂽVector(vec)�Ɋւ��ċ��ʕϊ�����
//  2D Transform�𐶐���������Transform�Ɠ��ꊷ����B�w��_
//  �ɉ����w�肳��Ȃ��ꍇ�A���_���g�p����B
MGTransf& MGTransf::set_reflect_2D(const MGVector& vec,
                              const MGPosition& origin)
{
    MGMatrix mat; mat.set_reflect_2D(vec);
	return set_matrix(mat,origin);
}

//  �w��_(origin)�𒆐S�ɔ����v���Ɏw��p�x(angle)��]������
//  2D Transf���쐬���āC���g��Transf�ƒu��������B
MGTransf& MGTransf::set_rotate_2D(double angle, const MGPosition& origin)
{
    MGMatrix mat; mat.set_rotate_2D(angle);
	return set_matrix(mat,origin);
}

//  �w��_��ʂ�A�^����ꂽVector�ɐ����ȕ��ʂɊւ��ċ��ʕϊ�����
// 3D Transf�𐶐���������Transf�Ɠ��ꊷ����B�w��_�ȗ����͌��_
//  �Ƃ���B
MGTransf& MGTransf::set_reflect_3D( 
         const MGVector & vec,
         const MGPosition & origin)
{
    MGMatrix mat; mat.set_reflect_3D(vec);
	return set_matrix(mat,origin);
}

//  �w��_��ʂ�w�����Vector���������̉��Ɏw��p�x��]����
//  Transf�𐶐���������Transf�Ɠ��ꊷ����B�w��_
//  �ɉ����w�肳��Ȃ��ꍇ�A���_���g�p����B
MGTransf& MGTransf::set_rotate_3D( 
         const MGVector & vec,
         double angle,                                      
         const MGPosition & origin)
{ 
    MGMatrix mat; mat.set_rotate_3D(vec,angle);
	return set_matrix(mat,origin);
}

//  �^����ꂽ�g��A�k����Transf��������Transf��
//  ���ꊷ����B
MGTransf& MGTransf::set_scale(double scale)
{
     m_affine.set_scale(scale);
     m_translation.clear();
	 return *this;
}

//  �e�������ňقȂ�g��A�k����Transf��������Transf
//  �Ɠ��ꊷ����B
MGTransf& MGTransf::set_diff_scale(double* scale)
{ 
     m_affine.set_diff_scale(scale);
     m_translation.clear();
	 return *this;
}

//Set up this transf from OpenGL matrix.
MGTransf& MGTransf::set_glMatrix(const double glMat[16]){
	m_affine.set_glMatrix(glMat);
	m_translation.resize(3);
	m_translation(0)=glMat[12];
	m_translation(1)=glMat[13];
	m_translation(2)=glMat[14];
	return *this;
}

//  �Q��

//Reference (i,j)-th element.
double MGTransf::ref(size_t i, size_t j) const{
	size_t dim=sdim();
	if(i<dim && j<dim)       return m_affine.ref(i,j);
	else if(i==dim && j<dim) return m_translation.ref(j);
	else if(i!=j)            return 0.;
	else                     return 1.;
}

//Change the space dimension of this MGTransf to the new sdim.
void MGTransf::resize(size_t sdim){
	m_affine.resize(sdim);
	m_translation.resize(sdim);
}

//Access to (i,j)-th element.
double& MGTransf::operator() (size_t i, size_t j){
	assert(i<=sdim() && j<sdim());
	if(i<sdim()) return m_affine(i,j);
	else         return m_translation(j);
}

//
// ���Z�q�̑��d��`
//

//Functor to apply this transform to the object
//Function's return is the reference to object.
MGObject& MGTransf::operator()(MGObject& object){
	return object*=(*this);
}
MGObject* MGTransf::operator()(MGObject* object){
	(*object)*=(*this);
	return object;
}

// ���g��Transf��vector�̉��Z���s���I�u�W�F�N�g�𐶐�����B
//Translation of Transf.
MGTransf MGTransf::operator+ (const MGVector& vec) const{
	return MGTransf(m_affine,m_translation+vec);
}

// ���g��Transf��vector�̉��Z���s�����g��Transf�Ƃ���B
//Translation of Transf.
MGTransf& MGTransf::operator+= (const MGVector& vec){
	m_translation+=vec;
	return *this;
}

// ���g��Transf��vector�̌��Z���s���I�u�W�F�N�g�𐶐�����B
//Translation of Transf.
MGTransf MGTransf::operator- (const MGVector& vec) const{
	return MGTransf(m_affine,m_translation-vec);
}

// ���g��Transf��vector�̌��Z���s�����g��Transf�Ƃ���B
//Translation of Transf.
MGTransf& MGTransf::operator-= (const MGVector& vec){
	m_translation-=vec;
	return *this;
}

// ���g��Transf��scale�̏�Z���s���I�u�W�F�N�g�𐶐�����B
//Scaling of the transf.
MGTransf MGTransf::operator* (double scale) const{
         MGTransf tr = *this;
         tr *= scale;
         return tr;
}

// ���g��Transf��scale�̏�Z���s���I�u�W�F�N�g�𐶐�����B
//Scaling of the transf.
MGTransf operator* (double scale, const MGTransf& tr){
	return tr*scale;
}

//���g��Transf�Ɨ^����ꂽMatrix�̏�Z���s���I�u�W�F�N�g�𐶐�����B
MGTransf MGTransf::operator* (const MGMatrix& mat) const{
         MGTransf tran = *this;
         tran *= mat;
         return tran;
}

//���g��Transf�Ɨ^����ꂽTransf�̏�Z���s���I�u�W�F�N�g�𐶐�����B
MGTransf MGTransf::operator* (const MGTransf& tran1) const{
         MGTransf tran2 = *this;
         tran2 *= tran1;
         return tran2;
}

 // �x�N�g���ƃg�����X�t�H�[���̏�Z���s��Transf�𐶐�
 //Genarate a Transf by multiplication of a vector and transf.
MGTransf operator* (const MGVector& v, const MGTransf& tr){
	return MGTransf(tr.affine(),v*tr.affine()+tr.translation());
}

//���g��Transf��scale�̏�Z���s�����g��Transf�Ƃ���B
//Scaling of the transf.
MGTransf& MGTransf::operator*= (double scale){
         m_affine *= scale;
         m_translation *= scale;
         return *this;
}

//���g��Transf�Ɨ^����ꂽMatrix�̏�Z���s���I�u�W�F�N�g�𐶐�����B
MGTransf& MGTransf::operator*= (const MGMatrix& mat){
         m_affine *= mat;
         m_translation *= mat;
         return *this;
}

//  ���g��Transf�Ɨ^����ꂽTransf�̏�Z���s��
//  ���g��Transf�Ƃ���B
MGTransf& MGTransf::operator*= (const MGTransf& tran1){
         m_affine *= tran1.m_affine;
         m_translation *= tran1.m_affine;
         m_translation += tran1.m_translation;
         return *this;
}

//  Boolean ���Z
//  ���g��Transf�Ɨ^����ꂽTransf�����������ǂ���
//  ��r���s���B
bool MGTransf::operator== (const MGTransf& tran1) const{
         return (m_affine==tran1.m_affine &&
		    m_translation==tran1.m_translation );
}          

bool MGTransf::operator!= (const MGTransf& tran1) const{
         return !((*this)==tran1);
}    
