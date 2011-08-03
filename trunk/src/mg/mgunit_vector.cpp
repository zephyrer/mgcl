/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Unit_vector.h"
#include "mg/Tolerance.h"
#include "mg/Default.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGUnit_vector.cc
// Implementation of MGUnit_vector

//
//Constructor
// 
// void コンストラクタ
MGUnit_vector::MGUnit_vector(size_t sdim):MGVector(0.0, 0.0, 1.0){
	assert(sdim>0);
	if(sdim!=3){
		resize(sdim);
		size_t sdimm1=sdim-1;
		for(size_t i=2; i<sdimm1; i++) m_element[i]=0.;
		m_element[sdimm1]=1.;
	}
	m_length=1.;	//m_element[3]=1.;
}

// ベクトルを指定してその単位ベクトルを生成する
MGUnit_vector::MGUnit_vector(const MGVector& vec):MGVector(vec.sdim()){
	size_t dim=sdim();
	if(!dim) (*this)=MGUnit_vector();
	else{
		size_t i;
		double length = vec.len();

		// 自身のベクトルが零ベクトルの時はデフォルトベクトルを生成
		if(MGMZero(length)){
			for(i=0; i<dim-1; i++) m_element[i]=0.;
			m_element[dim-1]=1.;
		}
		// 零ベクトル以外は与えられたベクトルの成分を長さで割る
		else {
			for(i=0; i<dim; i++) m_element[i]=vec.ref(i)/length;
		}
		m_length=1.;	//m_element[3]=1.;
	}
}

// Compute orthonormal system of the unit vestor, given sub(sv) vectors.
// (*this, v1, v2) organizes orthonormal system of 3D.
void MGUnit_vector::orthonormal(
	const MGVector& sv,
	MGUnit_vector& v1,
	MGUnit_vector& v2
)const{
	MGVector v(3,(*this)*sv);
	if(MGMZero(v.len())){
		if(*this==mgZ_UVEC) v1=MGVector(1.,0.);
		else if(*this==-mgZ_UVEC) v1=MGVector(0.,1.);
		else{
			double dx=fabs(ref(0)), dy=fabs(ref(1)), dz=fabs(ref(2));
			if(MGRZero(dz)) v1=MGVector(-ref(1),ref(0));
			else if(dx>dz){
				if(dx>dy) v1=MGVector(0.,1.,0.);
				else      v1=MGVector(0.,0.,1.);
			}else if(dy>dz) v1=MGVector(0.,0.,1.);
		}
		// Normalize m and n.
		v2=(*this)*v1; v1=v2*(*this);	
	}else{
		v2=v; v1=v2*(*this);
	}
}

//
// 演算子定義
//

// Assignment From Vector.
MGUnit_vector& MGUnit_vector::operator= (const MGVector& vec2){
	if(vec2.is_unit_vector()) MGVector::operator= (vec2);
	else MGVector::operator= (vec2.normalize());
	return *this;
}
	
//Update vector data by array of double.
//The result is unit of the updated vector.
MGUnit_vector& MGUnit_vector::operator=(const double* array){
	MGVector vec(sdim(),array);
	return *this=vec;
}

//Update the unit vector by adding vec2. The result is unit of the vector
//of two vector addition.
MGUnit_vector& MGUnit_vector::operator+= (const MGVector& vec2){
	return *this=(*this+vec2);
}

// 単項マイナス。自身の単位ベクトルを反転したオブジェクトを生成する
MGUnit_vector MGUnit_vector::operator- () const{
	size_t dim=sdim();
	MGUnit_vector temp(*this);
	for(size_t i=0; i<dim; i++) temp.m_element[i]=-m_element[i];
	return temp;
}

//Update the unit vector by subtractiong vec2. The result is unit of
//the vector of two vector subtraction.
MGUnit_vector& MGUnit_vector::operator-= (const MGVector& vec2){
	return *this=(*this-vec2);
}
 
//Update own vector by vector product output, changes to 3D vector.
//The result is unit of two vector product.
MGUnit_vector& MGUnit_vector::operator*= (const MGVector& vec2){
	return *this=(*this*vec2);
}
