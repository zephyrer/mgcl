/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGUnit_vector_HH_
#define _MGUnit_vector_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/Vector.h"

// MGUnit_vector.h 
// Header for class MGUnit_vector

/// Define a unit vector, is a MGVector.
class MGCLASS MGUnit_vector: public MGVector {

public:

/// String stream function
MGDECL friend std::ostream& operator<<(std::ostream&, const MGUnit_vector& );

//
// ////// Constructor.コンストラクタ /////////
//

/// Void constructor.
MGUnit_vector(
	 size_t sdim=3	///< space dimension
);

/// Unit vector from general vector.
/// ベクトルを指定してその単位ベクトルを生成する.
MGUnit_vector(const MGVector& v);

//
// //////// Operator overload.演算子定義 ////////
//

/// Assignment From Vector.
MGUnit_vector& operator= (const MGVector& vec2) ;

///
///Update vector data by array of double.
///The result is unit of the updated vector.
///
MGUnit_vector& operator= (const double*);

///Update the unit vector by adding vec2.
//The result is unit of the vector of two vector addition.
/// 自身のベクトルに与えられたベクトルを加算して自身のベクトルとする.
MGUnit_vector& operator+= (const MGVector& vec2);

///Unary minus. Negate the vector.
///単項マイナス。自身の単位ベクトルを反転したオブジェクトを生成
MGUnit_vector operator- () const;

///Update the unit vector by subtractiong vec2.
///The result is unit of the vector of two vector subtraction.
/// 自身のベクトルから与えられたベクトルを減算し自身のベクトルとする
MGUnit_vector& operator-= (const MGVector& vec2);

///Update own vector by vector product output.
//Cchanges to 3D vector. The result is a unit one of two vector product.
MGUnit_vector& operator*= (const MGVector& vec2);

//////// Member function. ////////

///Compute orthonormal system, given sub(sv) vectors.
///(*this, v1, v2) organizes orthonormal system of 3D.
///If sv.orthogonal(*this), v1=sv.normalize().
void orthonormal(const MGVector& sv
	, MGUnit_vector& v1, MGUnit_vector& v2) const;
	
private:

//プライベートメンバ関数: MGVector に定義されているが、MGUnit_vector 
//の利用者に利用できなくするために、private で定義する。
//The following private functions are defined to prohibit their use
//of inheritting MGVector.

// スカラーの乗算を行い自身のベクトルとする.
MGUnit_vector& operator *= ( double );

// スカラー除算を行い自身のベクトルとする.
MGUnit_vector& operator /= ( double );

};

/** @} */ // end of BASE group
#endif
