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
// ////// Constructor.�R���X�g���N�^ /////////
//

/// Void constructor.
MGUnit_vector(
	 size_t sdim=3	///< space dimension
);

/// Unit vector from general vector.
/// �x�N�g�����w�肵�Ă��̒P�ʃx�N�g���𐶐�����.
MGUnit_vector(const MGVector& v);

//
// //////// Operator overload.���Z�q��` ////////
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
/// ���g�̃x�N�g���ɗ^����ꂽ�x�N�g�������Z���Ď��g�̃x�N�g���Ƃ���.
MGUnit_vector& operator+= (const MGVector& vec2);

///Unary minus. Negate the vector.
///�P���}�C�i�X�B���g�̒P�ʃx�N�g���𔽓]�����I�u�W�F�N�g�𐶐�
MGUnit_vector operator- () const;

///Update the unit vector by subtractiong vec2.
///The result is unit of the vector of two vector subtraction.
/// ���g�̃x�N�g������^����ꂽ�x�N�g�������Z�����g�̃x�N�g���Ƃ���
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

//�v���C�x�[�g�����o�֐�: MGVector �ɒ�`����Ă��邪�AMGUnit_vector 
//�̗��p�҂ɗ��p�ł��Ȃ����邽�߂ɁAprivate �Œ�`����B
//The following private functions are defined to prohibit their use
//of inheritting MGVector.

// �X�J���[�̏�Z���s�����g�̃x�N�g���Ƃ���.
MGUnit_vector& operator *= ( double );

// �X�J���[���Z���s�����g�̃x�N�g���Ƃ���.
MGUnit_vector& operator /= ( double );

};

/** @} */ // end of BASE group
#endif
