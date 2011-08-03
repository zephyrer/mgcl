/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTransf_HH_
#define _MGTransf_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/Default.h"
#include "mg/Vector.h"
#include "mg/Matrix.h"

//  Defines class MGTransf

//  Forward Declarations
class MGUnit_vector;
class MGPosition;
class MGIfstream;
class MGOfstream;
class MGIgesOfstream;
class MGObject;

///MGTransf represents a transformation of a space dimension.
///While MGMatrix is for the transformation about the origin(0,0,0),
///MGTransf is a general transformation.
///Transformation consists of a matrix and a vector. The matrix defines
///an affine transformation and the vector defines a translation.
///Let T be a MGTransf and A be an object. Then MGTransf transformation of the
///object A is defined as:  A*T. (Not T*A)
class MGCLASS MGTransf {

public:

/// ���g��Transf��scale�̏�Z���s���I�u�W�F�N�g�𐶐�����B
///Scaling of the transf.
MGDECL friend MGTransf operator* (double scale, const MGTransf&);

/// �x�N�g���ƃg�����X�t�H�[���̏�Z���s��Transf�𐶐�
///Genarate a Transf by multiplication of a vector and transf.
MGDECL friend MGTransf operator* (const MGVector& v, const MGTransf& tr);

///String stream function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGTransf&);

//////////// Constructor //////////

/// General Space dimension Constructor.
/// Void Constructor
explicit MGTransf(size_t dim=0);

///  �S�ẴR���|�[�l���g���w�肵��Transf�𐶐��B
///Construct Transf from a matrix and a vector.
///***** This is the fundamental constructor.*****
MGTransf(const MGMatrix&, const MGVector&);

/// �e�������œ����� Scaling ��Transf�����B
///Transf from same scaling factor for all coordinates and a translation.
MGTransf(double scale, const MGVector&);

///Construct Transf by copying old Transf, changing space dimension,
///and ordering of coordinates.
MGTransf (size_t sdim, const MGTransf& transf,
			size_t start1=0, size_t start2=0);

///////// 2D Constructor.////////////

/// �e������(x,y)�ňقȂ� Scaling ��2D Transf�𐶐��B
///2D Transf from different scaling factors for x and y coordinates.
///No translation.
MGTransf(double scalex, double scaley);

/// Construct 2D space Transf to transform data for 'origin' to be
///origin(0,0,0) and for 'unit' to be x-coordimate.
MGTransf(const MGUnit_vector& unit,  ///< unit vector to be x-coordinate.
         const MGPosition& origin	 ///<origin to be origin.
);

//////////// 3D Constructor.////////////

/// �e������( x,y,z )�ňقȂ� Scaling ��3D Transf�𐶐��B
///3D Transf from different scaling factors for x, y and z coordinates.
///No translation.
MGTransf(double scalex, double scaley, double scalez);

///  �^����ꂽ�_�����_�ɂ���悤�ɕ��s�ړ����A�^����ꂽ�Q�̒P��vector
///  ���e�X�AX���AY���ɂ���悤��]������ 3D Transf�𐶐�����B�����A
///  �Q�ڂ̒P��vector���P�ڂ̒P��vector�ƒ������Ȃ��ꍇ�́A�Q��
///  �x�N�g���̂����ɗ�vector���܂ޕ��ʓ��Œ�������悤�ϊ�����vector��
///  �g�p����B
///3D Transf to transform P to be origin(0,0,0), uvecx to be x-axis
///and uvecy to be y-axis.
///If uvecx and uvecy does not cross at right angle, uvecy will be
///transformed.
MGTransf(const MGUnit_vector& uvecx, ///<1st unit vector.
         const MGUnit_vector& uvecy, ///<2nd unit vector
         const MGPosition& P		///<origin
 );	

///Transf to transform the line segment (P0, P1) to the line segment(Q0, Q1)
///P0 is transformed to Q0 and P1 is transformed to Q1.
///Space dimension of each points can be any number more than 1.
MGTransf(
	const MGPosition& P0, const MGPosition& P1,
    const MGPosition& Q0, const MGPosition& Q1);

//////////// Operator overload. /////////

///Reference (i,j)-th element.
double operator() (size_t i, size_t j) const{return ref(i,j);};

///Access to (i,j)-th element.
double& operator() (size_t i, size_t j) ;

///Functor to apply this transform to the object
///Function's return is the reference to object.
MGObject& operator()(MGObject& object);
MGObject* operator()(MGObject* object);

/// ���g��Transf��vector�̉��Z���s���I�u�W�F�N�g�𐶐�����B
///Translation of Transf.
MGTransf operator+ (const MGVector&) const;

/// ���g��Transf��vector�̉��Z���s�����g��Transf�Ƃ���B
///Translation of Transf.
MGTransf& operator+= (const MGVector&);

/// ���g��Transf��vector�̌��Z���s���I�u�W�F�N�g�𐶐�����B
///Translation of Transf.
MGTransf operator- (const MGVector&) const;

/// ���g��Transf��vector�̌��Z���s�����g��Transf�Ƃ���B
///Translation of Transf.
MGTransf& operator-= (const MGVector&);

/// ���g��Transf��scale�̏�Z���s���I�u�W�F�N�g�𐶐�����B
///Scaling of the transf.
MGTransf operator* (double scale) const;

/// ���g��Transf��Matrix�̏�Z���s���I�u�W�F�N�g�𐶐�����B
///Multiplication of transf and matrix.
MGTransf operator* (const MGMatrix&) const;

///���g��Transf�Ɨ^����ꂽTransf�̏�Z���s���I�u�W�F�N�g�𐶐�����B
///Multiplication of two transfs.
MGTransf operator* (const MGTransf& ) const;

///���g��Transf��scale�̏�Z���s�����g��Transf�Ƃ���B
///Scaling of the transf.
MGTransf& operator*= (double scale);

/// ���g��Transf��Matrix�̏�Z���s�����g��Transf�Ƃ���B
///Multiplication of transf and matrix.
MGTransf& operator*= (const MGMatrix&);

///���g��Transf�Ɨ^����ꂽTransf�̏�Z���s�����g��Transf�Ƃ���B
///Multiplication of two transfs.
MGTransf& operator*= (const MGTransf&);

///Boolean operation.
///���g��Transf�Ɨ^����ꂽTransf�����������ǂ�����r���s���B
///Equal comparison of two transf.
bool operator== (const MGTransf&) const;
bool operator!= (const MGTransf&) const;

//////////// Member Function /////////

///  Transf�̃A�t�B��������ԋp����B
///Return affine matrix of the transformation.
const MGMatrix& affine() const{return m_affine;};

///Convert this transf matrix to OpenGL Matrix.
void convert_to_glMatrix(
	double glMat[16]	///OpenGL Matrix will be output.
)const;

///Test if this is null.
bool is_null()const{return m_affine.is_null();};

///PD124=Transformation matrix.
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Reference (i,j)-th element.
double ref(size_t i, size_t j) const;

///Change the space dimension of this MGTransf to the new sdim.
void resize(size_t sdim);

///  Return space dimension
size_t sdim() const{return m_affine.sdim();};

/// Construct a transf to do the transformation of matrix around input point 
/// instead of origin of matrix, and replace own transf with it.
MGTransf& set_matrix(const MGMatrix& mat, const MGPosition& point);
   
///Set this as a null transf.
void set_null();

///  �w��_��ʂ�A�^����ꂽ�x�N�g���Ɋւ��ċ��ʕϊ����� 2D Transf
///  �𐶐���������Transf�Ɠ��ꊷ����B�w��_
///  �ɉ����w�肳��Ȃ��ꍇ�A���_���g�p����B
///2D mirror reflection transf about a line whose direction
///is vec and passes through point P.
MGTransf& set_reflect_2D(	const MGVector& vec,
							const MGPosition& P= mgORIGIN_2D);

///  �w��_�𒆐S�ɔ����v���Ɏw��p�x��]����
///  2D Transf�𐶐���������Transf�Ɠ��ꊷ����B�w��_
///  �ɉ����w�肳��Ȃ��ꍇ�A���_���g�p����B
///Rotation 2D matrix around point P by angle.
MGTransf& set_rotate_2D(	double angle,
							const MGPosition& P= mgORIGIN_2D); 

///  �w��_��ʂ�A�^����ꂽvector�ɐ����ȕ��ʂɊւ��ċ��ʕϊ�����
///  3D��Transf�𐶐���������Transf�Ɠ��ꊷ����B�w��_�ȗ����͌��_
///  �Ƃ���B
///3D mirror reflection transf about a plane whose normal
///is vec and passes through point P.
MGTransf& set_reflect_3D(
	const MGVector&,
	const MGPosition& P= mgORIGIN);

///  �w��_��ʂ�w�����vector���������̉��Ɏw��p�x��]����
///  Transf�𐶐���������Transf�Ɠ��ꊷ����B�w��_
///  �ɉ����w�肳��Ȃ��ꍇ�A���_���g�p����B
///3D rotation matrix around the straight line whose direction is vec
///and that passes through point P.
MGTransf& set_rotate_3D(
	const MGVector& vec,
	double angle,                                      
	const MGPosition& P= mgORIGIN); 

///  �^����ꂽ�g��A�k����Transf��������Transf��
///  ���ꊷ����B
///Update to scaling transf of same scaling value for each coordinate.
///Not change space dimension.
///No translation.
MGTransf& set_scale(double scale);

///  �^����ꂽ�g��A�k����Transf��������Transf��
///  ���ꊷ����B Scales are input through scale[sdim()].
///Update to scaling trans of different scaling values for each coordinate.
///Not change space dimension.
///No translation.
MGTransf& set_diff_scale(double* scale);

///Set up this transf from OpenGL matrix.
MGTransf& set_glMatrix(const double glMat[16]);

///  Transf�̕��s�ړ�������ԋp����B
///Return translation part of the transf.
const MGVector& translation() const{return m_translation;};

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:

//////////// Member Data /////////
  MGMatrix m_affine;		///< Affine transformation part.
  MGVector m_translation;	///< Translation part.

};

/** @} */ // end of BASE group
#endif
