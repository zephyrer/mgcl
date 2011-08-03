/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGMatrix_HH_
#define _MGMatrix_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

//  MGMatrix.h
//  Defines Class MGMatrix
//
//  Forward Declarations
class MGVector;
class MGUnit_vector;
class MGTransf;
class MGIfstream;
class MGOfstream;

///MGMatrix is a matix of m by m, where m is the space dimension.
///MGMatrix provides transformation around the origin. General transformation
///is provided by MGTransf.
///Let M be a matrix and A be an object(including MGVector).
///Then matrix transformation of the object A is defined as:  A*M. (Not M*A)
class MGCLASS MGMatrix {

public:

///���g��Matrix�Ɨ^����ꂽscale�̏�Z���s���I�u�W�F�N�g�𐶐��B
///Scaling of the matrix.
MGDECL friend MGMatrix operator* (double scale, const MGMatrix& mat);

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGMatrix&);

//////////// Constructor ////////////

///void constructor:void �R���X�g���N�^
explicit MGMatrix(size_t sdim=0);

///  Construct 2D matrix from two vectors.
MGMatrix(const MGVector&, const MGVector&);

///  Construct 3D matrix from three vectors.
MGMatrix(const MGVector&, const MGVector&, const MGVector& );

///  �e�������œ����� Scaling �̂��߂�Matrix�𐶐��B
///Scaling matrix of same scaling value for each coordinate.
MGMatrix(size_t dim, double scale);

///Construct dim dimension matrix from the array of double values.
MGMatrix(
	size_t dim,			///<dimension of the matrix.
	const double* values,///<array of values[dim*dim].
	bool column_wise=true///<If column_wise=true, (*this)(i,j)=values[i+dim*j],
						///<else(row_wise), (*this)(i,j)=values[j+dim*i], for 0<=i,j<=dim-1.
);

///  �^����ꂽVector���Ɏw��̊p�x��]������Matrix���쐬����B
///�@Rotation matrix around vec. Space dimension of vec can be
/// any number, i.e. can be more than 3.
MGMatrix(const MGVector& vec, double angle);

///Construct a matrix to rotate in (axis1, axis2) plane by algle theta,
///where angle is defined as cosv=cos(angle), and sinv=sin(angle).
MGMatrix(size_t dim,	///<Space dimension(can be more than 3).
  size_t axis1,		///<axis number 1
  size_t axis2,		///<axis number 2
  double cosv, double sinv);///<cosv=cos(angle), sinv=sin(angel).
					///<That is cosv*cosv+sinv*sinv must be 1.

///Construct Matrix by copying old Matrix, changing space dimension and
///ordering of old coordinates.
MGMatrix(size_t dim, const MGMatrix&, size_t start1=0, size_t start2=0);

///Copy constructor.
MGMatrix(const MGMatrix& mat);

////////////// Destructor./////////
~MGMatrix(){if(m_matrix) delete[] m_matrix;};

//////////// Operator Overload. ////////////

/// Assignment
MGMatrix& operator=(const MGMatrix& mat);

///Reference (i,j)-th element.
double operator() (size_t i, size_t j) const{return ref(i,j);}

///Access to (i,j)-th element.
double& operator() (size_t i, size_t j) {
	assert(i<sdim() && j<sdim());
	return m_matrix[i+j*m_sdim];
}

///���g��Matrix�Ɨ^����ꂽscale�̏�Z���s���I�u�W�F�N�g�𐶐��B
///Scaling of the matrix.
MGMatrix operator* (double) const;

///���g��Matrix�Ɨ^����ꂽscale�̏�Z�����g��Matrix�Ƃ���B
///Scaling of the matrix.
MGMatrix& operator*= (double);

///���g��Matrix�Ɨ^����ꂽMatrix�̏�Z���s���I�u�W�F�N�g�𐶐��B
///Matrix and matrix multiplication.
MGMatrix operator* (const MGMatrix&) const;

///���g��Matrix�Ɨ^����ꂽMatrix�̏�Z�����g��Matrix�Ƃ���B
///Matrix and matrix multiplication.
MGMatrix& operator*= (const MGMatrix&);

///���g��Matrix�Ɨ^����ꂽTransf�̏�Z���s���I�u�W�F�N�g�𐶐��B
///Matrix and Transf multiplication.
MGTransf operator* (const MGTransf&) const;

///  Boolean Operation
///  ���g��Matrix�Ɨ^����ꂽMatrix�����������ǂ���
///  ��r���s���B
///Equation comparison.
bool operator== (const MGMatrix&) const;
bool operator!= (const MGMatrix&) const;

////////////Member Function ////////////

///Convert this transf matrix to OpenGL Matrix.
void convert_to_glMatrix(
	double glMat[16]	///<OpenGL Matrix will be output.
)const;

///  �s�񎮂̒l��ԋp�B
double determinant() const;

///Construct a matrix to transform a unit vector on an axis
/// to a vector 'uvec', and replace own matrix with it.
///Inverse matrix of to_axis.
/// axis can be any number(can be more than 2).
MGMatrix& from_axis(
  const MGUnit_vector& uvec,	///< Unit vector for an axis.
  size_t axis=0);				///< Axis kind 0:x, 1:y, 2:z, ...

///Test if this is null.
bool is_null()const{return m_sdim==0;};

///Reference to (i,j)-th element of the matarix.
double ref(size_t i, size_t j) const;

///Construct a mirror reflection matrix about a plane whose normal
///is vec and that passes through the origin, then
///replace own matrix with it.
MGMatrix& reflection(const MGVector& vec);

///Resize, i.e., change, the space dimension.
///Result matreix will contain garbages.
void resize(size_t nsdim);

///Return space dimension
size_t sdim() const{return m_sdim;}

/// Construct 2D space Matrix to transform for 'unit' to be x-coordimate, and
/// replace own Matrix with it.
/// 2D version of to_axis.
MGMatrix& set_x_axis
  (const MGUnit_vector& unit); ///unit vector to be x-coordinate

///  ���_��ʂ�A�w��Vector�Ɋւ��ċ��ʕϊ����� 2D Matrix���쐬���C
///  ������Matrix�Ɠ��ꊷ����B
///Mirror reflection 2D matrix about a vector through origin.
MGMatrix& set_reflect_2D(const MGVector& );

///  ���_�̉��Ɏw��̊p�x��]������2D Matrix���쐬��,
///  ������ Matrix�Ɠ��ꊷ����B
///Rotation 2D matrix around origin by angle(in radian).
MGMatrix& set_rotate_2D(double angle);

///Rotation 2D matrix around origin by an angle.
///The angle is given by cval as cos(angle) and sval as sin(angle).
MGMatrix& set_rotate_2D(double cval, double sval);

///Construct a 3D matrix to transform a vector 'uvec' to be one of the axises,
/// and replace own matrix. Inverse matrix of set_vector.
/// 3D version of to_axis.
MGMatrix& set_axis(
  const MGUnit_vector& uvec,	///< Unit vector to be an axis.
  size_t axis=0);				///< Axis kind 0:x, 1:y, 2:z.

///Construct a matrix to transform a unit vector on an axis
/// to a vector 'uvec', and replace own matrix with it.
///Inverse matrix of set_axis.
/// 3D version of from_axis.
MGMatrix& set_vector(
  const MGUnit_vector& uvec,	///< Unit vector for an axis.
  size_t axis=0);				///< Axis kind 0:x, 1:y, 2:z.

///  �^����ꂽ�Q�̒P�ʃx�N�g�����e�X�AX���AY���ɂ���悤���_�̎����
///  ��]������ 3D Matrix �𐶐����C���g��Matrix�Ɠ���ւ���B�����A
///  �Q�ڂ̒P�ʃx�N�g�����P�ڂ̒P�ʃx�N�g���ƒ������Ȃ��ꍇ�́A�Q�̃x
///  �N�g���̂����ɗ��x�N�g�����܂ޕ��ʓ��Œ�������悤�ϊ������x�N�g����
///  �g�p����B
///3D Matrix to transform uvecx to be x-axis and uvecy to be y-axis.
///If uvecx and uvecy does not cross at right angle, uvecy will be
///transformed.
MGMatrix& set_xy_axis(
	const MGUnit_vector& uvecx,	///<Unit vector 1 for x axis.
    const MGUnit_vector& uvecy);///<Unit vector 2 for y axis.

///  X���AY�����e�X�A�^����ꂽ�Q�̒P�ʃx�N�g���ɂ���悤���_�̎����
///  ��]������ 3D Matrix �𐶐����C���g��Matrix�Ɠ���ւ���B�����A
///  �Q�ڂ̒P�ʃx�N�g�����P�ڂ̒P�ʃx�N�g���ƒ������Ȃ��ꍇ�́A�Q�̃x
///  �N�g���̂����ɗ��x�N�g�����܂ޕ��ʓ��Œ�������悤�ϊ������x�N�g����
///  �g�p����B
///3D matrix to transform x and y axis to be uvecx and uvecy.
/// This is the inverse matrix of set_xy_axis().
///If uvecy does not cross uvecx at right angle, uvecy is transformed to do so.
MGMatrix& set_xy_vector(
	const MGUnit_vector& uvecx,	///<Unit vector 1 for x axis.
    const MGUnit_vector& uvecy);///<Unit vector 2 for y axis.

///  ���_��ʂ�A�w��Vector�ɐ����ȕ��ʂɊւ��ċ��ʕϊ����� 3D Matrix
///  ��������Matrix�Ɠ��ꊷ����B
///3D mirror reflection matrix about a plane whose normal
///is vec and passes through the origin.
MGMatrix& set_reflect_3D (const MGVector&  vec);

///�i���_����_�Ƃ���jVector V0 �� V1�ɕϊ�����Matrix���쐬���A
///  ���g��Matrix�Ɠ��ꊷ����B
///Construct the matrix to rotate and scale that transform vector V0 to V1,
///and replace this matrix with the matrix.
///Space dimension of V0 and V1 can be any number greater than 1.
MGMatrix& set_rotate(const MGVector& V0, const MGVector& V1);

///  �^����ꂽVector���Ɏw��̊p�x��]������3D Matrix���쐬��������
///  Matrix�Ɠ��ꊷ����B
///3D rotation matrix around vec.
MGMatrix& set_rotate_3D (const MGVector& vec, double angle);

///3D rotation matrix around vec.
///The angle is given by cval as cos(angle) and sval as sin(angle).
MGMatrix& set_rotate_3D (const MGVector& vec, double cval, double sval);

///  �e�������œ����� Scaling �̂��߂�Matrix�𐶐����A������
///  Matrix�Ɠ��ꊷ����
///Scaling matrix of same scaling value for each coordinate.
///Not change space dimension.
MGMatrix& set_scale (double);

///  �e�������ňقȂ� Scaling �̂��߂�Matrix�𐶐����A������
///  Matrix�Ɠ��ꊷ����B
///Scaling matrix of different scaling values for each coordinate.
///Not change space dimension.
MGMatrix& set_diff_scale (double*);

///Set up this matrix from OpenGL matrix.
MGMatrix& set_glMatrix(const double glMat[16]);

///Set this as a null matrix.
void set_null();

///  �]�uMatrix�𐶐��B
///Transpose matrix.
MGMatrix transpose() const;

///Construct a matrix to transform a vector 'uvec' to be one of the axises,
/// and replace own matrix. Inverse matrix of from_axis.
/// axis can be any number(can be more than 2).
MGMatrix& to_axis(
	const MGUnit_vector& uvec,	///< Unit vector to be an axis.
	size_t axis=0);				///< Axis kind 0:x, 1:y, 2:z, ...

///Dump Functions.
///Calculate dump size
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Obtain the scaling factor of this matrix.
double scale()const;

////////////Member Data/////////////////

private:

  ///Elements of matrix is stored in m_matrix.
	size_t m_sdim;		///<Matrix is sdim by sdim.
	double* m_matrix;	///<Elements data area of length sdim*sdim
  
/// Member Function
private:

/// Multiply two same dimension's matrix this and m2.
MGMatrix multiply(const MGMatrix& m2) const;

};
	
/// �}�g���b�N�X�ɂ��x�N�g���̕ϊ����s���I�u�W�F�N�g�𐶐�
///Matrix transformation of the vector.
MGVector operator* (const MGVector& v, const MGMatrix& m);

/// �}�g���b�N�X�ɂ��x�N�g���̕ϊ����s�����g�̃x�N�g���Ƃ���
///Matrix transformation of the vector.
MGVector& operator*= (MGVector& v, const MGMatrix& m);

/// �}�g���b�N�X�ɂ��x�N�g���̕ϊ����s�����g�̃x�N�g���Ƃ���
///Update own vector by matrix transformation.
///The result is unit of transformed vector.
MGUnit_vector& operator*= (MGUnit_vector& v,const MGMatrix&);

/** @} */ // end of BASE group
#endif
