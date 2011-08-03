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

///自身のMatrixと与えられたscaleの乗算を行いオブジェクトを生成。
///Scaling of the matrix.
MGDECL friend MGMatrix operator* (double scale, const MGMatrix& mat);

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGMatrix&);

//////////// Constructor ////////////

///void constructor:void コンストラクタ
explicit MGMatrix(size_t sdim=0);

///  Construct 2D matrix from two vectors.
MGMatrix(const MGVector&, const MGVector&);

///  Construct 3D matrix from three vectors.
MGMatrix(const MGVector&, const MGVector&, const MGVector& );

///  各軸方向で等しい Scaling のためのMatrixを生成。
///Scaling matrix of same scaling value for each coordinate.
MGMatrix(size_t dim, double scale);

///Construct dim dimension matrix from the array of double values.
MGMatrix(
	size_t dim,			///<dimension of the matrix.
	const double* values,///<array of values[dim*dim].
	bool column_wise=true///<If column_wise=true, (*this)(i,j)=values[i+dim*j],
						///<else(row_wise), (*this)(i,j)=values[j+dim*i], for 0<=i,j<=dim-1.
);

///  与えられたVector回りに指定の角度回転させるMatrixを作成する。
///　Rotation matrix around vec. Space dimension of vec can be
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

///自身のMatrixと与えられたscaleの乗算を行いオブジェクトを生成。
///Scaling of the matrix.
MGMatrix operator* (double) const;

///自身のMatrixと与えられたscaleの乗算し自身のMatrixとする。
///Scaling of the matrix.
MGMatrix& operator*= (double);

///自身のMatrixと与えられたMatrixの乗算を行いオブジェクトを生成。
///Matrix and matrix multiplication.
MGMatrix operator* (const MGMatrix&) const;

///自身のMatrixと与えられたMatrixの乗算し自身のMatrixとする。
///Matrix and matrix multiplication.
MGMatrix& operator*= (const MGMatrix&);

///自身のMatrixと与えられたTransfの乗算を行いオブジェクトを生成。
///Matrix and Transf multiplication.
MGTransf operator* (const MGTransf&) const;

///  Boolean Operation
///  自身のMatrixと与えられたMatrixが等しいかどうか
///  比較を行う。
///Equation comparison.
bool operator== (const MGMatrix&) const;
bool operator!= (const MGMatrix&) const;

////////////Member Function ////////////

///Convert this transf matrix to OpenGL Matrix.
void convert_to_glMatrix(
	double glMat[16]	///<OpenGL Matrix will be output.
)const;

///  行列式の値を返却。
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

///  原点を通り、指定Vectorに関して鏡面変換する 2D Matrixを作成し，
///  既存のMatrixと入れ換える。
///Mirror reflection 2D matrix about a vector through origin.
MGMatrix& set_reflect_2D(const MGVector& );

///  原点の回りに指定の角度回転させる2D Matrixを作成し,
///  既存の Matrixと入れ換える。
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

///  与えられた２つの単位ベクトルを各々、X軸、Y軸にするよう原点の周りに
///  回転させる 3D Matrix を生成し，自身のMatrixと入れ替える。もし、
///  ２つ目の単位ベクトルが１つ目の単位ベクトルと直交しない場合は、２つのベ
///  クトルのかわりに両ベクトルを含む平面内で直交するよう変換したベクトルを
///  使用する。
///3D Matrix to transform uvecx to be x-axis and uvecy to be y-axis.
///If uvecx and uvecy does not cross at right angle, uvecy will be
///transformed.
MGMatrix& set_xy_axis(
	const MGUnit_vector& uvecx,	///<Unit vector 1 for x axis.
    const MGUnit_vector& uvecy);///<Unit vector 2 for y axis.

///  X軸、Y軸を各々、与えられた２つの単位ベクトルにするよう原点の周りに
///  回転させる 3D Matrix を生成し，自身のMatrixと入れ替える。もし、
///  ２つ目の単位ベクトルが１つ目の単位ベクトルと直交しない場合は、２つのベ
///  クトルのかわりに両ベクトルを含む平面内で直交するよう変換したベクトルを
///  使用する。
///3D matrix to transform x and y axis to be uvecx and uvecy.
/// This is the inverse matrix of set_xy_axis().
///If uvecy does not cross uvecx at right angle, uvecy is transformed to do so.
MGMatrix& set_xy_vector(
	const MGUnit_vector& uvecx,	///<Unit vector 1 for x axis.
    const MGUnit_vector& uvecy);///<Unit vector 2 for y axis.

///  原点を通り、指定Vectorに垂直な平面に関して鏡面変換する 3D Matrix
///  を既存のMatrixと入れ換える。
///3D mirror reflection matrix about a plane whose normal
///is vec and passes through the origin.
MGMatrix& set_reflect_3D (const MGVector&  vec);

///（原点を基点とする）Vector V0 を V1に変換するMatrixを作成し、
///  自身のMatrixと入れ換える。
///Construct the matrix to rotate and scale that transform vector V0 to V1,
///and replace this matrix with the matrix.
///Space dimension of V0 and V1 can be any number greater than 1.
MGMatrix& set_rotate(const MGVector& V0, const MGVector& V1);

///  与えられたVector回りに指定の角度回転させる3D Matrixを作成し既存の
///  Matrixと入れ換える。
///3D rotation matrix around vec.
MGMatrix& set_rotate_3D (const MGVector& vec, double angle);

///3D rotation matrix around vec.
///The angle is given by cval as cos(angle) and sval as sin(angle).
MGMatrix& set_rotate_3D (const MGVector& vec, double cval, double sval);

///  各軸方向で等しい Scaling のためのMatrixを生成し、既存の
///  Matrixと入れ換える
///Scaling matrix of same scaling value for each coordinate.
///Not change space dimension.
MGMatrix& set_scale (double);

///  各軸方向で異なる Scaling のためのMatrixを生成し、既存の
///  Matrixと入れ換える。
///Scaling matrix of different scaling values for each coordinate.
///Not change space dimension.
MGMatrix& set_diff_scale (double*);

///Set up this matrix from OpenGL matrix.
MGMatrix& set_glMatrix(const double glMat[16]);

///Set this as a null matrix.
void set_null();

///  転置Matrixを生成。
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
	
/// マトリックスによるベクトルの変換を行いオブジェクトを生成
///Matrix transformation of the vector.
MGVector operator* (const MGVector& v, const MGMatrix& m);

/// マトリックスによるベクトルの変換を行い自身のベクトルとする
///Matrix transformation of the vector.
MGVector& operator*= (MGVector& v, const MGMatrix& m);

/// マトリックスによるベクトルの変換を行い自身のベクトルとする
///Update own vector by matrix transformation.
///The result is unit of transformed vector.
MGUnit_vector& operator*= (MGUnit_vector& v,const MGMatrix&);

/** @} */ // end of BASE group
#endif
