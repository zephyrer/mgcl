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

/// 自身のTransfとscaleの乗算を行いオブジェクトを生成する。
///Scaling of the transf.
MGDECL friend MGTransf operator* (double scale, const MGTransf&);

/// ベクトルとトランスフォームの乗算を行いTransfを生成
///Genarate a Transf by multiplication of a vector and transf.
MGDECL friend MGTransf operator* (const MGVector& v, const MGTransf& tr);

///String stream function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGTransf&);

//////////// Constructor //////////

/// General Space dimension Constructor.
/// Void Constructor
explicit MGTransf(size_t dim=0);

///  全てのコンポーネントを指定してTransfを生成。
///Construct Transf from a matrix and a vector.
///***** This is the fundamental constructor.*****
MGTransf(const MGMatrix&, const MGVector&);

/// 各軸方向で等しい Scaling のTransf生成。
///Transf from same scaling factor for all coordinates and a translation.
MGTransf(double scale, const MGVector&);

///Construct Transf by copying old Transf, changing space dimension,
///and ordering of coordinates.
MGTransf (size_t sdim, const MGTransf& transf,
			size_t start1=0, size_t start2=0);

///////// 2D Constructor.////////////

/// 各軸方向(x,y)で異なる Scaling の2D Transfを生成。
///2D Transf from different scaling factors for x and y coordinates.
///No translation.
MGTransf(double scalex, double scaley);

/// Construct 2D space Transf to transform data for 'origin' to be
///origin(0,0,0) and for 'unit' to be x-coordimate.
MGTransf(const MGUnit_vector& unit,  ///< unit vector to be x-coordinate.
         const MGPosition& origin	 ///<origin to be origin.
);

//////////// 3D Constructor.////////////

/// 各軸方向( x,y,z )で異なる Scaling の3D Transfを生成。
///3D Transf from different scaling factors for x, y and z coordinates.
///No translation.
MGTransf(double scalex, double scaley, double scalez);

///  与えられた点を原点にするように平行移動し、与えられた２つの単位vector
///  を各々、X軸、Y軸にするよう回転させる 3D Transfを生成する。もし、
///  ２つ目の単位vectorが１つ目の単位vectorと直交しない場合は、２つの
///  ベクトルのかわりに両vectorを含む平面内で直交するよう変換したvectorを
///  使用する。
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

/// 自身のTransfとvectorの加算を行いオブジェクトを生成する。
///Translation of Transf.
MGTransf operator+ (const MGVector&) const;

/// 自身のTransfとvectorの加算を行い自身のTransfとする。
///Translation of Transf.
MGTransf& operator+= (const MGVector&);

/// 自身のTransfとvectorの減算を行いオブジェクトを生成する。
///Translation of Transf.
MGTransf operator- (const MGVector&) const;

/// 自身のTransfとvectorの減算を行い自身のTransfとする。
///Translation of Transf.
MGTransf& operator-= (const MGVector&);

/// 自身のTransfとscaleの乗算を行いオブジェクトを生成する。
///Scaling of the transf.
MGTransf operator* (double scale) const;

/// 自身のTransfとMatrixの乗算を行いオブジェクトを生成する。
///Multiplication of transf and matrix.
MGTransf operator* (const MGMatrix&) const;

///自身のTransfと与えられたTransfの乗算を行いオブジェクトを生成する。
///Multiplication of two transfs.
MGTransf operator* (const MGTransf& ) const;

///自身のTransfとscaleの乗算を行い自身のTransfとする。
///Scaling of the transf.
MGTransf& operator*= (double scale);

/// 自身のTransfとMatrixの乗算を行い自身のTransfとする。
///Multiplication of transf and matrix.
MGTransf& operator*= (const MGMatrix&);

///自身のTransfと与えられたTransfの乗算を行い自身のTransfとする。
///Multiplication of two transfs.
MGTransf& operator*= (const MGTransf&);

///Boolean operation.
///自身のTransfと与えられたTransfが等しいかどうか比較を行う。
///Equal comparison of two transf.
bool operator== (const MGTransf&) const;
bool operator!= (const MGTransf&) const;

//////////// Member Function /////////

///  Transfのアフィン部分を返却する。
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

///  指定点を通り、与えられたベクトルに関して鏡面変換する 2D Transf
///  を生成し既存のTransfと入れ換える。指定点
///  に何も指定されない場合、原点を使用する。
///2D mirror reflection transf about a line whose direction
///is vec and passes through point P.
MGTransf& set_reflect_2D(	const MGVector& vec,
							const MGPosition& P= mgORIGIN_2D);

///  指定点を中心に反時計回りに指定角度回転する
///  2D Transfを生成し既存のTransfと入れ換える。指定点
///  に何も指定されない場合、原点を使用する。
///Rotation 2D matrix around point P by angle.
MGTransf& set_rotate_2D(	double angle,
							const MGPosition& P= mgORIGIN_2D); 

///  指定点を通り、与えられたvectorに垂直な平面に関して鏡面変換する
///  3DのTransfを生成し既存のTransfと入れ換える。指定点省略時は原点
///  とする。
///3D mirror reflection transf about a plane whose normal
///is vec and passes through point P.
MGTransf& set_reflect_3D(
	const MGVector&,
	const MGPosition& P= mgORIGIN);

///  指定点を通り指定方向vectorを持つ直線の回りに指定角度回転する
///  Transfを生成し既存のTransfと入れ換える。指定点
///  に何も指定されない場合、原点を使用する。
///3D rotation matrix around the straight line whose direction is vec
///and that passes through point P.
MGTransf& set_rotate_3D(
	const MGVector& vec,
	double angle,                                      
	const MGPosition& P= mgORIGIN); 

///  与えられた拡大、縮小のTransfを既存のTransfと
///  入れ換える。
///Update to scaling transf of same scaling value for each coordinate.
///Not change space dimension.
///No translation.
MGTransf& set_scale(double scale);

///  与えられた拡大、縮小のTransfを既存のTransfと
///  入れ換える。 Scales are input through scale[sdim()].
///Update to scaling trans of different scaling values for each coordinate.
///Not change space dimension.
///No translation.
MGTransf& set_diff_scale(double* scale);

///Set up this transf from OpenGL matrix.
MGTransf& set_glMatrix(const double glMat[16]);

///  Transfの平行移動部分を返却する。
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
