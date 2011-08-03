/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGVector_HH_
#define _MGVector_HH_
/** @addtogroup BASE
 *  @{
 */
#include <stddef.h>
#include <vector>
#include "mg/MGCL.h"

// MGVector.h
// Header for MGVector.

//Forward Declaration
class MGUnit_vector;
class MGPosition;
class MGIfstream;
class MGOfstream;
class MGIgesOfstream;

///Vector of a general n space dimension.
class MGCLASS MGVector {

public:

/// �x�N�g���̉��Z
///Addition of two vectors.
MGDECL friend MGVector operator+(const MGVector& vec1,const MGVector& vec2);

/// �x�N�g���̌��Z
///Subtraction of two vectors.
MGDECL friend MGVector operator-(const MGVector& vec1,const MGVector& vec2);

/// �x�N�g���̓���
///Inner product of two vectors.
MGDECL friend double operator%(const MGVector& vec1,const MGVector& vec2);

///�x�N�^�̊O��
///vector product of two vectors.
MGDECL friend MGVector operator*(const MGVector& vec1,const MGVector& vec2);

/// �X�J���[�̏�Z���s���I�u�W�F�N�g�𐶐�
///Scalar multiplication.
MGDECL friend MGVector operator*(const MGVector& vec1,double scale);

/// �x�N�g���̃X�J���[��Z���s���I�u�W�F�N�g�𐶐�
///Scalar multiplication.
MGDECL friend MGVector operator* (double, const MGVector&);

/// �X�J���[���Z���s���I�u�W�F�N�g�𐶐�
///Scalar division.
MGDECL friend MGVector operator/(const MGVector& vec1,double scale);

///Test if this vector is less than v2.
///Comparison depends on two vectors' length.
inline friend
bool operator<(const MGVector& v1,const MGVector& v2){return v1.len()<v2.len();};
inline friend
bool operator<=(const MGVector& v1,const MGVector& v2){return v1.len()<=v2.len();};
inline friend
bool operator>(const MGVector& v1,const MGVector& v2){return v1.len()>v2.len();};
inline friend
bool operator>=(const MGVector& v1,const MGVector& v2){return v1.len()>=v2.len();};

/// �^����ꂽ�x�N�g���̐����̒l���r���A�����ł���� TRUE ��ԋp
///Test if two vectors are equal.
friend bool operator==(const MGVector& v1,const MGVector& v2);

/// �^����ꂽ�x�N�g���̐����̒l���r���Anot equal�̂Ƃ� TRUE ��ԋp
///Test if two vectors are equal.
inline friend bool operator!=(const MGVector& v1,const MGVector& v2){return !(v1==v2);}

///String stream function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGVector&);

/// �R�̃x�N�g�����狁�߂���s��̍s�񎮂̒l��ԋp����B
///Determinant of 3 by 3 matrix of 3 vectors.
MGDECL friend double MGDeterminant(
		const MGVector& v1, const MGVector& v2, const MGVector& v3
);

////////////// Constructor. �R���X�g���N�^ ////////////

///Void constructor.
explicit MGVector(size_t sdim=0);

///Construct 2D vector by providing each element data.
MGVector(double x, double y);

///Construct 3D vector by providing each element data.
MGVector(double x, double y , double z);

///Construct 4D vector by providing each element data.
MGVector(double x, double y, double z, double w);

/// �����l�@v �ł��ׂẴG�������g�����������ăx�N�g���𐶐�����B
///Vector of same value for each coordinate element.
MGVector(size_t sdim, double v);

/// double �̔z��ŃG�������g�����������ăx�N�g���𐶐�����B
///Vector from array of double v[sdim].
///***** This is the fundamental constructor.*****
MGVector(size_t sdim, const double* v);

///Construct a vector from a class MGPosition.
///Vector from a position.
MGVector(const MGPosition&);

///Construct a vector from a difference of two vectors.
MGVector(const MGVector &,	///Destination point
		const MGVector &);	///Source point

///Construct Vector by copying old Vector, changing space dimension and
///ordering of old coordinates.
/// (*this)(start1+i)=vec2(start2+i).
MGVector(
	size_t sdim,			///Space dimension.
	const MGVector& vec2,	///Original vector.
	size_t start1=0,		///id of constructing vector that indicates
							///from where to store the elements of vec2.
	size_t start2=0);		///id of vec2.

///Construct from std::vector<double>
MGVector(const std::vector<double>& darrays);

///Construct from std::valarray<double>
///MGVector(const std::valarray<double>& darrays);

///Copy constructor
MGVector ( const MGVector& );

//////////// Destructor //////////
~MGVector(){if(m_sdim>3) delete[] m_element;}

//////////// Operator Oveload ////////////

///Assignment
MGVector& operator =(const MGVector &);

///Return i-th element of the vector.
double operator() (size_t i) const{return ref(i);}  

///Return i-th element of the vector.
double operator[] (size_t i) const{return ref(i);}  

///Access to i-th Inteval.
///This is left hand side value operator. If only regference is needed,
/// operator[] should be used. 
double& operator()(size_t i);

///Update vector data by array of double.
MGVector& operator=(const double*);

/// ���g�̃x�N�g���ɗ^����ꂽ�x�N�g�������Z���Ď��g�̃x�N�g���Ƃ��� 
///Addition of two vectors.
MGVector & operator+= (const MGVector&);

/// �P���}�C�i�X�B���g�̃x�N�g���𔽓]���A�I�u�W�F�N�g�𐶐�
///Unary minus. Negate all the elements of the vector.
MGVector operator- () const;

/// ���g�̃x�N�g������^����ꂽ�x�N�g�������Z�����g�̃x�N�g���Ƃ���
///Subtraction of two vectors.
MGVector & operator -= ( const MGVector & );

/// �X�J���[�̏�Z���s�����g�̃x�N�g���Ƃ���
///Scalar multiplication.
MGVector& operator*= (double scale);

///Update own vector by vector product output, changes to 3D vector.
MGVector& operator*= (const MGVector& vec2);

/// �X�J���[���Z���s�����g�̃x�N�g���Ƃ���
///Scalar division.
MGVector& operator/= (double scale);

//////////// Member Function ////////////

/// ���g�̃x�N�g���Ɨ^����ꂽ�x�N�g���̂Ȃ��p�x�� Radian �ŕԋp
///Compute angle in radian of two vectors.
/// 0<= angle <pai.
double angle(const MGVector&) const;

///Compute angle in radian of two vectors.
/// 0<= angle <pai.
double anglepai(const MGVector& v2)const{return angle(v2);};

///Compute the angle in radian that is measured from this to v2 around the normal N.
///The angle's range is 0<= angle <2*pai.
///Although N is assumed to be parallel to N2=(*this)*v2, N may not perpendicular
///to v1 and v2, in which case, the projected N to N2 is used to measure the angle.
///v1.angle2pai(v2,N)+v2.angle2pai(v1,N)=2*pai always holds.
double angle2pai(const MGVector& v2, const MGVector& N)const;

/// ���g�̃x�N�g���Ɨ^����ꂽ�x�N�g���̂Ȃ��p�x�� cos�� �ŕԋp����
/// ���g���^����ꂽ�x�N�g������x�N�g���̎��́Acos���� 1.0 �Ƃ���
///Compute angle in cosine of two vectors.
double cangle(const MGVector&) const;

///Clear all the element by the value init.
MGVector& clear(double init=0.0);

///Return the 1st address of the array of the vector double data.
const double* data()const{return m_element;};
double* data(){return m_element;};

/// Generate a vector by interpolating two vectors. Input scalar is a ratio
/// and when zero, output vector is a copy of *this.
/// New vector vnew=(1-t)*(*this)+t*vec2.
MGVector interpolate(double t, const MGVector& vec2) const;

/// Generate a vector by interpolating two vectors by rotation.
/// Input scalar t is a ratio and when t is zero,
/// output vector is a copy of *this and when t=1., 	output is vec2.
/// New vector vnew=a*(*this)+b*vec2, where
/// a=sin(theta2)/sin(theta), b=sin(theta1)/sin(theta). Here,
/// theta=angle of *this and vec2. theta1=t*theta, theta2=theta-theta1.
/// theta may be zero.
///When ratio is not null, ratio[0]=a and ratio[1]=b will be returned.
MGVector interpolate_by_rotate(
	double t, const MGVector& vec2,
	double* ratio=0
) const;

///Test if this and v2 are on a single straight line.
///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGVector& v2
)const{	return (*this).parallel(v2);}

///Test if this, v2, and v3 are on a single straight line.
///Function's return value is true if the three points are on a straight,
///false if not.
bool is_collinear(
	const MGVector& v2,
	const MGVector& v3
)const;

///Test if this is null.
bool is_null()const{return m_sdim==0;}

/// ���g�̃x�N�g�����P�ʃx�N�g���Ȃ� TRUE ��ԋp
///Test if the vector is unit.
bool is_unit_vector() const;

///Return true when the vector is a zero vector.
bool is_zero_vector() const;

/// �x�N�g���̒�����ԋp����B
///Return vector length.
double len() const;

///Negate the vector.
void negate(){operator*=(-1.);};

/// ��ʃx�N�g����P�ʃx�N�g�������I�u�W�F�N�g�𐶐�
///Generate unit vector from the vector.
MGUnit_vector normalize() const;

/// ���g�̃x�N�g���Ɨ^����ꂽ�x�N�g�����������ǂ����ԋp
/// �����̂΂��� TRUE
///Test if two vectors are orthogonal, i.e. cross at right angle.
bool orthogonal(const MGVector& ) const;

///PD123=Direction.
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Compute the vector that is orthogonal to vec2 and is closest to this.
///"closest" means that the angle of the two vectors is minimum and
///the two vector length are equal.
MGVector orthogonize(const MGVector& vec2)const;

/// ���g�̃x�N�g���Ɨ^����ꂽ�x�N�g�������s�Ȃ� TRUE ��ԋp
///Test if two vectors are parallel.
bool parallel(const MGVector& ) const;

/// ���g�̃x�N�g�����x�N�g��(v2)�Ɏˉe�����x�N�g�������߂�B
/// v2 �� ��x�N�g���̂Ƃ�(*this)���Ԃ�B
MGVector project(const MGVector& v2) const;

///Reference to i-th element.
double ref(size_t i) const{ 
	if(i<sdim()) return m_element[i];
	else         return 0.;
}

///Resize the vector, that is , change space dimension.
///When this is enlarged, the extra space will contain garbages.
void resize(size_t new_sdim);

/// ���g�̃x�N�g���Ɨ^����ꂽ�x�N�g���̂Ȃ��p�x�� sin�l��ԋp
///Compute angle in sine of two vectors.
/// sanlge>=0.
double sangle(const MGVector& ) const;

///Get the space dimension
size_t sdim() const { return m_sdim; };

///Set this as a null vector.
void set_null();

///Change this to a unit vector.
void set_unit();

///Store vec2 data into *this.
///Store length is vec2.len().
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	size_t i,				///Displacement of *this.
	const MGVector& vec2,	///Vector 2.
	size_t j=0);			///Displacement of vec2.

///Store vec2 data into *this.
///Storing will be done rap-around. That is, if id i or j reached to
///each sdim(), the id will be changed to 0.
void store_at(
	size_t i,				///Displacement of *this.
	const MGVector& vec2,	///Vector 2.
	size_t j,				///Displacement of vec2.
	size_t len);			///Length to store 

///swap two coordinates.
///swap coordinates (i) and (j).
void swap(size_t i, size_t j);

///Dump Functions.
///Calculate dump size
size_t dump_size() const;
///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Friend Function

protected:
/// Protected data member
	size_t m_sdim;
	double* m_element;
	double m_data[3];	///For vector data of space dimension less or equal to 3.
	mutable double m_length;///To hold vector length if computed.
						///When length not computed, negative value will be set.

	///Set data at element i.
	///This should be used with care, since m__length will not be set.
	///Maintenance of m_length should be done by the user, that is, m_length=-1
	///must set if updated.
	double& set(size_t i) {return m_element[i];}

friend class MGLBRep;
friend class MGSBRep;

};

/// V1���x�N�g��(v2)�Ɏˉe�����x�N�g�������߂�B
/// v2 �� ��x�N�g���̂Ƃ�V1���Ԃ�B
MGVector project(const MGVector& V1, const MGVector& V2);

///Compute the angel around the normal N in radian range[0., 2*pia).
///angle(v1,v2,N)+angle(v2,v1,N)=2*pai always holds.
inline double angle(const MGVector& V1, const MGVector& V2, const MGVector& N){
	return V1.angle2pai(V2,N);
};

/** @} */ // end of BASE group
#endif
