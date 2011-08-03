/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSPointSeq_HH_
#define _MGSPointSeq_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGSPointSeq.h
//

///Forward Declaration
class MGBox;
class MGNDDArray;
class MGVector;
class MGPosition;
class MGMatrix;
class MGTransf;
class MGBPointSeq;
class MGPlane;
class MGStraight;
class MGIfstream;
class MGOfstream;

///Defines Spoint seq of a space dimension and of a size.
///MGSPointSeq is designed to represent a control polygon network of  B-Spline surface.
///MGSPointSeq has its capacity and actual length used.
///The dimensioning is like MGSPointSeq spoint(i,j,k).
/// i has u direction length nu, j has v direction length nv, and k has space dimension sd.
///Then minimum area size is nu*nv*sd, and 0<=i<nu, 0<=j<nv, 0<=k<sd.
class MGCLASS MGSPointSeq {

public:

/// 与えられたスケーリングで曲線の変換を行いオブジェクトを生成する。
///Scaling.
MGDECL friend MGSPointSeq operator* (double scale, const MGSPointSeq&);

///SDtring stream Function
MGDECL friend std::ostream& operator<< (std::ostream& out, const MGSPointSeq& sp);

//////////// Constructor ////////////

///Of m_lenght=size(for u and v dierction),
///Dimensioning of MGSointSeq is MGSPointSeq(i,j,k).
/// i:u-direction, j:v-direction, k:space dimension.
explicit MGSPointSeq(
	size_t sizeu=0,	///<size of u-direction
	size_t sizev=0,	///<size of v-direction
	size_t dim=0);	///<space dimension

///Construct a MGSPointSeq by copying original MGSPointSeq.
///Can change the ordering of coordinates 
MGSPointSeq(
	size_t dim,						///< New Space Dimension.
	const MGSPointSeq& old_brep,	///< Origianl B-Rep.
	size_t start1=0,///< Destination start space dimension order to store.
	size_t start2=0	///< Source start space dimension order to retrieve.
);

///Copy constructor.
MGSPointSeq(const MGSPointSeq& rhs);

////////////// Destructor //////////////
~MGSPointSeq(){if(m_spoint) delete[] m_spoint;};

////////////// Operator Oveload //////////////
										  
///Assignment
MGSPointSeq& operator=(const MGSPointSeq&);

///Access to (i,j,k)th element
double& operator()(size_t i, size_t j, size_t k); 
double operator()(size_t i, size_t j, size_t k) const{return ref(i,j,k);};

///Extract (i,j,k) elements for 0<=k<sdim() as a vector.
MGVector operator()(size_t i, size_t j) const;

/// 曲線の平行移動を行いオブジェクトを生成する。
///Translation.
MGSPointSeq operator+ (const MGVector& ) const;

/// 与ベクトルだけ曲線を平行移動して自身とする。
///Translation.
MGSPointSeq& operator+= (const MGVector& );

/// 曲線の逆方向に平行移動を行いオブジェクトを生成する。
///Translation.
MGSPointSeq operator- (const MGVector& ) const;

/// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
///Translation.
MGSPointSeq& operator-= (const MGVector& );

/// 与えられたスケーリングで曲線の変換を行いオブジェクトを生成する。
///Scaling.
MGSPointSeq operator* (double) const;

/// 与えられたスケーリングで曲線の変換を行い自身の曲線とする。
///Scaling.
MGSPointSeq& operator*= (double);

/// 与えられた変換で曲線の変換を行いオブジェクトを生成する。
///Matrix transformation.
MGSPointSeq operator* (const MGMatrix& ) const;

/// 与えられた変換で曲線の変換を行い自身の曲線とする。
///Matrix transformation.
MGSPointSeq& operator*= (const MGMatrix& );

/// 与えられた変換で曲線のトランスフォームを行いオブジェクトを生成する。
///General transformation.
MGSPointSeq operator* (const MGTransf& ) const;

/// 与えられた変換で曲線のトランスフォームを行い自身とする。
///General transformation.
MGSPointSeq& operator*= (const MGTransf& );

///Add and subtract operation of two MGSPointSeq.
MGSPointSeq operator+ (const MGSPointSeq& sp2) const;
MGSPointSeq& operator+= (const MGSPointSeq& sp2);
MGSPointSeq operator- (const MGSPointSeq& sp2) const;
MGSPointSeq& operator-= (const MGSPointSeq& sp2);

///Compare two SpointSeq if they are equal.
///Return true if equal.
bool operator== (const MGSPointSeq& ) const;

///Return true if not equal.
bool operator!= (const MGSPointSeq& brep) const{return !(operator== (brep));}

////////////// Member Function //////////////

///compute an average plane of the point sequence.
///Function's return value is:
/// 1: Point seq is a point.		2: Point seq is on a line.
/// 3: Plane is output. 
int average_plane(
	MGPosition& center	///< center of point seq will be output.
	, MGPlane& plane	///< Plane will be output, when average_plane=3.
	, MGStraight& line	///< Straight line will be output            =2.
	, double& deviation	///< Maximum deviation from point, line, or plane will be output.
)const;					

///Compute minimum box including all the points.
MGBox box()const;

///Returns a pointer to the area.
const double* data(size_t i=0, size_t j=0, size_t k=0) const;

///Returns a pointer to the area.
double* data(size_t i=0, size_t j=0, size_t k=0);

///Transformation of own for rational(MGRLBRep) Control Polygon.
///Scaling.
MGSPointSeq& homogeneous_transform(double);
///Add the vector.
MGSPointSeq& homogeneous_transform(const MGVector&);
///Multiply the matrix.
MGSPointSeq& homogeneous_transform(const MGMatrix&);
///Multiply the transform.
MGSPointSeq& homogeneous_transform(const MGTransf&);

///Test if this is a null SPointSeq.
bool is_null()const{return sdim()==0;};

///Returns the actual size of Spoint seq.
void length(size_t& lengthu, size_t& lengthv) const
{lengthu=m_lengthu;lengthv=m_lengthv;}

///Returns the actual size of Spoint seq.
size_t length_u() const{return m_lengthu;}

///Returns the actual size of Spoint seq.
size_t length_v() const{return m_lengthv;}

///Generate data point abscissa(utau=u-direction and vtau=v-direction)
///from data point ordinates SPointSeq.
///SPointSeq(*this) must be homogeneous, i.e. if SPointSeq is positional
///data, all of them must be positional data, should not include
///derivative data.
void make_data_point(MGNDDArray& utau, MGNDDArray& vtau) const;

///Compute non_homogeneous coordonate data without w coordinate element,
///assumed that this is homogeneous coordinate data,
///i.e., maximum space dimension element is w(weight) coordinate.
///Result data does not include weight elements.
MGSPointSeq non_homogeneous() const;

///Return (i,j)-th element data.
/// When k>=sdim(), returns 0.0  .
double ref(size_t i, size_t j, size_t k=0) const;

///Update size. Update of sdim not allowed.
void reshape(
	  size_t sizeu		///<New size about u
	, size_t suzev		///<		   about v
	, size_t startu=0	///<To which place to store the original data
	, size_t startv=0	///< along u and v direction.
);

///Change the size of the array to
///m_lengthu=m_capacityu=lenu, m_lengthv=m_capacityv=lenv, m_sdim=dim.
///Result will contain garbages.
void resize(size_t lenu, size_t lenv,  size_t dim);

///Reverse the ordering of the points.
void reverse(
	int is_u				///<if true, u-drection. if not, v.
);

///Returns the space dimension.
size_t sdim() const {return m_sdim;}

///Set this as a null.
void set_null();

///Set the length of effective data.
void set_length(size_t lengthu,
				size_t lengthv );

///Returns the sizes along u and v direction.
void capacity(size_t& capau, size_t& capav) const{capau=m_capacityu; capav=m_capacityv;}

///Returns the size of u-direction.
size_t capacity_u() const {return m_capacityu;}

///Returns the size of v-direction.
size_t capacity_v() const {return m_capacityv;}

///Store vector data vector(from+r) to this(i,j,to+r) for 0<=r<sdim().
///When (form+r) or (to+r) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_at(
	size_t i, size_t j,	///<id of this which indicates the placement.
	const MGVector& vctr,	///<Input vector.
	size_t to=0,	///<Indicates to where of this in the space dimension id.
	size_t from=0	///<Indicates from where of vector in the space dimension.
);

///Store vector data vector(from+r) to this(i,j,to+r) for 0<=r<sdim().
///When (form+r) or (to+r) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_at(
	size_t i, size_t j,	///<id of this which indicates the placement.
	const MGVector& vctr,///<Input vector.
	size_t to,	///<Indicates to where of this in the space dimension id.
	size_t from,///<Indicates from where of vector in the space dimension.
	size_t len	///<Length of data to store.
);

///Store data[r] data to this(i,j,to+r) for 0<=r<sdim().
///When (to+r) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_at(
	size_t i, size_t j,	///<id of this which indicates the placement.
	const double* data,	///<Input data array.
	size_t to=0		///<Indicates to where of this in the space dimension id.
);

///Store BPointSeq along u at v's id j. That is, 
///data bp(i,from+r) to this(i,j,to+r) for i=0,..., length_u()-1, and  0<=r<sdim().
///When (form+r) or (to+r) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_BP_along_u_at(
	size_t j,	///<id of this which indicates the placement.
	const MGBPointSeq& bp,	///<Input vector.
	size_t to=0,	///<Indicates to where of this in the space dimension id.
	size_t from=0	///<Indicates from where of vector in the space dimension.
);

///Store BPointSeq along v at u's id i. That is, 
///data bp(j,from+r) to this(i,j,to+r) for j=0,..., length_v()-1, and  0<=r<sdim().
///When (form+r) or (to+r) reached to maximum space dimension id, next id
///becomes 0(form the start).
void store_BP_along_v_at(
	size_t i,	///<id of this which indicates the placement.
	const MGBPointSeq& bp,	///<Input vector.
	size_t to=0,	///<Indicates to where of this in the space dimension id.
	size_t from=0);	///<Indicates from where of vector in the space dimension.

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:

//////////////Member Data ////////////
	size_t m_capacityu;	///< maximum size to store Spoint seq of
	size_t m_capacityv;	///< space dimension m_sdim, for u and v direction.
	size_t m_sdim;		///< Space Dimension.
	size_t m_lengthu;	///< length of Spoint actual data stored about u.
	size_t m_lengthv;	///< length of Spoint actual data stored about v.
	double* m_spoint;	///< Area to store Spoint.

///Compute minimum box including all the points.
MGBox* compute_box()const;

friend class MGSBRep;
friend class MGRSBRep;

};

/** @} */ // end of BASE group
#endif
