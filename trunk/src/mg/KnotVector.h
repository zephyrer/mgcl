/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGKnotVector_HH_
#define _MGKnotVector_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/NDDArray.h"

// MGKnotVector.h
//
// Forward declaration
class MGNDDArray;
class MGKnotArray;
class MGIfstream;
class MGOfstream;

/// Defines Knot vector of B-Representation.
///MGKnotVector is to represent a knot vector of a B-Spline.
///It is an array of double precision folating data that is non decreasing.
///Let n be a B-representaiton dimension, i.e., number of points of the control polygon.
///And let k be the order( or (k-1) be degree) of the B-rep, then the length of
///the kont vector is (n+k).
class MGCLASS MGKnotVector :public MGNDDArray{

public:

MGDECL friend MGKnotVector operator*(double scale, const MGKnotVector& t);

///String stream Function
MGDECL friend std::ostream& operator<<(std::ostream& strm, const MGKnotVector& t);

////////////Constructor////////////

/// Construct garbage knot vector of specified size.
explicit MGKnotVector(
	unsigned int order=0,	///< order
	unsigned int bdim=0,	///< B-Rep dimension
	const double* data=0	///< all of the knot data sequence if data!=NULL.
);

/// Construct uniform knot vector with inital value init_value
/// and incremental of increment.
MGKnotVector(
	unsigned int order,		///< order
	unsigned int bdim,		///< B-Rep dimension.
	double init_value,		///< Initial value
	double increment=1.0
);

///Obtain knot vector from Data Point.
MGKnotVector(
	const MGNDDArray& dtp,	///< input data point
	unsigned int order	///< order
);

///Add knots into original knot vector.
///Maximum multiplicity of the knot is guaranteed to be the order.
MGKnotVector(
	const MGKnotVector& vec2,	///<Original Knot Vector
	const MGKnotArray& knots	///<Knots to add with multiplicity.
);

/// Construct by extracting sub interval of vec2.
MGKnotVector(
	size_t start_id,			///<Start id of vec2(from 0).
	size_t num,					///<new B-Rep dimension.
	const MGKnotVector& vec2	///<Original knot vector.
);

///Construct new order knot vector. 
MGKnotVector(
	const MGKnotVector& knotv,	///< input knot vector
	size_t order				///< new order
);

///Generate knot vector as a part of original knotv.
///New knot vector's parame range is from t1 to t2.
///Should hold t1>=knotv.param_s() && t2<=knotv.param_e().
///Knots between t1< <t2 are copied from knotv.
MGKnotVector(const MGKnotVector& knotv, double t1, double t2);

/// Copy Constructor
MGKnotVector(const MGKnotVector& vec2);

//Destructor
//
//	~MGKnotVector();	//We use default destructor.

//////////// Operator overload ////////////

///Assignment
MGKnotVector& operator=(const MGKnotVector& vec2);

///Addition and subtraction of real number.
///All of the elements will be added or subtracted.
MGKnotVector operator+ (double) const;
MGKnotVector& operator+= (double);
MGKnotVector operator- (double) const;
MGKnotVector& operator-= (double);

/// 単項マイナス。
///Unary minus. Reverse the ordering of elements by changing all of the
///signs.
MGKnotVector operator- () const;

///Scaling.
MGKnotVector operator* (double) const;
MGKnotVector& operator*= (double);

///Compare two KnotVector if they are equal.
bool operator== (const MGKnotVector& ) const;	 ///Return true if equal.
bool operator!= (const MGKnotVector& knotv) const ///Return true if not equal.
{return !(operator== (knotv));}

////////////Member Function////////////

size_t bdim() const{return length()-m_order;};	///Return BRep Dimension

///Update array length.
///Updated array is so generated that the original proportions of
///neighbors hold as much as possible. 
///If original data point has multiplicities and nnew>=length(),
///original data point parameters and the multiplicities are preserved.
MGKnotVector& change_knot_number(size_t nnew);

///Change order
///For area adjustment. New knot vector possesses the same values for
///the area of the old one as long as area permitts.
///New vector is generally garbage knot vector.
MGKnotVector& change_order(unsigned order);	

/// Change parameter range.
void change_range(double ts, double te);

///Divide every spans. Every spans are subdivided into num equal spans.
///Result knot vector's bdim() becomes approximately bdim()*num.
void divide_span(size_t num);

/// Function's return value id is the index of B-coefficients that
/// should be multiplied to.
/// coef[j] is for (id+j)-th B-Coefficients, 0<= j <=order-1.
/// Multiplication with coef should be done like:
////
/// data=0.;
/// for(j=0; j<order; j++) data+=coef[j]*B_coef[id+j].
////
///left indicates whether left continuous(left=true), or right continuous
///(left=false) evaluation.
int eval_coef(
	double x,		///< Parameter value to evaluate
	double* coef,	///< coef's to multiply are returned.
					///< array of length order.
	unsigned nderiv=0,///<order of derivative, =0 means position
	int left=0		///<Left continuous(left=true)
					///<or right continuous(left=false).
)const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(double t) const;

///Test if this is null.
bool is_null()const{return MGNDDArray::is_null();};

/// Find where tau is located in knot. Return id is to evalute,
/// i.e. (*this)(id) < (*this)(id+1).
/// Returned id is :
///      order()-1<= id < bdim().
///left indicates whether left continuous(left=true), or right continuous
///(left=false) location.
int	locate(double tau, int left) const;
int	locate(double tau) const;

///Locate where data of multiplicity of multi is after start and before
///bdim(). index is the starting point index of this found first
///after start, index>=start.
///Function's return value locate_multi is actual multiplicity at the
///index, i.e. locate_multi>=multi if found.
///If position of the multiplicity is not found before bdim(), index=bdim()
///(index of the param_e()) and locate_multi=0 will be returned.
///multi must be >=1 and start must be <=bdim().
size_t locate_multi(size_t start, size_t multi, size_t& index) const;
	
///ノットベクトルをショートスパンのないように足しあわせる
///ノットベクトルのパラメータ範囲は等しいものとする。
///エラーのとき元のノットベクトルを返す.
///mixing two knotvectors.
///Mixing is so done as that no too close points are not included.
MGKnotVector& mix_knot_vector(const MGKnotVector& knot);

///Return order
unsigned int order() const{return m_order;};

/// Return end parameter value.
double param_e() const{return MGNDDArray::operator() (bdim());};

/// Return start parameter value.
double param_s() const{return MGNDDArray::operator() (m_order-1);};

///Return tolerance allowed in knot vector parameter space.
double param_error()const;

///Compute parameter span length.
double param_span() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///Round the input parameter value t into 
///the parameter range of the knot_vector
double range(double t) const;

///Reverse the ordering of knots.
///Parameter range of the original is preserved.
void reverse();

///Obtain parameter value if this knot vector is reversed by reverse().
double reverse_param(double t) const;

///Set B-rep dimension.
///Only set the B-Rep Dimension. Result knot vector may be garbage.
MGKnotVector& set_bdim(size_t brdim);

///Set this as a null.
void set_null(){MGNDDArray::set_null();m_order=m_current=0;};

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Resize this so that this order is order and this b-rep dimension is bdim.
///Result of size_change will contain garbages.
void size_change(size_t order, size_t bdim);

///Restore Function
int restore(MGIfstream& );

/////////Member data/////////

private:
	unsigned int   m_order;		///< Order of B-Representation.

void new_knot_coef(const MGKnotVector&, ///< old knot vector.
	size_t j_of_t, ///< Index of old Knot vector to get coef.
	double* coef ///< coef's to mutiply are returned.
);

};

/** @} */ // end of BASE group
#endif
