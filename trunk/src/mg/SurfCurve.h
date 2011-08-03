/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSurfCurve_HH_
#define _MGSurfCurve_HH_

#include "mg/Default.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/TrimmedCurve.h"

//
//Define MGSurfCurve Class.

class MGInterval;
class MGBox;
class MGVector;
class MGUnit_vector;
class MGPosition;
class MGPosition_list;
class MGTransf;
class MGCParam_list;
class MGStraight;
class MGEllipse;
class MGRLBRep;
class MGRSBRep;
class MGCCisect_list;
class MGCSisect_list;
class MGIfstream;
class MGOfstream;

/** @addtogroup GEO
 *  @{
 */

///MGSurfCurve is a curve defined by a surface and its parameter space line
///represented by (u,v). Let the surface S(u,v), and the parameter space
///curve of 2D be f(t). Then MGSurfCurve is defined as S(f(t)).
///MGSurfCurve is a TEMPORAL curve, and does not have update functions.
class MGCLASS MGSurfCurve:public MGCurve{

public:

MGDECL friend MGSurfCurve operator+ (const MGVector& v, const MGSurfCurve& cv2);
MGDECL friend MGSurfCurve operator* (double scale, const MGSurfCurve& cv2);

////////Constructor/////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGSurfCurve():MGCurve(),m_surface(0){;};

///Copy constructor
MGSurfCurve(const MGSurfCurve& sc);

///Default constructor
MGSurfCurve(
	const MGSurface& srf, ///< Surface
	const MGCurve& crv	///< Original uv-curve
); 

///Default constructor
MGSurfCurve(
	const MGFSurface& srf, ///< Surface
	const MGCurve& crv	///< Original uv-curve
); 

////////////Operator overload(演算子多重定義)////////////

public:

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGSurfCurve& operator=(const MGGel& gel2);
MGSurfCurve& operator=(const MGSurfCurve& el2);

///Comparison of two curves.
bool operator==(const MGSurfCurve& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGSurfCurve& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function/////////

///Test if m_curve is MGCompositeCurve. If composite, return
///the pointer. If not, return null.
const MGCompositeCurve* base_composite()const;

///Returns B-Rep Dimension.
size_t bdim() const	{return m_curve.bdim();}

///Return minimum box that includes the curve of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox box_limitted(
	const MGInterval & ///< Parameter Range of the curve.
) const;

///Changing this object's space dimension.
MGSurfCurve& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0		///< Source order of this object.
); 

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,		///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original.
){	assert(false);}

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGSurfCurve* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
///This function returned MGLBRep*(order = 4)
MGCurve* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

/// Return parameter curve's pointer
const MGCurve* curve() const {return &m_curve;};

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector eval(
	double,				///< Parameter value.
	size_t nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
)const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

/// Return This object's typeID
long identify_type() const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(double t) const;

///Provide divide number of curve span for function intersect.
size_t intersect_dnum() const;

///Intersection of Curve
MGCCisect_list isect(const MGCurve& curve2) const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGRLBRep& curve2)const;
MGCCisect_list isect(const MGEllipse& curve2)const;
MGCCisect_list isect(const MGLBRep& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;
MGCCisect_list isect(const MGBSumCurve& curve2)const;

///Intersection with a Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;

///Access to i-th element of knot.
double knot(size_t i) const;
	
///Returns the knot vector.
const MGKnotVector& knot_vector() const{return m_curve.knot_vector();}

///Update this by limiting the parameter range of the curve.
/// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
MGSurfCurve& limit(const MGInterval& );

///Returns the order.
unsigned order() const{return m_curve.order();}

/// Return ending parameter value.
double param_e() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

///Return parameter range of the curve(パラメータ範囲を返す)
MGInterval param_range() const;

/// Return starting parameter value.
double param_s() const;
	
///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGSurfCurve* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;
MGPosition_list perps(const MGRLBRep& crv2)const;
MGPosition_list perps(const MGEllipse& crv2)const;
MGPosition_list perps(const MGLBRep& crv2)const;
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const;

///Round t into curve's parameter range.
/// 入力パラメータをパラメータ範囲でまるめて返却する。
double range(double t) const;

///Return space dimension
size_t sdim() const ;

/// Return surface's pointer
const MGSurface* surface() const {return m_surface;};

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,		///<Sweep Direction.
	double start_dist,				///<distance to start edge.
	double end_dist) const;			///<distance to end edge.

///Unlimit parameter range of the curve(limitをはずす)
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
///(終点方向にlimitをはずす)
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
///(始点方向にlimitをはずす)
MGCurve& unlimit_start();

///Return curve type(曲線のタイプを返す)
MGCURVE_TYPE type() const;

/// Output function.
std::ostream& out(std::ostream&) const;

///IGES output function. PD126(approximated as a NURBS line).
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

protected:

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
///For SurfCurve this is not allowed to use.
std::auto_ptr<MGCurve> oneD(
	const double g[4]	///<Plane expression(a,b,c,d) where ax+by+cz=d.
)const{assert(false); std::auto_ptr<MGCurve> tmpPtr; return tmpPtr;}

///メンバデータを読み出す関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void WriteMembers(MGOfstream& buf) const;

std::string whoami()const{return "SurfCurve";};

private:

////////////Member Data//////////

	const MGSurface* m_surface;
	MGTrimmedCurve m_curve;

///Return minimum box that includes whole of the curve.
///曲線部分を囲むボックスを返す。
MGBox* compute_box() const;

//////////Following operators are declared to prohibit the use,
//////////since these cannot be provided as member functions.

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
MGSurfCurve& coordinate_exchange(size_t i, size_t j);

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list isect_noCompo(const MGCurve& curve2)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list isect_noCompo(const MGStraight& curve2)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list isect_noCompo(const MGSurfCurve& curve2)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisect_list isect_noCompo(const MGSurface& surf)const;

///isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisect_list isect_noCompo(const MGPlane& surf)const;

///isect of each elements of this m_curve,
///if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void isect_of_each(
	const MGCurve& curve2,	///<The isect objective curve.
	MGCCisect_list& list	///<Obtained isect will be appended.
)const;

///isect of each elements of this m_curve,
///if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void isect_of_each(
	const MGSurface& surf,	///<The isect objective surface.
	MGCSisect_list& list	///<Obtained isect will be appended.
)const;

///Negate the curve direction(曲線の方向を反転する)
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;
MGPosition negate_param(const MGPosition& t)const;

///perps of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list perps_noCompo(const MGCurve& curve2)const;

///perpendicular points of each elements of this m_curve,
///if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void perps_of_each(
	const MGCurve& curve2,	///<The perps objective curve.
	MGPosition_list& list	///<Obtained perpendicular points will be appended.
)const;

///Object transformation.
MGSurfCurve& operator+=(const MGVector& v);
MGSurfCurve& operator-=(const MGVector& v);
MGSurfCurve& operator*=(double scale);
MGSurfCurve& operator*=(const MGMatrix& mat);
MGSurfCurve& operator*=(const MGTransf& tr);

///Transformation object construction(NOT ALLOWED TO USE).
MGSurfCurve operator+ (const MGVector& v) const;
MGSurfCurve operator- (const MGVector& v) const;
MGSurfCurve operator* (double scale) const;
MGSurfCurve operator* (const MGMatrix& mat) const;
MGSurfCurve operator* (const MGTransf& tr) const;

friend class MGLBRep;
friend class MGStraight;
friend class MGPlane;
};

/** @} */ // end of GEO group
#endif
