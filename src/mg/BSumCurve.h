/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGBSumCurve_HH_
#define _MGBSumCurve_HH_

#include "mg/MGCL.h"
#include "mg/Curve.h"

class MGStraight;
class MGRLBRep;
class MGEllipse;
class MGLBRep;
class MGSurfCurve;
class MGBSumCurve;

/** @addtogroup GEO
 *  @{
 */

///Define MGBSumCurve Class(Boolean sum curve of three curves).
///MGBSumCurve is a curve to support the class MGBSumSurf, is boolean sum of
///three curves, m_g1, m_g2, and m_g12 whose parameter ranges are exactly the same.
///MGBSumCurve(t) is defined as m_g1(t)+m_g2(t)-m_g12(t).
class MGCLASS MGBSumCurve:public MGCurve{

MGDECL friend MGBSumCurve operator+ (const MGVector& v, const MGBSumCurve& lb);
MGDECL friend MGBSumCurve operator* (double scale, const MGBSumCurve&);

public:
/////////////Constructor///////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGBSumCurve();

///Copy constructor.
MGBSumCurve(const MGBSumCurve& curve);

///constructor of three curves.
///The ownership of g1, g2, and g12 will be transfered to this MGBSumCurve.
MGBSumCurve(MGCurve* g1, MGCurve* g2, MGCurve* g12);

///constructor of three curves.
MGBSumCurve(const MGCurve& g1, const MGCurve& g2, const MGCurve& g12);

//////////// Virtual Destructor ////////////
~MGBSumCurve();

//////////// Operator overload(演算子多重定義) ////////////

///Assignment.
///When the leaf object of this and geo2 are not equal, this assignment
///does nothing.
MGBSumCurve& operator=(const MGGel& gel2);
MGBSumCurve& operator=(const MGBSumCurve& el2);

///Object transformation.
MGBSumCurve& operator+=(const MGVector& v);
MGBSumCurve& operator-=(const MGVector& v);
MGBSumCurve& operator*=(double scale);
MGBSumCurve& operator*=(const MGMatrix& mat);
MGBSumCurve& operator*=(const MGTransf& tr);

///Transformation object construction
MGBSumCurve operator+ (const MGVector& v) const;
MGBSumCurve operator- (const MGVector& v) const;
MGBSumCurve operator* (double scale) const;
MGBSumCurve operator* (const MGMatrix& mat) const;
MGBSumCurve operator* (const MGTransf& tr) const;

///Comparison of two curves.
bool operator==(const MGBSumCurve& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGBSumCurve& gel2)const;
bool operator<(const MGGel& gel2)const;

//////////// Member Function /////////////

///Returns B-Rep Dimension.
size_t bdim() const;

///Return minimum box that includes the curve of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox box_limitted(
	const MGInterval& rang///< Parameter Range of the curve.
)const;

///Changing this object's space dimension.
MGBSumCurve& change_dimension(
	size_t sdim,	///< new space dimension
	size_t start1=0,///< Destination order of new object.
	size_t start2=0 ///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
);

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
MGBSumCurve& coordinate_exchange(size_t i, size_t j);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGBSumCurve* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be an MGRLBRep.
///Otherwise,  the new curve will be an MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGBSumCurve* copy_change_dimension(
	size_t sdim,	///< new space dimension
	size_t start1=0,///< Destination order of new line.
	size_t start2=0 ///< Source order of this line.
)const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector eval(
	double t,		///< Parameter value.
	size_t nderiv=0,///< Order of Derivative.
	int left=0		///<Left continuous(left=true)
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
long identify_type()const{return MGBSUMCRV_TID;};

///Intersection of MGBSumCurve and curve.
///MGBSumCurve と Curve の交点を求める。
MGCCisect_list isect(const MGCurve&) const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;

///Intersection with a Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;

///Access to i-th element of knot
double knot(size_t i) const;
		
///Returns the knot vector of the curve.
const MGKnotVector& knot_vector() const;

///Update this by limiting the parameter range of the curve.
/// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
MGBSumCurve& limit(const MGInterval& rng);

///Negate the curve direction(曲線の方向を反転する)
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///Returns the order.
unsigned order() const;

/// Return ending parameter value.
double param_e() const{return m_g1->param_e();};

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

/// Return starting parameter value.
double param_s() const{return m_g1->param_s();};

///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to the newed object, and should be deleted
///by calling program, or memory leaked.
MGBSumCurve* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;

///Return space dimension
size_t sdim() const{return m_g1->sdim();};

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
)const;

///Return curve type(曲線のタイプを返す)
MGCURVE_TYPE type() const{return MGCURVE_BSUM;};

///Unlimit parameter range of the curve(limitをはずす)
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
///(終点方向にlimitをはずす)
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
///(始点方向にlimitをはずす)
MGCurve& unlimit_start();

///Output to IGES stream file.
///BSumCurve is approximated as MGLBRep and output as PD126.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

/// Output  function.
std::ostream& out(std::ostream& ostrm) const;

protected:

///メンバデータを読み出す関数
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf)const;

///Provide divide number of curve span for function intersect.
size_t intersect_dnum()const;

std::string whoami()const{return "BSumCurve";};

protected:

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
)const;

private:

	MGCurve* m_g1;
	MGCurve* m_g2;
	MGCurve* m_g12;

///Return minimum box that includes whole of the curve.
///曲線部分を囲むボックスを返す。
MGBox* compute_box() const;

};

/** @} */ // end of GEO group
#endif
