/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCurve_HH_
#define _MGCurve_HH_

#include <vector>
#include <memory>
#include "mg/MGCL.h"
#include "mg/Default.h"
#include "mg/Position.h"
#include "mg/Geometry.h"
#include "mg/FSurface.h"
#include "mg/Pvector.h"
#include "mg/KnotVector.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"

//
//Define MGCurve Class.

class MGInterval;
class MGBox;
class MGVector;
class MGUnit_vector;
class MGPosition_list;
class MGTransf;
class MGCParam_list;
class MGPoint;
class MGStraight;
class MGEllipse;
class MGLBRep;
class MGRLBRep;
class MGSurfCurve;
class MGBSumCurve;
class MGTrimmedCurve;
class MGCompositeCurve;
class MGSurface;
class MGPlane;
class MGSphere;
class MGCylinder;
class MGSBRep;
class MGRSBRep;
class MGBSumSurf;
class MGCCisect_list;
class MGCSisect_list;
class MGIfstream;
class MGOfstream;
class MGFace;
class MGShell;
class MGCFisect_vector;
class MGPPRep;
class MGCommonON;

/** @addtogroup GEO
 *  @{
 */

///MGCurve is an abstract class which represents a whole curve.
class MGCLASS MGCurve:public MGGeometry{

public:
	
///向きが同じB表現曲線リストを接続する(LBRep, RLBRep同士のみ)。join_crvlに接続した曲線リストが入る。
///戻り値は、引数の曲線リストの向きが違うとき、同じB表現同士でなかったときfalseが返る。
friend int join(MGPvector<MGCurve>& crvl, MGPvector<MGCurve>& join_crvl);

//////////////Constructor///////////////

///Void constructor(初期化なしでオブジェクトを作成する。)
MGCurve();

///Copy constructor.
MGCurve(const MGCurve& curve);

//////////// Virtual Destructor ////////////
virtual ~MGCurve();

//////////// Operator overload(演算子多重定義) ////////////

///Assignment.
///When the leaf object of this and geo2 are not equal, this assignment
///does nothing.
virtual MGCurve& operator=(const MGCurve& gel2){MGGeometry::operator=(gel2);return *this;};

////////////Logical operator overload/////////

///Object transformation.
virtual MGCurve& operator+=(const MGVector& v)=0;
virtual MGCurve& operator-=(const MGVector& v)=0;
virtual MGCurve& operator*=(double scale)=0;
virtual MGCurve& operator*=(const MGMatrix& mat)=0;
virtual MGCurve& operator*=(const MGTransf& tr)=0;

///Comparison
virtual bool operator==(const MGCompositeCurve& crv)const;
virtual bool operator==(const MGTrimmedCurve& crv)const;
virtual bool operator==(const MGGel& gel2)const=0;
virtual bool operator<(const MGGel& gel2)const=0;

//////////// Member Function ////////////

///Generate arrow data of the tangent at the parameter value t of the curve.
///data[0] is the origin, data[1] is top of the arrow,
///data[2], [3] are two bottoms of arrowhead.
void arrow(double t,MGPosition data[4])const;

///Returns B-Rep Dimension.
virtual size_t bdim() const=0;

///Return minimum box that includes the curve of parameter interval.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
virtual MGBox box_limitted(
	const MGInterval& ///< Parameter Range of the curve.
) const = 0;

/// Calculate dividing Knots number for the initial approximation
///of the curve, used for precise approximation.
///分割数を求めるパラメータ範囲
virtual int calc_div_num(
	const MGInterval& interval)
const{return this->offset_div_num(interval);};

///Obtain ceter coordinate of the geometry.
virtual MGPosition center() const;

///Obtain ceter parameter value of the geometry.
virtual MGPosition center_param() const;

///Changing this object's space dimension.
virtual MGCurve& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
)=0;

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
virtual void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
)=0;

///Construct new geometry object by copying to newed area.
///User must delete this copied object by "delete".
virtual MGCurve* clone()const=0;

///Compute the closest point parameter value of this curve from a point.
virtual double closest(const MGPosition& point) const;

///Compute the point on the curve which is the intersection
///with x=point[0], or y=point[1] and the nearest point from point.
///If no intersection found, the nearer point parameter value out of start or end
///will be returned.
virtual double closest2D(const MGPosition& point,double& dist) const;

///Compute the closest point parameter value pair of this curve and curve2.
///MGPosition P of the function return contains this and curve2's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
virtual MGPosition closest(const MGCurve& curve2) const;

///曲線がCn連続かどうか調べる
///LBRep以外はかならずtrueが返却される
bool cn_continuity(size_t n)const;

///関数名：common
///目的：与えられた曲線と自身の共通部分があるかどうか調べる。
///引数：
///		const MGCurve&			curve2,		(I/ )	与えられる曲線
///		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
///		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
///					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
///				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
///		MGCCisect_list&			isect		( /O)	交点
///戻り値：
///		3:交点も共通部分も求まった
///		2:交点のみが求まった
///		1:共通部分のみが求まった
///		0:交点も共通部分もなかった
///		-1:共通エッジの収束計算エラー
///		-2:共通エッジが４個以上求まった(のっていないと見なす)
///追記：
///	曲線が共通かどうかの誤差にはline_zero()、をパラメータ範囲の収束計算の
///	誤差には、パラメータ範囲*rc_zero()を使用した
virtual int common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan,
	MGCCisect_list& isect
) const;

///関数名：common
///目的：与えられた曲線と自身の共通部分があるかどうか調べる。
///引数：
///		const MGCurve&			curve2,		(I/ )	与えられる曲線
///		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
///		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
///					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
///				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
///戻り値：
///		共通部分の数:	共通部分が求まった
///		0:				共通部分がなかった
///		-1:				共通エッジの収束計算エラー
///		-2:				共通エッジが４個以上求まった(のっていないと見なす)
///追記：
///	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
///  パラメータ範囲*rc_zero()を使用した
virtual int common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
virtual MGCurve& coordinate_exchange(size_t i, size_t j)=0;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
virtual MGCurve* copy_as_nurbs() const=0;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
virtual MGCurve* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const=0;

///Construct new curve object by copying to newed area,
///and limitting the parameter range to prange.
///Returned is a newed object and must be deleted.
virtual MGCurve* copy_limitted(const MGInterval& prange) const;

///Return curvature at the given point.
///When the curve is 2D, curvature has sign. when 3D, curvature is
///always plus.
/// 与えられた点における曲線の曲率を返却する。
virtual double curvature( double ) const;

///Return curve pointer if this MGGel is an MGCurve, else return null.
MGCurve* curve(){return this;};
const MGCurve* curve()const{return this;};

///線積分を求める。
///Compute curvilinear integral of the 1st two coordinates.
///This integral can be used to compute area sorounded by the curve.
///Second form is from param_s() to param_e();
///curvilinear_integral from t1 to t2 can be obtained by
///Integral of (x*dy-y*dx) about t, where curve is expressed by
///f(t)=(x(t),y(t)), dx=dx/dt, and dy=dy/dt.
virtual double curvilinear_integral(double t1, double t2) const;
virtual double curvilinear_integral() const{
	return curvilinear_integral(param_s(), param_e());
};

///Compute mean length of 1st derivative vector.
virtual double deriv_length()const;

///Compute direction unit vector of the geometry.
MGUnit_vector direction(const MGPosition& param) const;

///Return tangent vector at the given point.
/// 与えられた点における曲線の接ベクトルを返す。
virtual MGUnit_vector direction(double) const;

///////display member function.
virtual void display_arrows()const;
virtual void display_break_points()const;
virtual void display_curvatures(
	double	scale,	///<scaling of the graph.
	int		density,///<densitiy of the graph.
	bool	use_radius///<true:radius display, false:curvature display.
)const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
virtual int divide_multi(
	MGPvector<MGCurve>& crv_list,	///<divided curves will be appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at,
						///<When multiplicity<=0, order()-1 is assumed,
						///<When multiplicity>=order(), order() is assumed.
) const;

///get the a divide number for offset, intersection, or others.
int divide_number() const{return offset_div_num(param_range());};

virtual void drawSE(
	double span_length,	///<Line segment span length.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve,
						///<Draw will be performed from t0 to t1.
)const;

virtual void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const{drawSE(span_length,param_s(),param_e());};

///Return end point(終点を返却する)
virtual MGPosition end_point() const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
virtual MGVector eval(
	double,				///< Parameter value.
	size_t nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
) const = 0;

///Compute position, 1st and 2nd derivatives.
/// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
virtual void eval_all(
	double,			///<Input parameter value(パラメータ値)
	MGPosition&,	///<Position(位置)
	MGVector&,		///<1st derivative(1次微分値)
	MGVector&		///<2nd derivative(2次微分値)
) const;

///Compute 1st derivative.
/// 曲線上の与えられたパラメータ値における一次微分値をかえす。
virtual MGVector eval_deriv(double) const;

///Evaluate deviations of two curves(this and curve2) at npoint discrete points.
///(1)Search the common curve spans which have the distance within tolerance.
///(2)Compute the nearest points from npoint discrete points of this to curve2.
///Let sti=sts[i], then
///sti[0] is this curve's parameter value s, and sti[1] is the parameter value t
///of curve2 which is the nearest point from the point s.
///If this and curve2 have the minimum distance more than tolerance,
///sts.size()==1 and sts[0] is the minimum distance points of this and curve2.
void eval_discrete_deviation(
	const MGCurve& curve2,
	std::vector<MGPosition>& sts,
	int npoint=20,		///<indicates how many discrete points be obtained.
	double tolerance=0.1///<tolerance to get two edge to compute deviation.
)const;

///Evaluate line data at data point tau.
virtual void eval_line(
	const MGNDDArray& tau,	///<Data points.
	MGBPointSeq& value		///<Values evaluated. value(i,.)=eval(tau[i]);
)const;

///Compute positional data.
/// 与えられたパラメータ値に相当する自身上の点を返す。
virtual MGPosition eval_position(double) const;

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector evaluate(
	const MGPosition& t,	///< Parameter value.
				///<t's space dimension is geometry's manifold dimension.
	const size_t* nderiv=0	///<Order of derivative of i-th parameter
				///<in nderiv[i].
				///<When nderiv=null, nderiv[i]=0 is assumed for all i.
) const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
virtual void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
)=0;

///Compute Frenet_frame, curvature and torsion in 3D space.
virtual void Frenet_frame2(
	double t,			///<Input parameter value(パラメータ値)
	MGVector& V2,///<2nd derivative at t.
	MGVector& T,///<Tangent
	MGVector& N,///<Principal Normal
	MGVector& B	///<Binormal
)const;

///Compute Frenet_frame, curvature and torsion in 3D space.
virtual void Frenet_frame(
	double t,			///<Input parameter value(パラメータ値)
	MGUnit_vector& T,	///<Tangent
	MGUnit_vector& N,	///<Principal Normal
	MGUnit_vector& B,	///<Binormal
	double& curvature,	///<Curvature is always >=0.
	double& torsion
)const;

///Extracts control points.
///Fucntion's return value is 
///true if control points was obtained, false if not.
virtual bool get_control_points(
	MGBPointSeq& cpoints	///<Control points will be output.
)const{return false;};

///Test if this curve has the same direction with curve2 at the point s(of this)
/// and t(of curve2).
///Function's return value is true if they have the same direction.
///"same direction" means their tangent vectors have the angle less than 90 degree.
bool has_same_direction_at(
	double s,
	const MGCurve& curve2,
	double t
)const;

/// Return This object's typeID
virtual long identify_type() const=0;

///Test if input parameter value is inside parameter range of the line.
virtual bool in_range(double t) const;

///Test if input parameter value is inside parameter range of the line.
bool in_range(const MGPosition& t) const;

///Curve to curve intersection.
/// Curve と Curve の交点を求める。
MGCCisect_list intersect_brute_force(const MGCurve&) const;

///Curve to curve intersection.
///***Caution***intersect can be used only for finite curve, i.e.
///parameter range of the computation is only from param_s() to param_e().
///For example, intersect cannot be applied to infinite straight line.
virtual MGCCisect_list intersect(const MGCurve&) const;

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
MGisects intersection(const MGObject& obj2)const;
MGisects intersection(const MGCurve& obj2)const;
MGisects intersection(const MGFSurface& obj2)const;
MGisects intersection(const MGSurface& obj2)const;
MGisects intersection(const MGFace& obj2)const;
MGisects intersection(const MGShell& obj2)const;

///Provide divide number of curve span for function intersect.
virtual size_t intersect_dnum()const=0;

///intersections with a plane.
MGCSisect_list intersect_with_plane(const MGPlane& surf)const;

///Intersection of Curve and other geometry.
virtual MGCCisect_list isect(const MGCurve& curve2)const=0;
virtual MGCCisect_list isect(const MGStraight& curve2)const=0;
virtual MGCCisect_list isect(const MGRLBRep& curve2)const;
virtual MGCCisect_list isect(const MGEllipse& curve2)const;
virtual MGCCisect_list isect(const MGLBRep& curve2)const;
virtual MGCCisect_list isect(const MGSurfCurve& curve2)const=0;
virtual MGCCisect_list isect(const MGBSumCurve& curve2)const;
MGCCisect_list isect(const MGTrimmedCurve& curve2)const;
MGCCisect_list isect(const MGCompositeCurve& curve2)const;

virtual MGCSisect_list isect(const MGSurface& surf)const=0;
virtual MGCSisect_list isect(const MGPlane& surf)const=0;
virtual MGCSisect_list isect(const MGSphere& surf)const=0;
virtual MGCSisect_list isect(const MGCylinder& surf)const=0;
virtual MGCSisect_list isect(const MGSBRep& surf)const=0;
virtual MGCSisect_list isect(const MGRSBRep& surf)const=0;
virtual MGCSisect_list isect(const MGBSumSurf& surf)const=0;

MGCSisect_list isect(const MGFSurface& fs)const{return fs.isect(*this);};
virtual MGCSisect_list isect(const MGFace&)const;

///Intersection of a shell and a curve.
MGCFisect_vector isect(const MGShell& shl) const;

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
MGCParam_list isect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
) const;
	
/// 曲線が閉曲線かどうかを返す
///Test if this is a closed curve;
bool is_closed()const{return start_point()==end_point();};

///Test if this cure is co-planar with the 2nd curve curve2.
///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
virtual bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if the input parameter t is the start point parameter or not.
virtual bool is_startpoint_parameter(double t)const;

///Test if the input parameter t is the start point parameter or not.
virtual bool is_endpoint_parameter(double t)const;

///Test if the vector from P to this->eval(t) is perpendicular to
///the tangent of this curve at t.
bool is_perpendicular(const MGPosition& P, double t)const;

///Test if this cure is linear or not, that is, is straight or not.
///MGStraight expression will be out to straight if this is linear or not.
///Function's return value is true if linear.
virtual bool is_linear(MGStraight& straight)const;

///Test if this cure is planar or not.
///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
virtual bool is_planar(MGPlane& plane)const;

///Access to i-th element of knot
virtual double knot(size_t i) const=0;
	
///Returns the knot vector of the curve.
virtual const MGKnotVector& knot_vector() const=0;

///Returns the knot vector of the curve.
MGKnotVector& knot_vector();

///Cmpute curve length of the interval.
///If t1 is greater than t2, return negative value.
/// 与えられたパラメータ値間の曲線の長さを返す。
/// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
virtual double length(double t1, double t2) const;

///Compute whole curve length. If the curve is infinite, return -1.
/// 自身の曲線が有界の場合、その曲線の距離を返却する。非有界の場
/// 合はー１を返却をする。
virtual double length() const {return length(param_s(), param_e());}

///Inverse function of length. Compute the point that is away from
///the point t by length len.
/// lengthの逆関数。指定パラメータtで示される点から指定距離len
/// 曲線上に沿って離れた点を示すパラメータ値を返す。
virtual double length_param( double t, double len) const;

///Update this by limiting the parameter range of the curve.
/// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
virtual MGCurve& limit(const MGInterval& rng) = 0;
MGCurve& limit(double t0, double t1);

///Return manifold dimension, i.e. 0:point, 1:curve, 2:surface.
size_t manifold_dimension() const{ return 1;};

///Negate the curve direction(曲線の方向を反転する)
virtual void negate() = 0;

///Obtain parameter value if this curve is negated by "negate()".
virtual double negate_param(double t)const= 0;
///virtual MGPosition negate_param(const MGPosition& t)const;

///Transform the coordinates of boundary of this geometry so that
///new coordinate of boundary is the same coordinate as the new one of
///this geometry after negate() of this geometry is done.
///That is, boundary coordinates are parameter world of this geometry.
void negate_transform(MGGeometry& boundary)const;

///一定オフセット関数
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線リストが返却される。
///Offset of costant deviation from a curve.
///If the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
virtual MGPvector<MGCurve> offset(
	double ofs_value,							///<オフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///可変オフセット関数
///オフセット量は空間次元1の線B表現で与えられる。
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線リストが返却される。
///Offset of variable deviation from a curve.
///If the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
virtual MGPvector<MGCurve> offset(
	const MGLBRep& ofs_value_lb,				///<空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///C2連続曲線の一定オフセット関数
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線が返却される。
///costant offset curve of C2 continuous curve. if the norm_vector is given, the positive offset direction
///decide to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
virtual MGLBRep offset_c2(
	double ofs_value,							///<オフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///C2連続曲線の可変オフセット関数
///オフセット量は空間次元1の線B表現で与えられる。
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
///は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線が返却される。
///valuable offset curveof C2 continuous curve. if the norm_vector is given, the positive offset direction
///decide to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
virtual MGLBRep offset_c2(
	const MGLBRep& ofs_value_lb,				///<空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector = mgNULL_VEC	///<法線ベクトル
) const;

///オフセットで使用する、あるパラメータ範囲の分割数を求める
///get the number of division for offset
virtual int offset_div_num(
	const MGInterval& interval	///<分割数を求めるパラメータ範囲
)const;

///Test if given point is on the curve or not. If yes, return parameter
///value of the curve. Even if not, return nearest point's parameter t.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
virtual bool on(
	const MGPosition& point,///<point(指定点)
	double& t	///<Parameter of the curve(パラメータ) will be returned.
) const;

///Test if given point is on the geometry or not. If yes, return parameter
///value of the geometry. Even if not, return nearest point's parameter.
/// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
/// なくても最近傍点のパラメータ値を返す。
/// Function's return value is >0 if the point is on the geometry,
/// and 0 if the point is not on the geometry.
bool on(
	const MGPosition& P,///<Point(指定点)
	MGPosition&	t		///<Parameter of the geometry(パラメータ)
) const;

///Returns the order.
virtual	unsigned order() const=0;

///Compute parameter value of given point.
/// 自身の上の指定点を表すパラメータ値を返す。
/// If input point is not on the curve, return the nearest point on the
/// curve.
virtual double param(
	const MGPosition &	///<Point(指定点)
) const;

/// Return ending parameter value.
virtual double param_e() const=0;

///Obtain parameter space error.
virtual double param_error() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
virtual double param_normalize(double t) const=0;

///Return parameter range of the curve(パラメータ範囲を返す)
virtual MGInterval param_range() const;

///Round the parameter t into this parameter range.
double param_round_into_range(double t)const;

///Return parameter range of the geometry(パラメータ範囲を返す)
MGBox parameter_range() const;

/// Return starting parameter value.
virtual double param_s() const=0;

/// Return starting or ending parameter value that is nearer to the param t.
double param_se(double t) const;

///Compute parameter span length 
virtual double param_span() const{return param_e()-param_s();};

///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
virtual MGCurve* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const=0;

///Return perpendicular point from a point P,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual int perp_guess(
	double t0, double t1,	///<parameter range of this,
							///<(t0>=t1) indicates no range specified.
	const MGPosition& P,	///<Point(指定点)
	double tg,				///<Guess parameter values of this curve.
	double& t				///Output parameter
) const;

///Return perpendicular points of two curves,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual int perp_guess(
	double s0, double s1,	///<parameter range of this.
					///<When s0>=s1, no limit for this parameter range.
	const MGCurve& curve2,	///<2nd curve.
	double t0, double t1,	///<parameter range of curve2.
					///<When t0>=t1, no limit for curve2 parameter range.
	double sg, double tg,	///<Guess parameter values of the two curves
	///<sg: this curve's parameter, tg:curve2's parameter.
	MGPosition& st		///<perpendicular points' parameter values
						///<will be output.
	///<st(0): this curve's parameter, st(1):curve2's parameter.
) const;

///Compute a foot point of the perpendicular line from point p to
///the curve. If more than one points are found, return nearest one.
/// 指定点からの自身への垂線の足とパラメータ値を返す。
/// Function's return value is if point is obtained(1) or not(0)
virtual int perp_point(
	const MGPosition& p,	///<Point(指定点)
	double& t,				///<Parameter of the curve(パラメータ値)
	const double* g=0		///<guess parameter value of line
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///**** NOTE 1 ****
///perpendiculars is general function of perps, used in perps.
///General users should use function perps, not perpendiculars, since
///perps is optimized for each curve type.
///**** NOTE 2 ****
///perpendiculars can not be used for infinite parameter range curve.
///param_s() and param_e() of both curves must return their finite
///parameter range.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perpendiculars(
	const MGCurve& crv		///<The second curve
) const;

///Compute all foot points of the perpendicular line from point to
///the curve.
/// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
/// すべて求める。
virtual MGCParam_list perps(
	const MGPosition& P		///<Point(指定点)
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=curve2's parameter value.
virtual MGPosition_list perps(const MGCurve& crv2)const=0;
virtual MGPosition_list perps(const MGStraight& crv2)const=0;
virtual MGPosition_list perps(const MGRLBRep& crv2)const;
virtual MGPosition_list perps(const MGEllipse& crv2)const;
virtual MGPosition_list perps(const MGLBRep& crv2)const;
virtual MGPosition_list perps(const MGSurfCurve& crv2)const;
virtual MGPosition_list perps(const MGBSumCurve& crv2)const;
MGPosition_list perps(const MGCompositeCurve& crv2)const;
MGPosition_list perps(const MGTrimmedCurve& crv2)const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
virtual MGPosition pick_closest(const MGStraight& sl)const;

///Approximate this curve by a polyline and output to lb2.
///The tolerance of the approximation is error.
void polygonize(
	double error,	///<tolerance allowed for the approximation
	MGLBRep& lb2	///<Obtained polyline will be output as an MGLBRep of order2.
)const;

///Round t into curve's parameter range.
/// 入力パラメータをパラメータ範囲でまるめて返却する。
virtual double range(double t) const;

///Round t into geometry's parameter range.
/// 入力パラメータをパラメータ範囲でまるめて返却する。
///t's space dimension is geometry's manifold dimension.
MGPosition range(const MGPosition& t) const;

///Approcimate the curve by MGLBRep.
///The parameter of the original will not be changed.
std::auto_ptr<MGLBRep> MGCurve::rebuild(
	size_t order,		///<order of the new MGLBRep, must be >=4.
	double tol		///<tolerance allowed for the approximation
)const;

///ノット削除関数(B表現曲線のみ)
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///Remove redundant knot, and reduce the b-rep dimension.
///The tolerance used is MGTolerance::line_zero().
virtual void remove_knot();

///Update curve by rotating around straight line.
/// 指定点を通る指定ベクトルを軸として回転させたものを自身とする。
virtual MGCurve& rotate_self(
	const MGVector& v,			///<Vector of the line to rotate around.
	double,						///<Angle of rotation.
	const MGPosition & = mgORIGIN	///<A point on the line to rotate around.
);

///Return space dimension
virtual size_t sdim() const =0 ;

///Return start point(始点を返却する)
virtual MGPosition start_point() const;

///Return sweep surface from crv.
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
virtual MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
) const =0;

///Return tangent point from a point P,
///given guess starting paramter tg.
virtual int tangent_guess(
	double t0, double t1,	///<parameter range of this.
	const MGPosition& P,	///<Point(指定点)
	double tg,				///<Guess parameter values of the two curves
	double& t				///<Output parameter
)const;

///Trim the end part of this curve at the parameter t.
///The new curve range is [start_of_original, t]
///t must be inside this parameter rage, else does nothing.
void trim_end(double t);

///Trim the start part of this curve at the parameter t.
///The new curve range is [t,end_of_original]
///t must be inside this parameter rage, else does nothing.
void trim_start(double t);

///Return curve type(曲線のタイプを返す)
virtual MGCURVE_TYPE type() const =0;

///Unlimit parameter range of the curve(limitをはずす)
virtual MGCurve& unlimit() =0;

///Unlimit parameter range of the curve to the end point direction
///(終点方向にlimitをはずす)
virtual MGCurve& unlimit_end() =0;

///Unlimit parameter range of the curve to the start point direction
///(始点方向にlimitをはずす)
virtual MGCurve& unlimit_start() =0;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

protected:

///メンバデータを読み出す関数
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
virtual void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
virtual void WriteMembers(MGOfstream& buf) const;

///virtual std::vector<MGInterval> clip(const MGBox&) const;

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
virtual MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

virtual std::string whoami()const{return "Curve";};

protected:

///Approximate this curve as a MGLBRep curve from knot_vector[is] to [ie].
///This is an internal program of MGLBRep constructor.
void approximate_as_LBRep(
	MGLBRep& lb,		///<Approximated obrep will be set.
	size_t order,		///<new order
	size_t is, size_t ie///<approximation parameter range, from lb.knot_vector()[is] to [ie].
)const;

///Obtain an extrapolated PP-Rep curve by the parameter value.
void extrapolated_pp(
	double tau,		///<The parameter value at the end of extended point,
					///<When tau<param_s(), extension will be done at the starting point,
					///<When tau>param_e(), extension will be done at the end point.
	double dk,     ///<Coefficient of how curvature should vary at the connecting point.
	MGPPRep& pp
)const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
virtual MGCCisect_list isect_withC1LB(const MGLBRep& curve2)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
virtual MGCCisect_list isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
virtual std::auto_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const=0;

///Perpendicular points with C1 conitnuity LBRep.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
virtual MGPosition_list perps_withC1LB(
   const MGLBRep& lbC1
)const;

///Perpendicular points with SurfCurve
///whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
virtual MGPosition_list perps_with_noCompoSC(const MGSurfCurve& curve2)const;

///Perpendicular points with straight.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perpsSl(
	const MGStraight& sl	///<The second curve
)const;

///Mark this as updated.
virtual void update_mark(){	MGGeometry::update_mark();}

private:

///thisのスタートパラメータを与え共通パラメータ範囲を1つだけ取得する
///sparamはスタートパラメータであり、検索終了時のパラメータを返す
///戻り値は	1:範囲が求まった(thisの最後まで)
///			0:範囲が求まらなかった(thisの最後まで)
///			-1:範囲が求まった(thisの途中まで、まだ共通範囲があるかもしれない)
int common_one_span(
	const MGCurve& curve2,
	MGCommonON SEon[2],///<data if curve2's start or end point is on this curve.
	MGCommonON &sparam, ///<Input start parameter of this curve to search starting point of the next
					///<common part.
					///<On return from common_one_span, the next startign parameter will be set.
	double span,	///<parameter span to advance the next on point check of this curve.
	double comSpan[4]
)const;

///this上のcurve2に乗っているパラメータonparamと乗っていないパラメータoffparam
///を与え、ちょうどはずれる点を求める。
///チェックポイントの移動量 < (パラメータ範囲*rc_zero())になれば終了。
///戻り値は　収束して求まったパラメータ値である。
void common_boundary(
	const MGCurve& curve2,
	const MGCommonON& onparam,
	const MGCommonON& offparam,
	double& param1, double& param2
)const;

///Return minimum box that includes whole of the curve.
///曲線部分を囲むボックスを返す。
virtual MGBox* compute_box() const=0;

///法線ベクトルが指定されているときの処理
int offset_norm_proc(
	const MGLBRep& ofs_value_lb,		///<オフセット量
	const MGVector& norm_vector,		///<法線ベクトル
	MGPvector<MGCurve>& ofs_crvl	///<オフセットカーブリスト
)const;

///法線ベクトルが指定されていないときの処理
int offset_proc(
	const MGLBRep& ofs_value_lb,	///<オフセット量
	MGPvector<MGCurve>& ofs_crvl	///<オフセットカーブリスト
)const;

///法線ベクトルが指定されているC2連続曲線のオフセット
int offset_norm_c2_proc(
	const MGLBRep& ofs_value_lb,	///<オフセット量
	const MGVector& norm_vector,	///<法線ベクトル
	MGLBRep& ofs_brep		///<オフセットカーブ
)const;

///法線ベクトルが指定されていないC2連続曲線のオフセット
int offset_c2_proc(
	const MGLBRep& ofs_value_lb,	///<オフセット量
	MGLBRep& ofs_brep,				///<オフセットカーブ
	MGUnit_vector& preN,			///<前回のノーマルベクトル
	int& freverse			///<向きを逆にしているというフラグ
)const;

///曲線を折れで分割する(オフセット量を示す曲線の折れついても分割する)
int divide_multi_ofs(
	const MGLBRep& ofs_value_lb,		///<オフセット量を示す曲線
	MGPvector<MGCurve>& brep_list		///<分割した曲線リスト
)const;

///2本のB表現曲線を接続する(同じ種類のとき)
///join LBRep or RLBRep respectively.
///virtual MGCurve* join(const MGCurve& crv1) const;

///曲線をオフセットするのに十分分割したノットベクトルを返却する
///オフセット量曲線も考慮に入れ、分割数の多い方にあわせている。
MGKnotVector offset_make_knotvector(const MGLBRep& ofs_value_lb)const;

friend class MGFace;
friend class MGFSurface;
friend class MGSurface;
friend class MGLBRep;
friend class MGTrimmedCurve;
friend class MGSurfCurve;
friend class MGCompositeCurve;
friend class MGBSumCurve;

};

///@cond

//The class for function object for mgGausp to compute the length() of the curve.
class MGCLASS MGCurveLengthDrive{
	const MGCurve* m_curve;
public:
	MGCurveLengthDrive(const MGCurve* curve):m_curve(curve){;};
	double operator()(double t) const;
};

//The class for function object for mgDefint to compute the length_param() of the curve.
class MGCLASS MGCurveLenParamDrive {
	const MGCurve* m_curve;
	double m_len, m_ts;
public:
	MGCurveLenParamDrive(const MGCurve* curve, double len, double ts)
		:m_curve(curve), m_len(len), m_ts(ts){;};
	double operator()(double t)const;
};

///@endcond

///Compute curvature in 3D space, ie, the value is not negative.
MGEXTERN double MG_Curvature(
		const MGVector& v1,		///<First derivative.
		const MGVector& v2	///<Second derivative.
);

///Compute torsion.
MGEXTERN double MG_Torsion(
		const MGVector& v1,		///<First derivative.
		const MGVector& v2,		///<Second derivative.
		const MGVector& v3	///<Third derivative.
);

///Generate arrow data from (root, vecx, vecy).
MGDECL void one_arrow(
	const MGPosition& root,	///<root of the arrow
	const MGVector& vecx,	///<the vector from the roo to the head of the arrrow
	const MGUnit_vector& vecy,///<vecy that is normal to the vector from root to head
	MGPosition& head,		///<head of the arrow will be returned.
	MGPosition& headtail1,	///<two tail of arrowhead line segments will be returned.
	MGPosition& headtail2
);

namespace MGCL{

/// @brief  Creates a curve that has weight.
/// @param  curve 曲線オブジェクト
///Returned object is a newed object. User must delete it.
MGRLBRep* convert_to_rational(const MGCurve& curve);

}

/** @} */ // end of GEO group
#endif
