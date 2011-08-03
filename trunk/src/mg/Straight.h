/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGStraight_HH_
#define _MGStraight_HH_

#include "mg/EReal.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Curve.h"

// MGStraight.h
// Definition of MGStraight class

class MGBox;
class MGInterval;
class MGVector;
class MGTransf;
class MGCCisect;
class MGCSisect;
class MGRLBRep;
class MGEllipse;
class MGLBRep;
class MGSurfCurve;
class MGBSumCurve;
class MGCompositeCurve;
class MGTrimmedCurve;
class MGKnotVector;
class MGPlane;
class MGSurface;
class MGCCisect_list;
class MGPosition_list;
class MGIfstream;
class MGOfstream;

/** @addtogroup GEO
 *  @{
 */

/// MGStraight is a curve of any space dimension, represent a straight line.
/// Parameterization of the MGStraight is as below using parameter t:
/// point(t) = m_root_point + t*m_direction.
///MGStraight can be a line segment that has start and end points.
///Let t be the parameter of the straight line(i.e. length from m_root_point),
///then straight lien f(t) can be expressed as:
///  f(t)=m_root_point+m_direction*t.
class MGCLASS MGStraight: public MGCurve{
public:

MGDECL friend MGStraight operator+ (const MGVector& v, const MGStraight& sl);
MGDECL friend MGStraight operator* (double scale, const MGStraight&);

//////////// Conctructor. コンストラクタ ////////////
	
///Void constructor
MGStraight();

///Copy constructor.
MGStraight(const MGStraight& sl);

///Straight specifying all the member data. All of the data are employed as 
///member data of this straight.
explicit MGStraight (
	const MGEReal& endparam,	///<end parameter value
	const MGEReal& sparam,		///<start parameter value
	const MGVector& direction,	///<Direction, which will be the direction of this.
	const MGPosition& Origin	///<Origin
);
explicit MGStraight(
	double endparam,			///<end parameter value
	double sparam,				///<start parameter value
	const MGVector& direction,	///<Direction, which will be the direction of this.
	const MGPosition& Origin	///<Origin
);

/// 直線のタイプ，方向ベクトル，始点を指定して直線を生成する。
///Straight from straight line type, direction vector, and an origin.
///This constrcutor converts input vec to a unit vector.
///If you do not like the conversion, use set_straight().
MGStraight(
	MGSTRAIGHT_TYPE type,				///<Type
		///<Straight line type(直線の種類)
		///<enum MGSTRAIGHT_TYPE {
		///<	MGSTRAIGHT_EMPTY 		//Empty. 空
		///<	,MGSTRAIGHT_SEGMENT		//Line segment. 線分
		///<	,MGSTRAIGHT_HALF_LIMIT	//Half unlimit. 半直線
		///<	,MGSTRAIGHT_UNLIMIT		//Unlimit line for both direction. 無限直線
		///<};
	const MGVector& vec,				///<Direction
	const MGPosition & = mgORIGIN		///<Origin
);

/// ２点から直線を生成する。
///MGSTRAIGHT_SEGMENT straight from two points.
///Start point is start and end point is end.
///Parameter value of the start point is set to be 0.
MGStraight(
	const MGPosition& end,		///<End point.
	const MGPosition& start		///<Start point.
);

///始終点の座標、パラメータ値から直線を生成する。
///MGSTRAIGHT_SEGMENT straight from two points.
///Start point is start and end point is end.
///In this version, can specify start and end parameter values.
MGStraight(
	const MGPosition& endP,		///<End point.
	const MGPosition& startP,	///<Start point.
	const double endT,			///<End Parameter.
	const double startT = 0.0	///<Start Parameter.
);

/// 単位方向ベクトル、終点のパラメータ値、始点を指定して直線を生成する。
///Straight line from direction vector, end point parameter, and an origin.
MGStraight(
	const MGUnit_vector &,		///<Unit direction vector
	double,						///<Parameter value of end point
	const MGPosition& = mgORIGIN///<Origin
);
	
///Construct the infinite straight line that is a perpendicular bisect
///of the two point P1 and P2 and that is normal to the vector N.
///The line's direction is N*(P2-P1).
///N is the normal of the plane P1, P2, and the constructed line lie on.
MGStraight(
	const MGPosition& P1,	///<point 1.
	const MGPosition& P2,	///<point 2.
	const MGVector& N
);

///<Construct Straight Line copying original line. Able to change
///space dimension and ordering of axis.
MGStraight(
	size_t dim,					///<New space dimension.
	const MGStraight& linne2,	///<Original line.
	size_t start1=0, 			///<Destination order of new line.
	size_t start2=0 			///<Source order of original line.
);

//////////// Destructor ////////////////
~MGStraight();

//////////// Operator overload. 演算子の多重定義 ////////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGStraight& operator=(const MGGel& gel2);
MGStraight& operator=(const MGStraight& gel2);

/// 直線の平行移動を行いオブジェクトを生成する。
///Transformation
MGStraight operator+ (const MGVector& vec)const;
MGStraight operator- (const MGVector& vec)const;
MGStraight operator* (double scale)const;
MGStraight operator* (const MGMatrix&)const;
MGStraight operator* (const MGTransf&)const;

///Object transformation.
MGStraight& operator+=(const MGVector& v);
MGStraight& operator-=(const MGVector& v);
MGStraight& operator*=(double scale);
MGStraight& operator*=(const MGMatrix& mat);
MGStraight& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGStraight& sl2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGStraight& gel2)const;
bool operator<(const MGGel& gel2)const;

//////////// Member function. メンバ関数 ////////////

///Returns B-Rep Dimension.
size_t bdim()const{return 2;};

/// 指定線分を囲むボックスを返却する。
///Minimum box that includes the line limitted by interval l.
MGBox box_limitted(const MGInterval& l)const;

///Changing this object's space dimension.
MGStraight& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,		///<Parameter value for the start of the original. 
	double t2		///<Parameter value for the end of the original. 
);

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
MGStraight& coordinate_exchange(size_t i, size_t j);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGStraight* clone()const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs()const;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGStraight* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

/// 曲率を返却する。直線の場合は０。
///Return curvature, i.e. return 0.
double curvature( double ) const{return 0.;};

///Compute curvilinear integral of the 1st two coordinates.
///This integral can be used to compute area sorounded by the curve.
///(線積分）を求める。
double curvilinear_integral(double t1, double t2)const;

#ifdef __sgi
	double curvilinear_integral()const
	{return MGCurve::curvilinear_integral();}
#endif

/// 直線上のポイントにおける接線の傾きを返却する。
/// 直線のときは一定。
///Return direction vector of the line at a parameter.
MGUnit_vector direction(double)const{return m_direction;};

/// 線分の方向ベクトルを返却する。
///Return direction vector of the line.
const MGVector& direction()const{return m_direction;};

///Return direction vector length of the straight line.
double direction_len()const{return m_direction.len();};

/// 自身と与えられた点との距離を返す。
///Return the ditance between a position.
double distance( const MGPosition& = mgORIGIN )const;

/// 自身と与えられた直線との距離を返す。
///Return distance between two striaght lines.
double distance( const MGStraight& )const;

void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1///<line density to draw surface in wire mode.
)const;

void drawSE(
	double span_length,	///<Line segment span length.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve,
						///<Draw will be performed from t0 to t1.
)const;

///Return end parameter value.
const MGEReal& end()const{ return m_endparam;};

/// 終点のパラメータ値を返却する。
///Return end point parameter, valid only when type=MGSTRAIGHT_SEGMENT.
double end_param()const{return param_e();};

/// 終点の座標値を返却する。
///Return end point coordinate, valid only when type=MGSTRAIGHT_SEGMENT.
MGPosition end_point()const;

/// Evaluate n'th derivative data. nderiv=0 means
/// positional data evaluation.
MGVector eval(
	double t,			///< Parameter value.
	size_t nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
)const;

/// パラメータ値を与えて位置、一次微分値、二次微分値を求める。
///Evaluate positional data, 1st derivative, and 2nd derivative
///at parameter t. 
	void eval_all (
	double t,		///<Parameter of the straight.
	MGPosition&,	///<Positional data will be returned.
	MGVector&,		///<1st derivative will be returned.
	MGVector&		///<2nd derivative will be returned.
) const;

/// 直線上の与えられたパラメータ値における一次微分値を求める。
///Evaluate 1st derivative at parameter t.
MGVector eval_deriv(double t)const;

/// 与えられたパラメータ値に相当する直線上の点を返却する。
///Evaluate positional data at parameter t.
MGPosition eval_position(double t)const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line,
					///<If start is true extend on the start point.
);

///Test if thie straight is a finite line.
bool finite()const{return m_sparam.finite() && m_endparam.finite();};

/// Return This object's typeID
long identify_type()const;

///Test if infinite_above(infinite to larger parameter direction).
bool infinite_above()const{return m_endparam.plus_infinite();};

///Test if infinite_below(infinite to smaller parameter direction).
bool infinite_below()const{return m_sparam.minus_infinite();};

///Test if input parameter value is inside parameter range of the line.
bool in_range(double t)const;

///Provide divide number of curve span for function intersect.
size_t intersect_dnum() const{return 2;}

///Test if this cure is co-planar with the 2nd curve curve2.
///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if the input parameter t is the start point parameter or not.
bool is_startpoint_parameter(double t)const;

///Test if the input parameter t is the start point parameter or not.
bool is_endpoint_parameter(double t)const;

///Test if this cure is linear or not, that is, is straight or not.
///MGStraight expression will be out to straight if this is linear or not.
///Function's return value is true if linear.
bool is_linear(MGStraight& straight)const;

///Test if this cure is planar or not.
///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

/// Straight と Curve の交点を求める。
///Intersection of straight and a curve.
MGCCisect_list isect(const MGCurve&)const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGRLBRep& curve2)const;
MGCCisect_list isect(const MGEllipse& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;
MGCCisect_list isect(const MGBSumCurve& curve2)const;

///Intersection of Straight and Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;
MGCSisect_list isect(const MGFace& f)const;

///Access to i-th element of knot.
///i=0, 1 and returns start or end parameter value of the straight.
double knot(size_t i) const;

///Returns the knot vector of the curve.
const MGKnotVector& knot_vector() const;
MGKnotVector& knot_vector();

/// パラメータ値が昇順であたえられた場合正の値で、降順の場合
/// 負の値でパラメータ間の直線の代数的距離を返却する。
///Return curve length from parameter t1 to t2. If t1>t2, the result is
///minus value.
double length(double t1, double t2) const;

/// 自身の直線が有界の場合、その直線の距離を返却する。
/// 非有界のとき－１を返却する。
///When type=MGSTRAIGHT_SEGMENT, return the whole length of straight
///line segment.
///Else, return -1.
double length() const;

/// パラメータtで示される点から指定距離lenはなれた点のパラメータ値を返す。
///Return the parameter of the line away from point t by length len
///along the straight.
double length_param(double t, double len) const;

/// 自身の直線に指定されたｌｉｍｉｔを付与する。
///Compute sub straight line limitted by an interval.
MGStraight& limit(const MGInterval& );

///Compute sub straight line limitted by an box.
///This box's coordinates consist of world coordinates.
MGStraight& limit(const MGBox& box);

///Compute nearest point on the line to the origin.
MGPosition nearest_to_origin() const;

/// 直線の方向を反転する(方向ベクトルを逆向きにする。）。
/// 始終点があるときは入れ換える。
///Negate the direction of the curve.
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///一定オフセット関数
///オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
///法線ベクトルがヌルの場合、空間次元のエレメント番号が最も大きいエレメントが1の単位ベクトルを正とする。
///トレランスはline_zero()を使用する。戻り値は、オフセット曲線リストが返却される。
///costant offset curve. if the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the MGUnit_vector() direction.
///line_zero() is used. return value is number of offset curve.
MGPvector<MGCurve> offset(
	double ofs_value,								///<オフセット量
	const MGVector& norm_vector = mgNULL_VEC) const;///<法線ベクトル

/// 点が直線上にあるかを試験する。直線上にあれば，その点のパラメーター値を，
///　直線上になくても最近傍点のパラメーター値を返す。
///Test if input point is on the line or not. Even if the point is not
///on the line, return nearest point's parameter value.
///Function's return value is true if input point is on the line,
/// and  false if the point is not on the line.
bool on(
	const MGPosition&,	///<Input a point. 指定点
	double&				///<Parameter value of the curve will be returned,
						///<パラメータ値
) const;

/// 直線が平面上にあるか調べる。（平面上ならばtrue）
///Test if the straight is on a plane or not.
///Function's return value is true if on the line,
/// and  false if not on the line.
bool on(
	const MGPlane&		///< Plane
)const;

///Returns the order.
unsigned order() const{return 2;};

/// 直線上の与えられたポイントにおけるパラメータ値を返す。
///Return parameter of the straight of input point.
/// If input point is not on the curve, return the nearest point's
///parameter on the curve.
double param (
	const MGPosition &
) const;

/// Return ending parameter value.
double param_e() const{return m_endparam.value();};

///Obtain parameter space error.
double param_error() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance. For straight, the knots are start and end points.
double param_normalize(double t) const;

///Return parameter range of the curve(パラメータ範囲を返す)
MGInterval param_range() const;

/// Return starting parameter value.
double param_s() const{return m_sparam.value();};
	
///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGStraight* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

/// 与点から直線への垂線の足のパラメータ値を返却する。
///Return the foot of the straight line that is perpendicular from a point.
///Function's return value is parameter value of this straight line,
///may not be in_range.
double perp_param(
	const MGPosition &		///< 与点
) const;

/// 与えられたポイントから直線への垂線の足のパラメータ値を返却する。
///Return the foot of the straight line that is perpendicular from a point.
///Function's return value is if point is obtained(1) or not(0).
int perp_point(
	const MGPosition &,		///<point
	double& d1,				///<parameter value will be returned.
	const double* d2=NULL	///<guess parameter value, dummy, not used.
) const;
	
/// 与ポイントから直線へ下ろした垂線の足の，直線のパラメータ値を
/// すべて求める。
///Return all the foots of the  straight lines that is perpendicular
///to the line. Actually only one point for straight.
MGCParam_list perps(
	const MGPosition& point	///< 与ポイント
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
MGPosition_list perps(const MGRLBRep& crv2)const;
MGPosition_list perps(const MGEllipse& crv2)const;
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const;

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///Round the input parameter value t into the parameter range of the line.
double range(double) const;

/// 自身と与えられた直線の関係を調べる。
///Return two straight line's relationship.
///MGPSREL_UNKNOWN,     ///Unknown. 不明
///MGPSREL_TORSION,     ///Torsion. ねじれ
///MGPSREL_ISECT,       ///Intersection. 交差
///MGPSREL_PARALLEL,    ///Parallel. 平行
///MGPSREL_COIN,        ///Coincidence. 包含
///MGPSREL_VIRTUAL_ISECT///Virtually intersection
///                     (Intersection at extended line).
///                     指定直線の延長上の点での交差
MGPSRELATION relation(
	const MGStraight &,		///<2nd straight line.
	MGCCisect &		///<Intersection point will be returned if exist.
) const;
	
/// 自身と与えられた平面の関係を調べる。
///Return relationship with a plane.
MGPSRELATION relation(
	const MGPlane &,	///<Plane.
	MGCSisect &		///<Intersection point will be returned if exist.
) const;

const MGPosition& root_point() const{return m_root_point;};

///Return space dimension
size_t sdim() const;

/// 直線のタイプ，方向ベクトル，始点を指定して直線を生成する。
///Straight from straight line type, direction vector, and an origin.
///Construct a straight and replce this with it.
///This fuction does not convert input vec to a unit vector.
///If you like the conversion, use MGStraight() constructor.
MGStraight& set_straight(
	MGSTRAIGHT_TYPE type,	///<Type
		///<Straight line type(直線の種類)
		///<enum MGSTRAIGHT_TYPE {
		///<	MGSTRAIGHT_EMPTY 		//Empty. 空
		///<	,MGSTRAIGHT_SEGMENT		//Line segment. 線分
		///<	,MGSTRAIGHT_HALF_LIMIT	//Half unlimit. 半直線
		///<	,MGSTRAIGHT_UNLIMIT		//Unlimit line for both direction. 無限直線
		///<};
	const MGVector& vec,				///<Direction
	const MGPosition & = mgORIGIN			///<Origin
);

///Return straight line's direction.
const MGVector& sl_direction()const{return m_direction;}

///Return start point parameter value.
const MGEReal& start() const{ return m_sparam;};

/// 自身の直線の始点を返却する。
///Return start(root) point of the straight.
MGPosition start_point() const;

/// 直線のタイプを返却する。
///Return the straight line type.
///Straight line type(直線の種類)
///enum MGSTRAIGHT_TYPE {
///	MGSTRAIGHT_EMPTY 		//Empty. 空
///	,MGSTRAIGHT_SEGMENT		//Line segment. 線分
///	,MGSTRAIGHT_HALF_LIMIT	//Half unlimit. 半直線
///	,MGSTRAIGHT_UNLIMIT		//Unlimit line for both direction. 無限直線
///};
MGSTRAIGHT_TYPE straight_type() const;

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

/// 曲線タイプを返却する。
///Return curve type, i.e. MGCURVE_STRAIGHT.
MGCURVE_TYPE type() const{return MGCURVE_STRAIGHT;};

/// 自身の直線からｌｉｍｉｔを取り除く。
///Unlimit the parameter range, i.e. change to infinite striahgt line
///for both direction.
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
///(始点方向にlimitをはずす)
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
///(終点方向にlimitをはずす)
MGCurve& unlimit_start();

/// Output function.
///Output to ostream. メンバデータを標準出力に出力する。
std::ostream& out(std::ostream &) const;

int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

std::string whoami()const{return "Straight";};

protected:

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
)const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const;

///メンバデータを読み込む関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void WriteMembers(MGOfstream& buf) const;

private:

////////////Member data. メンバデータ//////////
/// Parameterization of the MGStraight is as below using parameter t:
/// point(t) = m_root_point + t*m_direction.

	MGPosition	    m_root_point;///<Root point. 基点
	MGVector		m_direction;///<Direction vector. 直線の方向ベクトル
	MGEReal			m_sparam;	///<Start point's parameter
	MGEReal		    m_endparam;	///<End point's parameter value. 
	mutable MGKnotVector* m_knotV;///<When knot_vector() is invoked, the knot vector is set.

///Compute the box of the whole of the curve.
///Returned is a newed object pointer.
MGBox* compute_box() const;

///isect2D returns parameter values of this(t) and l2(s)
/// of the intersection point of both 2D straight lines.
/// This and l2 are treated as infinite lines.
///Function's return value is:
/// true if two lines are parallel(or one of the directin is zero)
/// false if intersection was obtained.
bool isect2D(const MGStraight& l2, double& t,double& s) const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list isect_withC1LB(const MGLBRep& curve2)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///Compute parallel range of two straight lines.
///Two straight line this and lines must be parallel.
///Function's return value MGPosition_list list is:
/// list.entries() is 0(when no paralle part) or 2(when there is a parallel
/// part). When list.entries() is 2, let their entries be P1 and P2.
/// Then from P1(0) to P2(0) is the range of this straight line.
/// From P1(1) to P2(1) is the range of line2 straight line.
MGPosition_list relation_parallel(const MGStraight& line2) const;

///Function to avoid m_direction.len()=zero.
///The line's m_direction is set as a unit vector
///and m_endparam is set to zero.
void save_length_zero();

MGPSRELATION relation_parallel(
	const MGStraight& s,
	MGCCisect& ip ) const;

friend class MGCylinder;

};

/** @} */ // end of GEO group
#endif
