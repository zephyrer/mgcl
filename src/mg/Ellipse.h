/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGEllipse_HH_
#define _MGEllipse_HH_

#include "mg/Curve.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"

// MGEllipse.h
// header for class MGEllipse
//

class MGInterval;
class MGMatrix;
class MGTransf;
class MGCCisect_list;
class MGPosition_list;
class MGLBRep;
class MGStraight;
class MGIfstream;
class MGOfstream;
class MGKnotVector;
class MGIgesDirectoryEntry;

/** @addtogroup GEO
 *  @{
 */

/// MGEllipse is a class to define an ellipse of 2D or 3D.
/// Ellipse is expressed as below using parameter t(radian):
/// Point(t) = m_center + m_m * cos(t) + m_n * sin(t),
///where t0 <= t <= t1, -2 pai<=t0<t1<=2 pai and t1-t0<=2 pai.
class MGCLASS MGEllipse : public MGCurve {

public:

///Object transformation.
MGDECL friend MGEllipse operator+ (const MGVector& v, const MGEllipse& e);
MGDECL friend MGEllipse operator* (double scale, const MGEllipse& el);

////////////Constructor(コンストラクタ)////////////

///Void constructor(初期化なしで楕円を生成する)
MGEllipse();

///Copy constructor.
///If to_radian is true,
///the parameter range will be changed to radian values.
MGEllipse(const MGEllipse& el, bool to_radian=false);

/// Change space dimension and ordering of coordinates.
MGEllipse (size_t sdim		///< new space dimension
	, const MGEllipse& ellip///< original ellipse data
	, size_t start1=0		///<start position coordinate of new ellipse
	, size_t start2=0		///<start position coordinate of original
);

///Ellipase from center, major axis vector, minor axis vector,
///and the parameter range. If minor axis vector is not normal to major,
///it is set to normal.
/// 中心、２つのベクトル、パラメータ範囲から楕円を生成する.
///******This is the fundamental constructor.******
MGEllipse (
	const MGPosition& center,	///<Center:中心点
	const MGVector& major_axis,	///<Major axis:第１ベクタ
	const MGVector& minor_axis,	///<Minor axis:第２ベクタ
	const MGInterval& prange	///<Parameter range in radian:角度の範囲(ラジアン)
);

///A whole circle form center, radius, and the normal.
/// 中心、半径、パラメータ範囲と法線ベクトルを指定して真円を作成する。
MGEllipse (
	const MGPosition& center,	///<Center:中心点
	double r,					///<Radius:半径
	const MGVector& normal=mgZ_UVEC	///<Normal:法線
);

///Arc from center, start point, angle d, and the normal.
///The angle can be negative value. When it is negative, the arc's
///direction is clockwise around v. When positive, anti-clockwise.
///The start point's angle is zero, and the the end point's angle is d.
/// 中心、始点、角度と法線を指定して円弧を作成する。
MGEllipse (
	const MGPosition& center	///<Center:中心点
	, const MGPosition& spoint	///<Start point:始点
	, double d= mgDBLPAI		///< Angle in radian.
	, const MGVector& v=mgZ_UVEC///<Normal:法線
);

///  矩形に内接する楕円を生成する
///  plane    楕円が乗る平面に平行な平面
///  corner1  楕円が内接する矩形のコーナーの座標
///  corner2  corner1 の対頂角の座標のヒント
/// 
///  もし corner2 が corner1 を通る plane に平行な平面に
///  乗っていない場合、corner2 は平行な平面に面直投影され、
///  それが計算に用いられる。
MGEllipse(
	const MGPlane&    plane, 
	const MGPosition& corner1, 
	MGPosition&       corner2
);

///An arc of radius r whose start point is start, and the end point is end.
///and that is normal to the vector N.
///The circle lies on the plane whose normal is N and that passes through
///start and end.
///radius r is able to have minus value, in which case the longer part of the
///arc out of the whole circle is constructed.
///The center of the circle C is:
///C=M+MGUnit_vector(sign(r)*N*(end-start))*sqrt(r*r-d*d),
///where M=(start+end)*.5, and d is the distance between start and M.
MGEllipse(
  double r,		///<radius
  const MGPosition& start,
  const MGPosition& end,
  const MGVector& N,
  bool whole_circle=false		///<true if the whole circle is to generate.
);

///An arc of radius r whose start point is start, and the end point is end.
///The circle lies on the plane that the three points start, end, and reference
///lie on.
///The center of the circle C is:
///C=M+MGUnit_vector(sign(r)*N*(end-start))*sqrt(r*r-d*d),
///where M=(start+end)*.5, d is the distance between start and M,
///and N=(reference-start)*(end-start).
///That is, if r>0., the position reference indicates on which side
///the C lies against the vector (end-start), C and reference lie
///in the same half plane that the vector (end-start) divides.
///If r<0., C and reference lie in the opposite side.
///When r>0., smaller part arc of the circle is selected,
///and when r<0., larger part of the circle is selected.
MGEllipse(
  double r,		///<radius
  const MGPosition& start,
  const MGPosition& end,
  const MGPosition& reference,
  bool whole_circle=false		///<true if the whole circle is to generate.
 );

///Arc of the circle that osculates to two straight lines that passes
///point P and whose directional vectors are V1 and V2.
///The starting point is the contact point with the 1st line, and end point
///is the contact point with 2nd line. If radius d is positive, the arc is
///the one of point P side, and if negative, the other side.
/// １点Pと２つのベクトル(V1, V2)からなる２直線に接する指定半径の
/// 円弧を作成する。
/// 第一ベクトルとの接点を始点、他方を終点とし、
/// d>0 のときこれらで分割された円弧の指定点P側を
/// d<0 のときPと反対側を作成する。
MGEllipse (
	const MGPosition &P		///<Start point of two straight lines:直線の始点
	, const MGVector & V1	///<Directional vector 1:直線の方向ベクトル１
	, const MGVector & V2	///<Directional vector 2:直線の方向ベクトル２
	, double d				///<Radius:半径
);

///Arc from start point, through point, and end point.
/// 始点、通過点、終点を指定して円弧を作成する。
MGEllipse (
	const MGPosition& start,	///<Start point:始点
	const MGPosition& through,	///<Through point:通過点
	const MGPosition& end,		///<End point:終点
	bool whole_circle=false		///<true if the whole circle is to generate.
);

///Arc from center, start point, end point.
MGEllipse(
	const MGPosition& center,	///<center of the circle.
	const MGPosition& start,	///<Start point:始点
	const MGPosition& end,		///<End point:終点
	const MGVector& N,			///<Normal of the plane the circle lies on.
	bool whole_circle=false		///<true if the whole circle is to generate.
);

///Construct the arc whose start point is start,
///whose end point is end, and whose tangent at the start point is dir_s.
MGEllipse(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& dir_s,
	bool whole_circle=false		///<true if the whole circle is to generate.
);

///目的：
///		基本線２本と半径からコーナーＲを作成する
///		初期パラメータはコーナーＲを作成したい側に設定し、基本線は同一平面上にあること
///戻り値：
///		0:			正常終了
///		-1:			基本線の曲率半径より半径Rが大きい
///		-2:			基本線同士の交点が求まらない
///Start point of the generated ellipse is either t1 side or t2,
///which depends on normal direction. The ellipse direction is defined from normal,
///If the arc length from t1 to t2(arond normal) is smaller than from t2 to t1(around normal),
///the start point is on t1 side. If the arc length from t1 to t2 is longer than from t2 to t1,
///the start point is on t2 side.
MGEllipse (
	const MGCurve& crv1,			///<I:基本線1
	const MGCurve& crv2,			///<I:基本線2
	const MGUnit_vector& normal,	///<I:基本線のノーマルベクトル
	double dRadius,					///<I:コーナーＲの半径
	double& t1,						///<I:基本線１の初期パラメータ(コーナーＲを作成する側を指定)
									///<O:コーナーＲと接する基本線１のパラメータ値
	double& t2,						///<I:基本線２の初期パラメータ(コーナーＲを作成する側を指定)
									///<O:コーナーＲと接する基本線２のパラメータ値
	int& rc							///<O:リターンコード
);

///目的：
///		基本線２本とＲ止まり点からコーナーＲを作成する
///		初期パラメータはコーナーＲを作成したい側に設定し、基本線は同一平面上にあること
///戻り値：
///		0:			正常終了
///		-1:			入力値が不正
///		-2:			crv2上の接点が求まらない(交点なし)
///		-3:			crv2上の接点が求まらない(収束しない)
///Start point of the generated ellipse is either t1 side or t2,
///which depends on normal direction. The ellipse direction is defined from normal,
///If the arc length from t1 to t2(arond normal) is smaller than from t2 to t1(around normal),
///the start point is on t1 side. If the arc length from t1 to t2 is longer than from t2 to t1,
///the start point is on t2 side.
MGEllipse (
	const MGCurve&			crv1,	///<I:基本線1
	const MGCurve&			crv2,	///<I:基本線2
	const MGUnit_vector&	normal,	///<I:基本線のノーマルベクトル
	double					t1,		///<I:基本線1上のＲ止まり点のパラメータ値
	double&					t2,		///<I:基本線2の初期パラメータ(コーナーＲを作成する側を指定)
									///<O:コーナーＲと接する基本線２のパラメータ値
	int&					rc	///<O:リターンコード
);

///目的：
///		基本線３本からコーナーＲを作成する
///		crv1, crv2, crv3がそれぞれ始点、通過点、終点となるようにする。
///		初期パラメータはコーナーＲが接する近辺に設定し、基本線は同一平面上にあること
///戻り値：
///		0:			正常終了
///		-1:			連立方程式が解けなかった
///		-2:			収束しなかった(パラメータ範囲を超えてしまった)
///		-3:			収束しなかった(無限ループにおちいった)
MGEllipse (
	const MGCurve&			crv1,	///<I:基本線1(始点)
	const MGCurve&			crv2,	///<I:基本線2(通過点)
	const MGCurve&			crv3,	///<I:基本線3(終点)
	const MGUnit_vector&	normal,	///<I:基本線のノーマルベクトル
	double&					t1,		///<I:基本線1の初期パラメータ(コーナーＲの始点近辺)
									///<O:コーナーＲの始点における基本線１のパラメータ値
	double&					t2,		///<I:基本線2の初期パラメータ(コーナーＲの通過点近辺)
									///<O:コーナーＲの通過点における基本線２のパラメータ値
	double&					t3,		///<I:基本線3の初期パラメータ(コーナーＲの終点近辺)
									///<O:コーナーＲの終点における基本線２のパラメータ値
	int&					rc	///<O:リターンコード
);

///目的：
///		基本線と円弧端点と半径からコーナーＲを作成する
///		初期パラメータはコーナーＲを作成したい側に設定し、基本線と円弧端点は同一平面上にあること
///戻り値：
///		0:			正常終了
///		-1:			入力値が不正
///		-2:			半径値が小さすぎる
///		-3:			半径値が大きすぎる
///		-4:			基本線をオフセットしたカーブと円弧端点から作成した円の交点が求まらなかった
MGEllipse (
	const MGCurve&			crv,	///<I:基本線
	const MGPosition&		pos,	///<I:円弧端点
	const MGUnit_vector&	normal,	///<I:基本線のノーマルベクトル
	double					dRadius,///<I:コーナーＲの半径
	double&					dParam,	///<I:基本線の初期パラメータ(コーナーＲを作成する側を指定)
									///<O:コーナーＲと接する基本線１のパラメータ値
	int&					rc	///<O:リターンコード
);

///目的：
///		基本線と基本線上Ｒ止まりと円弧端点からコーナーＲを作成する
///		初期パラメータはコーナーＲを作成したい側に設定し、基本線と円弧端点は同一平面上にあること
MGEllipse (
	const MGCurve& crv,		///<I:基本線
	const MGPosition& P2,	///<I:円弧端点
	const MGUnit_vector& normal,///<I:基本線のノーマルベクトル
	double t				///<I:R止まり点のパラメータ
);

///Conversion constructor from MGIgesDirectoryEntry object.
MGEllipse(
	const MGIgesDirectoryEntry& de	///<type number must be 100(Circular Arc).
);

//////////Destructor//////////
~MGEllipse();

////////////Operator overload(演算子多重定義)////////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGEllipse& operator=(const MGGel& gel2);
MGEllipse& operator=(const MGEllipse& el2);

///Transformation object construction
MGEllipse operator+ (const MGVector&) const;
MGEllipse operator- (const MGVector&) const;
MGEllipse operator* (double scale) const;
MGEllipse operator* (const MGMatrix&) const;
MGEllipse operator* (const MGTransf&) const;

///Object transformation.
MGEllipse& operator+=(const MGVector& v);
MGEllipse& operator-=(const MGVector& v);
MGEllipse& operator*=(double scale);
MGEllipse& operator*=(const MGMatrix& mat);
MGEllipse& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGEllipse& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGEllipse& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member function:メンバ関数////////////

///Returns B-Rep Dimension.
size_t bdim() const {return 2;};

///Return minimum box that includes the curve of parameter interval.
/// パラメータ値による楕円上のポイントによって境界が与えられた
/// 部分または、楕円の周りを含む最小のボックスを返す。
MGBox box_limitted( const MGInterval&) const;

///Return the center of the ellipse:楕円の中心座標を返却する。
MGPosition center()const{return m_center;};

///Changing this object's space dimension.
MGEllipse& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Chane the parameter range of this ellipse to radian.
void change_param_to_radian();

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
);	

///Test if true circle. If circle, return true.
///円か調べる。円のときtrue。
bool circle() const {return m_circle!=0;};

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGEllipse* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGEllipse* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
MGEllipse& coordinate_exchange(size_t i, size_t j);

///  焦点と通過点を与えて楕円を生成する
///  F0  焦点
///  F1  焦点
///  P   楕円上の点
/// 
///  F0 と F1 は相異なる焦点でなければ失敗する。
///  P が2焦点を結ぶ直線上にあるとき失敗する。
bool create_ellipse(
	const MGPosition& F0,
	const MGPosition& F1,
	const MGPosition& P
	);

///Compute curvilinear integral of the 1st two coordinates.
///This integral can be used to compute closed area sorounded by the curve.
///(線積分）を求める。
double curvilinear_integral(double t1, double t2) const; ///From t1 to t2;

#ifdef __sgi
	double curvilinear_integral() const
	{	return curvilinear_integral(param_s(), param_e());};
#endif

void drawSE(
	double span_length,	///<Line segment span length.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve,
						///<Draw will be performed from t0 to t1.
)const;

///Return ellipse type:楕円のタイプを返却する。
MGELLIPSE_TYPE ellipse_type()const;

/// Evaluate n'th derivative data. nderiv=0 means
/// positional data evaluation.
MGVector eval(
	double,				///< Parameter value.
	size_t nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
) const;

/// Evaluatefirst derivative data.
/// 楕円上の与えられたパラメータ値における一次微分値を返却する。
MGVector eval_deriv ( double ) const;

///Compute position, 1st and 2nd derivatives.
/// パラメータを与え、位置、一次微分値、二次微分値を求める。
void eval_all (
	double,
	MGPosition &,	///<位置
	MGVector &,		///<一次微分値
	MGVector &		///<二次微分値
) const;

///Evaluate ellipse data.
///Input parameter t must be in radian.
MGVector eval_in_radian(
		double t,		///<Parameter value in radian.
		size_t nderiv=0	///< Order of Derivative.
) const;

///Evaluate ellipse data.
///Input parameter degree must be in degree.
MGVector eval_in_degree(
		double degree,	///<Parameter value in degree.
		size_t nderiv=0	///< Order of Derivative.
) const;

///Compute positional data.
/// 与えられたパラメータ値に相当する楕円上の点を返却する。
///Input parameter t must be in radian.
MGPosition eval_position(double t)const{
	t=range(t);
	return eval_in_radian2(gp_to_radian(t));
};

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line,
					///<If start is true extend on the start point.
);

///Compute the radian parameter of the general parameter t;
double gp_to_radian(double t)const{
	if(param_range_is_in_radian())
		return t;
	return m_prange[0]+(t-m_gprange[0])*m_gprange[2];
}

/// Return This object's typeID
long identify_type() const;

///Test if angle t is in the parameter range.
///t is increased or decreased by 2*PAI if necessary in testing.
///Radian version. Input t must be of radian.
bool in_RelativeRange_of_radian(double t) const;

///Provide divide number of curve span for function intersect.
size_t intersect_dnum() const;

///Intersection of ellipse and curve.
/// Ellipse と Curve の交点を求める。
MGCCisect_list isect(const MGCurve&) const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGEllipse& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;
MGCCisect_list isect(const MGBSumCurve& curve2)const;

///Intersection of Ellipse and Surface

///Intersection of MGEllipse and Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;

///Test if this cure is linear or not, that is, is straight or not.
///MGStraight expression will be out to straight if this is linear or not.
///Function's return value is true if linear.
bool is_linear(MGStraight& straight)const;

///Test if this cure is planar or not.
///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///Test if this ellipse is whole ellipse or not.
///Function's return value is true if this is whole.
bool is_whole_ellipse()const;

///Access to i-th element of knot.
///i=0, 1 and returns start or end parameter value of the ellipse.
double knot(size_t i) const;

///Returns the knot vector of the curve.
const MGKnotVector& knot_vector() const;
MGKnotVector& knot_vector();

///Update this by limiting the parameter range of the curve.
/// 自身の楕円に指定された範囲のlimitを付加する。
MGEllipse& limit(const MGInterval& );

///Negate the curve direction(曲線の方向を反転する)
/// 楕円の方向を反転する。方向ベクトルを逆にする。範囲があるとき
/// は始終点を入れ換える。
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///Unlimit parameter range of the curve(limitをはずす)
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
///(終点方向にlimitをはずす)
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
///(始点方向にlimitをはずす)
MGCurve& unlimit_start();

///Cmpute curve length of the interval.
///If t1 is greater than t2, return negative value.
/// 与えられたパラメータ間の曲線に沿った距離を返却する。
/// パラメータ値が昇順で与えられたとき正値、降順のとき負値を
/// 返す。
double length(double t1, double t2) const;

///Compute whole curve length.
/// 楕円の全体の長さを返却する。
double length() const;

///Inverse function of length. Compute the point that is away from
///the point t by length len.
/// パラメータtで示される点から指定距離lenはなれた点のパラメータ
/// を返す。
double length_param(double t, double len) const;

///Return major axis:長軸を返却する。
const MGVector& major_axis() const {return m_m;};

///Return major axis length:長軸の長さを返却する。
double major_len() const {return m_m.len();};

///Return minor axis:短軸を返却する
const MGVector& minor_axis() const {return m_n;};

///Return major axis length:長軸の長さを返却する。
double minor_len() const {return m_n.len();};

///Return normal:楕円のある平面の法線ベクトルを返却する。
const MGUnit_vector& normal() const {return m_normal;};

///Test if given point is on the curve or not. If yes, return parameter
///value of the curve. Even if not, return nearest point's parameter.
/// Function's return value is true(>0) if the point is on the curve,
/// and false(0) if the point is not on the curve.
/// 点が楕円上にあるか調べる。楕円上であれば，そのパラメータ値を，
/// そうでなくても最近傍点のパラメータ値を返す。
bool on(
	const MGPosition&,	///<A point:与ポイント
	double& d1			///<Parameter value will be returned:パラメータ
) const;

///Returns the order.
unsigned order() const{return 2;};

/// Return ending parameter value.
double param_e() const{
	if(m_gprange) return m_gprange[1];
	else return m_prange[1];
};

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance. For ellipse, the knots are start and end points.
double param_normalize(double t) const;

/// Return starting parameter value.
double param_s() const{
	if(m_gprange) return m_gprange[0];
	else return m_prange[0];
};

///Test if this ellipase 's parameter range is expressed in radian.
bool param_range_is_in_radian()const{return m_gprange==0;};

///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by the calling program, or memory leaked.
MGEllipse* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Compute a foot point of the perpendicular line from point p to
///the curve. If more than one points are found, return nearest one.
/// Function's return value is if point is obtained(>0) or not(0)
/// 与ポイントから楕円への垂線の足、そのポイントでの楕円の
///パラメータ値を返却する。
int perp_point (
	const MGPosition &p,///<A point:与ポイント
	double& d,			///<Parameter value will be returned:パラメータ
	const double*g=NULL	///< guess parameter value of the foot parameter
) const;			///< value of the ellipse.
	
///Compute all foot points of the perpendicular line from point p to
///the curve.
/// 与ポイントから楕円へ下ろした垂線の足の，楕円のパラメータ値を
/// すべて求める。
MGCParam_list perps(
	const MGPosition &	///<A point:与ポイント
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;

///Compute general parameter t from the radian parameter angle
double radian_to_gp(double angle)const{
	if(!m_gprange) return angle;
	return m_gprange[0]+(angle-m_prange[0])/m_gprange[2];
}

///Return the radius of the circle. This is valid only when circle() is true.
double radius()const{return m_r;};

///Return space dimension
size_t sdim() const ;

///Replace this ellipse with the arc whose start point is start,
///whose end point is end, and whose tangent at the start point is dir_s.
void set_arc(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& dir_s
);

///Return sweep surface from crv.
///Output is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
) const;	
///Return curve type(曲線のタイプを返す)
MGCURVE_TYPE type() const {return MGCURVE_ELLIPSE;};

///Debug function:デバッグ関数

///Output function.
///Output to stream file:メンバデータを標準出力に出力する。
std::ostream& out(std::ostream& ostrm) const;

///Output to IGES stream file.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

protected:

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const;

///メンバデータを読み出す関数
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

std::string whoami()const{return "Ellipse";};

private:

///Member data(メンバデータ)
	MGPosition	m_center;	///<Center(中心点)
	MGUnit_vector m_normal;	///<Normal of plane the ellipse lies on
							///<(楕円のある平面の法線ベクトル)*
	MGVector	m_m;		///<Major axis(長軸)
	MGVector	m_n;		///<Minor axis(短軸)
	double		m_r;	///< sqrt((a*a+b*b)/2) where a=m_m.len(), b=m_n.len()
	double	m_prange[2];///<Parameter range in radian from start to end point.
			///<-mgDBLPAI<=m_prange[.]<=mgDBLPAI and m_prange[1]-m_prange[0]<=mgDBLPAI.
	int	m_circle;	///<True when ellipse is a circle(m_m.len()==m_n.len()).
	double* m_gprange;	///<When m_gprange!=0, includes the general parameter range:
			///<m_gprange[0]=start parameter, m_gprange[1]=end, and
			///<m_gprange[2]=(m_prange[1]-m_prange[0])/(m_gprange[1]-m_gprange[0]).
	mutable MGKnotVector* m_knotV;
			///<When knot_vector() is invoked, the knot vector will be set.

///Compute if m_normal is the same as x, y, or z-axis or not.
///Return value is: 0,1,2: m_normal is the same as x, y, or z-axis, each.
///                -1, -2, -3: m_normal is not the same but nearest to
///                            x, y, or z-axis each.
int axis() const;

///Compute whole box of the curve. Retured is a pointer of a newed MGBox.
MGBox* compute_box() const;			///Whole of the curve.

///Compute positional data.
/// 与えられたパラメータ値に相当する楕円上の点を返却する。
///Input parameter t must be in radian.
///eval_in_radian2 does not round the parameter t into this parameter range.
MGVector eval_in_radian2(
		double t,		///<Parameter value in radian.
		size_t nderiv=0	///< Order of Derivative.
) const;

///copy all the ellipse specific data into this from ellipse2.
void copy_ellipse_data(const MGEllipse& ellipse2);

///Test if input angle is in the radian parameter range.
///When m_prange[0]<=angle and angle<=m_prangle[1], return true
///(tolerance included).
bool in_radian_range(double angle) const;

///Compute intesection points of  2D ellipse whose center is origin and
///a straight line of 2D.
///Return value of isect2D is number of intersetion points: 0, 1, or 2.
int isect2d(
	const MGPosition& sp,	///<Start point of straight line.
	const MGVector& dir,	///<Direction unit vector of the straight line
	double t[2],	///<Parameter of intersection point of the straight.
	double angle[2],	///< angles of ellipse in radian.
	int& tangen			///< Return if isect is tangent point(1) or not(0).
) const;

///Compute angles of elllipse points that are perpendicular to a point.
///This function is 2-D version of perp_point.
///Function's return value is number of points obtained.
size_t perp2d(double p,		///<x-coordinate of the point.
	double q,			///<y-coordinate of the point.
	double theta[4] ///<angles of the ellipse.
) const;

///Compute how distant param p1 is from param p2 in radian.
///p1 and p2 must be the values in radian.
double param_length(double p1, double p2) const;

///Round angle into curve's parameter range.
///angle is increased or decreased by 2*PAI if necessary in rounding.
///Input p and function's return value are both in radian value.
double RelativeRange_in_radian(double angle) const;

///Normalize parameter range [a0, a1] in radian of the ellipse,
///and set the parameter range in m_prange.
void set_param(
	double a0, 	///<Input paramter range in radian.
	double a1
);

///Compute m_normal, m_r, and m_circle from m_m and m_n.
void set_normal_r_c();

friend class MGCylinder;
friend class MGSphere;

};

/** @} */ // end of GEO group
#endif
