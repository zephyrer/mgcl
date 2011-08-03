/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSphere_HH_
#define _MGSphere_HH_
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/CSisect_list.h"
#include "mg/Surface.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"

// MGSphere.h
// Header for class MGSphere

class MGTransf;
class MGStraight;
class MGSBRep;
class MGIfstream;
class MGOfstream;
class MGIgesOfstream;

/** @addtogroup GEO
 *  @{
 */

///MGSphere is a Sphere in 3D space. Sphere f(u,v) is expressed
/// by two ellipses EL1(m_ellipseu) and EL2(m_ellipsev) as:
///f(u,v) = C+(M*cos(u)+N*sin(u))*cos(v)+B*sin(v), or
///f(u,v) = C+EL1(u)*cos(v)+B*sin(v),
///where EL1=M*cos(u)+N*sin(u), and EL2=C+N*cos(v)+B*sin(v).
///Here M is the major axis of EL1, N is the minor axis of EL1, N is
///also the major axis of EL2, and B=(unit vector of (M*N))*(N.len()),
///which is the minor axis of EL2. (M,N,B) make a orthonormal system.
///v is the angle with M axis in the (B,M) plane, and u is the angle with M in the (M,N) plane.
///v=0 parameter line makes the ellipse C+EL1, and u=pai/2 parameter line
///makes the ellipse EL2.
///MGSphereクラスは３次元空間における球を表すクラスである。
class MGCLASS MGSphere :public MGSurface{

public:

///translation
MGDECL friend MGSphere operator+ (const MGVector& v, const MGSphere& cyl);

///Scaling of the Sphere by a double.
MGDECL friend MGSphere operator* (double scale, const MGSphere& cyl);

///////////////Constructor コンストラクタ//////////////

///Void constructor 初期化なしで柱面を生成する。
MGSphere(void);

//Copy constructor.
//MGSphere(const MGSphere& cyl);

/// Construct a  whole sphere from the center and the radius.
MGSphere(
	const MGPosition& cntr,	///< Sphere center.
	double radius); 		///< Sphere radius.

/// Construct a  whole sphere from the center and the radius.
///Let MGUnit_vector N(B*M), M2(N*B). Then (M2,N,B) makes a orthonormal system,
///and this sphere is parameterized as:
///F(u,v)=cntr+radis*cos(v)(M*cos(u)+N*sin(u))+radis*sin(v)*B.
MGSphere(
	const MGPosition& cntr,	///< Sphere center.
	double radius,			///< Sphere radius.
	const MGUnit_vector& B,	///<axis
	const MGVector& M ///<reference direciotn that is approximately perpendiculat to B.
);

/// Construct a Sphere by changing this space dimension or
/// ordering the coordinates.
MGSphere(
	size_t dim,				///< New space dimension.
	const MGSphere& cyl,	///< Original Sphere.
	size_t start1=0, 		///< Destination order of new Surface.
	size_t start2=0 		///< Source order of original Surface.
);

///Whole sphere(u parameter range is from 0 to 2 pai) around minor axis of
///the input ellipse. The parameter range of the ellipse must be within the range
///from -pai/2 to pai/2, and the range becomes the v parameter range of the sphere.
///The input ellipse makes u=const(u=pai/2) v-paramarter line.
MGSphere(
	const MGEllipse& ellipse	///<ellispe  of the Sphere
);

///Sphere(u parameter range is urange) around minor axis of the input ellipse.
///The parameter range of the ellipse must be in the range from -pai/2 to pai/2,
///and the range becomes the v parameter range of the sphere.
///The input ellipse makes u=const(u=pai/2) v-paramarter line.
MGSphere(
	const MGEllipse& ellipse,	///<ellispe  of the Sphere
	MGInterval urange	///<parameter range along v in radian.
);

////////////Destructor////////////////
///~MGSphere();

////////////Operator overload 演算子の多重定義///////////////

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGSphere& operator=(const MGGel& gel2);
MGSphere& operator=(const MGSphere& gel2);

///Translation of the Sphere
MGSphere operator+ (const MGVector& ) const;

///Translation of the Sphere
MGSphere operator- (const MGVector& ) const;

///柱面のスケーリングを行い，柱面を作成する。
///Scaling of the Sphere by a double.
MGSphere operator* (double) const;

///Transformation of the Sphere by a matrix.
MGSphere operator* (const MGMatrix& ) const;

///Transformation of the Sphere by a MGTransf.
MGSphere operator* (const MGTransf& ) const;

///Object transformation.
MGSphere& operator+=(const MGVector& v);
MGSphere& operator-=(const MGVector& v);
MGSphere& operator*=(double scale);
MGSphere& operator*=(const MGMatrix& mat);
MGSphere& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGSphere& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGSphere& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGSphere& gel2)const{return !(gel2==(*this));};

///PD196=Spherical surface(parameterized).
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

//////////Debug function デバッグ関数////////////
/// Output function.
///Output to ostream メンバデータを標準出力に出力する。
std::ostream& out(std::ostream &) const;

//////////////Member function メンバ関数////////////

/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
///Box that includes limitted Sphere by box.
MGBox box_limitted(
	const MGBox& uvrange	///< Parameter Range of the surface.
) const;

///Changing this object's space dimension.
MGSphere& change_dimension(
	size_t sdim,		///< new space dimension
	size_t start1=0,	///< Destination order of new object.
	size_t start2=0 	///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
MGSphere& change_range(
	int is_u,		///<if true, (t1,t2) are u-value. if not, v.
	double t1,		///<Parameter value for the start of original. 
	double t2		///<Parameter value for the end of original. 
);

///Compute the closest point parameter value (u,v) of this surface
///from a point.
MGPosition closest(const MGPosition& point) const;

///Compute the closest point on a perimeter of the surface. The point is returned
///as the parameter value (u,v) of this surface.
MGPosition closest_on_perimeter(const MGPosition& point)const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGSphere* clone() const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGSphere* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Ask if this sphere has the degenerate point at v=min.
bool degenerate_at_v0()const;

///Ask if this sphere has the degenerate point at v=max.
bool degenerate_at_v1()const;

/// 与えられた点との距離を返却する。
///Return the distace between Sphere and the point.
double distance(const MGPosition& point) const;

///Return the constituent ellipses of the Sphere.
const MGEllipse& ellipseu() const{return m_ellipseu;};
const MGEllipse& ellipsev() const{return m_ellipsev;};
MGEllipse& ellipseu(){return m_ellipseu;};
MGEllipse& ellipsev(){return m_ellipsev;};

///Evaluate surface data.
MGVector eval(
	double u, double v		///< Parameter value of the surface.
	, size_t ndu=0			///< Order of derivative along u.
	, size_t ndv=0			///< Order of derivative along v.
) const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv	///< Parameter value of the surface.
	, size_t ndu=0			///< Order of derivative along u.
	, size_t ndv=0			///< Order of derivative along v.
)const{return eval(uv.ref(0),uv.ref(1),ndu,ndv);}

/// Exchange parameter u and v.
MGSurface& exchange_uv();

///Modify the original Surface by extrapolating the specified perimeter.
///The extrapolation is C2 continuous if the order >=4.
///The extrapolation is done so that extrapolating length is "length"
///at the position of the parameter value "param" of the perimeter.
MGSphere& extend(
	int perimeter,	///<perimeter number of the Surface.
					///< =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	///< parameter value of above perimeter.
	double length,	///<chord length to extend at the parameter param of the perimeter.
	double dk=0.  ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point.
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

/// Return This object's typeID
long identify_type() const;

bool in_range(double u, double v) const;

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
size_t intersect_dnum_u() const{ return m_ellipseu.intersect_dnum();};
size_t intersect_dnum_v() const{ return m_ellipsev.intersect_dnum();};

/// Surface と Curve の交点を求める。
///Compute intesection of Sphere and Curve.
MGCSisect_list isect(const MGCurve& curve)const;
MGCSisect_list isect(const MGStraight& line)const{ return isectSl(line);};
MGCSisect_list isect(const MGRLBRep& curve)const;
MGCSisect_list isect(const MGEllipse& curve)const;
MGCSisect_list isect(const MGLBRep& curve)const;
MGCSisect_list isect(const MGSurfCurve& curve)const;
MGCSisect_list isect(const MGBSumCurve& curve)const;

///Surface と Surface の交線を求める。
///Surface and Surface intersection.
MGSSisect_list isect(const MGSurface& srf2)const;
MGSSisect_list isect(const MGPlane& srf2)const;
MGSSisect_list isect(const MGSphere& srf2)const;
MGSSisect_list isect(const MGCylinder& srf2)const;
MGSSisect_list isect(const MGSBRep& srf2)const;
MGSSisect_list isect(const MGRSBRep& srf2)const;
MGSSisect_list isect(const MGBSumSurf& srf2)const;

///Return knot value of (infinite-minus, infinite-plus)
double knot_u(size_t i)const{return m_ellipseu.knot(i);};
double knot_v(size_t j)const{return m_ellipsev.knot(j);};

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const{return m_ellipseu.knot_vector();};
MGKnotVector& knot_vector_u(){return m_ellipseu.knot_vector();};

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const{return m_ellipsev.knot_vector();};
MGKnotVector& knot_vector_v(){return m_ellipsev.knot_vector();};

///Return the three axises of the sphere.
///MGSphere is a Sphere in 3D space. Sphere f(u,v) is expressed
/// by two ellipses EL1(m_ellipseu) and EL2(m_ellipsev) as:
///f(u,v) = C+(M*cos(u)+N*sin(u))*cos(v)+B*sin(v), or
///f(u,v) = C+EL1(u)*cos(v)+B*sin(v),
///where EL1=M*cos(u)+N*sin(u), and EL2=C+N*cos(v)+B*sin(v).
///Here M is the major axis of EL1, N is the minor axis of EL1, N is
///also the major axis of EL2, and B=(unit vector of (M*N))*(N.len()),
///which is the minor axis of EL2. (M,N,B) make a orthonormal system.
///v is the angle with M axis in the (B,M) plane, and u is the angle with M in the (M,N) plane.
///v=0 parameter line makes the ellipse C+EL1, and u=pai/2 parameter line
///makes the ellipse EL2.
const MGVector& M()const{return m_ellipseu.major_axis();};
const MGVector& N()const{return m_ellipseu.minor_axis();};
const MGVector& B()const{return m_ellipsev.minor_axis();};
const MGPosition& C()const{return m_ellipsev.m_center;};

/// 柱面を反転する。ノーマルを逆方向にする。
///Negate the normal of the Sphere.
void negate(int is_u);/// Negate along u-direction if is_u is ture,
					/// else along v-direction.

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const;

///C1連続曲面の一定オフセット関数
///オフセット方向は、ノーマル方向を正とする。トレランスはline_zero()を使用している。
///戻り値は、オフセット面のオートポインターが返却される。
///Surface offset. positive offset value is offset normal direction.
///the radius of curvature must be larger than offset value.
///line_zero() is used.
std::auto_ptr<MGSurface> offset_c1(
	double ofs_value,	///<オフセット量
	int& error	///<エラーコード 0:成功 -1:面におれがある
				///< -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
)const;

/// 指定点が面上にあるか調べる。（面上ならばtrue）
///Test if a point is on the Sphere. If on the Sphere, return true.
bool on(
	const MGPosition& point,///<A point to test 指定点
	MGPosition& puv			///<Parameter value of the Sphere will be returned.
)const;

///test if the surface normal is outgoing from the center or not.
///If the sphere normao is outgoing, retrun true.
bool outgoing()const;

///Obtain parameter space error.
double param_error() const;

/// Return ending parameter value.
double param_e_u()const{return m_ellipseu.param_e();};
double param_e_v()const{return m_ellipsev.param_e();};

/// パラメータ範囲を返す。
///Return parameter range of the Sphere(Infinite box).
MGBox param_range()const;

/// Return starting parameter value.
double param_s_u()const{return m_ellipseu.param_s();};
double param_s_v()const{return m_ellipsev.param_s();};

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
	MGCurve* parameter_curve(
	int is_u				///<Indicates x is u-value if is_u is true.
	, double x				///<Parameter value.
							///<The value is u or v according to is_u.
)const;

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
MGSphere* part(
	const MGBox& bx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
MGCurve* perimeter_curve(size_t i)const;

///Return how many perimeters this surface has.
size_t perimeter_num()const{return 4;};

/// 与えられた点にもっとも近い面上の点を返却する。パラメ
/// ータ値も返却する。
///Return the nearest point of the Sphere from P.
///Function's return value is always true.
int perp_point(
	const MGPosition& P,///< 与えられた点
	MGPosition& uv,		///<Parameter value of the Sphere will be output
	const MGPosition* uvguess=NULL	///< guess
)const;

///Return all(actually at most two) foots of perpendicular straight lines from P.
///When two points are output, the nearer point will be output first
///in MGPosition_list.
MGPosition_list perps(
	const MGPosition& P				///< Point of a space(指定点)
)const;

///Get the radius of the sphere
double radius()const{return m_ellipseu.radius();};

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///Round the input uv into parameter range of the Sphere, 
///return the same value as input.
MGPosition range(const MGPosition& uv) const;

///Return the space dimension.
size_t sdim() const;

///Return if this is sphere(i.e. the length of M, N, and B are all equal) or not.
bool sphere() const{return m_ellipseu.circle();};

///Get the center of the sphere
const MGPosition& sphere_center()const{return m_ellipsev.m_center;};

/// 曲面タイプを返却する。
///Return surface type of the Sphere.
MGSURFACE_TYPE type() const{return MGSURFACE_SPHERE;}

std::string whoami()const{return "Sphere";};

protected:

///Intersection of Surface and a straight line.
MGCSisect_list isectSl(
	const MGStraight& sl,
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
	) const;

///メンバデータを読み込む関数
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

private:

	MGEllipse  m_ellipseu;///<Ellipse of the Sphere(v=0 radian parameter line)
		///<The parameter range of m_ellipseu is within from 0 to 2*pai.
	MGEllipse  m_ellipsev;///<Ellipse of the Sphere(u=pai/2. radian parameter line)
		///<The parameter range of m_ellipsev is within from -pai/2 to pai/2.

///Compute this surface's box
void box_driver(MGBox& bx)const;

///Compute box of u=const parameter line
///The box will be or-ed to the input box bx.
///u must be the value in radian.
void box_uconst(MGBox& bx, double u)const;

///Compute box of v=const parameter line
///The box will be or-ed to the input box bx.
///v must be the value in radian.
void box_vconst(MGBox& bx, double v)const;

///Return minimum box that includes whole of the surface.
///Returned is a newed object pointer.
MGBox* compute_box() const;

///isect_area_length() returns initial area length for the intersection
///line.
size_t isect_area_length() const{return 10;};

///isect_dt computes incremental values of u and v direction for the intersection
///computation at parameter position (u,v).
void isect_dt(
	double u, double v, double& du, double& dv,
	double acuRatio=1.	///acuracy ratio.
) const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not,
							///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,				///<next u value will be output.
	double& v,				///<next v value will be output.
	size_t incr=1		///<Incremental valuse of B-coef's id.
) const;

///Return order of intersection line order of MGLBRep.
///The default is 4.
size_t isect_order() const{return 4;}

///Compute the intersection line of this and the plane pl.
MGSSisect_list intersectPl(const MGPlane& pl) const;

///オフセットするサンプルポイントの1パッチごとの分割数を求める
///全てのパッチ中の分割数で最大の値を返す
int offset_div_num() const{return 1;};

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& plane). This surf1D is used in isect for
///the argument of isect_startSphere, which will use surf1D to compute isect(plane).
///surf1D=0.(intersection with x=0. Sphere) is the intersection lines.
///This is not used for MGSphere, should not be invoked.
MGSBRep* surf1D(const MGPlane& plane)const{assert(false);return 0;};

///Compute (u,v), given a point P on the surface.
///Computed (u,v) may be outside real parameter range
///(sphere is regarded as a whole one).
void compute_uv(const MGPosition& P, double&u, double&v)const;

};

/** @} */ // end of GEO group
#endif
