/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPlane_HH_
#define _MGPlane_HH_
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Surface.h"

// MGPlane.h
// Header for class MGPlane
class MGBPointSeq;
class MGCCisect;
class MGCSisect;
class MGTransf;
class MGStraight;
class MGSBRep;
class MGIfstream;
class MGOfstream;
class MGFace;

/** @addtogroup GEO
 *  @{
 */

///MGPlane is infinite plane in 3D space.
///Using orthonormal two vector m_uderiv and m_vderiv in 3D space,
///plane function f(u,v)= m_root_point + u*m_uderiv + v*m_vderiv,
///where u and v are two parameter of surface representation. 
/// MGPlaneクラスは３次元空間における平面を表すクラスである。
/// MGPlaneクラスでは以下のようなパラメータ表現を使用します。
/// Point(u,v) = m_root_point + u * m_uderiv + v * m_vderiv
class MGCLASS MGPlane :public MGSurface{

public:

MGDECL friend MGPlane operator+ (const MGVector& v, const MGPlane& pl);
MGDECL friend MGPlane operator* (double scale, const MGPlane& pl);

///////////////Constructor コンストラクタ//////////////

///Void constructor 初期化なしで平面を生成する。
MGPlane(void);

///Copy constructor.
MGPlane(const MGPlane& pl);

/// Construct a plane by changing this space dimension or ordering the coordinates.
MGPlane(
	size_t dim,				///< New space dimension.
	const MGPlane& plane,	///< Original Plane.
	size_t start1=0, 		///< Destination order of new Surface.
	size_t start2=0		///< Source order of original Surface.
); 

/// Construct a plane from the coefficients of the plane equation.
///     a*x+b*y+c*z=d.
/// Coefficients a,b,c,d are provided by double array g[4].
MGPlane(
	const double g[4],		///<coefficients g[0]=a, b=g[1], c=g[2], d=g[3].
	const double* root_point=0///<When root_point!=0, root_point[.] are (x,y,z) values of the root point.
);

///Plane from normal of the plane and the distance from the origin(0,0,0).
/// 平面のノーマルと平面の原点からの距離を指定して面を生成する.
MGPlane(
	const MGUnit_vector& normal,	///<Normal of the plane.
	double d			///<distance from origin of the plane.
						///<When normal=(a,b,c,....), d=a*x+b*y+c*z+.... .
);

///Plane from a point on the plane and the normal.
/// 点と平面のノーマルを指定して面を生成する。
MGPlane(
	const MGUnit_vector& normal,
	const MGPosition& p
);

///Plane from a point and a straight line on the plane.
/// 直線と直線に乗らない点を指定して面を生成する。
MGPlane(
	const MGStraight& st,
	const MGPosition& point
);

///Plane from a point on the plane, u and v direction vector of the plane.
/// 点、ｕ方向、ｖ方向を指定して面を生成する。
///***** This is the fundamental constructor.*****
MGPlane(
	const MGVector& uderiv,
	const MGVector& vderiv,
	const MGPosition &origin
);

/// Construct a plane by interpolating two planes.
/// If two planes intersect, interpolation is rotation around the
/// intersection line.
///If two planes are parallel, interpolation is parallel move.
MGPlane(
	const MGPlane& plane1,
	const MGPlane& plane2,
	double t				///<Input ratio.
				///< When t=0, the plane will be plane1,
				///< When t=1, the plane will be plane2.
);

///Construct a plane from three points on the plane.
MGPlane(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
);

////////////Destructor////////////////
~MGPlane();

////////////Operator overload 演算子の多重定義///////////////

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGPlane& operator=(const MGPlane& gel2);
MGPlane& operator=(const MGGel& gel2);

///Transformation object construction
MGPlane operator+ (const MGVector& v) const;
MGPlane operator- (const MGVector& v) const;
MGPlane operator* (double scale) const;
MGPlane operator* (const MGMatrix& mat) const;
MGPlane operator* (const MGTransf& tr) const;

///Object transformation.
MGPlane& operator+=(const MGVector& v);
MGPlane& operator-=(const MGVector& v);
MGPlane& operator*=(double scale);
MGPlane& operator*=(const MGMatrix& mat);
MGPlane& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGPlane& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGPlane& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGPlane& gel2)const{return !(gel2==(*this));};

///Output to IGES stream file(=PD190).
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

//////////Debug function デバッグ関数////////////

///Output function.
///Output to ostream メンバデータを標準出力に出力する。
std::ostream& out(std::ostream &) const;

//////////////Member function メンバ関数////////////

///Gets parameters(a,b,c,d) of the plane expression a*x+b*y+c*z=d.
void abcd(double g[4]) const;
///g[.]=(a,b,c,d)

///Box that includes limitted plane by box.
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox box_limitted(
	const MGBox& uvrange	///< Parameter Range of the surface.
) const;

///Obtain ceter coordinate of the geometry.
MGPosition center() const{return m_root_point;};

///Changing this object's space dimension.
MGPlane& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
/////********** MGPlane does not accept change_range, does nothing.****///
MGPlane& change_range(
	int is_u,				///<if true, (t1,t2) are u-value. if not, v.
	double t1,				///<Parameter value for the start of original. 
	double t2				///<Parameter value for the end of original. 
){return *this;};

///Compute the closest point parameter value (u,v) of this surface
///from a point.
MGPosition closest(const MGPosition& point) const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGPlane* clone() const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGPlane* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

void display_arrows()const;

///Return the distace between plane and the origin(0,0,0).
/// 原点との距離を返却する。
double distance() const{return m_d;};

///Return the distace between plane and the point.
/// 与えられた点との距離を返却する。
double distance(const MGPosition& point) const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Make a display list without color of this gel.
///Return is the display list name.
size_t make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const{return MGObject::make_display_list_to_hilight(span_length,line_density);}

///Shade the object in world coordinates.
void shade(
	double span_length	///Line segment span length.
)const{drawWire(span_length);};

///Evaluate surface data.
MGVector eval(
	double u, double v		///< Parameter value of the surface.
	, size_t ndu=0			///< Order of derivative along u.
	, size_t ndv=0			///< Order of derivative along v.
)const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv	///< Parameter value of the surface.
	, size_t ndu=0			///< Order of derivative along u.
	, size_t ndv=0			///< Order of derivative along v.
)const{return eval(uv.ref(0),uv.ref(1),ndu,ndv);}

///Evaluate right continuous surface data.
///Evaluate all positional data, 1st and 2nd derivatives.
void eval_all(
	double u, double v,	///< Parameter value of the surface.
	MGPosition& f,		///< Positional data.
	MGVector&   fu,		///< df(u,v)/du
	MGVector&   fv,		///< df/dv
	MGVector&   fuv,	///< d**2f/(du*dv)
	MGVector&   fuu,	///< d**2f/(du**2)
	MGVector&   fvv		///< d**2f/(dv**2)
)const;

/// Exchange parameter u and v.
MGSurface& exchange_uv();

/// Return This object's typeID
long identify_type() const;

bool in_range(double u, double v) const{ return 1;};
bool in_range(const MGPosition& uv) const{ return 1;};

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
size_t intersect_dnum_u() const{ return 2;};
size_t intersect_dnum_v() const{ return 2;};

/// Surface と Curve の交点を求める。
///Compute intesection of Plane and Curve.
MGCSisect_list isect(const MGCurve& curve)const;
MGCSisect_list isect(const MGStraight& curve)const;
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

///isect_startHPL compute one intersection line of two surfaces, this and sf2,
/// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2).
/// isect_startHPL halts the computation when intersection
/// reached to a boundary of this or sf2, or reached to one of the points
/// in uvuv_list.
///The function's return value is:
/// =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///    =1: End point is a point on a perimeter of one of the surfaces.
///    =3: End point is one of boundary points in uvuv_list.
int isect_startHPL(
	const MGPosition& uvuv_startIn, ///<Starting point of the intersection line.
	MGPosition_list& uvuv_list,	///<isect_startHPL will halt when ip reached one of 
		///<the point in uvuv_list. isect_startHPL does not change uvuv_list(actually
		///<uvuv_list is const). uvuv's space dimension is at least 4,
		///<and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGSurface& sf2,	///<2nd surface.
	MGSSisect& ssi,			///<Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id
		///<When the end point of ip was one of the points of uvuv_list, that is, 
		///<when the function's return value was 3, uvuv_list's iterator
		///<of the point will be returned,
		///<When the end point was not a point of uvuv_list, end() of uvuv_list
		///<will be returned.
) const;

///Return knot value of (infinite-minus, infinite-plus)
double knot_u(size_t)const ;
double knot_v(size_t)const ;

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const;
MGKnotVector& knot_vector_u();

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const;
MGKnotVector& knot_vector_v();

///Negate the normal of the plane,平面を反転する。ノーマルを逆方向にする.
void negate(
		int is_u///< Negate along u-direction if is_u is ture, else along v-direction.
);

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const;

///Return the normal of the plane, 平面の法線を返却する.
MGVector normal(double u, double v) const{return m_normal;}
MGVector normal(const MGPosition& uv) const{return m_normal;}
const MGUnit_vector& normal() const{return m_normal;}

///Surface offset. positive offset value is offset normal direction.
///the radius of curvature is larger than offset value.line_zero() is used.
///C1連続曲面の一定オフセット関数
///オフセット方向は、ノーマル方向を正とする。トレランスはline_zero()を使用している。
///戻り値は、オフセット面のオートポインターが返却される。
std::auto_ptr<MGSurface> offset_c1(
	double ofs_value,	///<オフセット量.
	int& error			///<エラーコード 0:成功 -1:面におれがある
						///< -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー.
) const;

///Test if a point is on the plane. If on the plane, return true.
/// 指定点が面上にあるか調べる。（面上ならばtrue）.
bool on(const MGPosition& point) const;
bool on(
	const MGPosition& point,	///<A point to test 指定点.
	MGPosition& puv	///<Parameter value of point on the plane will be returned either
					///<point is on the plane or not.
) const;

///Test if a straight line is on the plane. Return true if on.
/// 直線が平面上にあるか調べる。（平面上ならばtrue）
bool on(const MGStraight&) const;

///Test if input (u,v) is parameter value on a perimeter of the surface.
///If u or v is on a perimeter, (u,v) will be updated to the perimeter value.
bool on_a_perimeter(
	double& u, double& v,///<Surface parameter (u,v).
	size_t& perim_num	///<if function returns true,
						///<the perimete rnumber is output.
)const{return false;};

///Test if input (u,v) is on the perimeter perim_num.
///If u or v is on a perimeter, true will be returned.
bool on_the_perimeter(
	size_t perim_num,	///<a perimete number is input.
	double u, double v	///<Surface parameter (u,v).
)const{return false;};

///Return plane parameter value of a point on the plane. 
///If input point is not on the plane, returned is
///the nearest point parameter of the plane.
MGPosition param(
	const MGPosition&	///< 指定点
) const;

///Obtain parameter space error.
double param_error() const;
double param_error_u() const;
double param_error_v() const;

/// Return ending parameter value.
double param_e_u() const{return mgInfiniteVal;}
double param_e_v() const{return mgInfiniteVal;}

/// パラメータ範囲を返す。
///Return parameter range of the plane(Infinite box).
MGBox param_range() const;

/// Return starting parameter value.
double param_s_u() const{return -mgInfiniteVal;}
double param_s_v() const{return -mgInfiniteVal;}

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u		///<Indicates x is u-value if is_u is true.
	, double x		///<Parameter value.
					///<The value is u or v according to is_u.
) const;

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
MGSurface* part(
	const MGBox& uvbx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
MGCurve* perimeter_curve(size_t i) const{return 0;};

///Return how many perimeters this surface has.
size_t perimeter_num() const{return 0;};

/// Construct perimeter (u,v) parameter position.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
/// t is perimeter parameter line's parameter value of u or v.
MGPosition perimeter_uv(unsigned i,double t) const;

/// 与えられた点にもっとも近い面上の点を返却する。パラメ
/// ータ値も返却する。
///Return the nearest point of the plane from P.
///Function's return value is always true.
int perp_point(
	const MGPosition& P,///< 与えられた点.
	MGPosition& uv,		///<Parameter value of the plane will be output.
	const MGPosition* uvguess=NULL	///< guess.
) const;

///Return all(actually one) foots of perpendicular straight lines from P.
MGPosition_list perps(
	const MGPosition& P	///< Point of a space(指定点).
) const;

///Test if the surface is planar or not.
///Returned is 0(false) if this is not planar, 1(true) if this planar.
int planar(
	MGPlane& plane,		///<Plane that might be closest to this.
						///<Plane is always output even if not planar.
	double& deviation	///<maximum deviation of this from the output plane.
) const;

///Test if part of the surface is planar or not within the tolerance tol.
///The part of the surface is input by the surface parameter range uvbox.
///Returned is 0(false) if this is not planar, 1(true) if planar.
///For plane, planar always returns true.
int planar(
	const MGBox& uvbox,///<This surface parameter range.
	double tol,	///<maximum deviation allowed to regard the sub surface as a plane.
	int* divideU=0///<Direction to subdivide will be output, if this was not planar,
				///<=1: u direction, =0: v direction.
)const;

///与えられた曲線を自身に面直またはベクトル投影して曲線リストを求める.
///投影曲線は面上のパラメータ曲線と3次元曲線としてそれぞれ順番に、
///vec_crv_uv, vec_crvに格納される。
///uv曲線のトレランスはrc_zero()を、3次元曲線はline_zero()をそれぞれ使用している。
///get perpendicular or vector projection curve list.
///uv projection curves are put into vec_crv_uv(rc_zero() is used),
///3d projection curves are put into vec_crv(line_zero() is used) respectively.
///引数：
///		const MGCurve&			crv,		(I/ )	given curve.
///		MGPvector<MGCurve>&	vec_crv_uv,		( /O)	uv projection curve.
///		MGPvector<MGCurve>&	vec_crv,		( /O)	3d projection curve.
///		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
///戻り値：
///		投影曲線の数:		投影曲線が求まった
///		0:			投影曲線が求まらなかった
///		-1:			内部処理エラー
///		-2:			収束処理エラー（収束しなかった）
///追記：引数vecが与えられないとき、面直投影する。
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int project(
	const MGCurve& crv,					///<given curve.
	MGPvector<MGCurve>& vec_crv_uv,		///<uv projection curve.
	MGPvector<MGCurve>& vec_crv,		///<3d projection curve.
	const MGVector& vec = mgNULL_VEC	///<projection vector,
						///<if vec = NULL then calculate perpendicular project.
) const;

int project(
	const MGStraight& sl,	///<given curve.
	MGPvector<MGCurve>& vec_crv_uv,	///<uv projection curve.
	MGPvector<MGCurve>& vec_crv,	///<3d projection curve.
	const MGVector& vec= mgNULL_VEC	///<projection vector,
						///<if vec = NULL then calculate perpendicular project.
)const;

///与えられた曲線を自身に面直またはベクトル投影して曲線リストを求める.
///投影曲線は3次元曲線としてvec_crvに格納される。
///3次元曲線のtoleranceはline_zero()を使用している。
///get perpendicular or vector projection curve list.
///3d projection curves are put into vec_crv(line_zero() is used).
///引数：
///		const MGCurve&			crv,		(I/ )	given curve.
///		MGPvector<MGCurve>&	vec_crv,		( /O)	3d projection curve.
///		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
///戻り値：
///		投影曲線の数:		投影曲線が求まった
///		0:			投影曲線が求まらなかった
///		-1:			内部処理エラー
///		-2:			収束処理エラー（収束しなかった）
///追記：引数vecが与えられないとき、面直投影する。
///Obtain the projected curve of a curve onto the surface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is general world coordinate curves('vec_crv')
int project(
	const MGCurve& crv,						///<given curve.
	MGPvector<MGCurve>& vec_crv,			///<3d projection curve.
	const MGVector& vec = mgNULL_VEC	///<projection vector,
							///<if vec = NULL then calculate perpendicular project.
) const;

/// 入力パラメータをパラメータ範囲でまるめて返却する.
///Round the input uv into parameter range of the plane, 
///return the same value as input.
MGPosition range(const MGPosition& uv) const{return uv;}

/// 与えられた平面との関係を返す.
///Relation of two planes.
MGPSRELATION relation(
	const MGPlane&,	///<Second plane
	MGStraight&		///<Intersection line of the two planes
					///<will be output if intersect.
) const;

/// 与えられた直線との関係を返す。
///Relation between plane and straight line.
MGPSRELATION relation(
	const MGStraight&,	///<Straight line
	MGCSisect&			///<Intersection point if intersect.
) const;

/// 平面のパラメータ表現の起点を返却する。
///Return root point of the plane.
const MGPosition& root_point() const{return m_root_point;}

///Return the space dimension.
size_t sdim() const;

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
MGPvector<MGCurve> skeleton(int density=1)const;

///Obtain all the parameter curves at knots of u and v knot vector.
MGPvector<MGCurve> skeleton_at_knots()const;

///split this fsurface at the parameter param.
void split(
	double param,///<parameter value of this fsurface. if is_u is true, param is u-value,
				///<else v-value.
	bool is_u,	///<indicates if param is u or v of the surface parameter (u,v).
	MGPvector<MGFSurface>& surfaces///<splitted surfaces will be output.
)const;

/// 曲面タイプを返却する。
///Return surface type of the plane.
MGSURFACE_TYPE type() const{return MGSURFACE_PLANE;}

/// 平面のパラメータ表現のｕ方向を返却する。
///Return u-direction vector of the plane.
const MGVector& u_deriv() const{return m_uderiv;}

/// 平面のパラメータ表現のｖ方向を返却する。
///Return v-direction vector of the plane.
const MGVector& v_deriv() const{return m_vderiv;}

/// 点を平面に投影した点の平面のパラメータ表現(u,v)を求める。
///Return uv parameter of the point projected from point p to the plane.
MGPosition uv(const MGPosition& p) const;

/// Vectorを平面に投影したVectorの平面のパラメータ表現(u,v)を求める。
///Return uv parameter of the point projected from point p to the plane.
///p is the end of the vector v originated from the root_point().
MGVector uv(const MGVector& v) const;

std::string whoami()const{return "Plane";};

protected:

///Test if the surface is flat or not within the parameter value rectangle of uvbox.
///Function's return value is:
///	true: if the surface is flat
///  false: if the surface is not falt.
///When this is not falt, the direction that indicates which direction the surface
///should be divided will be output.
///***** the flatness is tested only approximately. This is for exclusive use of
///planar().
bool flat(
	const MGBox& uvbox,
	double tol,		///<Tolerance allowed to regard flat
					///<(Allowed distance from a plane).
	int& direction,	///<   1: u-direction is more non flat,
					///<   0: v-direction is more non flat.
	MGPosition& P,	///<Position of the flat plane will be output.
	MGUnit_vector& N///<Normal of the flat plane will be output.
)const;

///Test if the surface is flat or not within the parameter value rectangle of uvbox.
///Function's return value is:
///	true: if the surface is flat
///  false: if the surface is not falt.
///When this is not falt, the direction that indicates which direction the surface
///should be divided will be output.
///This is the same as flat except that this does not have the arguments P, N.
///***** the flatness is tested only approximately. This is for exclusive use of
///tessellation.
bool flat_tess(
	double u0,double u1,///<u range from u0 to u1.
	double v0,double v1,///<v range from v0 to v1.
	double tol,		///<Tolerance allowed to regart flat
					///<(Allowed distance from a plane).
	bool& direction	///<   1: u-direction is more non flat,
					///<   0: v-direction is more non flat.
)const;

///Intersection of Surface and a straight line.
MGCSisect_list isectSl(
	const MGStraight& sl,
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
)const;

///メンバデータを読み込む関数
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

private:

	MGUnit_vector	m_normal;	///<Normal of the plane 平面の法線ベクトル.
	double		m_d;
					///<Distance from the origin(0,0,0) 平面の陰関数表現
					///< (m_d=ax+by+cz+....where m_normal=(a,b,c,....))
					///< すなわちm_dは原点と平面の距離.
	MGPosition	m_root_point;	///<A point on the plane,
								///<平面のパラメータ表現の基点.
	MGVector	m_uderiv;	///<U direction vector,
							///<平面のパラメータ表現のｕ方向.
	MGVector	m_vderiv;	///<V direction vector,
							///< 平面のパラメータ表現のｖ方向.
	mutable MGKnotVector* m_uknotV;
	mutable MGKnotVector* m_vknotV;
			///<When knot_vector_u,v() is invoked, the knot vector will be set,
			///<These two variables must be initialize to 0(null).

///Return minimum box that includes whole of the surface.
///Returned is a newed object pointer.
MGBox* compute_box() const;

///get the display arrow length vector.
void MGPlane::get_uv_display_vector(
	MGVector& u,
	MGVector& v
)const;

///isect_area_length() returns initial area length for the intersection
///line.
size_t isect_area_length() const{return 10;};

/// "isect_guess" computes one intersection point of surface and a curve,
/// given initail guess parameter values of surface and curve.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int isect_guess(
	const MGCurve& crv,		///<Curve.
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface,
	double ti,			///<Input initial guess parameter value of the line.
	MGPosition& uv,		///< Output parameter value obtained. 
	double& t			///< Output parameter value obtained. 
)const;

/// "isect_guess_straight" computes one intersection point of surface and
///a straight line, given initail guess parameter values of the surface and 
///the straight line.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int isect_guess_straight(
	const MGStraight& sl,	///<Straight line.
	double ti,			///<Initial guess parameter value of the straight.
	const MGPosition& uvi,	///<Input initial guess parameter value,
						///< of the i.p. of the surface. 
	double& t,			///<Straight parameter obtained.
	MGPosition& uvout		///<Surface parameter value obtained(u,v). 
) const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not,
							///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	size_t incr=1		///<Incremental valuse of B-coef's id.
) const;

///"isect_inner_dt" is a dedicated function of isect_startPt,
/// comutes adequate incremental parameter value(du,dv) and parameter line kind
///kdt(u=const or v=const).
void isect_inner_dt(
	size_t n,				///<num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,///<intersection point obtained last(of this).
	double& du, double& dv,	///<incremental length from previous to uvnow is input.
				///<New du or dv will be output according to kdt's return value.
	int& kdt,	///<Parameter kind used so far is input, will be output as:
				///<=1:parameter line kind(u=const), =0: v=const,
				///<=-1:should halt computation since incremental value is zero.
	double acuRatio=1.	///<Accurate ratio.
) const;

///Return order of intersection line order of MGLBRep.
///The default is 4.
size_t isect_order() const{return 4;}

///Compute intersections with MGLBRep lb that does not have C0 continuity in it.
MGCSisect_list isect_withC1LB(const MGLBRep& lb)const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisect_list isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///オフセットするサンプルポイントの1パッチごとの分割数を求める
///全てのパッチ中の分割数で最大の値を返す
int offset_div_num() const;

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& pl). This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
MGSBRep* surf1D(const MGPlane& pl)const{
	assert(false); return 0;
};

};

/** @} */ // end of GEO group
#endif
