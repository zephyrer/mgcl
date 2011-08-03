/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCylinder_HH_
#define _MGCylinder_HH_
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/CSisect_list.h"
#include "mg/Surface.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"

// MGCylinder.h
// Header for class MGCylinder

class MGTransf;
class MGStraight;
class MGSBRep;
class MGIfstream;
class MGOfstream;

/** @addtogroup GEO
 *  @{
 */

///MGCylinder is a Cylinder in 3D space.
///Cylinder is expressed by an ellipse and a straight line.
///Cylinder function  f(u,v) = m_ellipse(u) + m_axis(v),
///where u and v are two parameter of surface representation.
///Here, m_axis is a straight line whose root point is the origin.
/// MGCylinderクラスは３次元空間における円筒面を表すクラスである。
/// MGCylinderクラスでは以下のようなパラメータ表現を使用します。
/// f(u,v) = m_ellipse(u) + m_axis(v);
class MGCLASS MGCylinder :public MGSurface{

public:

MGDECL friend MGCylinder operator+ (const MGVector& v, const MGCylinder& cyl);

///////////////Constructor コンストラクタ//////////////

///Void constructor 初期化なしで柱面を生成する。
MGCylinder(void);

//Copy constructor.
//MGCylinder(const MGCylinder& cyl);

/// Construct a Cylinder by changing this space dimension or
/// ordering the coordinates.
MGCylinder(
	size_t dim,				///< New space dimension.
	const MGCylinder& cyl,	///< Original Cylinder.
	size_t start1=0, 		///< Destination order of new Surface.
	size_t start2=0 		///< Source order of original Surface.
);

/// Construct a cylinder of whole circle whose bottom center is bottom
///and top center is bottom+axis.
MGCylinder(
	const MGPosition& bottom,	///<Location on axis to position the cylinder
						///<defines zero v. 
	const MGVector axis,///<The axis vector for the cylinder. 
	double radius,		///<The radius of the cylinder.
	bool outgoing=true	///<Indicates if the surface normal is going to
						///<outside of the cylinder(true) or inside(false).
);

///Cylinder from an ellipse and the axis straight line.
///When axis'es start point is not ellipse's center, axis'esstart point
///will be set to the center.
MGCylinder(
	const MGEllipse& ellipse,	///<ellispe  of the cylinder
	const MGStraight& axis		///<axis of the cylinder.
			///<axis's root point will be neglected, always be set as
			///<the origin.
);

////////////Destructor////////////////
///~MGCylinder();

////////////Operator overload 演算子の多重定義///////////////

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGCylinder& operator=(const MGGel& gel2);
MGCylinder& operator=(const MGCylinder& gel2);

///Translation of the Cylinder
MGCylinder operator+ (const MGVector& ) const;

///Translation of the Cylinder
MGCylinder operator- (const MGVector& ) const;

///柱面のスケーリングを行い，柱面を作成する。
///Scaling of the Cylinder by a double.
MGCylinder operator* (double) const;

///Scaling of the Cylinder by a double.
MGDECL friend MGCylinder operator* (double scale, const MGCylinder& cyl);

///Transformation of the Cylinder by a matrix.
MGCylinder operator* (const MGMatrix& ) const;

///Transformation of the Cylinder by a MGTransf.
MGCylinder operator* (const MGTransf& ) const;

///Object transformation.
MGCylinder& operator+=(const MGVector& v);
MGCylinder& operator-=(const MGVector& v);
MGCylinder& operator*=(double scale);
MGCylinder& operator*=(const MGMatrix& mat);
MGCylinder& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGCylinder& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGCylinder& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGCylinder& gel2)const{return !(gel2==(*this));};

///Output to IGES stream file.
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
///Box that includes limitted Cylinder by box.
MGBox box_limitted(
	const MGBox& uvrange	///< Parameter Range of the surface.
) const;

///Obtain ceter coordinate of the geometry.
MGPosition center() const;

///Changing this object's space dimension.
MGCylinder& change_dimension(
	size_t sdim,		///< new space dimension
	size_t start1=0, 	///< Destination order of new object.
	size_t start2=0		///< Source order of this object.
); 

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
MGCylinder& change_range(
	int is_u,				///<if true, (t1,t2) are u-value. if not, v.
	double t1,				///<Parameter value for the start of original. 
	double t2				///<Parameter value for the end of original. 
);

///Compute the closest point parameter value (u,v) of this surface
///from a point.
MGPosition closest(const MGPosition& point) const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGCylinder* clone()const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGCylinder* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

/// 与えられた点との距離を返却する。
///Return the distace between Cylinder and the point.
double distance(const MGPosition& point) const;

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
) const{return eval(uv.ref(0),uv.ref(1),ndu,ndv);}

///Evaluate right continuous surface data.
///Evaluate all positional data, 1st and 2nd derivatives.
void eval_all(
	double u, double v,		///< Parameter value of the surface.
	MGPosition& f,			///< Positional data.
	MGVector&   fu,			///< df(u,v)/du
	MGVector&   fv,			///< df/dv
	MGVector&   fuv,		///< d**2f/(du*dv)
	MGVector&   fuu,		///< d**2f/(du**2)
	MGVector&   fvv			///< d**2f/(dv**2)
) const;

///Evaluate right continuous surface data.
///Evaluate all positional data, 1st and 2nd derivatives.
void eval_all(
	const MGPosition& uv,	///< Parameter value of the surface.
	MGPosition& f,			///< Positional data.
	MGVector&   fu,			///< df(u,v)/du
	MGVector&   fv,			///< df/dv
	MGVector&   fuv,		///< d**2f/(du*dv)
	MGVector&   fuu,		///< d**2f/(du**2)
	MGVector&   fvv			///< d**2f/(dv**2)
) const{ eval_all(uv(0),uv(1),f,fu,fv,fuv,fuu,fvv);}

/// Exchange parameter u and v.
MGSurface& exchange_uv();

///Modify the original Surface by extrapolating the specified perimeter.
///The extrapolation is C2 continuous if the order >=4.
///The extrapolation is done so that extrapolating length is "length"
///at the position of the parameter value "param" of the perimeter.
MGCylinder& extend(
	int perimeter,	///<perimeter number of the Surface.
					///< =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	///< parameter value of above perimeter.
	double length,	///<chord length to extend at the parameter param of the perimeter.
	double dk=0.  ///<Coefficient of how curvature should vary at
///    extrapolation start point. When dk=0, curvature keeps same, i.e.
///    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///    i.e. dK/dS=-K/length at extrapolation start point,
///    (S=parameter of arc length, K=Curvature at start point)
///    That is, when dk reaches to 1 from 0, curve changes to flat.
);

/// Return This object's typeID
long identify_type() const;

bool in_range(double u, double v) const;
bool in_range(const MGPosition& uv) const{return in_range(uv[0], uv[1]);};

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
size_t intersect_dnum_u() const{ return m_ellipse.intersect_dnum();};
size_t intersect_dnum_v() const{ return m_axis.intersect_dnum();};

/// Surface と Curve の交点を求める。
///Compute intesection of Cylinder and Curve.
MGCSisect_list isect(const MGCurve& curve)const;
MGCSisect_list isect(const MGStraight& line)const{return isectSl(line);};
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
double knot_u(size_t i)const{return m_ellipse.knot(i);};
double knot_v(size_t j)const{return m_axis.knot(j);};

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const{return m_ellipse.knot_vector();};
MGKnotVector& knot_vector_u(){return m_ellipse.knot_vector();};

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const{return m_axis.knot_vector();};
MGKnotVector& knot_vector_v(){return m_axis.knot_vector();};

/// 柱面を反転する。ノーマルを逆方向にする。
///Negate the normal of the Cylinder.
void negate(
	int is_u	///< Negate along u-direction if is_u is ture, else along v-direction.
);

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const;

///Return the normal of the Cylinder.
const MGEllipse& ellipse() const{return m_ellipse;};
const MGStraight& axis() const{return m_axis;};

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
///Test if a point is on the Cylinder. If on the Cylinder, return true.
bool on(
	const MGPosition& point,	///<A point to test 指定点
	MGPosition& puv				///<Parameter value of the Cylinder will be
								///<returned.
) const;

///Test if input (u,v) is parameter value on a perimeter of the surface.
///If u or v is on a perimeter, (u,v) will be updated to the perimeter value.
bool on_a_perimeter(
	double& u, double& v,	///<Surface parameter (u,v)
	size_t& perim_num	///<if function returns true,
						///<the perimete rnumber is output.
)const;

///Obtain parameter space error.
double param_error() const;

/// Return ending parameter value.
double param_e_u() const{return m_ellipse.param_e();};
double param_e_v() const{return m_axis.param_e();};

/// パラメータ範囲を返す。
///Return parameter range of the Cylinder(Infinite box).
MGBox param_range() const;

/// Return starting parameter value.
double param_s_u() const{return m_ellipse.param_s();};
double param_s_v() const{return m_axis.param_s();};

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u				///<Indicates x is u-value if is_u is true.
	, double x				///<Parameter value.
							///<The value is u or v according to is_u.
) const;

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
MGCylinder* part(
	const MGBox& bx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

/// i must be < perimeter_num().
///When perimeter_num()==0, this function is undefined.
MGCurve* perimeter_curve(size_t i) const;

///Return how many perimeters this surface has.
size_t perimeter_num() const{return 4;};

/// 与えられた点にもっとも近い面上の点を返却する。パラメ
/// ータ値も返却する。
///Return the nearest point of the Cylinder from P.
///Function's return value is always true.
int perp_point(
	const MGPosition& P,///< 与えられた点
	MGPosition& uv,		///<Parameter value of the Cylinder will be output
	const MGPosition* uvguess=NULL	///< guess
) const;

///Return all(actually one) foots of perpendicular straight lines from P.
MGPosition_list perps(
	const MGPosition& P				/// Point of a space(指定点)
) const;

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///Round the input uv into parameter range of the Cylinder, 
///return the same value as input.
MGPosition range(const MGPosition& uv) const;

///Return the space dimension.
size_t sdim() const;

/// 曲面タイプを返却する。
///Return surface type of the Cylinder.
MGSURFACE_TYPE type() const{return MGSURFACE_CYLINDER;}

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
	double tol,		///<Tolerance allowed to regart flat
					///<(Allowed distance from a Cylinder).
	int& direction,	///<   1: u-direction is more non flat.
					///<   0: v-direction is more non flat.
	MGPosition& P,	///<Position of the flat Cylinder will be output.
	MGUnit_vector& N///<Normal of the flat Cylinder will be output.
)const;

///This is the same as flat except that this does not have the arguments P, N.
///***** the flatness is tested only approximately. This is for exclusive use of
///tessellation.
bool flat_tess(
	double u0,double u1,	///<u range from u0 to u1.
	double v0,double v1,	///<v range from v0 to v1.
	double tol,		///<Tolerance allowed to regart flat
					///<(Allowed distance from a Cylinder).
	bool& direction,///<   1: u-direction is more non flat.
					///<   0: v-direction is more non flat.
	double max_edge_len
)const;

///Intersection of Surface and a straight line.
MGCSisect_list isectSl(
	const MGStraight& sl,
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
	) const;

///Obtain the v parameter value of the neareast point from the point P to the axis.
///Function's return value is if the parameter value v is in the range of this
///cylinder: -1: below the range, 0:in the range, 1:above the range.
int vrange(const MGPosition& P, double& v)const;

///メンバデータを読み込む関数
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

std::string whoami()const{return "Cylinder";};

private:

	MGEllipse	m_ellipse;	///<Ellipse of the cylinder
	MGStraight	m_axis;		///<Axis of the cylinder
	bool m_ortho;///<Indicates if normal of m_ellipse is parallel to m_axis,
				///<If parallel, m_ortho is true.

///Compute the axis point of the parameter v.
MGPosition axis_point(double v)const;

///Return minimum box that includes whole of the surface.
///Returned is a newed object pointer.
MGBox* compute_box() const;

///isect_area_length() returns initial area length for the intersection
///line.
size_t isect_area_length() const{return 10;};

///isect_direction() is used by isect_startPt() to define which constant
///parameter line should be used to compute intersection, and what
///incremental value be used for the parameter.
///Function's return value is direction to get next intersection(with dt).
///When =1: u=const direction, =0: v=const, =-1: cannot get intersection.
int isect_direction(
	const MGFSurface& sf2,	///<Second surface for the intersection.
	size_t m1,		///<id of uvuvS that indicates this surface's parameter
		///<position in uvuvS. (uvuvS(m1), uvuvS(m1+1))=(u,v) of this surface.
	MGPosition& uvuvS,///<start parameter (u,v) pair of this surface and sf2.
	double& du,	///<Incremental value of the parameter kind of kdt will be output.
	double& dv, ///<Right dt will be output according to the function's output =0,1.
	double acuRatio=1.	///<acuracy ratio.
)const;

///isect_dt computes incremental values of u and v direction for the intersection
///computation at parameter position (u,v).
void isect_dt(
	double u, double v, double& du, double& dv,
	double acuRatio=1.	///acuracy ratio.
) const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,///<last intersection point.
	int kdt,			///<Input if u=const v-parameter line or not.
						///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,			///<next u value will be output
	double& v,			///<next v value will be output
	size_t incr=1		///<Incremental valuse of B-coef's id.
)const;

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

///Compute the intersection line of this and the plane pl.
MGSSisect_list intersectPl(const MGPlane& pl) const;

///オフセットするサンプルポイントの1パッチごとの分割数を求める
///全てのパッチ中の分割数で最大の値を返す
int offset_div_num() const{return 1;};

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& plane). This surf1D is used in isect for
///the argument of isect_startCylinder, which will use surf1D to compute isect(plane).
///surf1D=0.(intersection with x=0. Cylinder) is the intersection lines.
///This is not used for MGCylinder, should not be invoked.
MGSBRep* surf1D(const MGPlane& plane)const{assert(false);return 0;};

};

/** @} */ // end of GEO group
#endif
