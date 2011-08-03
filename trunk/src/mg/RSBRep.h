/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGRSBRep_HH_
#define _MGRSBRep_HH_

#include "mg/Position.h"
#include "mg/SBRep.h"

// MGRSBRep.h
//

// Forward Declaration
class  MGSPointSeq;
class  MGKnotArray;
class  MGMatrix;
class  MGTransf;
class  MGStraight;
class  MGRLBRep;
class  MGPlane;
class  MGCSisect_list;
class  MGSSisect_list;
class  MGIfstream;
class  MGOfstream;

/** @addtogroup GEO
 *  @{
 */

/// Defines Surface B-Representation of rational form.
/// This NURBS is of homogeneous form, i.e., B-Coefficients have
/// weight included values. 
/// When usual(non-homogeneous) NURBS form is (xij, yij, zij, wij) ,
/// MGRSBRep form is (xij*wij, yij*wij, zij*wij, wij)
///				 for i=0,..., m-1, and j=0,..., n-1.
class MGCLASS MGRSBRep: public MGSurface {

public:

/// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
///Scaling.
MGDECL friend MGRSBRep operator* (double scale, const MGRSBRep& sb);

//////////// Constructor ////////////

///Default constructor(dummy surface brep).
MGRSBRep():MGSurface(){;};

///Construct MGRSBRep from the raw data.
MGRSBRep(
	const MGSPointSeq& bcoef,	
	///<Control Vertex of rational surface B-Rep that
	///<includes weight multiplied when homogeneous=true(1),
	///<and not includes when homogeneous =false.
	///<Mximum space dimension id of bcoef is for weight of the rational.
	const MGKnotVector& tu,		///<knot vector of u-direction
	const MGKnotVector& tv,		///<knot vector of v-direction
	int homogeneous=1
);

///Construct MGRSBRep from the raw data.
MGRSBRep(
	const MGSPointSeq& bcoef,	
	///<Control Vertex of rational surface B-Rep that does not includes weights.
	const MGSPointSeq& weights,	///<weights, weights(i,j,0) is for bcoef(i,j,.)
	const MGKnotVector& tu,		///<knot vector of u-direction
	const MGKnotVector& tv		///<knot vector of v-direction
);

///Construct surface of revolution, given planar MGRLBRep and
///rotation axis sl.
///Parameterization of the surface is:
///	u=const parameter line generates given rlb(when u=0.).
///  v=const parameter line generates a circle whose center is sl.
MGRSBRep(
	const MGRLBRep& rlb,	///<Planar MGRLBRep to rotate.
	const MGStraight& sl,	///<Rotation axis. This is treated as infinite
			///<one, even if it is not.
	double angle			///<Rotation angle in radian,
	///<-2*pai<=angle<=2*pai,
	///<If angle is positive, circle is anti-clockwise around direction Vector N
	///<of sl. If negative, circle is clockwise around N.
);

///Construct MGRSBRep by sweep NURBS and sweep length.
///The sweep surface is defined as:
///rlbrep(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGRSBRep(
	const MGRLBRep& rlbrep,		///<Sweep crv.
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
);

//**** 1. Approximation Constructor ****

///Approximate an original B-Rep by a new knot configuration.
///The new knot config must be inside the range of the original B-Rep
///parameter. However new knots may be coarse or fine.
MGRSBRep(
	const MGRSBRep& old,///<Original B-Rep.
	const MGKnotVector& ut,	///<knot vector of u-direction
	const MGKnotVector& vt,	///<knot vector of v-direction
	int &error				///<Error flag.
);

//**** 2.Conversion Constructor.****

/// Convert from Non ratoinal form to Rational form.
///When homogeneous==true(non zero), brep is homogeneous form MGSBRep.
///When homogeneous==false(zero), brep is ordinary MGSBRep and
///will be converted to MGRSBRep. That is, weight=1 elements will be
///added to the last space dimension element.
///***** This is the fundamental constructor when homogeneous==1. *****
explicit MGRSBRep(
	const MGSBRep& brep,///<Original SBRep. This can be ordinary SBRep, or 
	///<homogeneous form of MGRSBRep. When homogeneous form,
	///<the last space dimension elements are weights.
	int homogeneous=0	///<true(non zero): homogeneous form,
							///<false(zero):ordinary SBRep.
);

///Gets new B-Rep by adding knots to an original B-Rep.
MGRSBRep(
	const MGRSBRep& old,		///<Original B-Rep.
	const MGKnotArray& uknots,	///<Knots to add for u-direction
	const MGKnotArray& vknots	///<Knots to add for v-direction.
);

/*
///Gets new B-Rep by connecting two B-Rep to one.
///The parameter (which1,continuity,which2,opposite) can be obtained by
///public function continuity.
MGRSBRep(
	const MGRSBRep& brep1,	///B-Rep 1.
	int which1,				///which perimeter of brep1.
	int continuity,			///continuity. Must >=0.
	const MGRSBRep& brep2,	///B-Rep 2.
	int which2,				///which perimeter of brep2.
	int opposite			/// Input if parameter direction of which2
					/// is the same as which1 along common edge.
					/// If opposite is true, the direction is opposite.
);
*/

/// Gets new NURBS Surface by computing a part of the original.
///New one is exactly the same as the original except that it is partial.
///If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
///will be guaranteed. Here, n=bdim_u(), k=order_u(),
///t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
///About knot_v(j), the same.
/// Both u-range and v-range must be inside the range of old.
MGRSBRep(
	const MGBox& uvrange,	///<u and v parameter range.
	const MGRSBRep& old,	///<Original B-Rep.
	int multiple=0 ///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

/// Construct a Surface B-Rep by changing space dimension and ordering
///of coordinates.
MGRSBRep(
	size_t dim,				///< New space dimension.
	const MGRSBRep& sbrep,	///< Original Surface B-rep.
	size_t start1=0, 		///< Destination order of new Surface.
	size_t start2=0 		///< Source order of original Surface.
);

///リブ曲線列から面を作成する
///作成する面のノットベクトルはリブ曲線の向きをu,リブ列方向をvとする
///This constructor only generates MGSBRep even if curves are MGRLBRep of the same
///knot configuration. To avoid this, use createSurfaceFromRibs() that generates
///MGRSBRep when curves are MGRLBRep of the same knot configuration.
///Let v0=start parameter value, v1=terminate parameter value along v, then
///v=v0 const parameter line is curves[0], and v=v1 const parameter line is
///curves[n-1], where n=curves.size(). n must be greater or equal to 2.
///When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
///
///If MGRLBRep's in vecPtrRibRLBReps may have different knot configurations,
///use the global function createSurfaceFromRibs(declared in MGSBRep.h).
MGRSBRep(
	const std::vector<const MGRLBRep*>& vecPtrRibRLBReps,
	bool direction_adjustment=true//=true, curves[.] direction are adjusted to line
								//to the same direction.
);

///	MGRSBRep(const MGRSBRep&);  ///Copy constructor.
///  We can use default copy constructor.

///Destructor
///	~MGRSBRep();	///We can use default destructor.

//////////// Operator overload ////////////

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGRSBRep& operator=(const MGGel& gel2);
MGRSBRep& operator=(const MGRSBRep& gel2);

/// 曲面の平行移動を行いオブジェクトを生成する。
///Translation.
MGRSBRep operator+ (const MGVector& ) const;

/// 曲面の逆方向に平行移動を行いオブジェクトを生成する。
///Translation.
MGRSBRep operator- (const MGVector& ) const;

/// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
///Scaling.
MGRSBRep operator* (double) const;

/// 与えられた変換で曲面の変換を行いオブジェクトを生成する。
///Matrix transformation.
MGRSBRep operator* (const MGMatrix& ) const;

/// 与えられた変換で曲面のトランスフォームを行いオブジェクトを生成する。
///General transformation.
MGRSBRep operator* (const MGTransf& ) const;

///Object transformation.
MGRSBRep& operator+=(const MGVector& v);
MGRSBRep& operator-=(const MGVector& v);
MGRSBRep& operator*=(double scale);
MGRSBRep& operator*=(const MGMatrix& mat);
MGRSBRep& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool operator==(const MGRSBRep& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGRSBRep& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGRSBRep& gel2)const{return !(gel2==(*this));};
bool operator==(const MGSBRep& sb)const;

//////////// Member Function ////////////

/// 自身の曲面の全体の面積を返却する。
///Compute total surface area.
///	double area() const;

/// 与えられたパラメータ範囲の曲面の面積を返す。
///Compute surface area limitted by parameter range box.
///	double area(const MGBox& box) const;

///Returns B-Rep Dimension of u.
size_t bdim_u() const{return m_surface.bdim_u();}

///Returns B-Rep Dimension of v.
size_t bdim_v() const{return m_surface.bdim_v();}
	
/// 入力のパラメータ範囲の曲面部分を囲むボックスを返す。
///Compute minimum box that includes the surface.
MGBox box_limitted(const MGBox& bx) const;///Limited surface be the parameter box.

///Changing this object's space dimension.
MGRSBRep& change_dimension(
	size_t sdim,	///< new space dimension
	size_t start1=0,///< Destination order of new object.
	size_t start2=0 ///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
MGRSBRep& change_range(
	int is_u,		///<if true, (t1,t2) are u-value. if not, v.
	double t1,		///<Parameter value for the start of original. 
	double t2		///<Parameter value for the end of original. 
);

///Access to (i,j)th element of coef.
///Left-hand side version.
double& coef(size_t i, size_t j, size_t k)
{	assert(i<bdim_u() && j<bdim_v() && k<=sdim());
	return m_surface.coef(i,j,k);}

///Access to (i,j)th element of coef.
///(right-hand side version).
double coef(size_t i, size_t j, size_t k)const{return m_surface.coef(i,j,k);}

///Extract (i,j,k) elements for 0<=k<sdim() as a vector.
MGVector coef(size_t i, size_t j) const{return m_surface.coef(i,j);}

///Returns a pointer to the surface b-coef data.
const double* coef_data(size_t i=0, size_t j=0, size_t k=0) const
{return m_surface.coef_data(i,j,k);}

/*
///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two surfaces are discontinuous.
///  0: G0 continuity, i.e. two surfaces are connected,
///     but tangents are discontinuous
///  1: G1 continuity, i.e. two surfaces are connected,
///     and tangents are also continuous.
///  2: G2 continuity, i.e. two surfaces are connected,
///     and tangents and curvatures are also continuous.
int continuity(		/// Reuturn value is the continuity.
	const MGRSBRep& brep2,	/// Input second RSBRep
	int& which1,	/// Outputs which perimeter(which1) of this is
	int& which2,	/// connected to which(which2) of brep2.
					/// These are valid only when continuity>=0.
	int& opposite,	/// Outputs if parameter direction of which2
					/// is the same as which1 along common edge.
					/// If opposite is true, the direction is opposite.
	double& ratio	/// Ratio of 1st derivatives of the two surfaces will
					/// be returned.
			/// ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;
*/

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGRSBRep* clone() const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGRSBRep* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Display control polygons using mgGDL::MGDrawPointSeq(sp)
void display_control_polygon()const;

///uまたはv方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot(
    MGPvector<MGSurface>& srfl      ///<分割した曲面リスト
) const;

///Evaluate right continuous ndu'th and ndv'th derivative data.
///Function's return value is (d(ndu+ndv)f(u,v))/(du**ndu*dv**ndv).
/// ndu=0 and ndv=0 means positional data evaluation.
MGVector eval(
	double u, double v,	///< Parameter value of the surface.
	size_t ndu=0,	///< Order of Derivative along u.
	size_t ndv=0	///< Order of Derivative along v.
) const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv	///< Parameter value of the surface.
	, size_t ndu=0			///< Order of derivative along u.
	, size_t ndv=0			///< Order of derivative along v.
) const;

///Evaluate right continuous surface data.
///Evaluate all positional data and 1st and 2nd derivatives.
void eval_all(
	double u, double v,		///< Parameter value of the surface.
	MGPosition& f,			///< Positional data.
	MGVector&   fu,			///< df(u,v)/du
	MGVector&   fv,			///< df/dv
	MGVector&   fuv,		///< d**2f/(du*dv)
	MGVector&   fuu,		///< d**2f/(du**2)
	MGVector&   fvv			///< d**2f/(dv**2)
) const;

///Evaluate all of i and j'th derivative data for 0<=i<=ndu, 0<=j<=ndv.
/// Output. (d(i+j)f(u,v))/(du**i*dv**j) in deriv[r+j*dim+i*ndv*dim]
///for 0<=r<dim=sdim(), 0<=i<=nderiv and 0<=j<sdim(). 
void eval_all(
	double u, double v,		///<Parameter value to evaluate.
	size_t ndu, size_t ndv,	///<Order of Derivative along u and v direction.
	double* deriv	///<Output. (d(i+j)f(u,v))/(du**i*dv**j) in
					///<deriv[r+j*dim+i*(ndv+1)*dim] for 0<=r<dim=sdim().
					///<for 0<=i<=ndu and 0<=j<=ndv.
					///<deriv is an array of deriv[ndu+1][ndv+1][r].
) const;

/// Exchange parameter u and v.
MGSurface& exchange_uv(){m_surface.exchange_uv(); return *this;};

///Modify the original Surface by extrapolating the specified perimeter.
///The extrapolation is C2 continuous if the order >=4.
///The extrapolation is done so that extrapolating length is "length"
///at the position of the parameter value "param" of the perimeter.
MGRSBRep& extend(
	int perimeter,	///<perimeter number of the Surface,
					///< =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	///< parameter value of above perimeter.
	double length,	///<chord length to extend at the parameter param of the perimeter.
	double dk=0.  ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.,
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point,
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

///Return homogeneous Surface B-Representation of the rational B-Spline.
const MGSBRep& homogeneous() const {return m_surface;};

/// Return This object's typeID
long identify_type() const;

///Test if input parameter value is inside parameter range of the surface.
bool in_range(double u, double v) const{return m_surface.in_range(u,v);};
bool in_range(const MGPosition& uv) const{return m_surface.in_range(uv);};

/// Surface と Curve の交点を求める。
///Compute curve and surface intersection point(s)
MGCSisect_list isect(const MGCurve& curve)const;
MGCSisect_list isect(const MGStraight& sl)const{return isectSl(sl);};
MGCSisect_list isect(const MGRLBRep& curve)const;
MGCSisect_list isect(const MGEllipse& curve)const;
MGCSisect_list isect(const MGLBRep& curve)const;
MGCSisect_list isect(const MGSurfCurve& curve)const;
MGCSisect_list isect(const MGBSumCurve& curve)const;

///Surface と Surface の交線を求める。
///Surface and Surface intersection.
///Compute intersectio line(s) of two surface.
///Restriction:Currently if two surface do not have intersection on
///any of 4 perimeters, this function does not compute surface to surface
///intersection.
MGSSisect_list isect(const MGSurface& srf2)const;
MGSSisect_list isect(const MGPlane& srf2)const;
MGSSisect_list isect(const MGSphere& srf2)const;
MGSSisect_list isect(const MGCylinder& srf2)const;
MGSSisect_list isect(const MGSBRep& srf2)const;
MGSSisect_list isect(const MGRSBRep& srf2)const;
MGSSisect_list isect(const MGBSumSurf& srf2)const;

///Access to i-th element of u knot
///( left-hand side version)
double& knot_u(size_t i){return m_surface.knot_u(i);}

///Access to i-th element of u knot
///(right-hand side version)
double knot_u(size_t i) const{return m_surface.knot_u(i);}

double& knot_v(size_t i){return m_surface.knot_v(i);}

///Access to i-th element of v knot
///(right-hand side version)
double knot_v(size_t i) const{return m_surface.knot_v(i);}

///Returns a pointer to the u knot vector data.
const double* knot_data_u() const{return m_surface.knot_data_u();}

///Returns a pointer to the v knot vector data.
const double* knot_data_v() const{return m_surface.knot_data_v();}

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const
{return m_surface.knot_vector_u();}
MGKnotVector& knot_vector_u(){return m_surface.knot_vector_u();}

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const
{return m_surface.knot_vector_v();}
MGKnotVector& knot_vector_v(){return m_surface.knot_vector_v();}

///Compare two parameter values. If uv1 is less than uv2, return true.
///Comparison is done after prjected to i-th perimeter of the surface.
bool less_than(
	size_t i,	///<perimeter number.
	const MGPosition& uv1,
	const MGPosition& uv2
)const{return m_surface.less_than(i,uv1,uv2);};

/// 自身に指定したパラメータ範囲のlimitをつける。
///Update this by limitting the parameter range.
///uvrange is parameter value range of (umin, vmin) to (umax, vmax).
MGRSBRep& limit(const MGBox& uvrange){m_surface.limit(uvrange); return *this;}

///Change direction of the surface.
void negate(			
	int is_u		///< Negate along u-direction if is_u is ture,
					///< else along v-direction.
){m_surface.negate(is_u);}

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const
{ return m_surface.negate_param(uv,is_u);}

///Return non_homogeneous B-Coefficients with weights of
///the rational Surface B-Spline. This MGSPointSeq includes weights.
MGSPointSeq non_homogeneous_bcoef() const;

///Test if this is actually non_rational, i.e. , all of the weights are
///same values.
int non_rational() const;

///Returns the B-Rep order(u-direction).
unsigned order_u() const{return m_surface.knot_vector_u().order();}

///Returns the B-Rep order(v-direction).
unsigned order_v() const{return m_surface.knot_vector_v().order();}

/// Return ending parameter value.
MGPosition param_e() const{return m_surface.param_e();}
double param_e_u() const{return m_surface.param_e_u();}
double param_e_v() const{return m_surface.param_e_v();}

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u			///<Indicates x is u-value if is_u is true.
	, double x			///<Parameter value.
						///<The value is u or v according to is_u.
) const;

/// Compute parameter line.
MGRLBRep parameter_line(
	int is_u			///<Indicates x is u-value if is_u is true.
	, double x			///<Parameter value.
						///<The value is u or v according to is_u.
) const;

/// パラメータ範囲を返す。
///Return parameter range.
MGBox param_range() const;

/// Return starting parameter value.
MGPosition param_s() const{return m_surface.param_s();}
double param_s_u() const{return m_surface.param_s_u();}
double param_s_v() const{return m_surface.param_s_v();}

///Compute part of the surface limitted by the parameter range bx.
///bx(0) is the parameter (us,vs) and bx(1) is (ue,ve).
///That is u range is from us to ue , and so on.
MGRSBRep* part(
	const MGBox& bx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

///Retrieve perimeter i of this surface.
/// Compute perimeter Rational line B-Rep.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
MGRLBRep perimeter(size_t i) const;
MGCurve* perimeter_curve(size_t i) const;

///Return how many perimeters this surface has.
size_t perimeter_num() const{return 4;};

///Test if the RSBRep is planar or not.
///Returned is 0(false) if this is not planar, 1(true) if this planar.
int planar(
	MGPlane& plane,		///<Plane that might be closest to this.
						///<Plane is always output even if not planar.
	double& deviation	///<maximum deviation of this from the output plane.
) const;

///Test if part of the surface is planar or not within the tolerance tol.
///The part of the surface is input by the surface parameter range uvbox.
///Returned is 0(false) if this is not planar, 1(true) if planar.
int planar(
	const MGBox& uvbox,///<This surface parameter range.
	double tol,	///<maximum deviation allowed to regard the sub surface as a plane.
	int* divideU=0///<Direction to subdivide will be output, if this was not planar,
				///<=1: u direction, =0: v direction.
	) const;

/// 入力パラメータをパラメータ範囲でまるめて返却する。
///Round the input parameter value uv into 
///the parameter range of the surface.
MGPosition range(const MGPosition& uv)const{return m_surface.range(uv);}

///Change the B-Rep by decreasing B-Rep dimension by ndec. This is
///an approximation of the origimal B-Rep. Return value is error flag.
int reduce(
	int is_u,	///<if true, reduce b-rep dimension of u-direction.
	int ndec	///<Number of B-rep dimension to decrease .
){return m_surface.reduce(is_u,ndec);}

///Change an original B-Rep to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
MGRSBRep& refine(
	const MGKnotVector& uknot,	///< new knot of u-direction
	const MGKnotVector& vknot	///< new knot of v-direction
){m_surface.refine(uknot,vknot); return *this;}

///ノット削除関数
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///removal knot. line_zero tolerance is used.
void remove_knot();

///Returns the space dimension.
size_t sdim() const{return m_surface.sdim()-1;}

///Shrink this surface to the part limitted by the parameter range of uvbx.
///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
void shrink_to_knot(
	const MGBox& uvbx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
){m_surface.shrink_to_knot(uvbx,multiple);};

///Returns the B-coef's.
///Right hand side version.
const MGSPointSeq& surface_bcoef() const{return m_surface.surface_bcoef();}	

///Returns the B-coef's.
///Left hand side version.
MGSPointSeq& surface_bcoef(){return m_surface.surface_bcoef();}

///Compute surface integral of the 1st two coordinates.
///(面積分）を求める。
///This integral can be used to compute volume sorounded by the surface.
///double surface_integral(const MGBox&) const;

/// 曲面のタイプをを返す。
///Return the surface type.
MGSURFACE_TYPE type() const{return MGSURFACE_RSPLINE;}

/// ｌｉｍｉｔをはずす。
///Unlimit the parameter range. Return the same.
MGSurface& unlimit(){return *this;};

public:

///output to IGES, PD128
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
std::ostream& out(std::ostream&) const;

std::string whoami()const{return "RSBRep";};

protected:

///Intersection of Surface and a straight line.
MGCSisect_list isectSl(
	const MGStraight& sl,
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
)const;

///メンバデータを読み込む関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void WriteMembers(MGOfstream& buf) const;

private:

//////////////Member Data//////////////
	MGSBRep m_surface;		/// Maximum space dimension id is for weights.

///Obtain coefficient's space dimension.
///This function is used in isect_start etc.
size_t coef_sdim() const{return m_surface.sdim();};

///Return minimum box that includes whole of the surface.
///Returned is a newed object pointer.
MGBox* compute_box() const;

/*
///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two surfaces are discontinuous.
///  0: G0 continuity, i.e. two surfaces are connected,
///     but tangents are discontinuous
int continuity(		/// Reuturn value is the continuity.
	const MGRSBRep& brep2,	/// Input second RSBRep
	int is_u1,		/// Input if u-direction of this.
	int is_u2,		/// Input if u-direction of brep2.
	int opposite,	/// Input if parameter direction of which2 is equal or not.
	int& which1,	/// Outputs which perimeter(which1) of this is
	int& which2,	/// connected to which(which2) of brep2.
					/// These are valid only when continuity>=0.
	double& ratio	/// Ratio of 1st derivatives of the two surfaces will
					/// be returned.
			/// ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;
*/

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
size_t intersect_dnum_u() const;
size_t intersect_dnum_v() const;

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
) const;

///Return order of intersection line order of MGLBRep.
size_t isect_order() const;

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& pl). This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
MGSBRep* surf1D(
	const MGPlane& pl
)const;

///u方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot_u(
    MGPvector<MGRSBRep>& srfl) const;      ///<分割した曲面リスト

///v方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot_v(
    MGPvector<MGSurface>& srfl) const;      ///<分割した曲面リスト

};

///Compute binominal coefficients (i,j) for 0<=i<=m and 0<=j<=i
///in bc[(m+1)*i+j].
MGEXTERN void MGBinominal(size_t m, double* bc);

/** @} */ // end of GEO group
#endif
