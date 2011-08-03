/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFSurface_HH_
#define _MGFSurface_HH_

#include <deque>
#include "mg/Default.h"
#include "mg/Unit_vector.h"
#include "mg/Position_list.h"
#include "mg/Pvector.h"

class MGCompositeCurve;
class MGHHisect_vector;
class MGSSisect_list;
class MGCSisect_list;
class MGObject;
class MGCurve;
class MGFace;
class MGShell;

/** @addtogroup MGObjectRelated
 *  @{
 */

///Define MGFSurface Class.
///MGFSurface is an abstract class to provide the comman interfaces to
///MGFace and MGSurface.
class MGCLASS MGFSurface{

friend class MGFace;
friend class MGSurface;

public:
//////// Constructor ////////

///Null FSurface.
MGFSurface(){;}

///Copy constructor.
MGFSurface(const MGFSurface& fsurf){;};

//////////// Virtual Destructor ////////////
virtual ~MGFSurface(){;};

//////// operator overload////////

///Comparison operator.
bool operator<(const MGFSurface& f2)const;
bool operator>(const MGFSurface& f2)const{return f2<(*this);};

/////////Member function////////

///Generate arrow data of the tangent along u and v and the normal
///at the parameter value (u,v) of the FSurface.
///data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
///data[2], [3] are two bottoms of u-tangent arrowhead.
///data[0], [4], [5], [6] are the points of v-tangent arrow.
///data[0], [7], [8], [9] are the points of v-tangent arrow.
virtual void arrow(double u,double v, MGPosition data[10])const=0;
virtual void arrow(const MGPosition& uv, MGPosition data[10])const=0;

///Get the box of the object.
const MGBox& get_box() const;

///Return box of the parameter space of the FSurface.
///After trimmed one.
virtual const MGBox box_param2()const=0;

///Get the clone of this MGFSurface.
virtual MGFSurface* clone_fsurface()const=0;

///Get the clone of this as a MGFace.
///If this is MGSurface, it is converted to MGFace.
virtual MGFace* clone_as_face()const=0;

///Compute closest point from a point.
///Returned is the parameter value of the FSurface that is closest to point.
virtual MGPosition closest(const MGPosition& point)const=0;

///Compute closest point from a line to the boundary of the MGFSurface.
///Returned is the parameter value of the FSurface that is closest to point.
virtual MGPosition closest_on_boundary(const MGStraight& sl) const=0;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWireFS(
	double span_length,	///Line segment span length.
	int line_density=1	///line density to draw a surface in wire mode.
)const;

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void drawWireFS_to_highlight(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const;

///Evaluation
///Input parameter value is not checked if it is in_range() or not.
///Even if it is not in_range(), surFSurface evaluation will be executed.
virtual MGVector eval(
	double u, double v,	///<FSurface parameter value(u,v)
	 size_t ndu=0, size_t ndv=0///<Order of derivative.
)const=0;

virtual MGVector eval(
	const MGPosition& uv,	///<FSurface parameter value(u,v)
	size_t ndu=0, size_t ndv=0///<Order of derivative.
)const=0;

///Evaluate deviations of two faces(this and face2) at npoint
///discrete points.
///(1)Search the common edges which have the distance within tolerance.
///(2)Compute the nearest points from npoint discrete points of this to face2.
///Let uvuvi=uvuvs[i], then
///uvuvi[0], [1] are this face's parameter value(u1,v1), and uvuvi[2], [3] are
///parameter value(u2,v2) of face2 which is the nearest point from the point (u1, v1).
void eval_discrete_deviation(
	const MGFSurface& face2,
	std::vector<MGPosition>& uvuvs,
	int npoint=20,		///<indicates how many discrete points be obtained.
	double tolerance=0.1///<tolerance to get two edge to compute deviation.
)const;

///Obtain all the boundaries(i.e., outer boundary and all the inner boundaries)
MGPvector<MGCurve> get_all_boundaries(void)const;

///get face pointer if this is MGFace, else null will be returned.
virtual MGFace* get_face_pointer()=0;
virtual const MGFace* get_face_pointer()const=0;

///Get number of inner boundaries as the output of the function.
virtual size_t get_number_of_boundaries()const=0;

///get surface pointer. Null will never be returned if this is valid MGFSurface.
///That is, if this is MGFace, base surface will be returned.
virtual MGSurface* get_surface_pointer()=0;
virtual const MGSurface* get_surface_pointer()const=0;

///Test if this and 2nd object has common area about their box(),
///taking error into account.
virtual bool has_commonFS(const MGObject& obj2) const=0;

///Test if this FSurface has inner boundary loops or not.
///If this has one, return true.
bool hasInnerBoundaryLoop()const{ return number_of_inner_boundaries()>0;};

///Test if this FSurface has boundary loops or not in the specified box.
///If this has one, return true.
virtual bool hasLoop(const MGBox& uvbox)const{return false;};

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the FSurface. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
virtual MGPvector<MGCurve> inner_boundary(size_t i)const=0;

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the FSurface. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
virtual MGPvector<MGCurve> inner_boundary_param(size_t i)const=0;

///Test if parameter value (u,v) is in the range of the FSurface parameter.
virtual bool in_range(double u, double v)const=0;
virtual bool in_range(const MGPosition& uv) const=0;

///Test if (u,v) is inside the face.
///Function's return value is:
///  0:outside the face.
///  1:unknown.
///  2:inside the face, not on a boundary.
///  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
///  4:(u,v) is on the outer boundary.
///  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int in_range_with_on(double u, double v)const{return in_range_with_on(MGPosition(u,v));};
virtual int in_range_with_on(const MGPosition& uv)const=0;

///Intersection.
virtual MGHHisect_vector isect(const MGShell& shell2) const=0;
virtual MGSSisect_list isect(const MGFSurface& fsurf) const=0;
virtual MGSSisect_list isect(const MGFace& fsurf) const=0;
virtual MGSSisect_list isect(const MGSurface& fsurf) const=0;
virtual MGCSisect_list isect(const MGCurve& curv) const=0;

///Compute all the intersection points of this face's boundaries(out or inner)
///with face2, and vice versa.
///These intersection points are used to compute surface to surface
///intersection lines.
void intersect12Boundary(
	const MGFSurface& face2,		///<The second surface.
	MGPosition_list& uvuv_list	///<The intersection points will be output,
			///<One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is
			///<a parameter of this surface and (u2,v2) is a parameter of 
			///<surf.
)const;

///Compute intersection points of this face's boundary(outer and inners) with
///face2. If intersection points are found and the boundary is a loop,
///the point's edge pointer(of this) will be stored in a member uvuv of uvuvs.
///uvuv[7] is the edge pointer. If the boundary is not a loop(that is, a perimeter of
///Surfaces), uvuv.sdim()==7 and an edge pointer is not returned.
///When uvuv.sdim()==8, the edge pointer of uvuv[7] is accessed through union mgEdgeP.
///uvuvs[i] is i-th intersection points.
virtual size_t isect_boundary(
	const MGFSurface& face2,
	MGPosition_list& uvuvs,
	///<id1 and id2 are the ids of uvuv where this face's and f2's parameters
	///<are to be stored in a member of uvuvs.
	///<This face's (u,v) is stored in uvuv(id1) and (id1+1).
	///<f2's (u,v) is stored in uvuv(id2) and (id2+1).
	///<id2=0 if id1=2, and id2=2 if id1=0.
	size_t id1=0
)const=0;

///Compute intersection points between the boundary of iid-th inner boundary
///of this face and face2 to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
virtual size_t isect_incurves(
	const MGFSurface& face2,
	size_t iid,	///<Inner loop id of this face(from 0)
	MGPosition_list& uvuv_list,	///<intersection points will be appended,
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1	///<id of uvuv(a member of uvuv_list),
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter,
		///<id2=0 if id1=2, and id2=2 if id1=0.
) const=0;

///Compute intersection points of outer boundary curves of this face 
///with face2 to compute intersections.
///Function's return value is the number of ip's obtained(appended)
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
virtual size_t isect_outcurves(
	const MGFSurface& face2,
	MGPosition_list& uvuv_list,	///<intersection points will be appended,
		///<One member in the list is of sdim 7,
		///<and the last three elements are the ip direction vector.
	size_t id1			///<id of uvuv(a member of uvuv_list),
		///<uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter,
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const=0;

/// "isect_guess" computes one intersection point of surface and a curve,
/// given initail guess parameter values of surface and curve.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
virtual int isect_guess(
	const MGCurve& crv,		///<Curve.
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface. 
	double ti,			///<Input initial guess parameter value of the line.
	MGPosition& uv,		///< Output parameter value obtained. 
	double& t			///< Output parameter value obtained. 
)const;

/// "isect_guess" computes one intersection point of surface and a curve,
/// given initail guess parameter values of surface and curve.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
virtual int isect_guess(
	const MGStraight& sl,	///<Curve.
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface. 
	double ti,			///<Input initial guess parameter value of the line.
	MGPosition& uv,		///< Output parameter value obtained. 
	double& t			///< Output parameter value obtained. 
)const{return isect_guess_straight(sl,ti,uvi,t,uv);};

/// "isect_guess" computes one intersection point of surface and a curve,
/// given initail guess parameter values of surface and curve.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
virtual int isect_guess(
	const MGCompositeCurve& crv,	///<Curve
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface. 
	double ti,			///<Input initial guess parameter value of the line.
	MGPosition& uv,		///< Output parameter value obtained. 
	double& t			///< Output parameter value obtained. 
)const{return isect_guess_composite(crv,uvi,ti,uv,t);};

/// "isect_guess" computes one intersection point of surface and a curve,
/// given initail guess parameter values of surface and curve.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
virtual int isect_guess_composite(
	const MGCompositeCurve& crv,	///<Curve
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface. 
	double ti,			///<Input initial guess parameter value of the line.
	MGPosition& uv,		///< Output parameter value obtained. 
	double& t			///< Output parameter value obtained. 
)const;

/// "isect_guess_straight" computes one intersection point of surface and
///a straight line, given initail guess parameter values of the surface and 
///the straight line.
///Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
virtual int isect_guess_straight(
	const MGStraight& sl,	///<Straight line.
	double ti,			///<Initial guess parameter value of the straight.
	const MGPosition& uvi,	///<Input initial guess parameter value
						///< of the i.p. of the surface. 
	double& t,			///<Straight parameter obtained.
	MGPosition& uv		///<Surface parameter value obtained(u,v). 
)const;

///Access to i-th element of u knot
virtual double knot_u(size_t i) const{return 0.0;}

///Access to i-th element of v knot
virtual double knot_v(size_t i) const{return 0.0;}

///Returns the u knot vector.
virtual const MGKnotVector& knot_vector_u() const=0;
virtual MGKnotVector& knot_vector_u()=0;

///Returns the v knot vector.
virtual const MGKnotVector& knot_vector_v() const=0;
virtual MGKnotVector& knot_vector_v()=0;

///This is a newed MGFace or MGSurface object.
///If this is a MGFace, returns this pointer.
///If this is a MGSurface, construct a newed MGFace using this newed MGSurface,
///and returns the MGFace*.
virtual MGFace* make_face()=0;

///Make a display list without color of this gel.
///Return is the display list name.
size_t make_display_list_to_hilightFS(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///Negate the FSurface.
virtual void negateFS()=0;

///Compute normal vector(not unit) at uv.
virtual MGVector normal(const MGPosition& uv) const=0;

///Compute normal vector(not unit) at (u,v).
virtual MGVector normal(double u,double v) const=0;

///Get number of inner boundaries as the output of the function.
virtual size_t number_of_inner_boundaries()const{return 0;};

///Get the object point of this MGFSurface.
virtual const MGObject* object_pointer()const=0;
virtual MGObject* object_pointer()=0;

///Offset.
///distance is plus value if the direction is toward normal vector of the
///FSurface. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
virtual int offset_fs(double distance, MGPvector<MGFSurface>& vecOfsFSurface)const=0;

///Test if a point P is on the FSurface.
///Returned is true if the point P is on the FSurface.
///false(0) if P was not on the FSurface.
virtual bool on(const MGPosition& P,
		MGPosition& uv	///<Parameter value of the FSurface is returrned,
						///<Even if P is not on the FSurface, nearest point
						///<parameter value will be returned.
)const=0;

///Test if input (u,v) is parameter value on a perimeter of the base surface.
///If u or v is on a perimeter, they will be updated to the perimeter value.
virtual bool on_a_perimeter(
	double& u, double& v,	///<Surface parameter (u,v)
	size_t& perim_num	///<if function returns true,
						///<the perimete number is output.
)const=0;

/// Output virtual function.
virtual std::ostream& outFS(std::ostream& ostrm) const=0;

///Obtain outer_boundary curves(world coordinates representation) of the FSurface.
///Let the output of outer_boundary() be wcurves and of outer_boundary_param()
///be pcurves, then wcurves[i] corresponds to pcurves[i] one to one.
virtual MGPvector<MGCurve> outer_boundary()const=0;

///Obtain boundary curves(parameter space representation) of the FSurface.
///Let the output of boundary() be wcurves and of boundary_parameter()
///be pcurves, then wcurves[i] corresponds to  pcurves[i] one to one.
virtual MGPvector<MGCurve> outer_boundary_param()const=0;

///Obtain parameter value of the FSurface whose world coordinates are P.
virtual MGPosition param(const MGPosition& P)const=0;

/// Return ending parameter value.
virtual double param_e_u()const=0;
virtual double param_e_v()const=0;

/// Return starting parameter value of the base surface.
virtual double param_s_u()const=0;
virtual double param_s_v()const=0;

///Obtain parameter curves.
///In the case of MGSurface, parameter curve is only one. However, in the case
///of MGFace,  number of parameter curves are more than one.
virtual MGPvector<MGCurve> parameter_curves(
	int is_u,	///<True(!=0) if x is u-value.(i.e. obtain u=const line)
	double x	///<parameter value. u or v-value accordint to is_u.
)const=0;

///Obtain parameter space error.
double param_error() const;
double param_error_u() const;
double param_error_v() const;

/// パラメータ範囲を返す。
///Return parameter range.
virtual MGBox param_range() const=0;

///Return the foot of the perpendicular straight line from P.
///Computation is done from the guess parameter value.
///Function's return value is whether point is obtained(true) or not(false).
virtual bool perp_guess(
	const MGPosition& P,		///<Point
	const MGPosition& uvguess,	///< guess parameter value of the shell
	MGPosition& uv				///< Parameter value will be returned.
)const=0;

///Compute perpendicular points of a curve and the FSurface,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
virtual bool perp_guess(
	const MGCurve& curve,	///<curve.
	const MGPosition& uvguess,	///<Guess parameter value of the FSurface.
	double tguess,			///<Guess parameter value of the curve.
	MGPosition& uv,			///<perpendicular point's parameter values of the shell
	double& t				///<will be output.
)const=0;

///指定点から最も近い、垂線の足とパラメータ値を返す。
///Return the foot of the perpendicular straight line from p that is 
///nearest to point p.
/// Function's return value is whether point is obtained(1) or not(0)
virtual int perp_point(
	const MGPosition& p,///< 指定点(point)
	MGPosition& uv,		///<Parameter value of the surFSurface will be returned.
	const MGPosition* uvguess=0	///< guess parameter value of surFSurface
)const=0;

///Compute perpendicular points on the FSurface from a point P((x,y,z)).
///MGPosition uv in the MGPosition_list is:
///uv(0): u parameter, and uv(1): v parameter of the FSurface.
///Generally number of uv are more than one.
virtual MGPosition_list perps(const MGPosition& P) const=0;

///指定点から最も近い、垂線の足とパラメータ値を返す。
///Return the foot of the perpendicular straight line from p that is 
///nearest to point P.
///Function's return value is whether point is obtained(>0) or not(0)
virtual int perp_one(
	const MGPosition& P, ///< 指定点(point)
	MGPosition& uv 		///<Parameter value of the FSurface will be returned.
)const;

///与えられた曲線を自身に面直またはベクトル投影して曲線リストを求める。引数vecが与えられないとき、面直投影する。
///投影曲線は面上のパラメータ曲線と3次元曲線としてそれぞれ順番に、vec_crv_uv, vec_crvに格納される。
///uv曲線のトレランスはrc_zero()を、3次元曲線はline_zero()をそれぞれ使用している。
///戻り値：
///		投影曲線の数:		投影曲線が求まった
///		0:			投影曲線が求まらなかった
///		-1:			内部処理エラー
///		-2:			収束処理エラー（収束しなかった）
///Obtain the projected curve of a curve onto the FSurface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the FSurface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the FSurface(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
///(vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
virtual int project(
	const MGCurve& crv,						///<given curve.
	MGPvector<MGCurve>& vec_crv_uv,
		///<Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		///<Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec = mgNULL_VEC	///<projection vector,
							///<if vec = NULL then calculate perpendicular project.
) const;

///与えられた曲線を自身に面直またはベクトル投影して曲線リストを求める。引数vecが与えられないとき、面直投影する。
///投影曲線は3次元曲線としてvec_crvに格納される。
///uv曲線のトレランスはline_zero()を使用している。
///戻り値：
///		投影曲線の数:		投影曲線が求まった
///		0:			投影曲線が求まらなかった
///		-1:			内部処理エラー
///		-2:			収束処理エラー（収束しなかった）
///Obtain the projected curve of a curve onto the FSurface.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the FSurface if the vec is NULL.
///Output of 'project' is general world coordinate curves('vec_crv')
virtual int project(
	const MGCurve& crv,	///<given curve.
	MGPvector<MGCurve>& vec_crv,
		///<Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec = mgNULL_VEC
		///<if vec = NULL then calculate perpendicular project.
)const;

///面直に投影した点を返却する
///戻り値は、交点または面直点が求まったときは1、求まらなかったときは0を返却する
int project_normal(
	const MGPosition& pos,
	const MGPosition& uv_guess,	///<推量パラメータ
	MGPosition& uv
)const;

///Round the input parameter (u,v) of the FSurface to the nearest point of
///the FSurface parameter range.
virtual MGPosition range(const MGPosition& uv)const=0;

//Obtain main parameter lines of the FSurface without boundaries.
//inner_skeleton includes only inner parameter lines without boundaries.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
MGPvector<MGCurve> inner_skeleton(int density)const;

///Obtain boundary and main parameter lines of the FSurface.
///skeleton includes boundary() and inner parameter lines.
///density indicates how many inner parameter lines are necessary
///for both u and v directions.
virtual MGPvector<MGCurve> skeleton(int density=1)const;

///Obtain all the parameter curves at knots of u and v knot vector.
virtual MGPvector<MGCurve> skeleton_at_knots()const;

///split this fsurface at the parameter param.
virtual void split(
	double param,///<parameter value of this fsurface. if is_u is true, param is u-value,
				///<else v-value.
	bool is_u,	///<indicates if param is u or v of the surface parameter (u,v).
	MGPvector<MGFSurface>& surfaces///<splitted surfaces will be output.
)const=0;

///split this fsurface with splitters. splitters are 2D (u,v) surfaces's parameter curves.
void split(
	const std::vector<const MGCurve*>& splitters,	//splitter world curves.
	const MGVector&  dir,	//splitter projection direction.
							//If dir.is_null(), normal projection will be performed.
	MGPvector<MGFace>& faces//Result splitted face(s) will be appended.
			//If no splitting was performed, no faces will be appended.
)const;

///split this fsurface with splitters. splitters are 2D (u,v) surfaces's parameter curves.
void split(
	const MGPvector<MGCurve>& splitters,//splitter (u,v) curves.
	MGPvector<MGFace>& faces//Result splitted face(s) will be appended.
			//If no splitting was performed, no faces will be appended.
)const;

///Extract a sub surface with trimmers. trimmers are 3D curves and will be projected
///onto this surface tword the direction dir. If dir is null vector, surface normal
///prjection will be performed. Extraction is so performed that the smallest region
///enclosed by trimmers that includes the surface point uv will be extracted. 
void extract(
	const std::vector<const MGCurve*>& trimmers,	///<Trimmer curves
	const MGVector&  dir,	///<trimmers projection direction.
	const MGPosition& uv,	///<surface parameter (u,v) that indicates the region to extract.
							///<The smallest region that inclued uv will be extracted.
	std::auto_ptr<MGFace>& eface///<Result extracted face will be output.
)const;

///Trim this fsurface with trimmers. trimmers are 3D curves and will be projected
///onto this surface tword the direction dir. If dir is null vector, surface normal
///prjection will be performed. Trimming is so performed that the smallest region enclosed
///by trimmers that includes the surface point uv will be removed. 
void trim(
	const std::vector<const MGCurve*>& trimmers,	///<Trimmer curves
	const MGVector&  dir,	///<trimmers projection direction.
	const MGPosition& uv,	///<surface parameter (u,v) that indicates the region to remove,
							///<The smallest region that inclued uv will be removed.
	MGPvector<MGFace>& faces///<Result trimmed face(s) will be appended,
			///<If no trimming was performed, no faces will be appended.
)const;

///Compute unit normal vector at uv.
virtual MGUnit_vector unit_normal(const MGPosition& uv) const=0;

///Compute unit normal vector at (u,v).
virtual MGUnit_vector unit_normal(double u,double v) const=0;

protected:

///Obtain coefficient's space dimension.
///This function is used in isect_start etc.
virtual size_t coef_sdim() const=0;

///isect_area_length() returns initial area length for the intersection
///line.
virtual size_t isect_area_length() const=0;

///isect_direction() is used by isect_startPt() to define which constant
///parameter line should be used to compute intersection, and what
///incremental value be used for the parameter.
///Function's return value is direction to get next intersection(with dt).
///When =1: u=const direction, =0: v=const, =-1: cannot get intersection.
virtual int isect_direction(
	const MGFSurface& sf2,	///<Second surface for the intersection.
	size_t m1,		///<id of uvuvS that indicates this surface's parameter
		///<position in uvuvS. (uvuvS(m1), uvuvS(m1+1))=(u,v) of this surface.
	MGPosition& uvuvS,///<start parameter (u,v) pair of this surface and sf2.
	double& du,	///<Incremental value of the parameter kind of kdt will be output.
	double& dv, ///<Right dt will be output according to the function's output =0,1.
	double acuRatio=1.	///<acuracy ratio.
)const;

///isect_direction_with_direction() is used by isect_start() to define which constant
///parameter line should be use to compute intersection, and what
///incremental value be used for the parameter.
///Function's return value isect_direction() is 1 for u=const parameter
///line, and 0 for v=const parameter line.
int isect_direction_with_direction(
	double u, double v,		///<start parameter (u,v) of this surface.
	const MGVector& tangent,///<To indicate which direction isect line
							///<should march toward.
	double& du,				///<Incremental value sign of the parameter kind of
	double& dv		///<isect_direction_with_direction.
)const;

///isect_dt computes incremental values du and dv for the intersection
///computation at parameter position (u,v).
void isect_dt(
	double u, double v, double& du, double& dv,
	double acuRatio=1.	///<acuracy ratio.
) const;

///"isect_inner_dt" is a dedicated function of isect_startPt,
/// comutes adequate incremental parameter value(du,dv) and parameter line kind
///kdt(u=const or v=const).
virtual void isect_inner_dt(
	size_t n,				///<num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,///<intersection point obtained last(of this).
	double& du, double& dv,	///<incremental length from previous to uvnow is input.
				///<New du or dv will be output according to kdt's return value.
	int& kdt,	///<Parameter kind used so far is input, will be output as:
				///<=1:parameter line kind(u=const), =0: v=const,
				///<=-1:should halt computation since incremental value is zero.
	double acuRatio=1.	///<Accurate ratio.
)const;

///isect_dt_coef provides coef of how fine parameter increment should be,
///given num of intersection points computed so far.
double isect_dt_coef(size_t n) const;

///isect_div_id_max is maximum id of array of sect_div defined in
///isect_dt_coef. That is, isect_div_id_max+1 is the length of the array
///sect_div.
size_t isect_div_id_max()const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
virtual MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not.
							///< true:u=const, false:v=const.
	double du, double dv,	///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	size_t incr=0			///<Incremental valuse of B-coef's id.
)const=0;

///isect_start compute one intersection line of two surfaces, this and sf2,
/// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2)
/// and direction(stan). isect_start halts the computation when intersection
/// reached to a boundary of this or sf2, or reached to one of the points
/// in uvuv_list.
///The function's return value is:
/// =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///    =1: End point is a point on a perimeter of one of the surfaces.
///    =3: End point is one of boundary points in uvuv_list.
///    =4: End point is the starting point.
///    =7: isect_start halted the computation since intersection was lost
///     during the computation.
int isect_start(
	const MGPosition& uvuv_startIn, ///<Starting point of the intersection line.
	MGPosition_list& uvuv_list,	///<isect_start will halt when ip reached one of 
		///<the point in uvuv_list. isect_start does not change uvuv_list(actually
		///<uvuv_list is const.) uvuv's space dimension is at least 4,
		///<and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGFSurface& sf2,	///<2nd surface.
	MGSSisect& ssi,			///<Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id,
			///<When the end point of ip was one of the points of uvuv_list,
			///<uvuv_list's iterator of the point will be returned, that is, 
			///<when the function's return value was 3 or 5.
			///<When was not a point of uvuv_list, end() of uvuv_list will be
			///<returned.
	int& m1	///<id that indicates which surface was used as the main surface.
			///<m1=0: this surface, m1=2: sf2.
)const;

/// "isect_start_boundary" is a dedicated function of isect_start.
///"isect_start_boundary" computes one intersection point of two surfaces,
/// this surface's parameter line at uv+dt(according to kdt) and sb2,
/// given previous intersetion point(uv) and incremental value dt.
/// "isect_start_boundary" is used only when no intersection found at dt.
///Function's return value is 0: ip not found.
///                           2: ip found as intersection line end point.
int isect_start_boundary(
	const MGFSurface& sf2,	///<2nd surface b-rep.
	const MGPosition& uvuv_pre,	///<Starting parameter values of ip. 
	int kdt,		///<kdt=true: u=const parameter line,
					///<    else: v=const parameter line of this surface.
	double du, double dv,///<Incremental parameter length.
	size_t lid1,		///<id of parameter of this surface in uvuv_pre or uvuv_now.
	MGPosition& uvuv_now///<New parameter values of ip will be output.
)const;

///Compute the maximum difference between the intersection line and the two surfaces,
///this and sr2. The difference evaluation is done at the data point tau.
double isect_start_dif(
	const MGNDDArray& tau,	///<data points
	const MGLBRep& line,	///<the intersection line of this and sf2.
	const MGFSurface& sf2	///<second surface.
)const;

///isect_start_incr compute one intersection point of two surfaces,
/// this surface's parameter line at uv1+dt(according to kdt) and sf2,
/// given previous intersetion point(uv1,uv2) and incremental value dt.
///Here uv1 is uvuv_pre(lid1, lid1+1), and uv2 is other two values of uvuvpre.
/// isect_start_incr is a dedicated function of isect_start
/// and isect_start_boundary.
///Function's return value is true: if ip found,
///                           false: if ip not found.
int isect_start_incr(
	const MGFSurface& sf2,	///<2nd surface b-rep.
	const MGPosition& uvuv_pre,	///<Starting parameter values of ip. 
	int kdt,		///<kdt=true: u=const parameter line,
					///<    else: v=const parameter line of this surface.
	double du, double dv,///<Incremental parameter length.
	size_t lid1,		///<id of parameter of this surface in uvuv_pre or uvuv_now.
	MGPosition& uvuv_now///<New parameter values of ip will be output.
)const;

///isect_startPt compute an array of parameter value pairs of this surf and sf2
///for one intersection line of the two surfaces,
/// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2)
/// and direction(stan).
///isect_startPt is a dedicated function for isect_start.
///isect_startPt halts the computation when intersection
///reached to a boundary of this or sf2, or reached to one of the points
///in uvuv_list.
///The function's return value is:
///  =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///  =1: End point is a point on a perimeter of one of the surfaces.
///  =3: End point is one of boundary points in uvuv_list.
///  =4: End point is the starting point.
///  =7: isect_start halted the computation since intersection was lost
///     during the computation.
int isect_startPt(
	const MGPosition& uvuv_startIn, ///<Starting point of the intersection line.
	MGPosition_list& uvuv_list,	///<isect_startPt will halt when ip reached one of 
		///<the point in uvuv_list. isect_startPt does not change uvuv_list(actually
		///<uvuv_list is const.) uvuv's space dimension is at least 4,
		///<and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGFSurface& sf2,	///<2nd surface.
	double acuRatio,///<Accurate ratio, should be decreased by multiplyng .2
		///<(or a number less than 1.).
	MGBPointSeq& point,		///<Surface-surface intersection parameter values
		///<will be returned as:point(.,0) and point(.,1) for(u,v) of this surface
		///<point(.,2) and point(.,3) for(u,v) of this surface.
		///<point has the dimension of(.,7).
	MGPosition_list::iterator& uvuv_id,
			///<When the end point of ip was one of the points of uvuv_list,
			///<uvuv_list's iterator of the point will be returned, that is, 
			///<when the function's return value was 3 or 5,
			///<When was not a point of uvuv_list, end() of uvuv_list will be
			///<returned.
	int& m1	///<id that indicates which surface was used as the main surface,
			///<m1=0: this surface, m1=2: sf2.
)const;

///Compute an intersection line of this surface(a surface that is converted to
///1D SBRep surface(sf1d)) and a plane pl.
///sf1D is so converted that pl be x=0. plane. This surface is the original
/// surface and pl is the original plane.
///The function's return value is:
/// =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///    =1: End point is a point on a perimeter of one of the surfaces.
///    =3: End point is one of boundary points in uvuv_list.
///    =4: End point is the starting point.
///    =7: isect_start halted the computation since intersection was lost
///     during the computation.
int isect_startPlane(
	const MGPosition& uvuvS,///<Starting parameter value of the intersection.
	MGPosition_list& uvuv_list,
		///<isect_startPlane will halt when ip reached one of the point in
		///<uvuv_list. isect_startPlane does not change uvuv_list(actually
		///<uvuv_list is const.) uvuv's space dimension is at least 4,
		///<and the first 2 is (u,v) of this and the next 2 is (u,v) of pl. 
		///<When uvuv's space dimension is more than 4, it indicates that the uvuv
		///<is used to input approcimate tangent of the intersection.
	const MGPlane& pl,	///<Plane expression(2nd surface for the intersection).
	MGSSisect& ssi,		///<Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id
			///<When the end point of ip was one of the points of uvuv_list,
			///<uvuv_list's iterator of the point will be returned, that is, 
			///<when the function's return value was 3 or 5.
			///<When was not a point of uvuv_list, end() of uvuv_list will be
			///<returned.
)const;

///isect_startPlanePt compute an array of parameter value pairs of this surf and pl2
///for one intersection line of a surface and a plane,
/// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of pl2)
/// and direction in uvuv_startI(4-6) (optionally).
///isect_startPlanePt is a dedicated function for isect_startPlane.
///isect_startPlanePt halts the computation when intersection
///reached to a boundary of this, or reached to one of the points
///in uvuv_list.
///The function's return value is:
///  =0: Intersection was not obtained.
/// !=0: Intersection was obtained as follows:
///  =1: End point is a point on a perimeter of this surfaces.
///  =3: End point is one of boundary points in uvuv_list.
///  =4: End point is the starting point.
///  =7: isect_startPlanePt halted the computation since intersection was lost
///     during the computation.
int isect_startPlanePt(
	const MGPosition& uvuv_startIn, ///<Starting point of the intersection line.
	MGPosition_list& uvuv_list,	///<isect_startPlanePt will halt when ip reached one of 
		///<the point in uvuv_list. isect_startPlanePt does not change uvuv_list(actually
		///<uvuv_list is const.) uvuv's space dimension is at least 4,
		///<and the first 2 is (u,v) of this and the next 2 is (u,v) of pl2. 
	const MGPlane& pl2,	///<2nd surface(MGPlane).
	double acuRatio,///<Accurate ratio, should be decreased by multiplyng .2
		///<(or a number less than 1.).
	MGBPointSeq& point,		///<Surface-surface intersection parameter values
		///<will be returned as:point(.,0) and point(.,1) for(u,v) of this surface
		///<point(.,2) and point(.,3) for(u,v) of pl2.
		///<point will have the dimension of(.,7).
	MGPosition_list::iterator& uvuv_id
			///<When the end point of ip was one of the points of uvuv_list,
			///<uvuv_list's iterator of the point will be returned, that is, 
			///<when the function's return value was 3.
			///<When was not a point of uvuv_list, end() of uvuv_list will be
			///<returned.
)const;

///Compute the intersection lines of this surface and srf2(both are not planes).
MGSSisect_list isect_with_surf(
	MGPosition_list& uvuv_list,
	///<Let a member of uvuv_list be uvuv. Then, uvuv's space dimension is
	///<at least 4, and the first 2 is (u,v) of this and the next 2 is (u,v) of srf2. 
	///<When uvuv's space dimension is more than 4, it indicates that the uvuv
	///<is used to input approximate tangent of the intersection.
	const MGFSurface& srf2	///<2nd surface for the intersection.
)const;

///Compute the intersection lines of this(is not a plane) surface and a plane pl.
///sd1D that is converted to surf1D() about pl is necessary to input.
MGSSisect_list isect_with_plane(
	MGPosition_list& uvuv_list,
	///<Let a member of uvuv_list be uvuv. Then, uvuv's space dimension is
	///<at least 4, and the first 2 is (u,v) of this and the next 2 is (u,v) of pl. 
	///<When uvuv's space dimension is more than 4, it indicates that the uvuv
	///<is used to input approximate tangent of the intersection.
	const MGPlane& pl,	///<Plane expression(2nd surface for the intersection).
	const MGFSurface& fsrf2	///<Original FSurface before casting to plane.
)const;

///ベクトル投影は、カーブを折れで分割して行い、後で接続する
int projVector(
	const MGCurve& crv,
	MGPvector<MGCurve>& vec_crv_uv,
	MGPvector<MGCurve>& vec_crv,
	const MGVector& vec
) const;

private:

///Test in_range if this is a face, if not, do nothing.
///This is to accelerate the test of in_range in isect_guess().
///See isect_guess(MGCurve).
virtual bool in_range_face(const MGPosition& uv)const{return true;};

///uvuvE_is_a_midpoint() is used by isect_with_surf() to determine if
///uvuvE is a mid point at the intersectio and should be neglected.
bool uvuvE_is_a_midpoint(
	const MGFSurface& f2,
	int m1,		///<id of uvuvE that indicates this surface's parameter
		///<position in uvuvE. (uvuvE(m1), uvuvE(m1+1))=(u,v) of this surface.
	const MGSSisect& ssi, ///<ssi obtained so far.
	const MGPosition& uvuvE///<start parameter (u,v) pair of this surface and sf2.
)const;

///Obtain the projected curve of a curve onto the surface.
///project1span does not divide the curve at the c0 continuity, which is treated by
///project(). See projec().
///The direction of the projection is normal to the surface.
///Output of 'project1span' is curves of space dimension 5, which are (x,y,z,u,v).
///Here (x,y,z) is the world coordinates, and (u,v) is parameter of this surface.
void project1span(
	const MGCurve& crv,	///<The target curve to project.
	MGPvector<MGLBRep>& crv_xyzuv_vector///<Projected curve of(x,y,z,u,v) will be appended.						
)const;

///投影を実行して５次元曲線を取得する
///curveにおけるデータポイントと投影ベクトルから投影点列を作成し
///点列から５次元(xyz,uv)の曲線列を作成して返却する
void prj2OneCurve(
	const MGCurve& curve,	///<target curve to prject.
	std::deque<MGPosition>& ranges,///<start(ranges[0]) and end(ranges[1]) point parameter
							///<of the curve and the face.
							///<On return ranges will be so updated that processed rages[0] to [1]
							///<are pop_front().
	MGLBRep*& crvProjected	///<newed object will be returend when obtained.
							///<When not obtained, null will be returned.
)const;

///crv1上の投影可能なパラメータと投影不可能なパラメータを与えて、
///間にある投影可能な境界パラメータ値を求めて返却する
///チェックポイントの移動量 < (パラメータ範囲*rc_zero())になれば終了。
///Function's return value is MGPosition tuv, where
///tuv[0]=start point parameter of crv1,
///(tuv[1], tuv[2])=the corresponding surface(this surface's) parameter (u,v)
///of the start point of the range.
MGPosition proj2GetParamIter(
	const MGCurve& curve,///<target curve to project.
	const MGPosition& tuv,///<tuv[0] is the curve's parameter value that has a perp point onto the
				///<surface. and (tuv[1], tuv[2]) is the surface's parameter of the perp point.
				///<tuv[0] is eithere ts or te.
	double ts,	///<(ts, te) is curve's parameter range that indicates
	double te	///<between ts and te theremust be the boundary point.
				///<That is one of the following situations occurs:
				///<(1) there are no perp points at ts and there is a perp point at te,
				///<(2) there is a perp point at ts and there are no perp points at te,
)const;

///スタートパラメータを与え投影可能なパラメータ範囲を1つだけ取得する。
///戻り値は	1:範囲が求まった(thisの最後まで)
///			0:範囲が求まらなかった(thisの最後まで)
///			-1:範囲が求まった(thisの途中まで、まだ面直範囲があるかもしれない)
int prj2GetParamRange(
	const MGCurve& curve,
	int start_counter,	///<input start counter of the curve parameter incrementation.
	MGPosition range[2],
		///<range[0][0]=start point parameter of curve,
		///<range[0][1-2]=the corresponding surface(this surface's) parameter (u,v)
		///<of the start point of the range.
		///<Regarding to range[1] : the same.
	int& next_counter	///<Updated curve parameter incremental counter will be output.
)const;

///トレランスに応じた細分したデータポイントを求める
void prjGetDataPoints(
	const MGCurve& curve,
	const MGPosition& tuv0,///<parameter range of curve to get the data point,
	const MGPosition& tuv1,///<From tuv0 to tuv1.
	MGNDDArray& tau		///<data points will be output.
)const;

///ベクトル投影は、crvをvec方向にスイープした面と元の面との交線で求める
int projVectorProc(
	const MGCurve& crv,
	MGPvector<MGLBRep>& vec_crv,	///<５次元曲線が返る
	const MGVector& vec
) const;

///Get world position data at (tuv[1], tuv[2]) and store the coordinates
///and the parameter in bp.
void prj_store_bp(
	int i,
	double u, double v,
	MGBPointSeq& bp
)const;

};
/** @} */ // end of MGObjectRelated group
#endif
