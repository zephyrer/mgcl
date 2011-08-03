/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFace_HH_
#define _MGFace_HH_

#include <vector>
#include "mg/Default.h"
#include "mg/BPointSeq.h"
#include "mg/Unit_vector.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/FSurface.h"
#include "mg/Pvector.h"

#include "topo/Cell.h"
#include "topo/LEPoint.h"
#include "topo/FOuterCurve.h"

class MGSSisect_list;
class MGCSisect_list;
class MGPosition_list;
class MGLBRep;
class MGCellNB;
class MGLoop;
class MGVPCell_const;
class MGLEPoint;
class MGShell;
class MGFSurface;
class MGHHisect_vector;

//
//Define MGFace Class.

/** @addtogroup TOPO
 *  @{
 */

///MGFace is a trimmed surface.
///MGFace is an instance of MGCell, can be a constituent of MGShell.
///Many useful functions are provided in MGFSurface. See MGFSurface.
class MGCLASS MGFace: public MGCell, public MGFSurface{

public:

///Faceのスケーリングを行い，Faceを作成する。
///Scaling of the Face by a double.
MGDECL friend MGFace operator* (double s, const MGFace& face);

///////// Constructor /////////

///Null face.
MGFace():MGCell(){;}

///Copy constructor.
MGFace(const MGFace& face, bool copy_boundary=true, bool no_binder=false);

///Fundamental constructor.
///Construct a face from geometry of manifold dimension 2
///and the boundaries.
///The constructor takes the ownership of geo and MGBoundary in boundaries.
///boundaries must be loops.
MGFace(
	MGSurface* geo,
	std::vector<MGBoundary*>& boundaries,
	MGCell* binder
);

///Face of whole surface of no boundary.
MGFace(const MGSurface& surf)
:MGCell(surf),m_box_param(surf.param_range())
{;}

///Conversion constructor from MGFSurface to MGFace.
MGFace(const MGFSurface& surf);

///Face of whole surface of no boundary.
///Ownership of surf is transfered to the face.
///(that is surf must be a newed object.)
MGFace(MGSurface* surf)
:MGCell(surf),m_box_param(surf->param_range())
{;}

///Construct a face by copying boundaries(only parameter rep of the boundary)
///from argument boundaries.
///Second form is to input a newed surface. The constructor takes the ownership
///of the surf.
MGFace(const MGSurface& surf,
	const std::vector<MGBoundary*>& boundaries);
MGFace(MGSurface* surf,
	const std::vector<MGBoundary*>& boundaries);

///////// operator overload/////////

///Assignment.
///When the leaf object of this and cell2 are not equal, this assignment
///does nothing.
MGFace& operator=(const MGGel& gel2);
MGFace& operator=(const MGFace& gel2);

/// Faceに平行移動を行ないオブジェクトを生成する。
///Translation of the Face
MGFace operator+ (const MGVector& v) const;

/// Faceに逆方向の平行移動を行ないオブジェクトを生成する。
///Translation of the Face
MGFace operator- (const MGVector& v) const;

///Faceのスケーリングを行い，Faceを作成する。
///Scaling of the Face by a double.
MGFace operator* (double s) const;

/// 与えられた変換でFaceの変換を行い，Faceを作成する。
///Transformation of the Face by a matrix.
MGFace operator* (const MGMatrix& mat) const;

/// 与えられた変換によってトランスフォームをおこないFaceを生成する。
///Transformation of the Face by a MGTransf.
MGFace operator* (const MGTransf& tr) const;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGFace operator/ (double s) const{return (*this)*(1./s);};

///Comparison of two curves.
bool operator<(const MGFace& gel2)const{return is_less_than(gel2);};
bool operator<(const MGGel& gel2)const;

///PD144=MGFace. Output to PD144(Trimmed surface).
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
std::ostream& out(std::ostream& ostrm) const;

/// Output virtual function.
std::ostream& outFS(std::ostream& ostrm) const{return out(ostrm);};

/////////Member function/////////

///Add a new loop to this face as aboundary.
///When the old loops that are outside the nloop will be removed from this.
///nloop can be inner or outer.
void add_boundary(MGLoop* nloop);

///Generate arrow data of the tangent along u and v and the normal
///at the parameter value (u,v) of the surface.
///data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
///data[2], [3] are two bottoms of u-tangent arrowhead.
///data[0], [4], [5], [6] are the points of v-tangent arrow.
///data[0], [7], [8], [9] are the points of v-tangent arrow.
void arrow(double u,double v, MGPosition data[10])const;
void arrow(const MGPosition& uv, MGPosition data[10])const{
	arrow(uv[0],uv[1],data);
}

///Append new one boundary to boundary vectors.
///Returned is the number of boudaries after appending.
///bound must be a newed MGLoop, and the ownership is transfered to this.
///*** append_boundary does not check validity with other loops
///(e.x. already existed loops will be outside the new boudanry bound).
///If the validity check is necessary, use add_boudanry().
size_t append_boundary(MGBoundary* bound);

///Obtain binder face pointer.
///Null when this does not have binder.
MGFace* binder_face() const;

///Obtain all the boundary curves(world coordinates representation)
///of the face.
///That is, all of the outer boundaries and all of the inner boundaries.
MGPvector<MGCurve> face_boundaries()const;

///Return box of the parameter space of the face.
///After trimmed one.
const MGBox& box_param() const;

///Return box of the parameter space of the FSurface.
///After trimmed one.
const MGBox box_param2() const{return box_param();};

///Build a loop of this face, given a closed curve crv on this face. Although crv
///is generally a MGCompositeCurve, this may be not the case. Returned MGLoop is not
///added into this face as a boundary. User must add it after the direction is adjusted.
///That is, the output loop can be an outer or inner loop.
std::auto_ptr<MGLoop> build_loop(
	const MGCurve& crv	///<curve of world coordinates.
		///<Generally this is not on face and always is projectd onto the face.
)const;

///Make a clone of the cell.
///clone(), clone_without_boundaries() does not copy the binder cell relation.
MGFace* clone() const;
MGFace* clone_without_boundaries() const;

///Get the clone of this MGFSurface.
MGFSurface* clone_fsurface()const{return clone();};

///Get the clone of this as a MGFace.
///If this is MGSurface, it is converted to MGFace.
MGFace* clone_as_face()const{return clone();};

///Make a clone of this(this is a binder), and set binder and parameter cell
///relation between the new binder and  the parameter cell f.
MGFace* clone_binder(const MGCellBase& f) const;

///Compute closest point from a point.
///Returned is the parameter value of the face that is closest to point.
MGPosition closest(const MGPosition& point) const;
	
///Compute closest point from a line to the boundary of the MGFSurface.
///Returned is the parameter value of the FSurface that is closest to point.
MGPosition closest_on_boundary(const MGStraight& sl) const;

///////display member function.
void display_arrows()const;
void display_control_polygon()const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const{drawWireFS(span_length,line_density);};

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void draw3DVertex()const;

///Shade the object in world coordinates.
void shade(
	double span_length	///<Line segment span length.
)const;

///Test if directions of parameter curve and world curve of the face boundary
///is equal or not. This function can be used to test the pair of
///the output of outer_boundary() and outer_boundary_param(), or the pair of
///inner_boundary() and inner_boundary_param().
///Return is:
///true if equal direction, false if opposite direction.
bool equal_direction(
	const MGPvector<MGCurve>& wcurves,
		///<output of outer_boundary() or inner_boundary().
	const MGPvector<MGCurve>& pcurves,
		///<output of outer_boundary_param() or inner_boundary_param().
	size_t i
		///<id of the curve in wcurves and pcurves to test the direction.
)const;

///Erase i-th loop.
///void erase_boundary(size_t i);

///Evaluate.
///Input parameter value is not checked if it is in_range() or not.
///Even if it is not in_range(), surface evaluation will be executed.
MGVector eval(
	double u, double v,	///<Face parameter value(u,v)
	 size_t ndu=0///<Order of derivative.
	 , size_t ndv=0
) const;

MGVector eval(
	const MGPosition& uv,	///<Face parameter value(u,v)
	size_t ndu=0, size_t ndv=0///<Order of derivative.
) const;

///Extract all the loops of this face.
void extract_loops(std::vector<const MGLoop*>& loops)const;

///Extract sub  face that is bounded by networks loops.
///Extracted sub face is the smallest closed part of this face bounded by
///the networks that includes the parameter position uv(u,v).
void extract_sub_face(
	const MGPvector<MGLoop>& networks,///<(u,v) representation networks.
	const MGPosition& uv,
	std::auto_ptr<MGFace>& face///<Result extracted face will be output.
)const;

///Return MGFace pointer if this MGGel is an MGFace, else return null.
MGFace* face(){return this;};
const MGFace* face()const{return this;};

///Get the MGFSurface pointer if this is MGSurface or MGFace.
const MGFSurface* fsurface()const{return this;};
MGFSurface* fsurface(){return this;};

///Get inner_aboundary loops included in the input box.
std::vector<const MGLoop*> get_inner_boundary_loops(const MGBox& uvbox) const;

///get face pointer if this is MGFace, else null will be returned.
MGFace* get_face_pointer(){return this;};
const MGFace* get_face_pointer()const{return this;};

///get surface pointer. Null will never be returned if this is valid MGFSurface.
///That is, if this is MGFace, base surface will be returned.
MGSurface* get_surface_pointer(){return surface();};
const MGSurface* get_surface_pointer()const{return surface();};

///Get number of inner boundaries as the output of the function.
size_t get_number_of_boundaries()const{return number_of_boundaries();};

///Test if this and 2nd object has common area about their box(),
///taking error into account.
bool has_commonFS(const MGObject& obj2)const{return has_common(obj2);};

///Test if this face has boundary loops or not in the specified box.
///If this has one, return true.
bool hasLoop(const MGBox& uvbox)const;

///Test if this face has an inactive loop.
///If this has one, return true.
bool hasInactiveLoop()const;

///Test if this face has the outer boundary loop instead of perimeter boundary
///loops. If this has the outer boundary loop and has not perimeter boundary loops,
///return true.
bool hasOuterBoundaryLoop()const;

///Test if this face has perimeter boundary loops or not.
///If this has one, return true.
bool hasPerimeterBoundaryLoop() const;

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the face. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
MGPvector<MGCurve> inner_boundary(size_t i)const;

///Obtain i-th inner_boundary curves(world coordinates representation)
///of the face. Let the output of inner_boundary(i) be wcurves and
///of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
///to pcurves[j] one to one. Number of inner_boundary can be obtained
///by the function number_of_inner_boundary().
MGPvector<MGCurve> inner_boundary_param(size_t i)const;

///Return Object's type ID (TID)
long identify_type()const;

///Test if parameter value (u,v) is in the range of the face parameter.
bool in_range(double u, double v)const;
bool in_range(const MGPosition& uv)const;

///Test if (u,v) is inside the face.
///Function's return value is:
///  0:outside the face.
///  1:unknown.
///  2:inside the face, not on a boundary.
///  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
///  4:(u,v) is on the outer boundary.
///  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int in_range_with_on(const MGPosition& uv)const;

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

///Intersection.
MGCSisect_list isect(const MGCurve& curv) const;
MGSSisect_list isect(const MGFSurface& fsurf) const;
MGSSisect_list isect(const MGFace& fsurf) const;
MGSSisect_list isect(const MGSurface& fsurf) const;
MGHHisect_vector isect(const MGShell& shell2) const;

///Access to i-th element of u knot
double knot_u(size_t i) const{return surface()->knot_u(i);};

///Access to i-th element of v knot
double knot_v(size_t i) const{return surface()->knot_v(i);};

///Returns the u knot vector.
const MGKnotVector& knot_vector_u()const{return surface()->knot_vector_u();};
MGKnotVector& knot_vector_u(){return surface()->knot_vector_u();};

///Returns the v knot vector.
const MGKnotVector& knot_vector_v()const{return surface()->knot_vector_v();};
MGKnotVector& knot_vector_v(){return surface()->knot_vector_v();};

///Obtain i-th boundary loop of the face.
MGLoop* loop(size_t i);
MGLoop* loop(iterator i);
const MGLoop* loop(size_t i)const;
const MGLoop* loop(const_iterator i)const;

///Make a binder cell of this parameter cell.
///Returned is the binder pointer generated by new.
///The binder has no geometry, only has binder and parameter cell relationship.
MGCellNB* make_binder() const;

///This is a newed MGFace or MGSurface object.
///If this is a MGFace, returns this pointer.
///If this is a MGSurface, construct a newed MGFace using this newed MGSurface,
///and returns the MGFace*.
MGFace* make_face(){return this;};

///Make a display list without color of this gel.
///Return is the display list name.
size_t make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const{return make_display_list_to_hilightFS(span_length,line_density);}

///Make outer boundary if not existed.
void make_outer_boundary();

///Get manifold dimension.
unsigned manifold_dimension() const{return 2;};

///Obtain the i-th member partner face.
const MGFace* member_partner_face(size_t i)const;

///Negate the face.
void negate();
void negateFS(){negate();};

///Compute normal vector(not unit) at uv.
MGVector normal(const MGPosition& uv) const
{	return surface()->normal(uv.ref(0), uv.ref(1));};

///Compute normal vector(not unit) at (u,v).
MGVector normal(double u,double v) const
{	return surface()->normal(u, v);};

///Test if no outer boundary except the surface perimeters.
///That is, test if the following two conditions are satisfied:
///         1. no perimeter boundaries.
///         2. no outer boundary.
bool no_outer_boundaries()const;

///Get number of inner boundaries as the output of the function.
size_t number_of_inner_boundaries()const
{	size_t dummy; return number_of_inner_boundaries(dummy);}

///Get number of inner boundary loops.
///Returned i is the id of the first inner boundary loop if inner boundaries
///exist.
size_t number_of_inner_boundaries(size_t& i)const;

///Compute number of active loops.
size_t number_of_loops()const;

///Get number of perimeter boundary loop.
size_t number_of_perimeter_boundaries()const;

///Return MGObject pointer if this MGGel is an MGObject, else return null.
MGObject* object_pointer(){return object();};
const MGObject* object_pointer()const{return object();};

///Offset.
///distance is plus value if the direction is toward normal vector of the
///face. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
MGFace offset(double distance, int& error)const;

///Offset.
///distance is plus value if the direction is toward normal vector of the
///face. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
int offset(double distance, MGPvector<MGFace>& vecOfsFace)const;

///Offset.
///distance is plus value if the direction is toward normal vector of the
///FSurface. Minus if against the normal vector.
///エラーコード 0:成功 -1:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
int offset_fs(double distance, MGPvector<MGFSurface>& vecOfsFSurface)const;

///Test if a point P is on the face.
///Returned is true if the point P is on the face.
///false(0) if P was not on the face.
bool on(const MGPosition& P,
		MGPosition& uv	///<Parameter value of the face is returrned.
						///<Even if P is not on the face, nearest point
						///<parameter value will be returned.
)const;

///Test if input (u,v) is parameter value on a perimeter of the base surface.
///If u or v is on a perimeter, they will be updated to the perimeter value.
bool on_a_perimeter(
	double& u, double& v,		///<Surface parameter (u,v)
	size_t& perim_num	///<if function returns true,
						///<the perimete number is output.
)const;

///Obtain outer_boundary curves(world coordinates representation) of the face.
///Let the output of outer_boundary() be wcurves and of outer_boundary_param()
///be pcurves, then wcurves[i] corresponds to pcurves[i] one to one.
MGPvector<MGCurve> outer_boundary()const;

///Obtain boundary curves(parameter space representation) of the face.
///Let the output of boundary() be wcurves and of boundary_parameter()
///be pcurves, then wcurves[i] corresponds to  pcurves[i] one to one.
MGPvector<MGCurve> outer_boundary_param()const;

///Obtain parameter value of the face whose world coordinates are P.
MGPosition param(const MGPosition& P)const;

///Obtain parameter curves.
///In the case of surface, parameter curve is only one. However, in the case
///of face,  number of parameter curves are more than one.
MGPvector<MGCurve> parameter_curves(
	int is_u,		///<True(!=0) if x is u-value.(i.e. obtain u=const line)
	double x	///<parameter value. u or v-value accordint to is_u.
)const;

/// パラメータ範囲を返す。
///Return parameter range.
MGBox param_range() const{return box_param();};

/// Return ending parameter value.
double param_e_u()const;
double param_e_v()const;

/// Return starting parameter value of the base surface.
double param_s_u()const;
double param_s_v()const;

///Obtain parent shell that this face belongs to.
MGShell* parent_shell();
const MGShell* parent_shell()const;

///Obtain perimeter boundadary loop's curve representation.
///Returned are curves of perimeter boundaries, do not contain perimeter
///of the surface.
MGPvector<MGCurve> PBloop_curves() const;

///Return the foot of the perpendicular straight line from P.
///Computation is done from the guess parameter value.
///Function's return value is whether point is obtained(true) or not(false).
bool perp_guess(
	const MGPosition& P,		///<Point
	const MGPosition& uvguess,	///< guess parameter value of the shell
	MGPosition& uv				///< Parameter value will be returned.
) const;

///Compute perpendicular points of a curve and the face,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
bool perp_guess(
	const MGCurve& curve,	///<curve.
	const MGPosition& uvguess,	///<Guess parameter value of the face.
	double tguess,			///<Guess parameter value of the curve.
	MGPosition& uv,			///<perpendicular point's parameter values of the shell
	double& t				///<will be output.
) const;

///指定点から最も近い、垂線の足とパラメータ値を返す。
///Return the foot of the perpendicular straight line from p that is 
///nearest to point p.
/// Function's return value is whether point is obtained(1) or not(0)
int perp_point (
	const MGPosition& p,		///< 指定点(point)
	MGPosition& uv,		///<Parameter value of the surface will be returned.
	const MGPosition* uvguess=0	///< guess parameter value of surface
)const;

///Compute perpendicular points on the face from a point P((x,y,z)).
///MGPosition uv in the MGPosition_list is:
///uv(0): u parameter, and uv(1): v parameter of the face.
///Generally number of uv are more than one.
MGPosition_list perps(const MGPosition& P) const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Round the input parameter (u,v) of the face to the nearest point of
///the face parameter range.
MGPosition range(const MGPosition& uv) const;

///Remove inactive loops from this face.
void remove_inactive_loops();

///Get the shell pointer if this belongs to a shell.
MGShell* shell();
const MGShell* shell()const;

///Shrink the base surface of this face to the part limitted by the parameter range of uvbx.
///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
void shrink_base_surface_to_knot(
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

///Sort boundary occurreces in m_boundaries.
///Sorting is done according to operator< of MGBoundary.
///parameter space box will be set.
void sort_boundaries();

///split this fsurface at the parameter param.
void split(
	double param,///<parameter value of this fsurface. if is_u is true, param is u-value,
				///<else v-value.
	bool is_u,	///<indicates if param is u or v of the surface parameter (u,v).
	MGPvector<MGFSurface>& surfaces///<splitted surfaces will be output.
)const;

///Split the face giving networks loops. Splitting is done by finding the smallest
///closed areas out of networks.
void split(
	const MGPvector<MGLoop>& networks,
	MGPvector<MGFace>& faces///<Result trimmed face(s) will be appended.
)const;

///Get surface pointer.
MGSurface* surface();
const MGSurface* surface() const;

///Trim the face by the projection of a curve along a vector direction.
///If mgNULL_VEC is specified as direction, surface normal projection
///will be employed.
///crv has a direction. That is, when this face is divided into two faces,
///left part face of the crv is the face selected.
///
///When the projected curve on the face is connected to already existed
///boundary curve, no new boundary is generated, and inserted(connected)
///to the old boundary. However new projection is floating, that is, not
///connected to any old boundaries, or this boundary is the first boundary,
///new boundary is generated.
///Function's return value is error code:
///0= normal return
///(not error, this includes the case of inactive loop generation)
///2= tried to generate outer boundary loop inside perimeter boudary.
///3= tried to generate inner boundary loop that incudes active loop inside.
///4= tried to generate perimeter boudary loop that inactivates perimeter
///   boundary loops existed.
///5= tried to generate a loop outside the face.
int trim_projection(
	const MGCurve& crv,		///<curve of world coordinates.
		///<Generally this is not on face and always is projectd to the face.
	const MGVector& direction=mgNULL_VEC///<Projection directin vector.
);

///Trim the face giving parameter curve and world curve of a curve on the face.
///crv has a direction. That is, when this face is divided into two faces,
///left part face of the crv is the face selected.
///
///When the pcrv is connected to already existed boundary curve,
///no new boundary is generated, and inserted(connected) to the old boundary.
///However, new projection is floating, that is, not
///connected to any old boundaries, or this boundary is the first boundary,
///new boundary is generated.
///Function's return value is error code:
///0= normal return
///(not error, this includes the case of inactive loop generation).
///1= input pcrv includes a part that is outside surface parameter range.
///2= tried to generate outer boundary loop inside perimeter boudary.
///3= tried to generate inner boundary loop that incudes active loop inside.
///4= tried to generate perimeter boudary loop that inactivates perimeter
///   boundary loops existed.
///5= tried to generate a loop outside the face.
int trim(
	const MGCurve& pcrv	///<parameter(u,v) space curve of the face.
);

///Trim the face giving a loop new_loop that does not have the parent face.
///new_loop mus be parrameter representaion of this face and
///must not have intersections with the loops of this face except the end points
///of new_loop.
///
///Function's return value is error code:
///0= normal return
///(not error, this includes the case of inactive loop generation).
///2= tried to generate outer boundary loop inside perimeter boudary.
///3= tried to generate inner boundary loop that incudes active loop inside.
///4= tried to generate perimeter boudary loop that inactivates perimeter
///   boundary loops existed.
///5= tried to generate a loop outside the face.
int trim(const MGLoop& new_loop_in);

///Trim the face giving networks loops. Trimming is done by removing the smallest
///closed area out of networks that includes the parameter position uv(u,v).
void trim(
	const MGPvector<MGLoop>& networks,
	const MGPosition& uv,
	MGPvector<MGFace>& faces///<Result trimmed face(s) will be appended.
)const;

///Compute unit normal vector at uv.
MGUnit_vector unit_normal(const MGPosition& uv) const
{	return surface()->unit_normal(uv.ref(0), uv.ref(1));};

///Compute unit normal vector at (u,v).
MGUnit_vector unit_normal(double u,double v) const
{	return surface()->unit_normal(u, v);};

std::string whoami()const{return "Face";};

protected:

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

private:

	mutable MGBox m_box_param;///<Box of parameter space of the face.

///Generate a new MGCellBase pointer by newing the original MGCellBase.
///This is a proprietry routine of MGComplex copy.
///Copy all boundary data, (but does not copy own binder cell relation)
///and register boundary binder association of new and old into cmap.
MGFace* clone(MGCellMap& cmap)const;

///Test in_range if this is a face, if not, do nothing.
///This is to accelerate the test of in_range in isect_guess().
///See isect_guess(MGCurve).
bool in_range_face(const MGPosition& uv)const{return in_range(uv);};

///Obtain coefficient's space dimension.
///This function is used in isect_start etc.
size_t coef_sdim() const{return surface()->coef_sdim();};

///compute box of the cell in m_box.
///Currently this does not compute corrct box, compute m_extent box.
void compute_box() const;

///Compute parameter range box.
void compute_box_param() const;

///set box as null(to set the box as initial)
void set_box_as_null()const;

///Test if (u,v) is inside inner boundary. inside means not on the
///boundary and not included inside the face.
///If true is returned, the id of m_boundaies is returned.
///Function's return value is:
///  0:outside of all the inner loops(not on the loop)
///  1:unknown
///  2:inside an inner loop(not on the loop), and the loop id is returned in id.
/// otherwise:on the loop(size_t(MGEdge* of parameter edge))will be returned, 
///			 and the loop id is returned in id.
size_t inside_inner_boundary(const MGPosition& uv, size_t& id)const;

///Test if (u,v) is inside the outer boundary.
///Inside the outer boundary means that inside outer_boudary_param() or not.
///***Caution***
///(1)This must not be used for faces that do not have perimeter or outer boundary
///loop.
///(2)inside_outer_boundary does not check about the inner loops.
///
///Function's return value is:
///  0:outside the outer boundary(not on a loop)
///  1:unknown
///  2:inside the outer boundary(not on a loop)
/// otherwise:on the outer boundary loop
size_t inside_outer_boundary(
	const MGPosition& uv
)const;

///Intersection with MGFSurface.
MGSSisect_list intersect(const MGFSurface& face2)const;

///Compute intersection points of an inner parameter line of this face and f2.
///The intersection point is used to compute surface to surface intersection lines.
///Function's return value is at most one intersection point in uvuv_list.
///One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
///this surface and (u2,v2) is a parameter of surf.
MGPosition_list intersectInner(
	const MGFace& f2		///<The second face.
) const;

///Compute intersection points of an inner parameter line of this face and sf2.
///The intersection point is used to compute surface to surface intersection lines.
///Function's return value is at most one intersection point in uvuv_list.
///One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
///this surface and (u2,v2) is a parameter of surf.
MGPosition_list intersectInner(
	const MGSurface& sf2		///<The second surface.
) const;

///isect_area_length() returns initial area length for the intersection
///line.
size_t isect_area_length() const{return surface()->isect_area_length();};

///Compute intersection points of this face's boundary(outer and inners) with
///face2. If intersection points are found and the boundary is a loop,
///the point's edge pointer(of this) will be stored in a member uvuv of uvuvs.
///uvuv[7] is the edge pointer. If the boundary is not a loop(that is, a perimeter of
///Surfaces), uvuv.sdim()==7 and an edge pointer is not returned.
///When uvuv.sdim()==8, the edge pointer of uvuv[7] is accessed through union mgEdgeP.
///uvuvs[i] is i-th intersection points.
size_t isect_boundary(
	const MGFSurface& face2,
	MGPosition_list& uvuvs,
	///<id1 and id2 are the ids of uvuv where this face's and f2's parameters
	///<are to be stored in a member of uvuvs.
	///<This face's (u,v) is stored in uvuv(id1) and (id1+1).
	///<f2's (u,v) is stored in uvuv(id2) and (id2+1).
	///<id2=0 if id1=2, and id2=2 if id1=0.
	size_t id1=0
)const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not.
							///< true:u=const, false:v=const.
	double du, double dv,	///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	size_t incr=0			///<Incremental valuse of B-coef's id.
)const;

///Compute intersection points between the boundary of iid-th inner boundary
///of this face and face2 to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t isect_incurves(
	const MGFSurface& face2,
	size_t iid,	///<Inner loop id of this face(from 0)
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;

///Compute intersection points of outer boundary curves of this face 
///with face2 to compute intersections.
///Function's return value is the number of ip's obtained(appended)
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t isect_outcurves(
	const MGFSurface& face2,
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 7,
		///<and the last three elements are the ip direction vector.
	size_t id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;

///Compute intersection points between loop lp of this face and face2
///to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t isect_lpcurves(
	const MGFSurface& face2,		///<srf!=null, face2=null.
	const MGLoop& lp,				///<Loop id of this face.
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
)const;

///Compute intersection points between loop lpid of this face and face2
///to compute intersections of face with face2.
///Function's return value is the number of ip's obtained before appending
///into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t isect_lpcurves(
	const MGFSurface& face2,
	size_t lpid,				///<Loop id of this face.
	MGPosition_list& uvuv_list,	///<intersection points will be appended.
		///<One member in the list is of sdim 8,
		///<(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1			///<id of uvuv(a member of uvuv_list).
		///<uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		///<id2=0 if id1=2, and id2=2 if id1=0.
) const;

///Compute intersection lines, given end points of the i.l.
///isectEP does not execute has_common() with srf2. Users of isectEP are recommended
///to execute has_common().
MGSSisect_list isectEP(
	MGPosition_list& uvuv_list,	///<End points list of the intersection.
		///<On return, uvuv_list.size() will be 0.
	const MGFSurface& fsrf2		///<The second surface.
) const;

///Obtain parameter u(kcod=1), or v(kcod=0) of the intersection point of
///v=x(const) line(kcod=1), or u=x(const) line (kcod=0) with
///all the face boundaries.
std::vector<double> isect1D_with_boundaries(
	double x,							///<coordinate value of kcod.
	size_t kcod					///<Coordinate kind, =0:u, =1: v.
)const;

///Make this cell's binder cell's extent expression.
///Returned is a MGGeometry pointer generated by new.
///When this cell does not have star cell, null pointer will be returned.
///make_binder_extent() only makes the expression, and does nothing to
///the topology structure.
MGGeometry* make_binder_extent() const;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
void make_extent() const;

///Obtain outer boundary curve expression as the combination of
///loop pointers and perimeter id's.
std::vector<MGFOuterCurve> outer_curve() const;

///Obtain perimeter i's parameter range.
///Let rvec be std::vector<MGInterval> of the fucntion's output,
///then rvec[j] is j-th parameter range of i-th perimeter. 
std::vector<MGInterval> perimeter_param_range(int i) const;

///Dedicated function of range.
///Will check if point (u,v) is inside inner boundary or not.
///If inside an inner boundary, obtain the closest point to the boundary.
MGPosition range_check_inner_boundary(const MGPosition& uv) const;

///Remove parameter uv from uvs that is outside face parameter range.
void remove_outside_param(MGPosition_list& uvs)const;

///Proprietry routine for make_outer_boundary(), will append newed straight line
///edge to the end of lp if param (uv1 and uv2) are far away enough compared with error
///err.
void sl_edge_append(
	MGLoop*& lp,			///<lp that the edge should be append.
							///<If lp is null, new lp will be generated.
	size_t id,				///<perimeter num of the staight line edge.
	const MGPosition& uv2,	///<face parameter of the end point of sl.
	const MGPosition& uv1,	///<face parameter of the start point of sl.
	double err_sqr		///<Error square allowed to regard as same points.
);

///ボックス枠に囲まれる交点を持つUV曲線を生成する
void getTrimCrv(
	double uerror, double verror,///<u and v parameter error.
	const MGBox& box,	///<parameter (u,v) box of the surface.
	MGPvector<MGCurve>& vecCrv	///<paremter curves will be output
) const;

friend class MGShell;
friend class MGSurface;

};

///@cond

///mgEdgeP is used in output of  isect_lpcurves() to store Edge pointer in
///uvuv_list.
union MGCLASS mgEdgeP{const MGEdge* pointer; double doubleV;};

///@endcond

///Obtain the closest point from point uv to vector of curves.
///MGClosest_to_curves does not change wc_zero, and so calling program of
///MGClosest_to_curves should change it if necessary.
MGDECL MGPosition MGClosest_to_curves(
	const MGPosition& uv,				///<Point.
	const MGPvector<MGCurve>& curves	///<vector of curves.
);

/** @} */ // end of TOPO group
#endif
