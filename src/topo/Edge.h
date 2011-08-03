/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGEdge_HH_
#define _MGEdge_HH_

#include "mg/Box.h"
#include "mg/TrimmedCurve.h"
#include "topo/CellNB.h"

class MGStraight;
class MGSurface;
class MGPVertex;
class MGBVertex;
class MGLoop;
class MGFace;

//
//Define MGEdge Class.

/** @addtogroup TOPO
 *  @{
 */

///MGEdge is an instance of MGCellNB, represents a boundary element of 2D manifold.
///MGEdge constitues an MGLoop that is a boundary of MGFace.
///MGEdge can be a parameter cell or a binder cell. The coordinates of a parameter cell
///MGEdge is (u,v) surface parameter, and the ones of binder cell MGEdge is
///(x,y,z) of the world.
class MGCLASS MGEdge: public MGCellNB{

public:

///Edgeのスケーリングを行い，Edgeを作成する。
///Scaling of the Edge by a double.
MGDECL friend MGEdge operator* (double s, const MGEdge& e);

///////// Constructor /////////

///void constructor.
MGEdge();

///Copy constructor.
MGEdge(const MGEdge& e, bool copy_boundary=true, bool no_binder=false);

///Fundamental constructor.
///Construct an edge from geometry of manifold dimension 1.
///The constructor takes the ownership of geo and MGPVertex* in boundaries.
MGEdge(
	MGGeometry* geo,
	MGPVertex* boundaries[2],
	MGCellNB* binder
	);

///Make an edge of a boundary that has active start and
///end vertex if the curve is not infinite straight line.
///The second form that input MGCurve* takes the ownership of the crv
///into the MGEdge, must not delete the object and the object must be
///newed one.
MGEdge(const MGCurve& crv);
explicit MGEdge(MGCurve* crv);

///Make an edge of a boundary(MGBoundary1D that has active start and
///end vertex).
///range is the parameter range of crv.
///The second form that input MGCurve* takes the ownership of the crv
///into the MGEdge, must not delete the object and the object must be
///newed one.
MGEdge(const MGCurve& crv, const MGInterval& range);
MGEdge(MGCurve* crv, const MGInterval& range);

///Make an edge with a binder of a boundary
///(MGBoundary1D that has active start and end vertex).
MGEdge(
	const MGSurface&surf,///<Parent surface of which this edge makes a boundary
	const MGCurve& pcrv, ///<Parameter curve of the surface surf.
	const MGInterval& prange,///<param range of pcrv.
	const MGCurve& wcrv	///<World coordinate curve of the surface surf.
						///<wcrv will be trimmed by prange of pcrv.
	);

////////////Destructor////////////
~MGEdge();

///////// operator overload/////////

///Assignment.
///When the leaf object of this and cell2 are not equal, this assignment
///does nothing.
///does not change binder and partner relation,
///does not change parent complex.
MGEdge& operator=(const MGGel& gel2);
MGEdge& operator=(const MGEdge& gel2);

/// Edge に平行移動を行ないオブジェクトを生成する。
///Translation of the Edge
MGEdge operator+ (const MGVector& v) const;

/// Edgeに逆方向の平行移動を行ないオブジェクトを生成する。
///Translation of the Edge
MGEdge operator- (const MGVector& v) const;

///Edgeのスケーリングを行い，Edgeを作成する。
///Scaling of the Edge by a double.
MGEdge operator* (double s) const;

/// 与えられた変換でEdgeの変換を行い，Edgeを作成する。
///Transformation of the Edge by a matrix.
MGEdge operator* (const MGMatrix& mat) const;

/// 与えられた変換によってトランスフォームをおこないEdgeを生成する。
///Transformation of the Edge by a MGTransf.
MGEdge operator* (const MGTransf& tr) const;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGEdge operator/ (double s) const{return (*this)*(1./s);};

///Comparison of two objects.
bool operator==(const MGEdge& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGEdge& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGEdge& gel2)const{return !(gel2==(*this));};

///Object transformation.
MGEdge& operator+=(const MGVector& v);
MGEdge& operator-=(const MGVector& v);
MGEdge& operator*=(double scale);
MGEdge& operator*=(const MGMatrix& mat);
MGEdge& operator*=(const MGTransf& tr);

std::ostream& out(std::ostream& ostrm) const;

/////////Member Function/////////

///Test if active at start or end.
bool active_end() const{ return m_vertex[1]!=0;}
bool active_start() const{ return m_vertex[0]!=0;}

///Get after edge in the loop sequence.
///The aft_edge is the first neighbour edge.
const MGEdge* aft_edge(bool at_end=true, size_t* vertexID=0)const;
MGEdge* aft_edge(bool at_end=true, size_t* vertexID=0);

///Obtain binder edge pointer.
///Null when this does not have binder.
MGEdge* binder_edge() const;

///Obtain the box of the cell.
const MGBox& box() const;

///Obtain the center parameter value of this cell.
MGPosition center_param() const;

///Make a clone of the cell.
///clone() does not copy the binder cell relation.
MGEdge* clone() const;
///clone_without_boundaries() does not copy the binder cell relation.
MGEdge* clone_without_boundaries() const;

///Make a clone of this(this is a binder), and set binder and parameter cell
///relation between the new binder and  the parameter cell e.
MGEdge* clone_binder(const MGCellBase& e) const;

///Compute the continuities between this edge(edge1) and the edge2.
///This edge and edge2 must be parameter edges of each face.
///In distance, tangent, and normal, the following output will be set:
///distance[0-6] as:
///	[0] edge1's curve parameter that has the maximum distance with edge2.
///	[1] edge2's curve parameter that has the maximum distance with edge1.
///  [2] the evaluated maximum distance between edge1 and edge2 at distance[0] and [1]
///	[3] edge1's curve parameter that has the minimum distance with edge2.
///	[4] edge2's curve parameter that has the minimum distance with edge1.
///  [5] the evaluated minimum distance between edge1 and edge2 at distance[3] and [4]
///	[6] mean distance between edge1 and edge2.
///tangent[0-3] as:
///	[0] edge1's curve parameter that has the maximum tangent difference with edge2.
///	[1] edge2's curve parameter that has the maximum tangent difference with edge1.
///  [2] the evaluated maximum tangent difference between edge1 and edge2 at tangent[0] and [1].
///	[3] mean tangent difference between edge1 and edge2.
///normal[0-3] as:
///	[0] edge1's curve parameter that has the maximum normal difference with edge2.
///	[1] edge2's curve parameter that has the maximum normal difference with edge1.
///  [2] the evaluated maximum normal difference between edge1 and edge2 at normal[0] and [1].
///	[3] mean normal difference between edge1 and edge2.
void compute_continuity(
	const MGEdge& edge2,
	double diatance[7],
	double tangent[4],
	double normal[4]
)const;

///Connect this edge  to cell2(is a MGEdge). Both edges are parameter edges of faces.
///This cell is a pcell of a boundary of a higher manifold dimension's cell A,
///and cell2 is also is a pcell of a boundary of another cell B.
///That is, this cell is a part of a boundary of cell A,
///and cell2 is a part of a boundary of cell B.
///If cell A's manifold dimension is n, then this cell's manifold dimension is n-1.
///B's manifold dimension is same as A's. and cell2's manifold dimension is n-1.
void connect(MGCellBase& cell2);
void connect(MGEdge& cell2);

///Connect the start(id1=0) or end(id1=1) of this to the start(id2=0) or
/// the end(id2=1) of e2.
///If both edges of this and e2 are members of a complex, they must be the same.
///e2 must be a newed object, and the owneship is transfered to the system.
void connect_at_id(size_t id1, MGEdge* e2, size_t id2);

///Return curve pointer of this edge.
///Null when this does not have geometry.
///The expression is of parameter space of face.
MGCurve* base_curve();
const MGCurve* base_curve() const;

///Return curve pointer cut by start and end parameter range.
///Output is newed curve object, must be deleted.
///Null when this does not have geometry.
///The expression is of parameter space of face if this is parameter edge of a face.
///curve_limitted() does not return MGTrimmedCurve, returns real curve.
MGCurve* curve_limitted() const;

///Disconnect the start(id=0) or end(id=1) neighbourhood relation.
///disconnect does not free membership of this edge from its parent complex.
void disconnect_at_id(size_t id);

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void draw3DVertex()const;

///Obtain the end point of the edge.
MGPosition end_point()const{return eval(param_e());};

///Evaluate the nderiv's derivative at parameter t.
///Evaluate of the curve's data.
MGVector eval(double t, size_t nderiv=0)const;

///Evaluation of the star curves of the edge at the point t.
///When nderi=0, get a position of the surface at the boundary point t.
///The star curve is SurfCurve(face's surface, edge's curve).
///(The star curve has the same world coordinate with the binder curve's, but
///their direction may be opposite. The star curve has always the same direction
///as the loop.)
MGVector eval_star(
	double t,		///<Parameter value of this parameter edge's curve.
	size_t nderi=0	///<Order of derivative.
)const;

///Test if SurfCurve of the edge has equal direction to binder edge's direction.
///Returned is true if eaual, false if not.
bool equal_direction_to_binder()const;

///Get the star face pointer.
const MGFace* face() const;
MGFace* face();

///Get the 1st partner edge of this edge.
const MGEdge* first_partner() const;

///Free neighbourhood relationship at the end of the edge.
void free_end_neighbourhood();

///Free neighbourhood relation at j-th boundary's i-th pcell of this cell.
///If start, j=0. If end, j=1. i must be always 0, since one boundary has
///only one cell.
void free_neighbourhood(size_t i, size_t j=0);

///Free neighbourhood relationship at the start of the edge.
void free_start_neighbourhood();

///Return Object's type ID (TID)
long identify_type()const;

///Test if this edge's start point(when start=true) and edge2 is connected
///and their directions are the same.
///When start=false, this edge's end point is tested.
bool is_connected_and_same_direction(
	bool start,
	const MGEdge& edge2
)const;

///Test if this is a free edge.
///Free edges are ones that do not have partner edges.
bool is_free()const{ return number_of_partners()==0;};

///Connect this and e2.
///If start==true, start of this edge to end of e2;
///If start==false, end of this edge to start of e2;
///e2 must be a newed object, and the ownership is transfered to the system.
void join(bool start, MGEdge* e2);

///Return parent loop pointer.
const MGLoop* loop() const;
MGLoop* loop();

///Make a binder cell of this parameter cell.
///Returned is the binder pointer generated by new.
///The binder has no geometry, only has binder and parameter cell relationship.
MGCellNB* make_binder() const;

///Make a binder associated with the world curve rep.
///Returned is the binder edge pointer.
///If the parameter edge had already the binder,
///make_binder_with_curve only returns the pointer.
///*** This edge must be a member of a loop that is a boundary of a face.
MGEdge* make_binder_with_curve()const;

///Obtain manifold dimension.
unsigned manifold_dimension() const{return 1;};

///Obtain the i-th member partner edge. This must be a binder edge.
const MGEdge* member_partner_edge(size_t i)const;

///Compute the mid point of this edge.
///Mid point is the point of the paramete mid=(param_s()+param_e())*.5
MGPosition mid_point()const;

///Negate the direction of the cell.
void negate();

///Obtain all the neighbours.
///The neighbours do not contain this cell except when this cell is
///connected to this cell itself(closed cell).
std::vector<const MGCellNB*> neighbours() const;

///Test if the edge is a part of a surface perimeter.
bool on_surface_perimeter() const{return surface_perimeter()>=0;};
bool on_surface_perimeter(const MGFace& f) const{return surface_perimeter(f)>=0;};
bool on_surface_perimeter(const MGSurface& sf) const{return surface_perimeter(sf)>=0;};

///Return parameter space error of the cell.
double parameter_error()const;

///Obtain the parameter of the binder edge's curve that represent
///the same point as sp. sp is a parameter value of this parameter edge.
///Let S() is the star(surface) of this edge, and fp() is the curve of this cell
///which is a boundary of S(). And fb() is the binder curve of this edge.
///Then S(fp(sp))=fb(param_bcell(sp)).
///This is a parameter edge and have the binder, and the parameter sp is a parameter
///of this cell's curve. If this does not have a binder, return -1.
double param_bcell(double tp, const double* guess=0)const;

///This must be a parameter edge.
///Obtain the parameter of this parameter edge's curve that represent the same
///point as the binder edge's paramter tb.
///Let S() is the star(surface) of this edge, and fp() is the curve of this cell
///which is a boundary of S(). And fb() is the binder curve.
///Then S(fp(param_pcell(tb)))=fb(tb).
///This edge must have the binder edge, and the parameter tb is the parameter
///of the binder edge's curve. If this does not have a binder, return -1.
double param_pcell(double tb, const double* guess=0)const;

///Obtain end parameter value of the edge.
double param_e()const;

///Obtain start parameter value of the edge.
double param_s()const;

///Obtain partner edges.
///Partners represent same world's(same cell's parameter) coordinates.
///Parameter edges' partners are parameter edges.
///Binder edges' partners are binder edges.
///The partners do not include this edge except when star cell is
///connected to the star cell itself(closed only by the star cell).
std::vector<const MGEdge*> partner_edges() const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Approximate the parameter edge by a polyline and replace this edge
///expression by the polyline. Polyline approximation is so done that
///the correspoinding binder edge can be appximated by the polyline connecting
///each binder edge's point that corresponds to the each this edge's point.
///(1) This must be a parameter cell edge.
///(2) This edge must be a member of a loop which is a boundary of a face.
///(3) If this edge did not have a binder edge, polygonize generates the binder edge.
///(The tolerance used to generate the binder is MGTolerance::line_zero(),
/// not input error.)
///Input error is tolerance allowed between the polygon and the original curve.
void polygonize(double error);

///Get previous edge in the loop sequence.
///The pre_edge is the first neighbour edge.
const MGEdge* pre_edge(bool at_start=true) const;
MGEdge* pre_edge(bool at_start=true);

///Get parameter range of the edge.
MGInterval range()const;

///Set binder cell edge to this parameter cell.
///This curve's coordinates are of parameter space of a face. And input crv's
///coordinates are world coordinate of the face.
///range is the parameter range of wcrv.
///Parameter range of the wcrv is from start to end of the wcrv when no range
///is specified.
///Function return value is the binder's pointer generated.
MGEdge* set_binder_edge(const MGCurve& wcrv)const;
MGEdge* set_binder_edge(const MGCurve& wcrv, const MGInterval& range)const;

///These forms give the ownership of wcrv to the edge.
///That is, wcrv must be newed one and users must not delete it.
///Others are same as above "set_binder_edge(const MGCurve& wcrv)" form.
///Function return value is the binder's pointer generated.
MGEdge* set_binder_edge(MGCurve* wcrv)const;
MGEdge* set_binder_edge(MGCurve* wcrv, const MGInterval& range)const;

///Set start point(boundary) data.
void set_end(double t);	///Parameter value of the start point.

///Set start point(boundary) data.
void set_start(double t);///Parameter value of the start point.

///Set binder relation to m_vertex[i].
///i is 0 for the start of the edge, and is 1 for the end.
void set_i_th_binder(size_t i, MGBVertex& binder)const;

///Obtain star surface.
///Star cell of this must be a face. If not, return null.
///If does not have star surface, returns null.
const MGSurface* star_surface()const;

///Obtain the end point of the edge.
MGPosition start_point()const{return eval(param_s());};

///Get the perimeter number where this edge is on.
///If this is not on any perimeter, -1 will be returned.
int surface_perimeter() const;
int surface_perimeter(const MGSurface& sf) const;
int surface_perimeter(const MGFace& face) const;

///Trim the edge at parameter t.
///When start=true, trim start, and the result is from t to end.
///When start=false, trim end, and the result is from start to t.
void trim(double t, bool start);

///Get trimmed curve representation of the edge.
MGTrimmedCurve trimmed_curve() const;

///Get the vertex at the start or end.
const MGPVertex* vertex(size_t id)const{return m_vertex[id];};
const MGPVertex* vertex_start()const{return m_vertex[0];};
const MGPVertex* vertex_end()const{return m_vertex[1];};
MGPVertex* vertex(size_t id){return m_vertex[id];};
MGPVertex* vertex_start(){return m_vertex[0];};
MGPVertex* vertex_end(){return m_vertex[1];};

///Return world curve pointer of this edge. That is, curve pointer
///of this edge's binder edge.
///May be null when no binder, or the binder does not have an extent.
MGCurve* world_curve();
const MGCurve* world_curve() const;

std::string whoami()const{return "Edge";};

protected:

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

private:

	MGPVertex* m_vertex[2];	///<[0]:for start, [1] for end point.
							///<When m_vertex[.]=0, the side has no boundary.
							///<MGPVertex must be a newed object.
	mutable MGBox m_box;	///<Box of the edge.
	mutable double m_perror;///<parameter error allowed for the parameter space
							///<of the edge.
	mutable int m_equal_to_binder;
		///<Flag if this curve's direction is equal to the binder or not.
		///<=0:undefined, =1: equal, =-1:opposite direction.

///Transform the boundary binders.
void bn_binder_tr(const MGVector& v);
void bn_binder_tr(double s);
void bn_binder_tr(const MGMatrix& mat);
void bn_binder_tr(const MGTransf& tr);

///Set the box data as null.
void set_box_as_null() const;

///Generate a new MGCellBase pointer by newing the original MGCellBase.
///This is a proprietry routine of MGComplex copy.
///Copy all boundary data, (but does not copy own binder cell relation)
///and register boundary binder association of new and old into cmap.
MGEdge* clone(MGCellMap& cmap) const;

void compute_box() const;

///Compute continuity, given the evaluation interval and the division number.
///This is a parameter edge that has the star face.
void compute_continuity2(
	const MGInterval& span,///<this edge's binder edge's parameter span
	int npoint,				///<division number of this edge's interval sspan
	const MGEdge& edge2,	///<the second parameter edge that has the star face.
	double distance[7],		///<evaluated data will be set in distance, tangent, and normal.
	double tangent[4],		///<See compute_continuity
	double normal[4]
)const;

///Copy boundary data of cell2 into this.
void copy_all_boundaries(const MGCellBase& cellin);

///Copy all boundaries of cell into this, and binders association
///of the boundaries in the cmap.
///Binder cells of cell will be registered in cmap.
void copy_all_boundaries(const MGCellBase& cellin, MGCellMap& cmap);

///Copy m_box data of cell2 into this.
void copy_box(const MGCellBase& cellin) const;

///Copy m_perror data of cell2 into this.
void copy_perror(const MGCellBase& cellin) const;

///Get boundary biders of all the boundaries.
///Binders will be appended to cvec.
void get_all_boundary_binders(std::vector<MGCellNB*>& cvec) const;

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

///Negate the boundary.
void negate_boundary();

friend class MGFace;
friend class MGLoop;

};

/** @} */ // end of TOPO group
#endif
