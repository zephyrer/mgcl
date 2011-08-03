/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPVertex_HH_
#define _MGPVertex_HH_

#include "mg/Default.h"
#include "topo/CellBase.h"

class MGBox;
class MGOfstream;
class MGIfstream;
class MGCellNB;
class MGBVertex;
class MGEdge;

/** @addtogroup TOPO
 *  @{
 */

//
//Define MGPVertex Class.

///MGPVertex is a parameter cell of the manifold dimension 0.
///MGPVertex is a boundary of an Edge(start or end) to hold edge's parameter data.
///Logically PVertex has characters of both a parameter cell and a boundary.
///MGPVertex cannot be a binder cell, and so does not have partner cells.
///MGPvertex can be only partners of a binder.
///This is the reason that MGPVertex is a subclass of MGCellBase.
///Since MGPVertex's manifold dimension is 0, MGPVertex does not have
///boundaries.
class MGCLASS MGPVertex: public MGCellBase{

public:

///////// Constructor /////////

///Void constructor.
MGPVertex():m_edge(0){;};

///Copy constructor.
///Parent edge will be cleared.
MGPVertex(const MGPVertex& v):m_edge(0), m_t(v.m_t){;};

///Copy constructor with parent edge.
///e is set as the parent edge.
MGPVertex(const MGPVertex& v, MGEdge* e):m_edge(e), m_t(v.m_t){;};

///Fundamental constructor.
///Construct from the parameter value t of edge.
///e is set as the parent edge.
explicit MGPVertex(double t, MGEdge* e=0):m_edge(e), m_t(t){;};

///////// operator overload/////////

///Assignment.
///does not change binder and partner relation.
///edge pointer will be cleared.
MGPVertex& operator=(const MGGel& gel2);
MGPVertex& operator=(const MGPVertex& gel2);

///Object transformation.
MGPVertex& operator+=(const MGVector& v);
MGPVertex& operator-=(const MGVector& v);
MGPVertex& operator*=(double scale);
MGPVertex& operator*=(const MGMatrix& mat);
MGPVertex& operator*=(const MGTransf& tr);

///Comparison of two objects.
bool operator<(const MGPVertex& gel2)const;
bool operator<(const MGGel& gel2)const;

/////////Member Function/////////

///Get binder.
MGBVertex* binder_vertex()const;

///Obtain the box into which the topology is included.
const MGBox& box() const{return mgNULL_BOX;};

///Make a clone of the cell.
///clone(), clone_without_boundaries() does not copy the binder cell relation.
MGPVertex* clone() const;
MGCellBase* clone_without_boundaries() const{
	return static_cast<MGCellBase*>(clone());};

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	double span_length,	///<Line segment span length.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void draw3DVertex(
	) const;

///Return Object's type ID (TID)
long identify_type()const;

///Ask if this is binder cell.
bool is_bcell() const{return false;};

///Test if this is the start vertex or end verstex on the edge.
bool is_start_vertex()const;

///Make a binder cell of this parameter cell.
///This is a parameter cell and the binder will be newed.
///Returned is the binder pointer generated.
///The binder has no geometry, only has binder and parameter cell relationship.
MGCellNB* make_binder() const;

///Get manifold dimension.
unsigned manifold_dimension() const{return 0;};

///Negate the direction of the cell.
void negate(){;};

/// Output function.
std::ostream& out(std::ostream&) const;

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
///This will be never invoked.
MGPosition pick_closest(const MGStraight& sl)const{return MGDefault::null_position();};

///Set the parameter value.
void set_t(double t){m_t=t;};

///Set the edge pointer.
void set_edge(MGEdge* e){m_edge=e;};

///Obtain star cells.
const MGCellNB* star() const;
MGCellNB* star();

///Return the edge pointer.
const MGEdge* edge() const{return m_edge;};

///Return the edge pointer.
MGEdge* edge(){return m_edge;};

///Return the parameter value.
double t() const{return m_t;};
std::string whoami()const{return "PVertex";};

protected:

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

private:

	MGEdge* m_edge;	///<The edge whose boundary is this PVertex.
	double m_t;		///<Parameter value of the curve of an edge.

///Copy boundary data of cell2 into this.
void copy_all_boundaries(const MGCellBase& cell2){;};

///Copy all boundaries of cell into this, and binders association
///of the boundaries in the cmap.
///Binder cells of cell will be registered in cmap.
void copy_all_boundaries(const MGCellBase& cell2, MGCellMap& cmap){;};
void copy_box(const MGCellBase& cb)const{;};
void copy_perror(const MGCellBase& cb)const{;};

///Generate a new MGCellBase pointer by newing the original MGCellBase.
///This is a proprietry routine of MGComplex copy.
///Copy all boundary data, (but does not copy own binder cell relation)
///and register boundary binder association data into cmap.
///MGPVertex does not have boundary. So this function is the same as
///clone().
MGPVertex* clone(MGCellMap& cmap) const{
	return clone();
};

///Make the binder cell's extent expression of this parameter cell.
///Returned is a MGGeometry pointer generated by new.
///When this cell does not have star cell, null pointer will be returned.
MGGeometry* make_binder_extent() const;

};

/** @} */ // end of TOPO group
#endif
