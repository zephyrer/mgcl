/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGBVertex_HH_
#define _MGBVertex_HH_

#include "mg/Default.h"
#include "mg/Position.h"
#include "topo/CellNB.h"

class MGBox;
class MGOfstream;
class MGIfstream;
class MGPVertex;

//
//Define MGBVertex Class.

/** @addtogroup TOPO
 *  @{
 */

///MGBVertex is 0 manifold dimension binder cell, is an point.
///MGBVertex is a binder cell of MGPVertex, and has manifold dimension 0.
///MGBVertex is not used as a parameter cell(boundary of an edge).
///Since MGBVertex's manifold dimension is 0, MGBVertex does not have
///boundaries.
class MGCLASS MGBVertex:public MGCellNB{

public:
///////// Constructor /////////

///Void constructor.
MGBVertex(){;};

///Copy constructor.
///MGBVertex(const MGBVertex&v);

///Fundamental constructor.
///Construct a BVertex from geometry of manifold dimension 0
///(MGPoint*, may be null).
///The constructor takes the ownership of geo.
MGBVertex(MGGeometry* geo);

///Construct from MGPosition data.
MGBVertex(const MGPosition& V);

~MGBVertex();

///////// operator overload/////////

///Assignment.
///When the leaf object of this and cell2 are not equal, this assignment
///does nothing.
MGBVertex& operator=(const MGGel& gel2);
MGBVertex& operator=(const MGBVertex& gel2);

///Comparison of two objects.
bool operator<(const MGBVertex& gel2)const;
bool operator<(const MGGel& gel2)const;

/////////Member Function/////////

///Obtain the box into which the topology is included.
const MGBox& box() const{return mgNULL_BOX;};

///Obtain the center parameter value of this cell.
MGPosition center_param() const{return MGPosition();};

///Make a clone of the cell.
///clone(), clone_without_boundaries() does not copy the binder cell relation.
MGBVertex* clone() const{return new MGBVertex(*this);};
MGBVertex* clone_without_boundaries() const{return new MGBVertex(*this);};

///Make a clone of this(this is a binder), and set binder and member
///partner relation between the new binder and  the cell c.
MGBVertex* clone_binder(const MGCellBase& c) const;

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

///Free neighbourhood relation at j-th boundary's i-th pcell of this cell.
void free_neighbourhood(size_t i, size_t j=0){;};

///Return Object's type ID (TID)
long identify_type()const;

///Make a binder cell of this cell.
///Returned is the binder pointer newed.
///The binder has no geometry, only has binder and partner member relationship.
MGCellNB* make_binder() const;

///Get manifold dimension.
unsigned manifold_dimension() const{return 0;};

///Obtain the i-th member partner PVertex.
const MGPVertex* member_partner_vertex(size_t i)const;

///Obtain all the neighbours.
///The neighbours do not contain this cell except when this cell is
///connected to this cell itself(closed cell).
///A vertex has no boundaries, and has no neighbours.
std::vector<const MGCellNB*> neighbours() const{
	return std::vector<const MGCellNB*>();
}

/// Output function.
std::ostream& out(std::ostream&) const;

///Return parameter space error of the cell.
double parameter_error()const{return 0.;};

///Compute the parameter value of the closest point from the straight to
///this object.
///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
///This will be never invoked.
MGPosition pick_closest(const MGStraight& sl)const;

///Return extent data(i.e. MGPoint*), which may be null if this has no extent.
const MGPoint* point() const;

///Get the position data of this vertex.
///Returns MGPoint data if this has extent. Otherwise obtain from partner's
///star edge by evaluating the edge's data.
MGPosition position() const;

protected:

///Transform the boundary binders.
///Since a vertex has no boundary, no process is done.
void bn_binder_tr(const MGVector& v){;};
void bn_binder_tr(double s){;};
void bn_binder_tr(const MGMatrix& mat){;};
void bn_binder_tr(const MGTransf& tr){;};

///Set the box data as null.
void set_box_as_null() const{;};

///Generate a new MGCellBase pointer by newing the original MGCellBase.
///This is a proprietry routine of MGComplex copy.
///Copy all boundary data, (but does not copy own binder cell relation)
///and register boundary binder association of new and old into cmap.
MGBVertex* clone(MGCellMap& cmap) const{return (MGBVertex*)MGCellNB::clone(cmap);};

///compute box of the cell in m_box.
void compute_box() const{;};

///Copy boundary data of cell2 into this.
void copy_all_boundaries(const MGCellBase& cell2){;};

///Copy all boundaries of cell into this, and binders association
///of the boundaries in the cmap.
///Binder cells of cell will be registered in cmap.
void copy_all_boundaries(const MGCellBase& cell2, MGCellMap& cmap){;};

///Copy m_box data of cell2 into this.
void copy_box(const MGCellBase& cell2) const{;};

///Copy m_perror data of cell2 into this.
void copy_perror(const MGCellBase& cell2) const{;};

///Get boundary biders of all the boundaries.
///Binders will be appended to cvec.
void get_all_boundary_binders(std::vector<MGCellNB*>& cvec) const{;};

///Negate the boundary.
void negate_boundary(){;};

std::string whoami()const{return "BVertex";};

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

private:
	
///Make the binder cell's extent expression of this parameter cell.
///Returned is a MGGeometry pointer generated by new.
///When this cell does not have star cell, null pointer will be returned.
MGGeometry* make_binder_extent() const;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
void make_extent() const;

};

/** @} */ // end of TOPO group
#endif
