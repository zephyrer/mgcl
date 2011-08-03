/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCell_HH_
#define _MGCell_HH_

#include <vector>
#include "mg/Box.h"
#include "mg/Position.h"
#include "topo/CellNB.h"
class MGComplex;
class MGBoundary;
class MGGeometry;
class MGCellMap;

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGBoundary*>;
#pragma warning( pop )
#endif

//
//Define MGCell Class.

/** @addtogroup TOPO
 *  @{
 */

///MGCell is a general cell that has bound.
///MGCell's additional data to CellNB are Boundary, box, and perror.
///There are two types of cells. One is parameter cell(pcell) and
///the other is binder cell(bcell). They are exclusive, that is, if
///a cell is a parameter cell, the cell cannot be binder cell and
///vice versa.
///Parameter cell is a constituent of a complex.
///Binder cell is a binder of parameter cells. Plural cells are connected
///through a binder.
///MGCell is an abstrct class.
class MGCLASS MGCell:public MGCellNB{

public:

typedef std::vector<MGBoundary*>::iterator boundaryItr;
typedef std::vector<MGBoundary*>::const_iterator const_boundaryItr;
typedef std::vector<MGBoundary*>::reverse_iterator boundaryRItr;
typedef std::vector<MGBoundary*>::const_reverse_iterator const_boundaryRItr;
typedef boundaryItr iterator;
typedef const_boundaryItr const_iterator;
typedef boundaryRItr reverse_iterator;
typedef const_boundaryRItr const_reverse_iterator;

///////Constructor/////////

///Void constructor. Constructor of pcell.
MGCell();

///Copy constructor. Result cell is not a member of any complex.
///Binders and boundaries of cell will not be copied.
///Copy of boudaries can be done by copy_all_boundaries().
MGCell(const MGCell& cell);

///MGCell of whole geometry(no boundary), under parent.
///Constructor of pcell.
///The second form that input MGGeometry* takes the ownership of the geo
///into the MGCell, must not delete the object and the object must be
///newed one.
MGCell(const MGGeometry& geo);
explicit MGCell(MGGeometry* geo);

///Construct a parameter cell from all the necessary data,
///geo, vector of boundaries, and the binder. The newly constructed parameter
///cell will be a partner member of the binder 'binder'.
///Constructor takes the ownership of goe and MGBoundary in boundaries.
MGCell(
	MGGeometry* geo,
	std::vector<MGBoundary*>& boundaries,
	MGCell* binder
);

///Parameter Cell with boundaries.
///Only pcells(parameter representation) in boundaries are copied.
///Binders(world coordinate representation) in boundaries are discarded.
MGCell(const MGGeometry& geo, 
		const std::vector<MGBoundary*>& boundaries);
MGCell(MGGeometry* geo, 
		const std::vector<MGBoundary*>& boundaries);

////////////Virtual Destructor////////////
virtual ~MGCell();

///Assignment.
///does not change binder and partner relation,
///does not change parent complex.
virtual MGCell& operator=(const MGCell& gel2);

///comparison
virtual bool operator<(const MGCell& gel2)const;

/////////////Member Function///////////////

///Append new one boundary to boundary vectors.
///Returned is the number of boudaries after appending.
virtual size_t append_boundary(MGBoundary* bound);

///Obtain i-th boundary pointer.
MGBoundary* boundary(size_t i) {return m_boundaries[i];};
const MGBoundary* boundary(size_t i) const {return m_boundaries[i];};

///Obtain boundaries of this cell.
const std::vector<MGBoundary*>& boundaries() const{return m_boundaries;};

///Obtain i-th boundary's j-th pcell direction
///(direction of boundary measured by this cell's
///coordinate along the boundary).
///The direction is represented by the center of the boundary.
MGVector boundary_direction(size_t i, size_t j) const;

///Obtain iterator of m_boundaries.
const_boundaryItr boundaryIterator(const MGBoundary* bnd) const;
boundaryItr boundaryIterator(MGBoundary* bnd);

///Obtain the box of the cell.
const MGBox& box() const;

///Obtain the center parameter value of this cell.
MGPosition center_param() const;

///Make a clone of the cell.
///clone(), clone_without_boundaries() does not copy the binder cell relation.
virtual MGCell* clone() const=0;
virtual MGCell* clone_without_boundaries() const=0;

///Make a clone of this(this is a binder), and set binder and parameter cell
///relation between the new binder and  the parameter cell pcell.
virtual MGCell* clone_binder(const MGCellBase& c) const=0;

///Connect i1-th boundary's j1-th pcell of this to i2-th boundary's
///j2-th pcell of cell2.
///**** This connect can be applied to any manifold dimension's cell.
void connect(size_t i1, size_t j1,
	MGCell* cell2, size_t i2, size_t j2);

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void draw3DVertex()const=0;

///Erase i-th boundary.
///erase_boundary remove from this cell's bounary and destruct the boundary.
void erase_boundary(iterator i);
void erase_boundary(size_t i);

///erase_boundary removes from this cell's bounary and destruct the boundary.
void erase_boundary(MGBoundary* bnd);

///Free specified boundary(bound) from a member of parent cell's boundaries.
///Return MGBoundary if freed normally.
///If bound was not a member of the boundaries, return 0.
///Only free, does not destruct the boundary.
virtual MGBoundary* free_boundary(const MGBoundary* bound);

///Free neighbourhood relation at j-th boundary's i-th pcell of this cell.
void free_neighbourhood(size_t i, size_t j=0);

///Return Object's type ID (TID)
virtual long identify_type()const=0;

///Make a binder cell of this parameter cell.
///Returned is the binder pointer generated by new.
///The binder has no geometry, only has binder and parameter cell relationship.
virtual MGCellNB* make_binder() const=0;

///Obtain manifold dimension.
virtual unsigned manifold_dimension() const=0;

///Obtain all the neighbours.
///The neighbours do not contain this cell except when this cell is
///connected to this cell itself(closed cell).
std::vector<const MGCellNB*> neighbours() const;

///Return neighbours at the j-th boundary's i-th pcell.
///The neighbours do not contain this cell except the case that this cell is
///connected to this cell itself(closed cell) at cell i of boundary j.
std::vector<const MGCellNB*> neighbours(size_t i, size_t j=0) const;

///Return number of boundaries.
size_t number_of_boundaries() const{ return m_boundaries.size();};

///Return parameter space error of the cell.
double parameter_error()const;

///Prepend new one boundary to boundary vectors.
///Returned is the number of boudaries after prepending.
virtual size_t prepend_boundary(MGBoundary* bound);

///Sort boundary occurreces in m_boundaries.
///Sorting is done according to operator< of MGBoundary.
virtual void sort_boundaries();

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

protected:
	mutable MGBox m_box;	///<Box of this cell. =null box when 0D cell.
		///<Initially this is null,
		///<and will be computed by compute_box() when necessary.
	std::vector<MGBoundary*> m_boundaries;
		///<vector of boundaries who bound this cell.
	mutable double m_perror;	///<Error allowed for the parameter space
								///<of the cell.

///Generate a new MGCellBase pointer by newing the original MGCellBase.
///This is a proprietry routine of MGComplex copy.
///Copy all boundary data, (but does not copy own binder cell relation)
///and register boundary binder association of new and old into cmap.
virtual MGCell* clone(MGCellMap& cmap)const=0;

///compute box of the cell in m_box.
virtual void compute_box() const;

///set box as null(to set the box as initial)
virtual void set_box_as_null()const{m_box.set_null();};

///Copy all boundaries into this.
void copy_all_boundaries(const MGCellBase& cell);

///Copy all boundaries of cell into this, and binders association
///of the boundaries in the cmap.
///Binder cells of cell will be registered in cmap.
void copy_all_boundaries(const MGCellBase& cell, MGCellMap& cmap);

///Copy m_box data of cell2 into this.
void copy_box(const MGCellBase& cell2) const;

///Copy m_perror data of cell2 into this.
void copy_perror(const MGCellBase& cell2) const;

///Assignment.
///does not change binder and partner relation,
///does not change parent complex.
MGCell& set_cell(const MGCell& cell);

virtual std::string whoami()const{return "Cell";};

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

private:

///Transform the boundary binders.
void bn_binder_tr(const MGVector& v);
void bn_binder_tr(double s);
void bn_binder_tr(const MGMatrix& mat);
void bn_binder_tr(const MGTransf& tr);

///Get boundary biders of all the boundaries.
///Binders will be appended to cvec.
void get_all_boundary_binders(std::vector<MGCellNB*>& cvec) const;

///Make this cell's binder cell's extent expression.
///Returned is a MGGeometry pointer generated by new.
///When this cell does not have star cell, null pointer will be returned.
///make_binder_extent() only makes the expression, and does nothing to
///the topology structure.
virtual MGGeometry* make_binder_extent() const=0;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
virtual void make_extent() const=0;

///Negate the boundaries.
void negate_boundary();

friend class MGComplex;
friend class MGBoundary;

};

/** @} */ // end of TOPO group
#endif
