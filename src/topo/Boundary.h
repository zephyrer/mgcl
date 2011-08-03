/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGBoundarynD_HH_
#define _MGBoundarynD_HH_

#include <list>
#include <vector>
#include "topo/Complex.h"

class MGUnit_vector;
class MGCellNB;
class MGLoop;
class MGCell;

//
//Define MGBoundary Class.

/** @addtogroup TOPO
 *  @{
 */

///MGBoundary is a boundary of more than 1 manifold dimension. That is,
///a boundary of a face, a volume, or other general manifold dimension's
///cells. Cells stored in parent class complex constitute this boundary.
///This boundary cell's geometry's space dimension is
///the same as the parent cell's manifold dimension, since
/// boundary cell is parameter space world.
///For edges whose manifold_dimension() is 1, this class is not used.
///Instead, Edges conatain the parameter cell MGPVertex.
///(MGPVertex is a parameter cell and also a boundary.)
///MGBoundary is an abstract class.
class MGCLASS MGBoundary: public MGComplex{

public:

/////////Constructor/////////

///Void constructor.
MGBoundary();

///Constructor of one parameter cell
explicit MGBoundary(MGCellNB* pcell);

///Constructor from list of member pcells.
explicit MGBoundary(
	std::list<MGCellNB*>& pcells); ///Boundary data pcells that constitue complex.

///Copy constructor.
MGBoundary(const MGBoundary& boundary);///original boundary.

///Virtual Destructor
virtual ~MGBoundary();

/////////operator overload/////////

///Assignment.
///When the leaf object of this and comp2 are not equal, this assignment
///does nothing.
virtual MGBoundary& operator=(const MGBoundary& gel2);

///Object transformation.
virtual MGBoundary& operator+=(const MGVector& v){MGComplex::operator+=(v);return *this;};
virtual MGBoundary& operator-=(const MGVector& v){MGComplex::operator-=(v);return *this;};
virtual MGBoundary& operator*=(double scale){MGComplex::operator*=(scale);return *this;};
virtual MGBoundary& operator*=(const MGMatrix& mat){MGComplex::operator*=(mat);return *this;};
virtual MGBoundary& operator*=(const MGTransf& tr){MGComplex::operator*=(tr);return *this;};

/////////Member Function/////////

///Test if this is an active boundary.
virtual bool active() const=0;

///Return the box of this boundary.
const MGBox& box() const{return MGComplex::box();};

///Make a clone.
///Returned is pointer of newed object, must be deleted.
///When parent is specified, clone's parent is set to the parent.
virtual MGBoundary* clone(MGCell& parent) const=0;
virtual MGBoundary* clone()const=0;

///Make a clone that has not binders.
virtual MGBoundary* clone_without_binders(MGCell& parent) const=0;
virtual MGBoundary* clone_without_binders() const=0;

///Test if this is closed boundary.
virtual bool closed() const=0;

///Obtain the direction of star cell at the i-th pcell of this boundary.
///Star cell's direction, not boundary's direction.
MGUnit_vector direction_star(size_t i) const;

///Test if PCells exist in this boundary.
///If no pcells are included, empty() returns true.
bool empty();

///Test if this boundary's star cell's direction is equal to bounda2's
///star cell's direction along this boundary's i and bound2's boundary j-th
///parameter cell.
///Not testing boundary's direction, but star cell's direction.
bool equal_direction(size_t i, const MGBoundary& bound2, size_t j) const;

///Return Object's type ID (TID)
virtual long identify_type()const=0;

///Get manifold dimension.
virtual unsigned manifold_dimension() const=0;

///Reverse the direction of the boundary.
///(Coordinate transformation is not performed.)
virtual void negate();

///Negate the boundary according to the parent cell negation.
///That is,
///1. Transform the coordinates of the bondary cell.
///(This transfromation depends on how the parent cell is transformed
///when negate() is invoked. So, the member cells of this boundary
///are transformed by negate_transoform of the parent cell.)
///2. Reverse the direction of the parameter cells(negate each cell).
///3. Reverse the ordering of the parameter cells.
//(*****Does not negate the binders*****)
virtual void negate_as_boundary(const MGCellNB* parent=0);

///Obtain how many parameter cells are included in the boundary.
size_t number_of_pcells() const;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

///Set binder relation to i-th parameter cell.
///***The binder must be newed object and the owenership is transfered
///to boundary(precisely, to parameter cell that corresponds to the binder).
void set_binder(size_t i, MGCellNB& binder)const;

///Set parent cell.
///Returned is the conventional parent cell attached to
///before execution of this set_parent.
MGCell* set_parent(MGCell& new_parent) const;

///Get the star cell.
const MGCellNB* star() const;
MGCellNB* star();

virtual std::string whoami()const{return "Boundary";};

//////////////////PROTECTED MEMBER///////////////////////////////
protected:

	mutable MGCell* m_parent_cell;///<Cell that has this boundary as a boundary.

///Copy constructor.
///Binder cells of the pcells in boundary will be registered in cmap.
MGBoundary(
	const MGBoundary& boundary,///<original boundary.
	MGCellMap& cmap		 ///<cellmap to register binder association.
);

///Make a clone.
///The forms that have cmap as an argumetnt is to register binder association.
///Binder cells of the pcells in this boundary will be registered in cmap.
///Returned is pointer of newed object, must be deleted.
///When parent is specified, clone's parent is set to the parent.
virtual MGBoundary* clone(MGCell& parent, MGCellMap& cmap) const=0;
virtual MGBoundary* clone(MGCellMap& cmap) const=0;

///Connect i-th pcell of this boundary to j-th pcell of boud2.
///Returned is the pointer of complex of the parent pcell of this boundary.
///The parent pcell of this must be a memeber of a complex.
///If both of this and bound2 are a member of a complex, they must be
///the same.
void connect_bound(size_t i, MGBoundary* bound2, size_t j);

///Copy boundary.
///This boundary data is cleared and bnd's boundary is copied into this.
virtual void copy_boundary(const MGBoundary& bnd);

///Copy boundary data, but does not copy the binders.
///This boundary data is cleared and bnd's boundary is copied into this.
virtual void copy_boundary_without_binders(const MGBoundary& bnd);

///Disconnect i-th pcell of this boundary from its partnership relation.
///disconnect does not free membership of the parent cell
///from its parent complex.
void disconnect(size_t i);

///Free all binders of this boundary.
///That is, if num n of partners of a binder of a pcell of this boundary
///is one, free from the parent complex, and does not free from the pcell.
///If n is more than one, free from the pcell, and does not free from
///the parent complex.
void free_binders();

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///When the leaf object of this and comp2 are not equal, this assignment
///does nothing.
MGBoundary& set_boundary(const MGBoundary& gel2);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

private:

friend MGCell;

};

/** @} */ // end of TOPO group
#endif
