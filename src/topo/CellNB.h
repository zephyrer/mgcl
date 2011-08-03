/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCellNB_HH_
#define _MGCellNB_HH_

#include <vector>
#include "topo/CellBase.h"
#include "mg/Position.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<const MGCellBase*>;
#pragma warning( pop )
#endif

//
//Defines MGCellNB Class.

class MGBox;
class MGComplex;
class MGGeometry;
class MGCellBase;
class MGBVertex;
class MGCellMap;
class mgTLTessellate;

/** @addtogroup TOPO
 *  @{
 */

///CellNB is a cell without boundaries(No Boundaries).
///MGCellNB ia an abstract class and the super class of
///MGBVertex, MGEdge, and MGCell.
///There are two types of cells. One is parameter cell(pcell) and
///the other is binder cell(bcell). They are exclusive, that is, if
///a cell is a parameter cell, the cell cannot be binder cell and
///vice versa.
///Parameter cell is a constituent of a complex and MGCellNB's are stored
///in MGComplex.
///Binder cell is a binder of parameter cells. Plural cells are connected
///through a binder.
class MGCLASS MGCellNB:public MGCellBase{

public:

typedef std::vector<MGCellBase*>::iterator partnerItr;
typedef std::vector<const MGCellBase*>::const_iterator const_partnerItr;
typedef std::vector<MGCellBase*>::reverse_iterator partnerRItr;
typedef std::vector<const MGCellBase*>::const_reverse_iterator const_partnerRItr;

enum CELL_KIND{
	UNKNOWN=0,
	PCELL=1,
	BCELL=2
};

///////Constructor/////////

///Void constructor. Constructor of pcell.
MGCellNB();

///Copy constructor. Result cell is not a member of any complex.
///Partners of cell will not be copied.
MGCellNB(const MGCellNB& cell);

///CellNB of whole geometry(no boundary).
///The second form that input MGGeometry* takes the ownership of the geo
///into the MGCellNB, must not delete the object and the object must be
///newed one.
MGCellNB(const MGGeometry& geo);
explicit MGCellNB(MGGeometry* geo);

///This constructor takes the ownership of geo.
///Construct a parameter cell whose binder cell is 'binder'.
MGCellNB(
	MGGeometry* geo,
	MGCellNB* binder	///<Can be null.
);

////////////Virtual Destructor////////////
virtual ~MGCellNB();

///////operator overload//////

///Assignment.
///When the leaf object of this and cb2 are not equal, this assignment
///does nothing.
virtual MGCellNB& operator=(const MGCellNB& gel2);

///Object transformation.
virtual MGCellNB& operator+=(const MGVector& v);
virtual MGCellNB& operator-=(const MGVector& v);
virtual MGCellNB& operator*=(double scale);
virtual MGCellNB& operator*=(const MGMatrix& mat);
virtual MGCellNB& operator*=(const MGTransf& tr);

/////////////Member Function///////////////

///Add partner to this binder cell.
///Can change pcell(of no binder) to binder cell.
///This must be newed binder since this pointer will be registered
///in the partner's binder pointer.
void add_partner(const MGCellBase& partner);

///Obtain the box of the cell.
const MGBox& box() const=0;

///Obtain the center of this cell.
MGPosition center() const;

///Obtain the center parameter value of this cell.
virtual MGPosition center_param() const=0;

///Make a clone of the cell.
///clone() does not copy the binder cell relation.
virtual MGCellNB* clone() const=0;

///Make a clone of the cell without boundaries.
///clone_without_boundaries() does not copy the binder cell relation.
virtual MGCellNB* clone_without_boundaries() const=0;

///Make a clone of this(this is a binder), and set binder and member
///partner relation between the new binder and  the cell c.
virtual MGCellNB* clone_binder(const MGCellBase& c) const=0;

///Obtain the direction of the cell.
virtual MGUnit_vector direction() const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
virtual void draw3DVertex()const=0;

///Get extent geometry, may be null if this does not have extent.
const MGGeometry* extent() const{ return m_extent;};
MGGeometry* extent() {return m_extent;};

///Free(but does not delete) the extent geometry.
///Freed extent is returned as the function's return value.
MGGeometry* free_extent();

///Free from membership of the parent complex.
///free_from_parent() does not maintain the box of the complex this cell
///belonged to. And so, users of free_from_parent() must do it.
MGComplex* free_from_parent();

///Free neighbourhood relation at j-th boundary's i-th pcell of this cell.
virtual void free_neighbourhood(size_t i, size_t j=0)=0;

///Free specified partner(cellin).
void free_partner(const MGCellBase* cellin) const;

///Return Object's type ID (TID)
virtual long identify_type()const=0;

///Ask if this is binder cell.
bool is_bcell() const{return m_partners.size()>=1;};

///Make a binder cell of this parameter cell.
///Returned is the binder pointer generated by new.
///The binder has no geometry, only has binder and parameter cell relationship.
virtual MGCellNB* make_binder() const=0;

///Obtain manifold dimension.
virtual unsigned manifold_dimension() const=0;

///Obtain the i-th member partner. This must be a binder cell.
const MGCellBase* member_partner(size_t i)const{return m_partners[i];};

///Obtain member partners. This must be a binder cell.
const std::vector<const MGCellBase*>& member_partners()const{return m_partners;};
std::vector<const MGCellBase*>& member_partners(){return m_partners;};

///Negate the direction of the cell.
virtual void negate();

///Return nummber of partners stored in m_partners.
size_t number_of_partner_members() const{return m_partners.size();};

///Obtain all the neighbours.
///The neighbours do not contain this cell except when this cell is
///connected to this cell itself(closed cell).
virtual std::vector<const MGCellNB*> neighbours() const=0;

///Return parameter space error of the cell.
virtual double parameter_error()const=0;

///Obtain parent complex.
const MGComplex* parent_complex() const{return m_parent_complex;};
MGComplex* parent_complex() {return m_parent_complex;};

///Set extent of this cell.
virtual void set_extent(MGGeometry* extent=0);

///Obtain star cells.
const MGCellNB* star() const;
MGCellNB* star();

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;
virtual std::string whoami()const{return "CellNB";};

protected:

	mutable MGComplex* m_parent_complex;
		///<World of this cell.
		///<=Null, or m_parent_complex's parent boundary is null
		///<this cell is ordinary world cell, not parameter world.
	MGGeometry* m_extent;	///<Geometry
	mutable std::vector<const MGCellBase*> m_partners;
		///<vector of partner cells who share this bcell(for bcell).

///check if boundary's binder transformation is necessary or not.
bool bn_binder_tr_necessary()const{	return !parent_complex();}

///Generate a new MGCellBase pointer by newing the original MGCellBase.
///This is a proprietry routine of MGComplex copy.
///Copy all boundary data, (but does not copy own binder cell relation)
///and register boundary binder association of new and old into cmap.
virtual MGCellNB* clone(MGCellMap& cmap) const;

///compute box of the cell in m_box.
virtual void compute_box() const=0;

///Cell comparison.
bool is_less_than(const MGCellNB& cell2)const;

///Read Object's member data.
virtual void ReadMembers(MGIfstream& buf);

///set member datas.
///virtual void set_members(
///	MGGeometry* geo);

///Write Object's Member Data
virtual void WriteMembers(MGOfstream& buf) const;

///Assignment.
///does not change binder and partner relation,
///does not change parent complex.
MGCellNB& set_cellnb(const MGCellNB& cell2);

private:

///Transform the boundary binders.
virtual void bn_binder_tr(const MGVector& v)=0;
virtual void bn_binder_tr(double s)=0;
virtual void bn_binder_tr(const MGMatrix& mat)=0;
virtual void bn_binder_tr(const MGTransf& tr)=0;

///Set the box data as null.
virtual void set_box_as_null() const=0;

///Copy boundary data of cell2 into this.
virtual void copy_all_boundaries(const MGCellBase& cell2)=0;

///Copy all boundaries of cell into this, and binders association
///of the boundaries in the cmap.
///Binder cells of cell will be registered in cmap.
virtual void copy_all_boundaries(const MGCellBase& cell2, MGCellMap& cmap)=0;

///Copy m_box data of cell2 into this.
virtual void copy_box(const MGCellBase& cell2) const=0;

///Copy m_perror data of cell2 into this.
virtual void copy_perror(const MGCellBase& cell2) const=0;

///Free specified partner(const_partnerItr).
void free_partner(const_partnerItr itr) const;

///Get boundary biders of all the boundaries.
///Binders will be appended to cvec.
virtual void get_all_boundary_binders(std::vector<MGCellNB*>& cvec) const=0;

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

///Merge two bcells.
///bcell2 will be destructed.
virtual void merge_bcell(MGCellNB* bcell2);

///Negate the boundary.
virtual void negate_boundary()=0;

friend class MGCellBase;
friend class MGComplex;
friend class MGEdge;
friend class mgTLDataVector;
friend class mgTLTessellate;

};

/** @} */ // end of TOPO group
#endif
