/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Unit_vector.h"
#include "mg/Geometry.h"
#include "topo/Cell.h"
#include "topo/Boundary.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//
//Define MGBoundary Class.

//Constructor

//Void constructor.
MGBoundary::MGBoundary():MGComplex(), m_parent_cell(0){;}

//Constructor of one parameter cell
MGBoundary::MGBoundary(MGCellNB* pcell)
:MGComplex(pcell), m_parent_cell(0){;}

//Constructor from list of member pcells.
MGBoundary::MGBoundary(
	std::list<MGCellNB*>& pcells) //Boundary data pcells that constitue complex.
:MGComplex(pcells), m_parent_cell(0){;}

//Copy constructor.
MGBoundary::MGBoundary(
	const MGBoundary& boundary)
:MGComplex(boundary), m_parent_cell(0){	//original boundary.
	;
}

//Copy constructor with mapping.
MGBoundary::MGBoundary(
	const MGBoundary& boundary,	//original boundary.
	MGCellMap& cmap)				//cellmap to register binder association.
:MGComplex(boundary,cmap), m_parent_cell(0){
	;
}

//Virtual Destructor
MGBoundary::~MGBoundary(){
	if(m_parent_cell){
		MGCell::boundaryItr itr=m_parent_cell->boundaryIterator(this);
		if(itr!=m_parent_cell->m_boundaries.end()){
			free_binders();
			m_parent_cell->m_boundaries.erase(itr);
		}
		m_parent_cell=0;
	}
}

////////operator overload/////////
MGBoundary& MGBoundary::operator=(const MGBoundary& gel2){
	return set_boundary(gel2);
}

//When the leaf object of this and comp2 are not equal, this assignment
//does nothing.
MGBoundary& MGBoundary::set_boundary(const MGBoundary& gel2){
	MGComplex::operator=(gel2);
	if(m_parent_cell){
		m_parent_cell->free_boundary(this);
		m_parent_cell=0;
	}
	return *this;
}

//////////////Member Function/////////////

//Connect i-th pcell of this boundary to j-th pcell of bound2.
//Returned is the pointer of complex of the parent pcell of this boundary,
//may be null.
//The parent pcell of this may be a member of a complex or not.
//If bound2's parent cell is a member of a complex, it must be
//the same as this parent cell's(that is, in that case,
// this parent cell must be a member of a complex.)
void MGBoundary::connect_bound(size_t i, MGBoundary* bound2, size_t j){
	MGCellNB* cell1=pcelli(i);
	MGCellNB* cell2=bound2->pcelli(j);
	cell1->connect(*cell2);
}

//Copy boundary.
//This boundary data is cleared and bnd's boundary is copied into this.
void MGBoundary::copy_boundary(const MGBoundary& bnd){
	MGComplex::operator=(bnd);
}

//Copy boundary data, but does not copy the binders.
//This boundary data is cleared and bnd's boundary is copied into this.
void MGBoundary::copy_boundary_without_binders(const MGBoundary& bnd){
	erase_all_elements();
	m_box=bnd.m_box;
	copy_without_binders(bnd);
}

//Obtain the direction of star cell at i-th pcell of this boundary.
//Star cell's direction, not boundary's direction.
MGUnit_vector MGBoundary::direction_star(size_t i) const{
	assert(star());

	const MGCellNB* pcell=pcelli(i);
	MGPosition param=pcell->center();
	return star()->extent()->direction(param);
}

//Disconnect i-th pcell of this boundary from its partnership relation.
//disconnect does not free membership of the parent cell
//from its parent complex.
void MGBoundary::disconnect(size_t i){
	MGCellNB* bcell=binder(i);
	if(bcell){
		if(bcell->number_of_partners()) bcell->free_partner(pcelli(i));
		else delete bcell;
			//Delete the bcell
			//since the bcell is proprietry to this boudary's parent cell.
	}
}

//Test if PCell exists in this boundary.
bool MGBoundary::empty(){
	return !pcell_exist();
}

//Test if this boundary's star cell's direction is equal to bounda2's
//star cell's direction along this boundary i and bound2's boundary j.
//Not testing boundary's direction, but star cell's direction.
bool MGBoundary::equal_direction(
	size_t i, const MGBoundary& bound2, size_t j) const{
	MGUnit_vector dir1=direction_star(i);
	MGUnit_vector dir2=bound2.direction_star(j);
	return dir1%dir2>0.;
}

//Free all binders of this boundary.
//That is, if num n of partners of a binder of a pcell of this boundary
//is one, free from the parent complex, and does not free from the pcell.
//If n is more than one, free from the pcell, and does not free from
//the parent complex.
void MGBoundary::free_binders(){
	const_pcellItr cell=pcell_begin(), celle=pcell_end();;
	for(; cell!=celle; cell++){
		const MGCellNB* pcell=*cell;
		if(pcell){
			MGCellNB* bndr=pcell->binder();
			if(bndr){
				if(bndr->number_of_partner_members()>1)
					bndr->free_partner(pcell);
				else
					bndr->free_from_parent();
			}
		}
	}
}

//Reverse the direction of the boundary.
//(Coordinate transformation is not performed.)
void MGBoundary::negate(){
	pcellItr i,cs, ce;
	i=cs=pcell_begin(); ce=pcell_end();
	//Negate each parameter cell.
	for(; i!=ce; i++)
		(*i)->negate();

	//Reverse the ordering of the parameter cells.
	std::reverse(cs,ce);
}

//Negate the boundary according to the parent cell negation.
//That is,
//1. Transform the coordinates of the bondary cell.
//(This transfromation depends on how the parent cell is transformed
//when negate() is invoked. So, the member cells of this boundary
//are transformed by negate_transoform of the parent cell.)
//2. Reverse the direction of the parameter cells(negate each cell).
//3. Reverse the ordering of the parameter cells.
//(*****Does not negate the binders*****)
void MGBoundary::negate_as_boundary(const MGCellNB* parent){
	const MGCellNB* parent_cell=parent;
	if(!parent){
		assert(star()); //This must be a boudary of a cell when parent not specified.
		parent_cell=star();
	}
	const MGGeometry* parent_geo=parent_cell->extent();

	pcellItr cs, ce;
	cs=pcell_begin(); ce=pcell_end();
	for(; cs!=ce; cs++){
	//1. Transform the coordinates of each cell.
		parent_geo->negate_transform(*((*cs)->extent()));
	//2. Negate each parameter cell.
		(*cs)->negate();//std::cout<<(**cs);//////////********
	}
	//3. Reverse the ordering of the parameter cells.
	std::reverse(pcell_begin(), ce);
	m_box.set_null();
}

//count number of pcells of the boundary.
size_t MGBoundary::number_of_pcells() const{
	return MGComplex::number_of_pcells();
}

// Output virtual function.
std::ostream& MGBoundary::out(std::ostream& outpt) const{
	outpt<<"Boundary::m_parent_cell="<<m_parent_cell<<",";
	MGComplex::out(outpt);
	return outpt;
}

//Set binder relation to i-th parameter cell.
void MGBoundary::set_binder(size_t i, MGCellNB& binder)const{
	const_pcellItr pcellp=pcell_begin();
	std::advance(pcellp,i);
	(*pcellp)->set_binder(binder);
}

//Set parent cell.
//Returned is the conventional parent cell attached to
//before execution of this set_parent.
MGCell* MGBoundary::set_parent(MGCell& new_parent)const{
	MGCell* parent=m_parent_cell;
	if(parent){
		if(parent==&new_parent) return parent;
		new_parent.free_boundary(this);
	}
	m_parent_cell=&new_parent;
	return parent;
}

//Get the star cell.
const MGCellNB* MGBoundary::star() const{return m_parent_cell;};
MGCellNB* MGBoundary::star(){return m_parent_cell;};
