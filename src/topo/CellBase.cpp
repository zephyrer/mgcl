/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "topo/CellBase.h"
#include "topo/CellNB.h"
#include "topo/CellMap.h"
#include "topo/Complex.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Define MGCellBase Class.
//A CellBase includes only a binder cell pointer if exists.

//There are two types of cells. One is parameter cell(pcell) and
//the other is binder cell(bcell). They are exclusive, that is, if
//a cell is a parameter cell, the cell cannot be binder cell and
//vice versa.
//Parameter cell is a constituent of a complex.
//Binder cell is a binder of parameter cells. Plural cells are connected
//through a binder.
//MGCellBase ia an abstrct class.

//////Constructor///////

//////////Virtual Destructor//////////

MGCellBase::~MGCellBase(){free_partnership();}

//////operator overload/////

MGCellBase& MGCellBase::operator=(const MGCellBase& gel2){
	return set_cellbase(gel2);
}

//Assignment.
//When the leaf object of this and topo2 are not equal, this assignment
//does nothing.
MGCellBase& MGCellBase::set_cellbase(const MGCellBase& cb2){
	if(this==&cb2)
		return *this;

	free_partnership();
	MGTopology::operator =(cb2);
	return *this;
}

///////////Member Function/////////////

//Generate new MGCellNB by newing the original MGCellNB.
//Binder cells of this cell's boundaries will be registered in cmap.
//The binder of this cell itself will be registered in cmap_parent.
MGCellBase* MGCellBase::clone_pcell_with_bcell(
	MGCellMap& cmap,	//cmap is a map for this boundary binders.
	MGCellMap* cmap_parent//this is a map for this cell's binder association.
)const{
	//1. Copy cell.
	MGCellBase* new_cell=clone(cmap);

	//2. copy binder if exist.
	if(m_binder){
		MGCellNB* new_binder;
		if(cmap_parent){
			MGCellMap::map::iterator cmitr=cmap_parent->m_map.find(m_binder);
			if(cmitr!=cmap_parent->m_map.end()){
				//If the binder is already registered.
				new_binder=(*cmitr).second;
					//Binders are always MGCellNB.
				new_cell->set_binder(*new_binder);
			}else{
				//If the binder is not registered, generate a binder
				//and register its association.
				new_binder=m_binder->clone_binder(*new_cell);
				cmap_parent->m_map.insert(MGCellMap::map::value_type(m_binder,new_binder));
			}
		}else new_binder=m_binder->clone_binder(*new_cell);
	}

	return new_cell;
}

//Connect this cell to cell2.
//This cell is a pcell of a boundary of a higher manifold dimension's cell A,
//and cell2 is also is a pcell of a boundary of another cell B.
//That is, this cell is a part of a boundary of cell A,
//and cell2 is a part of a boundary of cell B.
//If cell A's manifold dimension is n, then this cell's manifold dimension is n-1.
//B's manifold dimension is same as A's. and cell2's manifold dimension is n-1.
void MGCellBase::connect(MGCellBase& cell2){
	assert(manifold_dimension()>=1 &&
		manifold_dimension()==cell2.manifold_dimension());

	MGCellNB *star_cell1=star();
	MGCellNB *star_cell2=cell2.star();

	MGComplex* complex=star_cell1->parent_complex();
	MGComplex* complex2=star_cell2->parent_complex();
	MGCellNB* bcell1=binder();
	MGCellNB* bcell2=cell2.binder();
	if(bcell1 && (bcell1==bcell2))
		return;//If already connected.

	assert(!complex || !complex2 || complex==complex2);
	//If both has complexes, they must be the same.

	if(complex&&!complex2){
		if(star_cell2->is_bcell())
			complex->append_bcell(star_cell2);
		else
			complex->append_pcell(star_cell2);
	}else if(!complex && complex2){
		if(star_cell1->is_bcell())
			complex2->append_bcell(star_cell1);
		else
			complex2->append_pcell(star_cell1);
		complex=complex2;
	}
	//Now both of star_cell1 and 2 are members of complex, if complex exists.

	if(bcell1 && bcell2){//If both had binder merge them.
		bcell1->merge_bcell(bcell2);
	}else if(bcell1){//If bcell1 had binder and bcell2 had not.
		cell2.set_binder(*bcell1);//Set binder relation to bound2.
	}else if(bcell2){//If bcell2 had binder and bcell1 had not.
		set_binder(*bcell2);//Set binder relation to bound1.
		bcell1=bcell2;
	}else{	//If both had not binders, generate one.
		bcell1=make_binder();
		cell2.set_binder(*bcell1);//Set binder relation to bound2.
	}
	if(complex)
		complex->append_bcell(bcell1);
}

//Free partnership relation of this cell.
//If the binder is prprietry use of this cell, it will be deleted.
void MGCellBase::free_partnership(){
	MGCellNB* bcell=binder();
	if(bcell){
		if(bcell->number_of_partner_members()>1)
			bcell->free_partner(this);
		else
			delete bcell;
		//Delete the bcell
		//since the bcell is proprietry to this boudary's parent cell.
		m_binder=0;
	}
}

//Make a binder associated with the extent geometry rep.
//Returned is the binder MGCellNB pointer.
//If this already had the binder with the extent,
//make_binder_with_curve() only returns the pointer.
//If this had the binder without extent, make_binder_with_extent() will
//generate the extent expression.
const MGCellNB* MGCellBase::make_binder_with_extent()const{
	assert(star());
	assert(star()->extent());

	MGGeometry* geo;
	MGCellNB* bnder=binder();
	if(bnder){
		if(!bnder->extent()){
			//If binder did not have curve representation, make it.
			geo=make_binder_extent();
			bnder->set_extent(geo);
		}
	}else{	//Generate the binder.
		bnder=make_binder();
		geo=make_binder_extent();
		bnder->set_extent(geo);
	}
	return bnder;
}

//Return number of partners.
//Partners do not inclue own cell.
size_t MGCellBase::number_of_partners() const{
	size_t n=0;
	if(m_binder)
		n=m_binder->m_partners.size()-1;
	//decrease 1 since Partners do not inclue own cell.
	return n;
}

//Obtain partner cells.
//Partners represent same world's(same cell's parameter) coordinates.
//Parameter cell's partners are parameter cells.
//Binder cell's partners are binder cells.
//The partners do not include this pcell except when star cell is
//connected to the star cell itself(closed only by the star cell).
//Let prtnrs[.] is the function's output, then prtners[0] is
//the netxt partner of this, and prtnrs[last] is the previous one of this
//in the partners array of the binder.
std::vector<const MGCellBase*> MGCellBase::partners() const{
	const MGCellNB* bcell=binder();
	if(bcell){//If have binder, get bcell's partner.
		size_t m=number_of_partners();
		std::vector<const MGCellBase*> prtnrs(m);
		size_t i,j=0; const MGCellBase* prtnr;
		//Find this in bcell->m_partners[.].
		for(i=0; i<=m; i++){
			prtnr=bcell->m_partners[i];
			if(prtnr==this)
				break;
		}
		size_t mp1=m+1;
		for(j=0; j<m; j++){
			size_t id=++i%mp1;
			prtnrs[j]=bcell->m_partners[id];
		}
		return prtnrs;
	}else
		return std::vector<const MGCellBase*>();
}

//Set binder cell relation to this cell.
//***The binder must be newed object and the owenership is transfered
//to this cell.
void MGCellBase::set_binder(MGCellNB& binder)const{
	if(m_binder!=&binder)
		binder.add_partner(*this);
}

// Output virtual function.
std::ostream& MGCellBase::out(std::ostream& ostrm) const{
	ostrm<<","<<manifold_dimension()<<"D";
	if(is_bcell())
		ostrm<<"BCell=";
	else
		ostrm<<"PCell=";
	ostrm<<this;
	if(m_binder)
		ostrm<<",binder="<<m_binder;
	MGTopology::out(ostrm);
	return ostrm;
}
