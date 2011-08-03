/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Geometry.h"
#include "mg/Point.h"
#include "topo/CellNB.h"
#include "topo/CellMap.h"
#include "topo/Complex.h"
#include "topo/Boundary.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCellNB Class.
//CellNB is a cell without boundaries(No Boundaries).

//There are two types of cells. One is parameter cell(pcell) and
//the other is binder cell(bcell). They are exclusive, that is, if
//a cell is a parameter cell, the cell cannot be binder cell and
//vice versa.
//Parameter cell is a constituent of a complex.
//Binder cell is a binder of parameter cells. Plural cells are connected
//through a binder.
//MGCellNB ia an abstract class.

//////Constructor///////

//Void constructor. Constructor of pcell.
MGCellNB::MGCellNB():m_extent(0), m_parent_complex(0){;}

//Copy constructor. Result cell is not a member of any complex.
//Partners of cell will not be copied.
MGCellNB::MGCellNB(const MGCellNB& cell)
:MGCellBase(cell), m_extent(0), m_parent_complex(0){
	if(cell.m_extent)
		m_extent=cell.m_extent->clone();
}

//CellNB of whole geometry(no boundary), under parent.
//The second form that input MGGeometry* takes the ownership of the geo
//into the MGCellNB, must not delete the object and the object must be
//newed one.
MGCellNB::MGCellNB(const MGGeometry& geo)
:m_extent(geo.clone()),m_parent_complex(0){;}
MGCellNB::MGCellNB(MGGeometry* geo)
:m_extent(geo),m_parent_complex(0){;}

//This constructor takes the ownership of geo.
//Construct a parameter cell whose binder cell is 'binder'.
MGCellNB::MGCellNB(
	MGGeometry* geo,
	MGCellNB* binder)
:m_extent(geo),m_parent_complex(0){
	if(binder) binder->add_partner(*this);
}

//////////Virtual Destructor//////////
MGCellNB::~MGCellNB(){
	free_from_parent();//Free this cell from parent complex.
//When this destructor is invoked, all of the boundaries must be deleted.
	const_partnerItr itr=m_partners.begin();
	const_partnerItr itre=m_partners.end();
	if(itr!=itre){
		//When this is a binder cell.
		for(; itr!=itre; itr++)
			(*itr)->m_binder=0;
		m_partners.clear();
	}
	delete m_extent;
}

//////operator overload/////
MGCellNB& MGCellNB::operator=(const MGCellNB& gel2){
	return set_cellnb(gel2);
}

//Assignment.
//does not change binder and partner relation,
//does not change parent complex.
MGCellNB& MGCellNB::set_cellnb(const MGCellNB& cell2){
	MGCellBase::operator=(cell2);
	if(m_extent)
		delete m_extent;
	if(cell2.m_extent)
		m_extent=cell2.m_extent->clone();
	else
		m_extent=0;
	return *this;
}

// CellNBに平行移動を行ない自身のCellNBとする。
//Translation of the cell
MGCellNB& MGCellNB::operator+= (const MGVector& v){
	if(m_extent){
		*m_extent +=v; set_box_as_null();
		bn_binder_tr(v); //change binder.
	}
	return *this;
}

// CellNBに逆方向の平行移動を行ない自身のCellNBとする。
//Translation of the cell
MGCellNB& MGCellNB::operator-= (const MGVector& v){
	if(m_extent){
		*m_extent -=v; set_box_as_null();
		bn_binder_tr(-v); //change binder.
	}
	return *this;
}

//CellNBのスケーリングを行い自身のCellNBとする。
//Scaling of the cell by a double.
MGCellNB& MGCellNB::operator*= (double s){
	if(m_extent){
		*m_extent *=s; set_box_as_null();
		bn_binder_tr(s); //change binder.
	}
	return *this;
}

// 与えられた変換でCellNBの変換を行い自身のCellNBとする。
//Transformation of the cell by a matrix.
MGCellNB& MGCellNB::operator*= (const MGMatrix& mat){
	if(m_extent){
		*m_extent *=mat; set_box_as_null();
		bn_binder_tr(mat); //change binder.
	}
	return *this;
}

// 与えられた変換によってトランスフォームをおこない自身のCellNBにする。
//Transformation of the cell by a MGTransf.
MGCellNB& MGCellNB::operator*= (const MGTransf& tr){
	if(m_extent){
		*m_extent *=tr; set_box_as_null();
		bn_binder_tr(tr); //change binder.
	}
	return *this;
}

//Cell comparison.
bool MGCellNB::is_less_than(const MGCellNB& cell2)const{
	if(this==&cell2)
		return false;

	const MGComplex* comp1=parent_complex();
	if(!comp1)
		return true;
	const MGComplex* comp2=cell2.parent_complex();
	if(!comp2)
		return false;

	if(comp1!=comp2)
		return (*comp1)<(*comp2);

	//Now this and cell2 have the same parent complex.
	//Comparison is done by the appearance order of this or cell2 in the complex.
	MGComplex::const_pcellItr i=comp1->pcell_begin(), ie=comp1->pcell_end();
	for(;i!=ie; i++){
		if((*i)==this) return true;
		if((*i)==&cell2) return false;
	}
	return true;//This statement will never be executed.
}

///////////Member Function/////////////

bool MGCellBaseComparison(const MGCellBase* cell1, const MGCellBase* cell2){
   return *cell1 < *cell2;
}

//Add partner to this binder cell.
//Can change pcell(of no binder) to binder cell.
//This must be newed binder since this pointer will be registered
//in the partner's binder pointer.
void MGCellNB::add_partner(const MGCellBase& partner){
	const MGCellNB* binder=partner.binder();
	if(binder==this)
		return;

	assert(m_partners.size()==0 ||
		m_partners[0]->manifold_dimension()==partner.manifold_dimension());

	if(binder)
		binder->free_partner(&partner);
	partner.m_binder=this;
	size_t npartners=m_partners.size();
	if(npartners<=1)
		m_partners.push_back(&partner);
	else{
		std::vector<const MGCellBase*>::iterator
			i=std::lower_bound(m_partners.begin(),m_partners.end(),&partner,MGCellBaseComparison);
		m_partners.insert(i,&partner);
	}

	//When this binder does not belong to a complex,
	//add this to the complex of the star cell of the partner.
	if(!m_parent_complex){
		const MGCellNB* starcell=partner.star();
		if(starcell){
			const MGComplex* complex=starcell->parent_complex();
			if(complex)
				complex->append_bcell(this);
		}
	}
}

//Obtain the center of this cell.
MGPosition MGCellNB::center() const{
	MGPosition cntr;
	if(m_extent)
		cntr=m_extent->evaluate(center_param());
	return cntr;
}

//Generate a new MGCellBase pointer by newing the original MGCellBase.
//This is a proprietry routine of MGComplex copy.
//Copy all boundary data, (but does not copy own binder cell relation.)
//and register new and old boundary binders association into cmap.
MGCellNB* MGCellNB::clone(MGCellMap& cmap) const{
	MGCellNB* cell=clone_without_boundaries();
	if(manifold_dimension()>=1){
		cell->copy_all_boundaries(*this,cmap);
		cell->copy_box(*this);
		cell->copy_perror(*this);
	}
	return cell;
}

//Obtain the direction of the cell.
MGUnit_vector MGCellNB::direction() const{
	MGPosition param=center_param();
	return extent()->direction(param);
}

//Free(but does not delete) the extent geometry.
//Freed extent is returned as the function's return value.
MGGeometry* MGCellNB::free_extent(){
	MGGeometry* extent=m_extent;
	m_extent=0;
	return extent;
}

//Free from membership of the parent complex.
//free_from_parent() does not maintain the box of the complex this cell
//belonged to. And so, users of free_from_parent() must do it.
MGComplex* MGCellNB::free_from_parent(){
	MGComplex* parent=m_parent_complex;
	if(!parent) return parent;

	//1. free binder of this cell.
	MGCellNB* bnder=binder();
	if(bnder){
		if(bnder->m_partners.size()>1)
			bnder->free_partner(this);
		else
			bnder->free_from_parent();
	}

	//2. free this cell from parent complex.
	MGComplex::cellItr itri, itrend;

	itri=parent->m_pcells.begin(); itrend=parent->m_pcells.end();
	for(; itri!=itrend; itri++){
		if(this==(*itri)) {
			parent->m_pcells.erase(itri);
			m_parent_complex=0;
			return parent;
		}
	}

	itri=parent->m_bcells.begin(); itrend=parent->m_bcells.end();
	for(; itri!=itrend; itri++){
		if(this==(*itri)) {
			parent->m_bcells.erase(itri);
			m_parent_complex=0;
			return parent;
		}
	}
	
	return parent;
}

//Free specified partner(const_partnerItr).
void MGCellNB::free_partner(const_partnerItr itrin) const{
	size_t i=itrin-m_partners.begin();
	const_partnerItr itr_next=itrin; itr_next++;
	const_partnerItr itre=m_partners.end();
	while(itr_next!=itre)
		m_partners[i++]=(*itr_next++);
	m_partners.pop_back();
}

//Free specified partner(cellin).
void MGCellNB::free_partner(const MGCellBase* cellin) const{
	assert(is_bcell());	//This must be bcel.
	const_partnerItr itr=m_partners.begin();
	const_partnerItr itre=m_partners.end();
	for(; itr!=itre; itr++){
		if(cellin==*itr){
			free_partner(itr);
			cellin->m_binder=0;
			break;
		}
	}
}

//Merge two bcells.
//bcell2 will be destructed.
void MGCellNB::merge_bcell(MGCellNB* bcell2){
	assert(is_bcell() && bcell2->is_bcell());
	assert(manifold_dimension()==bcell2->manifold_dimension());

	MGComplex *complex1, *complex2;
	complex1=parent_complex(); complex2=bcell2->parent_complex();
	assert(complex1 || complex2);//one of the two must have complex.
	assert(((complex1 && complex2) && complex1==complex2)
		   || !(complex1&&complex2));
	//If both has complexes, they must be the same.

	//Merge extent.
	if(!m_extent){
		m_extent=bcell2->m_extent;//copy extent.
		copy_box(*bcell2);
		bcell2->m_extent=0;
	}

	//Merge own binder.
	const_partnerItr itr=bcell2->m_partners.begin();
	const_partnerItr itre=bcell2->m_partners.end();
	for(; itr!=itre; itr++){
		//Change binder from bcell2 to this.
		add_partner(**itr);
	}
	delete bcell2;
}

//Negate the direction of the cell.
void MGCellNB::negate(){
	if(m_extent){
		//1. Negate each boudary.
		negate_boundary();

		//2. Negate own extent.
		m_extent->negate();
		set_box_as_null();
	}
}

//Set extent of this cell.
void MGCellNB::set_extent(MGGeometry* extent){
	if(m_extent)
		delete m_extent;
	m_extent=extent; set_box_as_null();
	MGComplex* parent;
	if(parent=parent_complex())
		parent->m_box.set_null();
}

//Obtain star cells.
const MGCellNB* MGCellNB::star() const{
	const MGCellNB* astar=0;
	if(is_bcell()){
		astar=m_partners[0]->star();
		if(astar) astar=astar->star();
	}else{//Case that this is parameter cell.
		if(m_parent_complex){
			const MGBoundary* bound=
				dynamic_cast<const MGBoundary*>(m_parent_complex);
			if(bound) astar=bound->star();
		}
	}
	return astar;
}
MGCellNB* MGCellNB::star(){
	const MGCellNB* cthis=this;
	return const_cast<MGCellNB*>(cthis->star());
}

// Output virtual function.
std::ostream& MGCellNB::out(std::ostream& ostrm) const{
	if(m_parent_complex)
		ostrm<<",parent="<<m_parent_complex;
	ostrm<<",extent="<<m_extent<<"::"<<std::endl;
	if(m_extent)
		ostrm<<*m_extent;
	ostrm<<",CellNB::partners="<<m_partners.size();
	MGCellNB::const_partnerItr ps=m_partners.begin(),	pe=m_partners.end();
	if(ps!=pe){
		ostrm<<"::"; ostrm<<(*ps++);
		for(; ps!=pe; ps++) ostrm<<","<<(*ps);
	}
	MGCellBase::out(ostrm);
	return ostrm;
}
