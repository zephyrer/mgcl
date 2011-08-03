/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Geometry.h"
#include "mg/Tolerance.h"
#include "topo/Boundary.h"
#include "topo/Cell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCell Class.
//MGCell is a general cell that has bound.
//MGCell's additional data to CellNB are Boundary, box, and perror.
//MGCell ia an abstrct class and the super class of MGFace, Solid, and
//other general cells of more than 1 manifold dimension.

//There are two types of cells. One is parameter cell(pcell) and
//the other is binder cell(bcell). They are exclusive, that is, if
//a cell is a parameter cell, the cell cannot be binder cell and
//vice versa.
//Parameter cell is a constituent of a complex.
//Binder cell is a binder of parameter cells. Plural cells are connected
//through a binder.

//////Constructor///////

//Void constructor. Constructor of pcell.
MGCell::MGCell():m_perror(-1.){;}

//Copy constructor. Result cell is not a member of any complex.
//Binders and boundaries of cell will not be copied.
//Copy of boudaries can be done by copy_all_boundaries().
MGCell::MGCell(const MGCell& cell)
:MGCellNB(cell), m_perror(cell.m_perror){;}

//MGCell of whole geometry(no boundary), under parent.
//Constructor of pcell.
//The second form that input MGGeometry* takes the ownership of the geo
//into the MGCell, must not delete the object and the object must be
//newed one.
MGCell::MGCell(const MGGeometry& geo)
:MGCellNB(geo), m_perror(-1.){;}
MGCell::MGCell(MGGeometry* geo)
:MGCellNB(geo), m_perror(-1.){;}

//Construct a parameter cell from all the necessary data,
//geo, vector of boundaries, and the binder. The newly constructed parameter
//cell will be a partner member of the binder 'binder'.
//Constructor takes the ownership of goe and MGBoundary in boundaries.
MGCell::MGCell(
	MGGeometry* geo,
	std::vector<MGBoundary*>& boundaries,
	MGCell* binder)
:MGCellNB(geo, binder), m_perror(-1.), m_boundaries(boundaries){;}

//Parameter cell with boundaries.
//Only pcells(parameter representation) in boundaries are copied.
//Binders(world coordinate representation) in boundaries are discarded.
MGCell::MGCell(const MGGeometry& geo, 
		const std::vector<MGBoundary*>& boundaries)
:MGCellNB(geo), m_perror(-1.){
	const_boundaryItr i=boundaries.begin(), iend=boundaries.end();
	while(i!=iend){
		MGBoundary* bnd=(*i++)->clone_without_binders(*this);
		append_boundary(bnd);
	}
}
MGCell::MGCell(MGGeometry* geo, 
		const std::vector<MGBoundary*>& boundaries)
:MGCellNB(geo), m_perror(-1.){
	const_boundaryItr i=boundaries.begin(), iend=boundaries.end();
	while(i!=iend){
		MGBoundary* bnd=(*i++)->clone_without_binders(*this);
		append_boundary(bnd);
	}
}

//////////Virtual Destructor//////////
MGCell::~MGCell(){
	while(!m_boundaries.empty()) delete m_boundaries.back();
}

MGCell& MGCell::operator=(const MGCell& gel2){
	return set_cell(gel2);
}

//Assignment.
//does not change binder and partner relation,
//does not change parent complex.
MGCell& MGCell::set_cell(const MGCell& cell){
	while(!m_boundaries.empty())
		delete m_boundaries.back();

	MGCellBase::operator=(cell);
	m_box=cell.m_box;
	m_perror=cell.m_perror;
	copy_all_boundaries(cell);
	return *this;
}

//Comparison operator.
bool MGCell::operator<(const MGCell& e2)const{
	return is_less_than(e2);
}

///////////Member Function/////////////

//Append new one boundary to boundary vectors.
//Returned is the number of boudaries after appending.
size_t MGCell::append_boundary(MGBoundary* bound){
	if(bound){
		m_boundaries.push_back(bound);
		bound->set_parent(*this);
	}
	set_box_as_null();
	return m_boundaries.size();
}

//Transform the boundary binders.
void MGCell::bn_binder_tr(const MGVector& v){
	if(bn_binder_tr_necessary()){
		boundaryItr bs,be;
		bs=m_boundaries.begin(); be=m_boundaries.end();
		for(; bs!=be; bs++) (**bs).binderTr(v);
	}
}
void MGCell::bn_binder_tr(double s){
	if(bn_binder_tr_necessary()){
		boundaryItr bs,be;
		bs=m_boundaries.begin(); be=m_boundaries.end();
		for(; bs!=be; bs++) (**bs).binderTr(s);
	}
}
void MGCell::bn_binder_tr(const MGMatrix& mat){
	if(bn_binder_tr_necessary()){
		boundaryItr bs,be;
		bs=m_boundaries.begin(); be=m_boundaries.end();
		for(; bs!=be; bs++) (**bs).binderTr(mat);
	}
}
void MGCell::bn_binder_tr(const MGTransf& tr){
	if(bn_binder_tr_necessary()){
		boundaryItr bs,be;
		bs=m_boundaries.begin(); be=m_boundaries.end();
		for(; bs!=be; bs++) (**bs).binderTr(tr);
	}
}

//Obtain i-th boundary's j-th pcell direction
//(direction of boundary measured by this cell's
//coordinate along the boundary).
//The direction is represented by the center of the boundary.
MGVector MGCell::boundary_direction(size_t i, size_t j) const{
	assert(i<number_of_boundaries());
	assert(j<m_boundaries[i]->number_of_pcells());
	assert(manifold_dimension()>1);

	size_t nderi1[]={1,0};
	size_t nderi2[]={0,1};

	unsigned n=manifold_dimension();
	const MGCellNB* pc=m_boundaries[i]->pcelli(j);
	MGPosition bound_center=pc->center();
	MGPosition bound_center_param=pc->center_param();
	MGVector dt=pc->extent()->evaluate(bound_center_param, nderi1);
	MGVector dS1=extent()->evaluate(bound_center, nderi1);
	MGVector dS2=extent()->evaluate(bound_center, nderi2);
	return dS1*dt(0)+dS2*dt(1);
}

//Obtain iterator of m_boundaries.
MGCell::const_boundaryItr MGCell::boundaryIterator(const MGBoundary* bnd) const{
	const_boundaryItr itr=m_boundaries.begin(), itre=m_boundaries.end();
	for(; itr!=itre; itr++){ if((*itr)==bnd) break;}
	return itr;
}
MGCell::boundaryItr MGCell::boundaryIterator(MGBoundary* bnd){
	boundaryItr itr=m_boundaries.begin(), itre=m_boundaries.end();
	for(; itr!=itre; itr++) {if((*itr)==bnd) break;}
	return itr;
}

//Obtain the box of the cell.
const MGBox& MGCell::box() const{
	if(m_box.is_null()) compute_box();
	return m_box;	
}

//Obtain the center parameter value of this cell.
MGPosition MGCell::center_param() const{
	MGPosition cntr;
	const_boundaryItr bb=m_boundaries.begin(), be=m_boundaries.end();
	if(bb==be) cntr=m_extent->center_param();//If does not have boundary.
	else{
		//This manifold dimension is always greater than 1.
		const MGBoundary* bnd=dynamic_cast<const MGBoundary*>(*bb);
		//Barycenter of the vertices of the first boundary is the center.
		cntr=bnd->center();
	}
	return cntr;
}

//compute box of the cell in m_box.
void MGCell::compute_box() const{if(m_extent) m_box=m_extent->box();}

//Connect i1-th boundary's j1-th pcell of this to i2-th boundary's
//j2-th pcell of cell2.
//**** This connect can be applied to any manifold dimension's cell.
void MGCell::connect(size_t i1, size_t j1,
					MGCell* cell2, size_t i2, size_t j2){
	assert(i1<number_of_boundaries() &&
		j1<m_boundaries[i1]->number_of_pcells());
	assert(i2<cell2->number_of_boundaries() &&
		j1<cell2->m_boundaries[i1]->number_of_pcells());
	assert(!is_pcell() || (is_pcell()&&cell2->is_pcell()));

	MGBoundary *bound1=m_boundaries[i1], *bound2=cell2->m_boundaries[i2];
	bound1->connect_bound(j1,bound2,j2);
}

//Copy all boundaries of cell into this, and binders association
//of the boundaries in the cmap.
//Binder cells of cell will be registered in cmap.
void MGCell::copy_all_boundaries(const MGCellBase& cell,MGCellMap& cmap){
	assert(dynamic_cast<const MGCell*>(&cell));

	const MGCell* cb=dynamic_cast<const MGCell*>(&cell);
	const_boundaryItr bs,be;
	bs=cb->m_boundaries.begin(); be=cb->m_boundaries.end();
	for(; bs!=be; bs++){
		m_boundaries.push_back((*bs)->clone(*this,cmap));
	}
	set_box_as_null();
}

//Copy all boundaries into this.
void MGCell::copy_all_boundaries(const MGCellBase& cellin){
	assert(dynamic_cast<const MGCell*>(&cellin));

	const MGCell* cell=dynamic_cast<const MGCell*>(&cellin);
	const_boundaryItr bs,be;
	bs=cell->m_boundaries.begin(); be=cell->m_boundaries.end();
	for(; bs!=be; bs++){
		m_boundaries.push_back((*bs)->clone(*this));
	}
	set_box_as_null();
}

//Copy m_box data of cell2 into this.
void MGCell::copy_box(const MGCellBase& cell2) const{
	assert(dynamic_cast<const MGCell*>(&cell2));
	const MGCell* c=dynamic_cast<const MGCell*>(&cell2);
	m_box=c->m_box;
}

//Copy m_perror data of cell2 into this.
void MGCell::copy_perror(const MGCellBase& cell2) const{
	assert(dynamic_cast<const MGCell*>(&cell2));
	const MGCell* c=dynamic_cast<const MGCell*>(&cell2);
	m_perror=c->m_perror;
}

//Erase i-th boundary.
//erase_boundary remove from this cell's bounary and destruct the boundary.
void MGCell::erase_boundary(size_t i){
	assert(i<number_of_boundaries());
	iterator itr=m_boundaries.begin();
	std::advance(itr,i);
	erase_boundary(itr);
}
void MGCell::erase_boundary(iterator i){
	delete *i;
	set_box_as_null();
}

//Erase i-th boundary.
//erase_boundary remove from this cell's bounary and destruct the boundary.
void MGCell::erase_boundary(MGBoundary* bnd){
	boundaryItr itr=m_boundaries.begin();
	boundaryItr itre=m_boundaries.end();
	for(; itr!=itre; itr++){
		if(bnd==(*itr)){
			delete *itr;
			break;
		}
	}
	set_box_as_null();
}

//Free specified boundary(bound) from a member of parent cell's boundaries.
//Return MGBoundary if freed normally.
//If bound was not a member of the boundaries, return 0.
//Only free, does not destruct the boundary.
MGBoundary* MGCell::free_boundary(const MGBoundary* bound){
	boundaryItr itr=m_boundaries.begin();
	boundaryItr itre=m_boundaries.end();
	for(; itr!=itre; itr++){
		if(bound==*itr){
			MGBoundary* bn=*itr;
			m_boundaries.erase(itr);
			set_box_as_null();
			return bn;
		}
	}
	return 0;
}

//Free neighbourhood relation at j-th boundary's i-th pcell of this cell.
void MGCell::free_neighbourhood(size_t i, size_t j){
	assert(j<number_of_boundaries() && i<m_boundaries[j]->number_of_pcells());
	//Disconnect cell.
	m_boundaries[j]->disconnect(i);
}

//Get boundary biders of all the boundaries.
//Binders will be appended to cvec.
void MGCell::get_all_boundary_binders(std::vector<MGCellNB*>& cvec) const{
	const_boundaryItr itr=m_boundaries.begin();
	const_boundaryItr itre=m_boundaries.end();
	std::vector<MGCellNB*> bndbdrs;
	for(;itr!=itre; itr++){	//Loop for all the boundaries.
		MGComplex::const_pcellItr ps=(*itr)->pcell_begin(), pe=(*itr)->pcell_end();
		for(; ps!=pe; ps++){
			MGCellNB* bndr=(*ps)->binder();
			if(bndr) cvec.push_back(bndr);
		}
	}
}

class MGCLASS MGBoundaryCompare{
public:
	bool operator()(const MGBoundary* b1, const MGBoundary* b2) const
	{ return (*b1)<(*b2);};
};

//Negate the boundaries.
void MGCell::negate_boundary(){
	//Negate each boudary.
	boundaryItr bitr=m_boundaries.begin(), bitrend=m_boundaries.end();
	for(; bitr!=bitrend; bitr++)
		(*bitr)->negate_as_boundary(this);
	//Sort the boundaries.
	std::sort(m_boundaries.begin(),m_boundaries.end(),MGBoundaryCompare());
}

//Obtain all the neighbours.
//The neighbours do not contain this cell except when this cell is
//connected to this cell itself(closed cell).
std::vector<const MGCellNB*> MGCell::neighbours() const{
	const_boundaryItr itr=m_boundaries.begin();
	const_boundaryItr itre=m_boundaries.end();
	std::vector<const MGCellNB*> nbrs;
	for(;itr!=itre; itr++){	//Loop for all the boundaries.
		size_t n=(*itr)->number_of_pcells();
		for(size_t i=0; i<n; i++){	//Loop for all cells in a boundary.
			const MGCellNB* vc=(*itr)->pcelli(i);
			std::vector<const MGCellBase*> prtnrs=vc->partners();
			std::vector<const MGCellBase*>::iterator
				prtnr_now=prtnrs.begin(), prtnr_end=prtnrs.end();
			for(;prtnr_now!=prtnr_end; prtnr_now++){
				const MGCellNB* nbr=(*prtnr_now)->star();
				if(nbr){
					if(std::find(nbrs.begin(),nbrs.end(),nbr)==nbrs.end())
						nbrs.push_back(nbr);
				}
			}
		}
	}
	return nbrs;
}

//Return neighbours at the j-th boundary's i-th pcell.
//The neighbours do not contain this cell except the case that this cell is
//connected to this cell itself(closed cell) at cell i of boundary j.
std::vector<const MGCellNB*> MGCell::neighbours(size_t i, size_t j) const{
	assert(j<number_of_boundaries());
	assert(i<(boundary(j)->number_of_pcells()));

	const MGCellNB* vc=boundary(j)->pcelli(i);
	std::vector<const MGCellBase*> prtnrs=vc->partners();
	std::vector<const MGCellBase*>::iterator
				prtnr_now=prtnrs.begin(), prtnr_end=prtnrs.end();
	std::vector<const MGCellNB*> nbrs;
	for(;prtnr_now!=prtnr_end; prtnr_now++){
		const MGCellNB* nbr=(*prtnr_now)->star();
		if(nbr) nbrs.push_back(nbr);
	}
	return nbrs;
}

//Return parameter space error of the cell.
double MGCell::parameter_error()const{
	if(m_perror<=0.){
		if(m_extent) m_perror=m_extent->parameter_error();
		else return MGTolerance::wc_zero();
	}
	return m_perror;
}

//Prepend new one boundary to boundary vectors.
//Returned is the number of boudaries after prepending.
size_t MGCell::prepend_boundary(MGBoundary* bound){
	m_boundaries.insert(m_boundaries.begin(),bound);
	bound->set_parent(*this);
	set_box_as_null();
	return m_boundaries.size();
}

//Sort boundary occurreces in m_boundaries.
//Sorting is done according to operator< of MGBoundary.
void MGCell::sort_boundaries(){
	std::sort(m_boundaries.begin(),m_boundaries.end(),MGBoundaryCompare());
}

// Output virtual function.
std::ostream& MGCell::out(std::ostream& ostrm) const{
	ostrm<<",box="<<m_box;
	ostrm<<",boundaries="<<m_boundaries.size();
	const_boundaryItr bs=m_boundaries.begin(),be=m_boundaries.end();
	size_t n=0;
	if(bs!=be){
		ostrm<<"::"<<std::endl;
		ostrm<<(**bs).whoami()<<n++<<"="<<(**bs++);
		for(; bs!=be; bs++){
			ostrm<<std::endl;
			ostrm<<n++<<"="<<(**bs);
		}
	}
	MGCellNB::out(ostrm);
	return ostrm;
}
