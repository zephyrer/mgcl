/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Geometry.h"
#include "mg/Point.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "topo/PVertex.h"
#include "topo/Complex.h"
#include "topo/CellNB.h"
#include "topo/Boundary.h"
#include "topo/CellMap.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implement MGComplex Class.

//Constructor

//Void constructor (under a boundary).
MGComplex::MGComplex(){;}

//Construct of one cell.
//The second form takes the ownership of the cell, must be newed object.
MGComplex::MGComplex(const MGCellNB& cell){
	MGCellNB* cnb=cell.clone();
	m_pcells.push_back(cnb);
}
MGComplex::MGComplex(MGCellNB* cell):m_pcells(1,cell){;}

//Fundamental constructor.
//Construct from a list of pcells.
//This constructor takes the ownership of all pcells in pcells.
MGComplex::MGComplex(std::list<MGCellNB*>& pcells):m_pcells(pcells){;}

//Copy constructor.
MGComplex::MGComplex(
	const MGComplex& complex)	//Original complex.
:MGTopology(complex),m_box(complex.m_box){
	copy_all_elements(complex);
}

MGComplex::MGComplex(
	const MGComplex& complex,
	MGCellMap& cmap)		//cellmap to register binder association.
:MGTopology(complex),m_box(complex.m_box){
	copy_all_elements(complex,cmap);
}

//Virtual Destructor
MGComplex::~MGComplex(){
	erase_all_elements();
}

///////operator overload///////

// Complexに平行移動を行ないオブジェクトを生成する。
//Translation of the Complex
MGComplex MGComplex::operator+ (const MGVector& v) const{
	MGComplex ncomp(*this);
	ncomp+=v;
	return ncomp;
}

// Complexに平行移動を行ない自身のComplexとする。
//Translation of the Complex
MGComplex& MGComplex::operator+= (const MGVector& v){
	pcellItr cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++) **cs +=v;

	ce=m_bcells.end();
	for(cs=m_bcells.begin(); cs!=ce; cs++) **cs +=v;

	m_box +=v;
	return *this;
}

// Complexに逆方向の平行移動を行ないオブジェクトを生成する。
//Translation of the Complex
MGComplex MGComplex::operator- (const MGVector& v) const{
	MGComplex ncomp(*this);
	ncomp-=v;
	return ncomp;
}

// Complexに逆方向の平行移動を行ない自身のComplexとする。
//Translation of the Complex
MGComplex& MGComplex::operator-= (const MGVector& v){
	return operator+=(-v);
}

//Complexのスケーリングを行い，Complexを作成する。
//Scaling of the Complex by a double.
MGComplex MGComplex::operator* (double s) const{
	MGComplex ncomp(*this);
	ncomp*=s;
	return ncomp;
}

//Complexのスケーリングを行い，Complexを作成する。
//Scaling of the Complex by a double.
MGComplex operator* (double s, const MGComplex& complex){
	MGComplex ncomp(complex);
	ncomp*=s;
	return ncomp;
}

//Complexのスケーリングを行い自身のComplexとする。
//Scaling of the Complex by a double.
MGComplex& MGComplex::operator*= (double s){
	pcellItr cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++) **cs *=s;

	ce=m_bcells.end();
	for(cs=m_bcells.begin(); cs!=ce; cs++) **cs *=s;

	m_box.set_null();
	return *this;
}

// 与えられた変換でComplexの変換を行い，Complexを作成する。
//Transformation of the Complex by a matrix.
MGComplex MGComplex::operator* (const MGMatrix& mat) const{
	MGComplex ncomp(*this);
	ncomp*=mat;
	return ncomp;
}

// 与えられた変換でComplexの変換を行い自身のComplexとする。
//Transformation of the Complex by a matrix.
MGComplex& MGComplex::operator*= (const MGMatrix& mat){
	pcellItr cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++) **cs *=mat;

	ce=m_bcells.end();
	for(cs=m_bcells.begin(); cs!=ce; cs++) **cs *=mat;

	m_box.set_null();
	return *this;
}

// 与えられた変換によってトランスフォームをおこないComplexを生成する。
//Transformation of the Complex by a MGTransf.
MGComplex MGComplex::operator* (const MGTransf& tr) const{
	MGComplex ncomp(*this);
	ncomp*=tr;
	return ncomp;
}

// 与えられた変換によってトランスフォームをおこない自身のComplexにする。
//Transformation of the Complex by a MGTransf.
MGComplex& MGComplex::operator*= (const MGTransf& tr){
	pcellItr cs, ce, csave;
	ce=m_pcells.end();
	for(cs=m_pcells.begin(); cs!=ce; cs++) **cs *=tr;

	ce=m_bcells.end();
	for(cs=m_bcells.begin(); cs!=ce; cs++) **cs *=tr;

	m_box.set_null();
	return *this;
}

MGComplex& MGComplex::operator=(const MGComplex& gel2){
	set_complex(gel2);
	return *this;
}
MGComplex& MGComplex::operator=(const MGGel& gel2){
	const MGComplex* gel2_is_this=dynamic_cast<const MGComplex*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Assighment.
//When the leaf object of this and topo2 are not equal, this assignment
//does nothing.
MGComplex& MGComplex::set_complex(const MGComplex& comp2){
	MGTopology::operator =(comp2);
	erase_all_elements();
	m_box=comp2.m_box;
	copy_all_elements(comp2);
	return *this;
}

//Cell comparison.
bool MGComplex::operator<(const MGComplex& cell2)const{
	if(this==&cell2)
		return false;

	const MGCellNB* c1=star();
	const MGCellNB* c2=cell2.star();
	return number_of_pcells()<cell2.number_of_pcells();
}
bool MGComplex::operator<(const MGGel& gel2)const{
	const MGComplex* gel2_is_this=dynamic_cast<const MGComplex*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//Member Function

//Append a binder cell to the end of bcell sequence.
MGComplex::cellItr MGComplex::append_bcell(MGCellNB* cell) const{
	MGComplex* comptemp=const_cast<MGComplex*>(this);
	//We can do const_cast, since bcell does not affect
	//the shape of the comlex.
	return comptemp->add_cell(cell,false,m_bcells.end());
}

//Insert a Cell(binder or parameter) before the position loc.
//loc is the iterator of m_pcells or m_bcells according to pcell.
MGComplex::cellItr MGComplex::add_cell(
	MGCellNB* cell,	//Cell to insert
	bool pcell,		//indicates if input cell is parameter or binder cell.
	MGComplex::cellItr loc
					//Iterator that indicates insert position. Insert cell befoer 
					//the iterator loc.
){
	MGComplex::cellItr newloc=loc;
	MGComplex* compold=cell->parent_complex();
	if(compold!=this){
		//1. Free from the old parent(if exist).
		if(compold){
			if(pcell){
				cell->free_from_parent();
				compold->m_box.set_null();
			}else
				cell->free_from_parent();
		}

		//2. Add own cell.
		container_type* list;
		if(pcell)
			list=&m_pcells;
		else
			list=&m_bcells;
		newloc=list->insert(loc, cell);
		cell->m_parent_complex=this;
		if(pcell){
			m_box.set_null();
			MGCellNB* bcell=cell->binder();
			if(bcell){//If cell is pcell and had a binder, then add the binder
					  //to the complex of parent cell of this complex.	
				MGCellNB* starCell=star();
				if(starCell){
					MGComplex* PComp=starCell->parent_complex();
					if(PComp) PComp->add_cell(bcell,false,PComp->m_bcells.end());
				}
			}
		}

		//3. Append binders of the boundaries.
		if(cell->manifold_dimension()>=1){
			//Add binders of the boundaries.
			std::vector<MGCellNB*> cvec;
			cell->get_all_boundary_binders(cvec);
			size_t n=cvec.size();
			for(size_t i=0; i<n; i++)
				if(cvec[i]&&(cvec[i]->parent_complex()!=this))
					add_cell(cvec[i],false,m_bcells.end());
		}
	}
	return newloc;
}

//Append a PCell to the end of pcell sequence.
MGComplex::cellItr MGComplex::append_pcell(MGCellNB* cell){
	assert(!cell || (cell&&cell->extent()));//Extent must be attached.
	return add_cell(cell,true, m_pcells.end());
}

//Cehck if bcell exist.
bool MGComplex::bcell_exist()const{
	return m_bcells.begin()!=m_bcells.end();
}

//Obtain i-the pcell. MGCellNB version.
const MGCellNB* MGComplex::bcelli(size_t i) const{
	const_bcellItr bcellp=bcell_begin();
	std::advance(bcellp,i);
	return *bcellp;
}
MGCellNB* MGComplex::bcelli(size_t i){
	pcellItr bcellp=pcell_begin();
	std::advance(bcellp,i);
	return *bcellp;
}

//Obtain i-th pcell iterator.
MGComplex::const_bcellItr MGComplex::bcellIterator(size_t i) const{
	const_bcellItr bcellp=bcell_begin();
	std::advance(bcellp,i);
	return bcellp;
}
MGComplex::bcellItr MGComplex::bcellIterator(size_t i){
	pcellItr bcellp=bcell_begin();
	std::advance(bcellp,i);
	return bcellp;
}

//Obtain the binder of i-the pcell, may be null.
MGCellNB* MGComplex::binder(size_t i) const{
	return pcelli(i)->binder();
}

//Obtain binders of all the pcells of the boundary.
//i-th binder of the fucntion's returned value is the binder
//of i-th pcell of the boundary, and may be null.
std::vector<MGCellNB*> MGComplex::binders() const{
	size_t n=number_of_pcells();
	std::vector<MGCellNB*> bdrs(n);
	const_pcellItr cell=pcell_begin();
	for(size_t i=0; i<n; i++){
		const MGCellNB* pcell=*cell++;
		if(pcell) bdrs[i]=pcell->binder();
		else bdrs[i]=0;
	}
	return bdrs;
}

//Boundaryに平行移動を行ない自身のBoundaryとする。
//Translation of the boundary
void MGComplex::binderTr(const MGVector& v){
	std::vector<MGCellNB*> bndrs=binders();
	std::vector<MGCellNB*>::iterator bs=bndrs.begin(), be=bndrs.end();
	while(bs!=be){
		MGCellNB* bndr=*bs++;
		if(bndr) *bndr+=v;
	}
}

//Boundaryのスケーリングを行い自身のBoundaryとする。
//Scaling of the boundary by a double.
void MGComplex::binderTr(double s){
	std::vector<MGCellNB*> bndrs=binders();
	std::vector<MGCellNB*>::iterator bs=bndrs.begin(), be=bndrs.end();
	while(bs!=be){
		MGCellNB* bndr=*bs++;
		if(bndr) *bndr *=s;
	}
}

// 与えられた変換でBoundaryの変換を行い自身のBoundaryとする。
//Transformation of the boundary by a matrix.
void MGComplex::binderTr(const MGMatrix& mat){
	std::vector<MGCellNB*> bndrs=binders();
	std::vector<MGCellNB*>::iterator bs=bndrs.begin(), be=bndrs.end();
	while(bs!=be){
		MGCellNB* bndr=*bs++;
		if(bndr) *bndr *=mat;
	}
}

// 与えられた変換によってトランスフォームをおこない自身のBoundaryにする。
//Transformation of the boundary by a MGTransf.
void MGComplex::binderTr(const MGTransf& tr){
	std::vector<MGCellNB*> bndrs=binders();
	std::vector<MGCellNB*>::iterator bs=bndrs.begin(), be=bndrs.end();
	while(bs!=be){
		MGCellNB* bndr=*bs++;
		if(bndr) *bndr *=tr;
	}
}
//Return the box of this cell.
const MGBox& MGComplex::box() const{
	if(m_box.is_null()) compute_box();
	return m_box;
}

//Compute barycenter of all the vertex(binder cell of 0D manifold
//dimension).
MGPosition MGComplex::center() const{
	MGPosition cntr;
	size_t num=0;
	const_bcellItr binder=m_bcells.begin(), bce=m_bcells.end();
	for(;binder!=bce;binder++){

	if((*binder)->manifold_dimension()==0){//When vertex.
		const MGGeometry* geo=(*binder)->extent();
		if(geo){//If binder had geometry, we use it.
			const MGPoint* pnt=dynamic_cast<const MGPoint*>(geo);
			cntr+=pnt->position(); num++;
		}else{	//If binder did not have geometry, we use partner's.
			const MGCellBase* cb=(*binder)->m_partners[0];
			const MGPVertex* pv=dynamic_cast<const MGPVertex*>(cb);
			const MGGeometry* geo1d=cb->star()->extent();
			const MGCurve* crv=dynamic_cast<const MGCurve*>(geo1d);
			if(crv){cntr+=crv->eval(pv->t()); num++;}
		}
	}

	}
	if(num) cntr/=double(num);
	return cntr;
}

//Construct new object by copying to newed area.
//User must delete this copied object by "delete".
MGComplex* MGComplex::clone() const{return new MGComplex(*this);}

//Copy all the pcells of comp into this.
void MGComplex::copy_all_elements(const MGComplex& comp){
	MGCellMap cmap;
	const_pcellItr el=comp.pcell_begin(),el_end=comp.pcell_end();
	for(; el!=el_end; el++){
		MGCellBase* cb=(*el)->clone_pcell_with_bcell(cmap);
		append_pcell(static_cast<MGCellNB*>(cb));//clone of CellNB is CellNB.
	}
	m_box.set_null();
}

//Copy all the pcells of comp into this.
void MGComplex::copy_all_elements(
	const MGComplex& comp,
	MGCellMap& cmap_parent)	//cellmap to register binder association.							
{
	MGCellMap cmap;
	const_pcellItr el=comp.pcell_begin(),el_end=comp.pcell_end();
	for(; el!=el_end; el++){
		MGCellBase* cb=(*el)->clone_pcell_with_bcell(cmap,&cmap_parent);
		append_pcell(static_cast<MGCellNB*>(cb));//clone of CellNB is CellNB.
	}
	m_box.set_null();
}
//Copy all pcells of comp into this,
//but does not copy binders of the pcells.
void MGComplex::copy_without_binders(const MGComplex& comp){
	MGCellMap cmap;
	const_pcellItr el=comp.pcell_begin(),el_end=comp.pcell_end();
	for(; el!=el_end; el++){
		MGCellBase* cb=(*el)->clone(cmap);
		append_pcell(static_cast<MGCellNB*>(cb));//clone of CellNB is CellNB.
	}
	m_box.set_null();
}

//Erase all elements.
//erase_all_elements does free and destruct all the elements.
//erase_all_elements does not maintain m_box.
void MGComplex::erase_all_elements(){
	pcellItr cs, ce, csave;
	cs=m_pcells.begin(); ce=m_pcells.end();
	while(cs!=ce){
		csave=cs; csave++;
		delete *cs;
		cs=csave;
	}

	bcellItr bs, be, bsave;
	bs=m_bcells.begin(); be=m_bcells.end();
	while(bs!=be){
		bsave=bs; bsave++;
		delete *bs;
		bs=bsave;
	}

	m_bcells.clear();
	m_pcells.clear();
}

//Erase bcell element if it is registered in bcell sequence of this complex.
//Destruct the bcell if found in this complex.
void MGComplex::erase_binder(const MGCellNB* bcell){
	bcellItr itr=m_bcells.begin(), itre=m_bcells.end();
	for(; itr!=itre; itr++){
		if(bcell==(*itr)){
			delete *itr;
			return;
		}
	}
}

//Erase first pcell element, including the binder cell
//that will be free when the pcell is erased.
void MGComplex::erase_first_pcell(){
	if(pcell_exist())
		erase_pcell(m_pcells.front());
}

//Erase last pcell element, including the binder cell
//that will be free when the pcell is erased.
void MGComplex::erase_last_pcell(){
	if(pcell_exist())
		erase_pcell(m_pcells.back());
}

//Erase the pcell element, including binder cells
//that will be free when the last pcell is erased.
void MGComplex::erase_pcell(MGCellNB* pcell){
	delete pcell;
	m_box.set_null();
}

//Get fisrt pcell pointer.
const MGCellNB* MGComplex::first_pcell() const{
	return m_pcells.front();
}
MGCellNB* MGComplex::first_pcell(){
	return m_pcells.front();
}

//Test if this complex includes the MGCellNB cell as a contituent.
//Returns true if cell is included in this complex.
bool MGComplex::includes(const MGCellNB* cell)const{
	const_pcellItr ps=m_pcells.begin(), pe=m_pcells.end();
	for(; ps!=pe; ps++){	
		if((*ps)==cell)
			return true;
	}
	return false;
}

//Get last pcell pointer.
const MGCellNB* MGComplex::last_pcell() const{
	return m_pcells.back();
}
MGCellNB* MGComplex::last_pcell(){
	return m_pcells.back();
}

//Get manifold dimension.
unsigned MGComplex::manifold_dimension() const{return 0;}

// Output virtual function.
std::ostream& MGComplex::out(std::ostream& ostrm) const{
	ostrm<<"Complex="<<this;
	MGComplex::const_pcellItr ps=m_pcells.begin();
	MGComplex::const_pcellItr pe=m_pcells.end();
	ostrm<<"::"<<"PCells="<<m_pcells.size()<<"::";
	size_t n=0;
	ostrm<<std::endl;
	for(;ps!=pe;ps++){
		ostrm<<(**ps).whoami()<<n++<<"="<<(**ps);
	}
	MGComplex::const_bcellItr bs=m_bcells.begin(), be=m_bcells.end();
	size_t m=m_bcells.size();
	ostrm<<std::endl;
	if(m){
		ostrm<<"BCells="<<m<<"::";
	}
	n=0;
	for(;bs!=be;bs++) {
		ostrm<<std::endl;
		const MGCellNB& bsn=**bs;
		ostrm<<bsn.whoami()<<n++<<"="<<bsn;
	}
	MGTopology::out(ostrm);
	return ostrm;
}

//Cehck if pcell exist.
bool MGComplex::pcell_exist()const{
	return m_pcells.begin()!=m_pcells.end();
}

//Obtain i-the pcell. MGCellNB version.
const MGCellNB* MGComplex::pcelli(size_t i) const{
	const_pcellItr pcellp=pcell_begin();
	std::advance(pcellp,i);
	return *pcellp;
}
MGCellNB* MGComplex::pcelli(size_t i){
	pcellItr pcellp=pcell_begin();
	std::advance(pcellp,i);
	return *pcellp;
}

//Obtain i-th pcell iterator.
MGComplex::const_pcellItr MGComplex::pcellIterator(size_t i) const{
	const_pcellItr pcellp=pcell_begin();
	std::advance(pcellp,i);
	return pcellp;
}
MGComplex::pcellItr MGComplex::pcellIterator(size_t i){
	pcellItr pcellp=pcell_begin();
	std::advance(pcellp,i);
	return pcellp;
}

//Obtain pcells that constitute the boundary.
//Let pcellvec[.] be pcells' return value and bindervec[.] be
//binders's return value. Then pcellvec[i] corresponds to bindervec[i].
//bindervec[i] is binder cell of i-th pcell element of the boundary.
std::vector<MGCellNB*> MGComplex::pcells(){
	size_t n=number_of_pcells();
	std::vector<MGCellNB*> pcv(n);
	const_pcellItr pcellp, pcellend;
	pcellp=pcell_begin();
	pcellend=pcell_end();
	size_t i=0;
	for(;pcellp!=pcellend; pcellp++) pcv[i++]=*pcellp;
	return pcv;
}
std::vector<const MGCellNB*> MGComplex::pcells() const{
	size_t n=number_of_pcells();
	std::vector<const MGCellNB*> pcv(n);
	const_pcellItr pcellp, pcellend;
	pcellp=pcell_begin();
	pcellend=pcell_end();
	size_t i=0;
	for(;pcellp!=pcellend; pcellp++) pcv[i++]=(*pcellp);
	return pcv;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGComplex::pick_closest(const MGStraight& sl)const{
	const_pcellItr ps=m_pcells.begin(), pe=m_pcells.end();
	MGPosition prm;
	if(ps==pe) return prm;
	prm=(**ps).pick_closest(sl);
	double t=sl.closest((**ps).extent()->evaluate(prm));
	for(ps++; ps!=pe; ps++){
		MGPosition prm2=(**ps).pick_closest(sl);
		double t2=sl.closest((**ps).extent()->evaluate(prm2));		
		if(t2>t){
			t=t2; prm=prm2;
		}
	}
	return prm;
}

//Prepend a PCell to the end of pcell sequence.
MGComplex::cellItr MGComplex::prepend_pcell(MGCellNB* cell){
	assert(!cell || (cell&&cell->extent()));//Extent must be attached.
	return add_cell(cell,true, m_pcells.begin());
}

//Compute the box from the scratch.
void MGComplex::compute_box() const{
	m_box.set_null();
	const_pcellItr ps=m_pcells.begin(), pe=m_pcells.end();
	for(;ps!=pe;ps++) m_box |=(*ps)->box();
}

//Get star cell pointer if this complex is a boundary of a cell.
//Else, null will be returned.
const MGCellNB* MGComplex::star() const{
	const MGBoundary* bndry=dynamic_cast<const MGBoundary*>(this);
	if(bndry) return bndry->star();	//If this is a boundary.
	else return 0;
}
MGCellNB* MGComplex::star(){
	MGBoundary* bndry=dynamic_cast<MGBoundary*>(this);
	if(bndry) return bndry->star();	//If this is a boundary.
	else return 0;
}
