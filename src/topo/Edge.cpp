/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/Surface.h"
#include "mg/SurfCurve.h"
#include "mg/Tolerance.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Loop.h"
#include "topo/Edge.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//Implements MGEdge Class.
//MGEdge is a 1D minifold instance of MGCell.

/////// Constructor ///////

//void constructor.
MGEdge::MGEdge():m_perror(0.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
}

//Copy constructor.
MGEdge::MGEdge(const MGEdge& e, bool copy_boundary, bool no_binder)
:MGCellNB(e), m_perror(-1.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
	if(copy_boundary) copy_all_boundaries(e);
	if(no_binder)
		return;

	if(is_pcell()){//If pcell, copy binder edge.
		const MGEdge* bindr=e.binder_edge();
		if(bindr){
			MGEdge* eb=new MGEdge(*bindr);	//Copy binder.
			MGSurfCurve* scrv=dynamic_cast<MGSurfCurve*>(extent());
			if(scrv) eb->set_extent();
			set_binder(*eb);
		}
	}
}

//Fundamental constructor.
//Construct an edge from geometry of manifold dimension 1
//and the boundaries.
//The constructor takes the ownership of geo and boundaries.
MGEdge::MGEdge(
	MGGeometry* geo,
	MGPVertex* boundaries[2],
	MGCellNB* binder)
:MGCellNB(geo, binder), m_perror(-1.),m_equal_to_binder(0){
	for(size_t i=0; i<2; i++){
		m_vertex[i]=boundaries[i];
		if(m_vertex[i]) m_vertex[i]->set_edge(this);
	}
	assert(geo->manifold_dimension()==1);
}

//Make an edge of a boundary(MGBoundary1D that has active start and
//end vertex if the curve is not infinite straight line).
//The second form that input MGCurve* takes the ownership of the crv
//into the MGEdge, must not delete the object and the object must be
//newed one.
MGEdge::MGEdge(const MGCurve& crv): m_perror(-1.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
	MGInterval prange=crv.param_range();
	*this=MGEdge(crv,prange);
}
MGEdge::MGEdge(MGCurve* crv):MGCellNB(crv), m_perror(-1.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
	if(!crv) return;
	MGInterval prange=crv->param_range();
	if(prange.finite_below())
		set_start(prange.low_point());
	if(prange.finite_above())
		set_end(prange.high_point());
}

//Make an edge of a boundary(MGBoundary1D) that has active start and
//end vertex.
//range is the parameter range of crv.
//The second form that input MGCurve* takes the ownership of the crv
//into the MGEdge, must not delete the object and the object must be
//newed one.
MGEdge::MGEdge(const MGCurve& crv, const MGInterval& range)
:MGCellNB(crv), m_perror(-1.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
	MGInterval prange=crv.param_range();
	prange&=range;
	if(prange.finite_below())
		set_start(prange.low_point());
	if(prange.finite_above())
		set_end(prange.high_point());
}
MGEdge::MGEdge(MGCurve* crv, const MGInterval& range)
:MGCellNB(crv), m_perror(-1.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
	MGInterval prange=crv->param_range();
	prange&=range;
	if(prange.finite_below())
		set_start(prange.low_point());
	if(prange.finite_above())
		set_end(prange.high_point());
}

//Make an edge with a binder of a boundary
//(MGBoundary1D that has active start and end vertex).
MGEdge::MGEdge(
	const MGSurface&surf,//Parent surface of which this edge makes a boundary
	const MGCurve& pcrv, //Parameter curve of the surface surf.
	const MGInterval& prange,//param range of pcrv.
	const MGCurve& wcrv  //World coordinate curve of the surface surf.
						//wcrv will be trimmed by prange of pcrv.
):m_perror(-1.),m_equal_to_binder(0){
	m_vertex[0]=m_vertex[1]=0;
	*this=MGEdge(pcrv,prange);
	double t1=param_s(), t2=param_e(), s1,s2;
	wcrv.on(surf.eval(pcrv.eval(t1)),s1);
	wcrv.on(surf.eval(pcrv.eval(t2)),s2);
	MGEdge* binder_edge=new MGEdge(wcrv,MGInterval(MGInterval(s1),s2));
	set_binder(*binder_edge);
}

//////////Destructor//////////
MGEdge::~MGEdge(){
	if(m_vertex[0])
		delete m_vertex[0];
	if(m_vertex[1])
		delete m_vertex[1];
}

/////// operator overload///////

//Assignment.
//does not change binder and partner relation,
//does not change parent complex.
MGEdge& MGEdge::operator=(const MGEdge& e){
	if(this==&e)
		return *this;

	MGCellNB::operator=(e);
	for(size_t i=0; i<2; i++){
		if(e.m_vertex[i]){
			if(m_vertex[i]) m_vertex[i]->set_t(e.m_vertex[i]->t());
			else{
				m_vertex[i]=new MGPVertex(*(e.m_vertex[i]),this);
			}
		}else{
			delete m_vertex[i]; m_vertex[i]=0;
		}
	}
	m_box=e.m_box;
	m_perror=e.m_perror;
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator=(const MGGel& gel2){
	const MGEdge* gel2_is_this=dynamic_cast<const MGEdge*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGEdge::operator==(const MGEdge& e2)const{
	return this==&e2;
}
bool MGEdge::operator==(const MGGel& gel2)const{
	const MGEdge* gel2_is_this=dynamic_cast<const MGEdge*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGEdge::operator<(const MGEdge& gel2)const{
	if(loop()==gel2.loop())
		return MGCellNB::operator<(gel2);

	const MGCurve* crv1=curve_limitted();
	if(!crv1)
		return true;
	const MGCurve* crv2=gel2.curve_limitted();
	if(!crv2)
		return false;
	return (*crv1)<(*crv2);
}
bool MGEdge::operator<(const MGGel& gel2)const{
	const MGEdge* gel2_is_this=dynamic_cast<const MGEdge*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

// Edge に平行移動を行ないオブジェクトを生成する。
//Translation of the Edge
MGEdge MGEdge::operator+ (const MGVector& v) const{
	MGEdge nEdge(*this);
	nEdge+=v;
	return nEdge;
}

// Edgeに逆方向の平行移動を行ないオブジェクトを生成する。
//Translation of the Edge
MGEdge MGEdge::operator- (const MGVector& v) const{
	MGEdge nEdge(*this);
	nEdge-=v;
	return nEdge;
}

//Edgeのスケーリングを行い，Edgeを作成する。
//Scaling of the Edge by a double.
MGEdge MGEdge::operator* (double s) const{
	MGEdge nEdge(*this);
	nEdge*=s;
	return nEdge;
}

//Edgeのスケーリングを行い，Edgeを作成する。
//Scaling of the Edge by a double.
MGEdge operator* (double s, const MGEdge& Edge){
	MGEdge nEdge(Edge);
	nEdge*=s;
	return nEdge;
}

// 与えられた変換でEdgeの変換を行い，Edgeを作成する。
//Transformation of the Edge by a matrix.
MGEdge MGEdge::operator* (const MGMatrix& mat) const{
	MGEdge nEdge(*this);
	nEdge*=mat;
	return nEdge;
}

// 与えられた変換によってトランスフォームをおこないEdgeを生成する。
//Transformation of the Edge by a MGTransf.
MGEdge MGEdge::operator* (const MGTransf& tr) const{
	MGEdge nEdge(*this);
	nEdge*=tr;
	return nEdge;
}

//Object transformation.
MGEdge& MGEdge::operator+=(const MGVector& v){
	MGCellNB::operator+=(v);
	set_box_as_null();
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator-=(const MGVector& v){
	MGCellNB::operator-=(v);
	set_box_as_null();
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator*=(double scale){
	MGCellNB::operator*=(scale);
	set_box_as_null();
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator*=(const MGMatrix& mat){
	MGCellNB::operator*=(mat);
	set_box_as_null();
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator*=(const MGTransf& tr){
	MGCellNB::operator*=(tr);
	set_box_as_null();
	m_equal_to_binder=0;
	return *this;
}

///////Member Function///////

//Get after edge in the loop sequence.
//The aft_edge is the first neighbour edge.
const MGEdge* MGEdge::aft_edge(bool at_end,size_t* vertexID) const{
	const MGEdge* aft=0;
	size_t vid;
	if(at_end)
		vid=1;
	else
		vid=0;

	if(m_vertex[vid]){
		std::vector<const MGCellBase*> ptrs=m_vertex[vid]->partners();
		if(ptrs.size()){
			const MGPVertex* pv=static_cast<const MGPVertex*>(ptrs[0]);
			aft=static_cast<const MGEdge*>(pv->star());
			if(vertexID){
				*vertexID= (pv->is_start_vertex())? 0:1;
			}
		}
	}
	return aft;
}
MGEdge* MGEdge::aft_edge(bool at_end,size_t* vertexID){
	const MGEdge* cthis=this;
	return const_cast<MGEdge*>(cthis->aft_edge(at_end,vertexID));
}

//Obtain binder edge pointer.
MGEdge* MGEdge::binder_edge() const{
	MGCellNB* bcell=binder();
	return static_cast<MGEdge*>(bcell);
}

//Transform the boundary binders.
void MGEdge::bn_binder_tr(const MGVector& v){
	if(bn_binder_tr_necessary()){
		for(size_t i=0; i<2; i++){
			if(m_vertex[i]){
				MGCellNB* bndr=m_vertex[i]->binder();
				if(bndr)
					(*bndr)+=v;
			}
		}
	}
}
void MGEdge::bn_binder_tr(double s){
	if(bn_binder_tr_necessary()){
		for(size_t i=0; i<2; i++){
			if(m_vertex[i]){
				MGCellNB* bndr=m_vertex[i]->binder();
				if(bndr)
					(*bndr)*=s;
			}
		}
	}
}
void MGEdge::bn_binder_tr(const MGMatrix& mat){
	if(bn_binder_tr_necessary()){
		for(size_t i=0; i<2; i++){
			if(m_vertex[i]){
				MGCellNB* bndr=m_vertex[i]->binder();
				if(bndr)
					(*bndr)*=mat;
			}
		}
	}
}
void MGEdge::bn_binder_tr(const MGTransf& tr){
	if(bn_binder_tr_necessary()){
		for(size_t i=0; i<2; i++){
			if(m_vertex[i]){
				MGCellNB* bndr=m_vertex[i]->binder();
				if(bndr)
					(*bndr)*=tr;
			}
		}
	}
}

//Obtain the box of the cell.
const MGBox& MGEdge::box() const{
	if(m_box.is_null()) compute_box();
	return m_box;
}

//Set the box data as null.
void MGEdge::set_box_as_null() const{
	m_box.set_null();
}

//Obtain the center parameter value of this cell.
MGPosition MGEdge::center_param() const{
	double t=param_s();
	t+=param_e();
	t*=0.5;
	return MGPosition(1,&t);
}

//Make a clone of this edge.
//clone() does not copy the binder cell relation.
MGEdge* MGEdge::clone() const{
	return new MGEdge(*this,true,true);
}

//Generate a new MGCellBase pointer by newing the original MGCellBase.
//This is a proprietry routine of MGComplex copy.
//Copy all boundary data, (but does not copy own binder cell relation)
//and register boundary binder association of new and old into cmap.
MGEdge* MGEdge::clone(MGCellMap& cmap) const{
	MGEdge* cell=clone_without_boundaries();
	cell->copy_all_boundaries(*this,cmap);
	cell->copy_box(*this);
	cell->copy_perror(*this);
	return cell;
}

//Make a clone without boundaries of this edge.
//clone_without_boundaries() does not copy the binder cell relation.
MGEdge* MGEdge::clone_without_boundaries() const{
	return new MGEdge(*this,false,true);
}

//Make a clone of this(this is a binder), and set binder and parameter cell
//relation between the new binder and  the parameter cell e.
MGEdge* MGEdge::clone_binder(const MGCellBase& e) const{
	assert(dynamic_cast<const MGEdge*>(&e));

	MGEdge* bedge=new MGEdge(*this);
	if(dynamic_cast<MGSurfCurve*>(bedge->base_curve()))
		bedge->set_extent(0);
	e.set_binder(*bedge);
	return bedge;
}

void MGEdge::compute_box() const{
	const MGCurve* crv=base_curve();
	if(crv)
		m_box=crv->box_limitted(range());
	else
		m_box.set_null();
}

//Connect this edge  to cell2(is a MGEdge). Both edges are parameter edges of faces.
//This cell is a pcell of a boundary of a higher manifold dimension's cell A,
//and cell2 is also is a pcell of a boundary of another cell B.
//That is, this cell is a part of a boundary of cell A,
//and cell2 is a part of a boundary of cell B.
//If cell A's manifold dimension is n, then this cell's manifold dimension is n-1.
//B's manifold dimension is same as A's. and cell2's manifold dimension is n-1.
void MGEdge::connect(MGCellBase& cell2){
	MGEdge* edge2=dynamic_cast<MGEdge*>(&cell2);
	if(!edge2)
		return;

	m_equal_to_binder=edge2->m_equal_to_binder=0;
	MGCellBase::connect(cell2);
}
void MGEdge::connect(MGEdge& edge2){
	m_equal_to_binder=edge2.m_equal_to_binder=0;
	MGCellBase::connect(edge2);
}

//Connect the start(id1=0) or end(id1=1) of this to the start(id2=0) or
// the end(id2=1) of e2.
//If both edges of this and e2 are members of a complex, they must be
//the same.
void MGEdge::connect_at_id(size_t id1, MGEdge* e2, size_t id2){
	assert(id1<=1 && id2<=1);
	assert(m_vertex[id1] && e2->m_vertex[id2]);

	MGComplex* complex=parent_complex();
	MGComplex* complex2=e2->parent_complex();
	assert(!complex || !complex2 || complex==complex2);
	//If both has complexes, they must be the same.

	if(complex && !complex2){
		if(e2->is_bcell())
			complex->append_bcell(e2);
		else{
			if(id1==0)
				complex->prepend_pcell(e2);
			else
				complex->append_pcell(e2);
		}
	}else if(!complex && complex2){
		if(is_bcell())
			complex2->append_bcell(this);
		else {
			if(id1==0)
				complex2->prepend_pcell(this);
			else
				complex2->append_pcell(this);
		}
		complex=complex2;
	}
	//Now both of e1 and e2 are members of a complex, if a complex existed.

	MGBVertex* bcell1=m_vertex[id1]->binder_vertex();
	MGBVertex* bcell2=e2->m_vertex[id2]->binder_vertex();
	if(bcell1 && bcell2){//If both had binder merge them.
		 if(bcell1!=bcell2)
			 bcell1->merge_bcell(bcell2);
	}else if(bcell1){//If bcell1 had binder and bcell2 had not.
		e2->set_i_th_binder(id2,*bcell1);//Set binder relation to e2.
	}else if(bcell2){//If bcell2 had binder and bcell1 had not.
		set_i_th_binder(id1,*bcell2);//Set binder relation to this.
		bcell1=bcell2;
	}else{	//If both had not binders, generate one.
		bcell1=static_cast<MGBVertex*>(m_vertex[id1]->make_binder());
		e2->set_i_th_binder(id2,*bcell1);
	}
	if(complex)
		complex->append_bcell(bcell1);
}

//Copy boundary data of cellin into this.
void MGEdge::copy_all_boundaries(const MGCellBase& cellin){
	assert(dynamic_cast<const MGEdge*>(&cellin));

	const MGEdge* edge=dynamic_cast<const MGEdge*>(&cellin);
	for(size_t i=0; i<2; i++){
		const MGPVertex* pv=edge->m_vertex[i];
		if(m_vertex[i]){
			if(pv) *(m_vertex[i])=*pv;
			else{ delete m_vertex[i]; m_vertex[i]=0;}
		}else{
			if(edge->m_vertex[i])
				m_vertex[i]=pv->clone();
		}
		if(m_vertex[i]) m_vertex[i]->set_edge(this);
	}
	m_box.set_null();
}

//Copy all boundaries of cell into this, and binders association
//of the boundaries in the cmap.
//Binder cells of cellin will be registered in cmap.
void MGEdge::copy_all_boundaries(const MGCellBase& cellin, MGCellMap& cmap){
	assert(dynamic_cast<const MGEdge*>(&cellin));

	const MGEdge* edge=dynamic_cast<const MGEdge*>(&cellin);
	for(size_t i=0; i<2; i++){
		if(m_vertex[i]) delete m_vertex[i];
		if(edge->m_vertex[i]){
			MGCellBase* cb=edge->m_vertex[i]->clone_pcell_with_bcell(cmap, &cmap);
			m_vertex[i]=static_cast<MGPVertex*>(cb);
			m_vertex[i]->set_edge(this);
		}else m_vertex[i]=0;
	}
	m_box.set_null();
}

//Copy m_box data of cell2 into this.
void MGEdge::copy_box(const MGCellBase& cellin) const{
	assert(dynamic_cast<const MGEdge*>(&cellin));

	const MGEdge* edge=dynamic_cast<const MGEdge*>(&cellin);
	m_box=edge->m_box;
}

//Copy m_perror data of cell2 into this.
void MGEdge::copy_perror(const MGCellBase& cellin) const{
	assert(dynamic_cast<const MGEdge*>(&cellin));

	const MGEdge* edge=dynamic_cast<const MGEdge*>(&cellin);
	m_perror=edge->m_perror;
}

//Return curve pointer of this edge.
MGCurve* MGEdge::base_curve(){
	MGGeometry* geo=extent();
	return dynamic_cast<MGCurve*>(geo);
}
const MGCurve* MGEdge::base_curve() const{
	const MGGeometry* geo=extent();
	return dynamic_cast<const MGCurve*>(geo);
}

//Return curve pointer cut by start and end parameter range.
//Output is newed curve object, must be deleted.
MGCurve* MGEdge::curve_limitted() const{
	const MGCurve* crv=base_curve();
	if(crv) return crv->copy_limitted(range());
	else return 0;
}

//Disconnect the start(id=0) or end(id=1) partnership relation.
//disconnect does not free membership of the parent cell
//from its parent complex.
void MGEdge::disconnect_at_id(size_t id){
	assert(id<=1 && m_vertex[id]);

	MGBVertex* bcell=m_vertex[id]->binder_vertex();
	if(bcell){
		if(bcell->number_of_partner_members() > 1) bcell->free_partner(this);
		else delete bcell;
			//Delete the bcell since the bcell is proprietry use to this edge.
	}
}

//Test if SurfCurve of the edge has equal direction to binder edge's direction.
//Returned is true if eaual, false if not.
bool MGEdge::equal_direction_to_binder()const{
	if(m_equal_to_binder)
		return (m_equal_to_binder==1);

	const MGEdge* bedge=binder_edge();
	if(!bedge)
		return false;	//If does not have binder.
	const MGSurface* srf=star_surface();
	if(!srf)
		return false;		//If does not have star surface. 

	m_equal_to_binder=srf->equal_direction(trimmed_curve(), bedge->trimmed_curve());
	return (m_equal_to_binder==1);
}

//Evaluate the nderiv's derivative at parameter t.
//Evaluate of the curve's data.
MGVector MGEdge::eval(double t, size_t nderiv)const{
	assert(base_curve());
	return base_curve()->eval(t,nderiv);
}

//Evaluation of the star curves of the edge at the point t.
//When nderi=0, get a position of the surface at the boundary point t.
//The star curve is SurfCurve(face's surface, edge's curve).
//(The star curve has the same world coordinate with the binder curve's, but
//their direction may be opposite. The star curve has always the same direction
//as the loop.)
MGVector MGEdge::eval_star(
	double t,		//Parameter value of this parameter edge's curve.
	size_t nderi	//Order of derivative.
)const{
	assert(face()->surface());//This edge must be a boundary member of a face.
	const MGCurve* crv=base_curve();
	if(crv){
		MGSurfCurve starCurve(*(face()->surface()),*crv);
		return starCurve.eval(t, nderi);
	}
	else return MGVector();
}

//Get the star face pointer.
const MGFace* MGEdge::face() const{
	return dynamic_cast<const MGFace*>(star());
}
MGFace* MGEdge::face(){
	return dynamic_cast<MGFace*>(star());
}

//Get the 1st partner edge of this edge.
const MGEdge* MGEdge::first_partner() const{
	const MGEdge* pedge=0;
	std::vector<const MGCellBase*> cell=partners();
	if(cell.size()){
		pedge=dynamic_cast<const MGEdge*>(cell[0]);
	}
	return pedge;
}

//Free neighbourhood relationship at the end of the edge.
void MGEdge::free_end_neighbourhood(){
	if(m_vertex[1]) disconnect_at_id(1);
}

//Free neighbourhood relation at j-th boundary's i-th pcell of this cell.
//If start, j=0. If end, j=1. i must be always 0, since one boundary has
//only one cell.
void MGEdge::free_neighbourhood(size_t i, size_t j){
	assert(i==0 && j<=1);

	if(j==0) free_start_neighbourhood();
	else free_end_neighbourhood();
}

//Free neighbourhood relationship at the start of the edge.
void MGEdge::free_start_neighbourhood(){
	if(m_vertex[0]) disconnect_at_id(0);
}

//Get boundary biders of all the boundaries.
//Binders will be appended to cvec.
void MGEdge::get_all_boundary_binders(std::vector<MGCellNB*>& cvec) const{
	for(size_t i=0; i<2; i++){
		MGCellNB* bv=0;
		if(m_vertex[i]) bv=m_vertex[i]->binder();
		if(bv) cvec.push_back(bv);
	}
}

//Connect this and e2.
//If start==true, start of this edge to end of e2;
//If start==false, end of this edge to start of e2;
//e2 must be a newed object, and the ownership is transfered to the system.
void MGEdge::join(bool start, MGEdge* e2){
	size_t i=1,j=0;
	if(start){ i=0; j=1;}
	connect_at_id(i,e2,j);
}

//Return parent loop pointer.
const MGLoop* MGEdge::loop()const{
	const MGComplex*  cmp=parent_complex();
	return dynamic_cast<const MGLoop*>(cmp);
}
MGLoop* MGEdge::loop(){
	MGComplex* cmp=parent_complex();
	return dynamic_cast<MGLoop*>(cmp);
}

//Negate the boundary.
void MGEdge::negate_boundary(){
	MGPVertex* vertx=m_vertex[0];
	m_vertex[0]=m_vertex[1];
	m_vertex[1]=vertx;
	MGCurve* crv=base_curve();
	if(!crv)
		return;

	for(size_t i=0; i<2; i++){
		vertx=m_vertex[i];
		if(vertx){
			double param=vertx->t();
			vertx->set_t(crv->negate_param(param));
		}
	}
}

//Compute the mid point of this edge.
//Mid point is the point of the paramete mid=(param_s()+param_e())*.5
MGPosition MGEdge::mid_point()const{
	double t=(param_s()+param_e())*.5;
	return eval(t);
}

//Negate the direction of the cell.
void MGEdge::negate(){
	MGCellNB::negate();
	m_equal_to_binder*=-1;
}

//Obtain all the neighbours.
//The neighbours do not contain this cell except when this cell is
//connected to this cell itself(closed cell).
std::vector<const MGCellNB*> MGEdge::neighbours() const{
	std::vector<const MGCellNB*> nbrs;
	for(size_t i=0; i<2; i++){	//Loop for the two m_vertex[].
		std::vector<const MGCellBase*> prtnrs;
		if(m_vertex[i]) prtnrs=m_vertex[i]->partners();
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
	return nbrs;
}

//Get the perimeter number where this edge is on.
//If this is not on any perimeter, -1 will be returned.
int MGEdge::surface_perimeter() const{
	const MGFace* f=face(); if(!f) return -1;
	const MGSurface* sf=f->surface(); if(!sf) return -1;
	return surface_perimeter(*sf);
}
int MGEdge::surface_perimeter(const MGFace& face) const{
	const MGSurface* sf=face.surface(); if(!sf) return -1;
	return surface_perimeter(*sf);
}
int MGEdge::surface_perimeter(const MGSurface& sf) const{
	size_t pnum;
	if(sf.on_perimeter(*(base_curve()),pnum)) return pnum;
	return -1;
}

ostream& MGEdge::out(ostream& ostrm) const{
	ostrm<<"<<Edge="<<this<<"=";
	ostrm<<"m_vertex["<<m_vertex[0]<<","<<m_vertex[1]<<"],";
	if(m_vertex[0])
		ostrm<<endl<<",Vertex0="<<*(m_vertex[0]);
	if(m_vertex[1])
		ostrm<<endl<<",Vertex1="<<*(m_vertex[1]);
	ostrm<<endl<<",m_box="<<m_box<<",m_perror="<<m_perror<<",m_equal_to_binder=";
	ostrm<<m_equal_to_binder;
	MGCellNB::out(ostrm);
	ostrm<<"=Edge>>"<<endl;
	return ostrm;
}

//Obtain the parameter of the binder edge's curve that represent
//the same point as sp. sp is a parameter value of this parameter edge.
//Let S() is the star(surface) of this edge, and fp() is the curve of this cell
//which is a boundary of S(). And fb() is the binder curve of this edge.
//Then S(fp(sp))=fb(param_bcell(sp)).
//This is a parameter edge and have the binder, and the parameter sp is a parameter
//of this cell's curve. If this does not have a binder, return -1.
double MGEdge::param_bcell(double sp, const double* guess)const{
	const MGEdge* bdr=binder_edge();
	const MGSurface* srf=star_surface();
	if(!bdr || !srf) return -1.;

	const MGCurve& bcrv=*(bdr->base_curve());
	MGPosition P=srf->eval(eval(sp));
	double tb;
	if(guess){
		double tguess=*guess;
		tb=tguess;
		if(!bcrv.perp_guess(1.,0.,P,tguess,tb)) bcrv.on(P,tb);
	}else{
		bcrv.on(P,tb);// std::cout<<(P-bcrv.eval(tb)).len()<<std::endl;
	}
	return tb;
}

//This must be a parameter edge.
//Obtain the parameter of this parameter edge's curve that represent the same
//point as the binder edge's paramter tb.
//Let S() is the star(surface) of this edge, and fp() is the curve of this cell
//which is a boundary of S(). And fb() is the binder curve.
//Then S(fp(param_pcell(tb)))=fb(tb).
//This edge must have the binder edge, and the parameter tb is the parameter
//of the binder edge's curve. If this does not have a binder, return -1.
double MGEdge::param_pcell(double tb, const double* guess)const{
	const MGEdge* bdr=binder_edge();
	const MGCurve& pcrv=trimmed_curve();
	const MGSurface* srf=star_surface();
	if(!bdr || !srf) return -1.;

	return srf->param_of_pcurve(tb,bdr->trimmed_curve(),pcrv,guess);
}

//Obtain start or end parameter value of the edge.
double MGEdge::param_e()const	//End parameter value.
{
	if(active_end()) return m_vertex[1]->t();

	double t=mgInfiniteVal;
	const MGCurve* crv=base_curve();
	if(crv) t=crv->param_e();
	return t;
}
double MGEdge::param_s()const	//Start parameter value.
{
	if(active_start()) return m_vertex[0]->t();

	double t=-mgInfiniteVal;
	const MGCurve* crv=base_curve();
	if(crv) t=crv->param_s();
	return t;
}

//Return parameter space error of the cell.
double MGEdge::parameter_error()const{
	if(m_perror<=0.){
		if(m_extent) m_perror=m_extent->parameter_error();
		else return MGTolerance::wc_zero();
	}
	return m_perror;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGEdge::pick_closest(const MGStraight& sl)const{
	MGCurve* crv=curve_limitted();
	MGPosition prm=crv->pick_closest(sl);
	delete crv;
	return prm;
}

//Obtain partner edges.
//Partners represent same world's(same cell's parameter) coordinates.
//Parameter edges' partners are parameter edges.
//Binder edges' partners are binder edges.
//The partners do not include this edge except when star cell is
//connected to the star cell itself(closed only by the star cell).
std::vector<const MGEdge*> MGEdge::partner_edges() const{
	std::vector<const MGCellBase*> cells=partners();
	std::vector<const MGEdge*> pedges;
	size_t ncells=cells.size();
	for(size_t i=0; i<ncells; i++){
		const MGEdge* edgei=static_cast<const MGEdge*>(cells[i]);
		pedges.push_back(edgei);
	}
	return pedges;
}

//Get previous edge in the loop sequence.
//The pre_edge is the first neighbour edge.
const MGEdge* MGEdge::pre_edge(bool at_start) const{
	const MGEdge* pre=0;
	size_t vid;
	if(at_start)
		vid=0;
	else
		vid=1;
	if(m_vertex[vid]){
		std::vector<const MGCellBase*> ptrs=m_vertex[vid]->partners();
		size_t m=ptrs.size();
		if(m){
			pre=static_cast<const MGEdge*>(ptrs[m-1]->star());
		}
	}
	return pre;
}
MGEdge* MGEdge::pre_edge(bool at_start){
	const MGEdge* cthis=this;
	return const_cast<MGEdge*>(cthis->pre_edge(at_start));
}

//Get parameter range of the edge.
MGInterval MGEdge::range()const{
	MGEReal t0(MGINFINITE_MINUS), t1(MGINFINITE_PLUS);
	if(active_start()) t0=m_vertex[0]->t();
	if(active_end()) t1=m_vertex[1]->t();
	const MGCurve* crv=base_curve();
	if(MGREqual_base(t0,t1,crv->param_span())){
		double error=crv->param_error();error*=.5;
		double tmid=(t0.value()+t1.value())*.5;
		return MGInterval(tmid+error, tmid-error);//This is empty interval.
	}
	return MGInterval(t0,t1);
}

//Set binder cell edge to this parameter cell.
//This curve's coordinates are parameters of a face. And input wcrv's
//coordinates are world coordinate of the face.
//Parameter range of the curve is from start to end of the wcrv when no range
//is specified.
MGEdge* MGEdge::set_binder_edge(const MGCurve& wcrv)const{
	MGEdge* bndr=new MGEdge(wcrv);
	set_binder(*bndr);
	return bndr;
}
MGEdge* MGEdge::set_binder_edge(const MGCurve& wcrv, const MGInterval& range)
const{
	MGEdge* bndr=new MGEdge(wcrv,range);
	set_binder(*bndr);
	return bndr;
}
MGEdge* MGEdge::set_binder_edge(MGCurve* wcrv)const{
	MGEdge* bndr=new MGEdge(wcrv);
	set_binder(*bndr);
	return bndr;
}
MGEdge* MGEdge::set_binder_edge(MGCurve* wcrv, const MGInterval& range)const{
	MGEdge* bndr=new MGEdge(wcrv,range);
	set_binder(*bndr);
	return bndr;
}

//Set end point(boundary) data.
void MGEdge::set_end(
	double t	//Parameter value of the start point.
){
	if(m_vertex[1]){
		m_vertex[1]->free_partnership();
		m_vertex[1]->set_t(t);
	}else{
		m_vertex[1]=new MGPVertex(t,this);
	}
}

//Set start point(boundary) data.
void MGEdge::set_start(
	double t		//Parameter value of the start point.
){
	if(m_vertex[0]){
		m_vertex[0]->free_partnership();
		m_vertex[0]->set_t(t);
	}else{
		m_vertex[0]=new MGPVertex(t, this);
	}
}

//Set binder relation to m_vertex[i].
//i is the id of m_active, m_t, m_vertex.
void MGEdge::set_i_th_binder(size_t i, MGBVertex& binder)const{
	assert(i<=1 && m_vertex[i]);
	binder.add_partner(*(m_vertex[i]));
}

//Obtain star surface.
//Star cell of this must be a face. If not, return null.
//If does not have star surface, returns null.
const MGSurface* MGEdge::star_surface()const{
	const MGCellNB* str=star();
	if(str)
		return dynamic_cast<const MGSurface*>(str->extent());
	else return 0;
}

//Trim the edge at parameter t.
//When start=true, trim start, and the result is from t to end.
//When start=false, trim end, and the result is from start to t.
void MGEdge::trim(double t, bool start){
	//Set pcell parameter.
	if(start){
		free_start_neighbourhood(); set_start(t);
	}else{
		free_end_neighbourhood(); set_end(t);
	}

	//Set binder cell parameter.
	MGCurve* wcrv=world_curve();
	if(wcrv){
		//When binder curve exist, we change the parameter range.
		const MGSurface* srf=star_surface();
		if(srf){
			MGPosition uv(eval(t));
			MGPosition P(srf->eval(uv));
			double t_w; wcrv->on(P, t_w);
			MGEdge* ew=binder_edge();
			if(equal_direction_to_binder()){
				if(start) ew->set_start(t_w); else ew->set_end(t_w);
			}else{
				if(start) ew->set_end(t_w); else ew->set_start(t_w);
			}
		}else{
			//In this case, we cannot modidfy world_curve of this edge.
			//All we can do is to free(delete) the world_curve.
			free_partnership();
		}
	}
	m_box.set_null();
}

//Get trimmed curve representation of the edge.
MGTrimmedCurve MGEdge::trimmed_curve() const{
	return MGTrimmedCurve(*(base_curve()),range());
}

//Divide the edge into two parts and make a vertex at t(parameter value
//of the edge curve).
//MGEdge::make_vertex(double t);

//Return world curve pointer of this edge. That is, curve pointer
//of this edge's binder edge. May be null when no binder, or the binder
//does not have an extent.
MGCurve* MGEdge::world_curve(){
	MGCurve* crv=0;
	MGEdge* bndr=binder_edge();
	if(bndr) crv=bndr->base_curve();
	return crv;
}
const MGCurve* MGEdge::world_curve() const{
	const MGCurve* crv=0;
	MGEdge* bndr=binder_edge();
	if(bndr) crv=bndr->base_curve();
	return crv;
}
