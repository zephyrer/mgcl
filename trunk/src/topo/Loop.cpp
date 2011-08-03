/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/CompositeCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/SurfCurve.h"
#include "mg/FSurface.h"
#include "mg/Surface.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/Position_list.h"
#include "mg/Tolerance.h"
#include "topo/LEPoint.h"
#include "topo/LPoint.h"
#include "topo/Edge.h"
#include "topo/LCisect_vector.h"
#include "topo/LLisect_vector.h"
#include "topo/Loop.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGLoop Class.
//MGLoop is a boundary of a face, a boundary of 2D manifold cell.
//MGLoop always accepts parameter space curve and world space curve
//of a boundary curve, and constructs a boundary of a face from the
//two types of curves.

//Input curves direction indicate which part of the face will be target
//part after trimed by the boundary. In 2D space (u,v) of the parameter
//space, LEFT side of the parameter curve along the curve's direction
//is the target part of face.

///////Constructor////////

//Void constructor.
MGLoop::MGLoop()
:m_kind(UNDEFINED),m_area(0.),	m_perim_num_s(0), m_perim_num_e(0){;}

//Copy constructor.
MGLoop::MGLoop(const MGLoop& loop2):MGBoundary(loop2), m_kind(UNDEFINED){;}

//Fundamental constructor.
//Construct from boundary complex(i.e. MGLoop).
//This constructor takes the ownership of MGCell* in boundary.
MGLoop::MGLoop(
	std::list<MGCellNB*> boundaries)
		//Boundary data of the super class MGBoundary.
:MGBoundary(boundaries)
,m_kind(UNDEFINED),m_area(0.),m_perim_num_s(0), m_perim_num_e(0){;}

//Construct a loop of one edge.
MGLoop::MGLoop(MGEdge* edge)
:m_kind(UNDEFINED),m_area(0.),	m_perim_num_s(0), m_perim_num_e(0){
	append_pcell(edge);
}

//Construct a Loop of one edge of one curve cell.
//param_curve is parameter space representation of the face
//of which this loop will be a boundary, will make parameter cell.
//world_curve is world coordinate representation of the face
//of which this loop will be a boundary, will make binder cell.
//range1 is parameter range of the curve param_curve,
//range2 is parameter range of the curve world_curve.
//When range1,2 are not specified, the start and the end of the curve are
//treated as their range.
//***param_curve and world_curve must have the same direction.
MGLoop::MGLoop(const MGCurve& param_curve, const MGCurve& world_curve)
:m_kind(UNDEFINED),m_area(0.), m_perim_num_s(0), m_perim_num_e(0){
	std::auto_ptr<MGCurve> pcurve(param_curve.clone());
	std::auto_ptr<MGCurve> wcurve(world_curve.clone());
	make_loop(pcurve,wcurve);
}
MGLoop::MGLoop(
	const MGCurve& param_curve, const MGInterval& range1,
	const MGCurve& world_curve, const MGInterval& range2)
:m_kind(UNDEFINED),m_area(0.),	m_perim_num_s(0), m_perim_num_e(0){
	std::auto_ptr<MGCurve> pcurve(param_curve.clone());
	pcurve->limit(range1);
	std::auto_ptr<MGCurve> wcurve(world_curve.clone());
	wcurve->limit(range2);
	make_loop(pcurve,wcurve);
}
MGLoop::MGLoop(
	std::auto_ptr<MGCurve>& param_curve,
	std::auto_ptr<MGCurve>& world_curve)
:m_kind(UNDEFINED),m_area(0.), m_perim_num_s(0), m_perim_num_e(0){
	make_loop(param_curve,world_curve);
}

//Copy constructor with binder cell map.
MGLoop::MGLoop(
	const MGLoop& loop,		//original Loop.
	MGCellMap& cmap)	//cellmap to register binder association.
//Binder cells of the pcells in loop will be registered in cmap.
:MGBoundary(loop, cmap),m_area(loop.m_area),m_kind(loop.m_kind)
,m_perim_num_e(loop.m_perim_num_e), m_perim_num_s(loop.m_perim_num_s){;}

//Make loop
void MGLoop::make_loop(
	std::auto_ptr<MGCurve>& param_curve,//parameter curve
	std::auto_ptr<MGCurve>& world_curve	//World curve		
){
	MGCurve* pcurve=param_curve.release();
	MGCurve* wcurve=world_curve.release();
	MGCompositeCurve* pcurveC=dynamic_cast<MGCompositeCurve*>(pcurve);
	MGCompositeCurve* wcurveC=dynamic_cast<MGCompositeCurve*>(wcurve);
	size_t npcurveC;
	if(pcurveC){
		npcurveC=pcurveC->number_of_curves();
		if(npcurveC==1){
			//transform pcurveC to an ordinary curve.
			//(which will be stored in pcurve).
			pcurve=pcurveC->release_front();
			delete pcurveC;
			pcurveC=0;
			if(wcurveC){
				size_t nwcurveC=wcurveC->number_of_curves();
				if(nwcurveC==1){
					wcurve=wcurveC->release_front();
					delete wcurveC;
				}
			}
		}
	}

	if(!pcurveC){
	//Case that param_curve is not MGCmpositeCurve.
		//Generate loop complex.
		MGEdge* e=new MGEdge(pcurve);
		e->set_binder_edge(wcurve);
		append_pcell(e);
		return;
	}

	//Case that param_curve is MGCmpositeCurve.
	//Generate MGEdge for each component curve of pcurveC.
	bool use_wcurveC=false;//Flag to indicate that wcurveC[i] be used
		//for pcurveC[i]'s binder.
	if(wcurveC){
		if(npcurveC==wcurveC->number_of_curves())
			use_wcurveC=true;
	}
	for(size_t i=0; i<npcurveC; i++){
		MGEdge* e=new MGEdge(pcurveC->release_front());
		append(e);
		if(use_wcurveC){
			MGCurve* wcurvei=wcurveC->release_front();
			e->set_binder_edge(wcurvei);
		}
	}
	delete pcurveC;
	delete wcurveC;
}

///////Destructor///////

///////operator overload///////

//Assignment.
//When the leaf object of this and bnd2 are not equal, this assignment
//does nothing.
MGLoop& MGLoop::operator=(const MGLoop& loop2){
	if(this==&loop2)
		return *this;

	MGBoundary::operator=(loop2);
	m_kind=UNDEFINED;
	m_area=0.;
	return *this;
}
MGLoop& MGLoop::operator=(const MGGel& gel2){
	const MGLoop* gel2_is_this=dynamic_cast<const MGLoop*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Object transformation.
MGLoop& MGLoop::operator+=(const MGVector& v){
	MGBoundary::operator+=(v);
	m_kind=UNDEFINED;
	m_area=0.;
	return *this;
}
MGLoop& MGLoop::operator-=(const MGVector& v){
	MGBoundary::operator-=(v);
	m_kind=UNDEFINED;
	m_area=0.;
	return *this;
}
MGLoop& MGLoop::operator*=(double scale){
	MGBoundary::operator*=(scale);
	m_kind=UNDEFINED;
	m_area=0.;
	return *this;
}
MGLoop& MGLoop::operator*=(const MGMatrix& mat){
	MGBoundary::operator*=(mat);
	m_kind=UNDEFINED;
	m_area=0.;
	return *this;
}
MGLoop& MGLoop::operator*=(const MGTransf& tr){
	MGBoundary::operator*=(tr);
	m_kind=UNDEFINED;
	m_area=0.;
	return *this;
}

//This operator is to sort loops in the order:
//  1. Perimeter boundary.
//  2. Outer boundary.
//  3. Inner boundary.
//  4. Inactive loop.
bool MGLoop::operator<(const MGLoop& loop2)const{
	const MGFace* f1=face();
	const MGFace* f2=loop2.face();
	if(f1!=f2)
		return MGComplex::operator<(loop2);

	size_t i1,i2,j1,j2;
	bool is_peri1=both_end_on_perimeter(i1,i2);
	bool is_peri2=loop2.both_end_on_perimeter(j1,j2);
	if(is_peri1 && !is_peri2) return true;
	else if(!is_peri1 && is_peri2) return false;
	else if(is_peri1 && is_peri2){
		//Now both are on perimeter. i1,i2,j1,j2 are valid.
		if(i2<i1) return true;
		if(j2<j1) return false;
		if(i2<j1 || i2<j2 || i1<j1 || i1<j2) return true;
		if(i2>j1 || i2>j2 || i1>j1 || i1>j2) return false;

		//Now i1==i2==j1==j2.
		MGPosition p1=start_point(), p2=loop2.start_point();
		MGPosition q1=end_point(), q2=loop2.end_point();

		size_t is_v=i1%2; double t1,t2;
		t1=p1(is_v), t2=q1(is_v); if(i1>=2){ t1*=-1.; t2*=-1.;}
		if(t2<t1) return true;
		t1=p2(is_v), t2=q2(is_v); if(i1>=2){ t1*=-1.; t2*=-1.;}
		if(t2<t1) return false;

		const MGSurface* srf=surface();
		if(srf) return srf->less_than(i1,p1,p2);
		else return true;
	}else{
		//Now both are not on perimeter.
		bool cls1=closed(), cls2=loop2.closed();
		if(!cls1 && !cls2) return true;
		else if(cls1 && !cls2) return true;
		if(!cls1 && cls2) return false;
		else{
		//Now both are closed.
			if(is_outer_boundary()) return true;//Outer boundary.
			else if(loop2.is_outer_boundary()) return false;
					//this is inner and loop2 is outer.
			return this<(&loop2);
		}
	}
}
bool MGLoop::operator<(const MGGel& gel2)const{
	const MGLoop* gel2_is_this=dynamic_cast<const MGLoop*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}
bool MGLoop::operator<(const MGComplex& gel2)const{
	const MGLoop* gel2_is_this=dynamic_cast<const MGLoop*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

///////Member Function///////

//Test if this is active boundary.
bool MGLoop::active() const{
	get_kind();
	return m_kind>0;
}

//Append edge to the end of loop.
//append connects the edge's start to the end of the loop.
void MGLoop::append(MGEdge* e){
	if(edge_exist())
		last_edge()->join(false,e);
	else
		append_pcell(e);
	m_kind=UNDEFINED;
}

//Build one edge of srf from the curve wcrv on srf and common information
//pspan and peri_num, which are a perimeter peri_num's parameter spans(psapn).
//wcrv must not be a MGCompositeCurve.
//One edge is generated and append to this loop.
void MGLoop::append_edge_from_crv(
	const MGSurface& srf,
	const MGCurve& wcrv,//curve of world coordinates on this face that may
		//coincide to a perimeter of the surface of srf.
	double& tLast,//wcrv's parameter value to start is input and the end param value
		//of the last edge generated will be output.
	double terror,//wcrv's parameter space error.
	const std::vector<double>& pspan,//common parameter value of a perimeter peri_num.
	int peri_num,//pspan and peri_num are output of getPerimeterCommon(). Refer to it.
	bool orientation_is_opposite//orientation flag of wcrv to edge. True if opposite.
){
	std::auto_ptr<MGCompositeCurve> ccrv(new MGCompositeCurve);//temporal curve to join this loop.

	size_t ncommon=pspan.size()/4;
	for(size_t i=0; i<ncommon; i++){
		size_t i4=i*4;
		double t0=pspan[i4], t1=pspan[1+i4];
		double s0=pspan[2+i4], s1=pspan[3+i4];
		if(t0>tLast+terror){
			MGCurve* uvcrv=srf.get_parameterCurve(MGTrimmedCurve(wcrv,tLast,t0));
			ccrv->connect_to_end(uvcrv);//std::cout<<(*ccrv)<<std::endl;;//********
			tLast=t0;
		}
		int periN=peri_num;
		MGPosition uv0=srf.perimeter_uv(periN,s0),uv1=srf.perimeter_uv(periN,s1);
		ccrv->connect_to_end(new MGStraight(uv1,uv0));//std::cout<<(*edge)<<std::endl;;//********
		tLast=t1;
	}

	double tEnd=wcrv.param_e();
	if(tLast < (tEnd-terror)){
		MGCurve* uvcrv=srf.get_parameterCurve(MGTrimmedCurve(wcrv,tLast,tEnd));
		ccrv->connect_to_end(uvcrv);//std::cout<<(*edge)<<std::endl;;//********
	}
	if(orientation_is_opposite)
		ccrv->negate();

	MGCurve* crv;
	size_t ncrv=ccrv->number_of_curves();
	if(!ncrv){
		MGPosition uvS, uvE;
		srf.perp_one(wcrv.eval(tLast),uvS);
		srf.perp_one(wcrv.eval(tEnd),uvE);
		crv=new MGStraight(uvE,uvS);
	}else if(ncrv==1)
		crv=ccrv->release_front();
	else
		crv=ccrv.release();
	append(new MGEdge(crv));
	tLast=tEnd;
}

//Compute curvilinear integral of the parameter space of the area
//sorrounded by the loop.
double MGLoop::area()const{
	get_kind();
	if(m_kind==OUTER_LOOP || m_kind==INNER_LOOP)
		return m_area;
	return m_area=compute_area();
}

//Test if loop is a perimeter boundary loop or not.
//If yes, gets perimeter id.
//pid_s, e are valid only when both_end_on_perimeter() is true,
//contain perimeter number of the start or end of the loop.
bool MGLoop::both_end_on_perimeter(size_t& pid_s, size_t& pid_e,const MGFSurface* srf) const{
	get_kind(srf);
	pid_s=m_perim_num_s; pid_e=m_perim_num_e;
	return m_kind==PERIMITER_LOOP;
}

//Make a clone.
MGLoop* MGLoop::clone(MGCell& parent) const{
	MGLoop* lpp=clone();
	lpp->set_parent(parent);
	return lpp;
}
MGLoop* MGLoop::clone() const{
	MGLoop* lpp=new MGLoop(*this);
	return lpp;
}

//Make a clone.
//The forms that have cmap as an argumetnt is to register binder association.
//Binders of pcells in this loop will be registered in cmap.
//Returned is pointer of newed object, must be deleted.
//When parent is specified, clone's parent is set to the parent.
MGLoop* MGLoop::clone(MGCell& parent,MGCellMap& cmap) const{
	MGLoop* lpp=clone(cmap);
	lpp->set_parent(parent);
	return lpp;
}
MGLoop* MGLoop::clone(MGCellMap& cmap) const{
	MGLoop* lpp=new MGLoop(*this,cmap);
	return lpp;
}

//Make a clone that has not binders.
MGLoop* MGLoop::clone_without_binders(MGCell& parent) const{
	MGLoop* lpp=clone_without_binders();
	lpp->set_parent(parent);
	return lpp;
}
MGLoop* MGLoop::clone_without_binders() const{
	MGLoop* lpp=new MGLoop();
	lpp->copy_boundary_without_binders(*this);
	return lpp;
}

//Test if this is closed boundary.
bool MGLoop::closed() const{
	if(!edge_exist())
		return false;

	return (first_edge()->pre_edge()==last_edge());//When closed.
}

//Compute the closest point from the point P to this loop.
//Returned is the loop's point, and distance is the length of distance
//between P and MGLEPoint.
//All the coordinate values are of parameter space of the face's surface.
MGLEPoint MGLoop::closest(const MGPosition& P, double& distance) const{
	MGLEPoint lp;
	const_cellItr pcitr=pcell_begin(),pcitre=pcell_end(),pcisave;

	double t=0.,len=-1.;
	MGTolerance::push(); MGTolerance::set_wc_zero(error());

	double lent, ttemp;
	const MGEdge* e;
	for(; pcitr!=pcitre; pcitr++){
		lent=(*pcitr)->box().distance(P);
		if(len<0. || lent<len){
			e=static_cast<const MGEdge*>(*pcitr);
			MGTrimmedCurve crv(e->trimmed_curve());
			//std::cout<<crv<<std::endl;
			crv.on(P,ttemp);
			lent=(crv.eval(ttemp)-P).len();
			if(len<0. || lent<len){t=ttemp; len=lent; pcisave=pcitr;}
		}
	}
	distance=len;
	MGTolerance::pop();
	lp.set(pcisave,t);
	return lp;
}

//Compute curvilinear integral of the loop.
double MGLoop::compute_area()const{
	double area=0.;
	const_cellItr i=pcell_begin(),ie=pcell_end();
	for(; i!=ie; i++){
		const MGEdge* edg=edge_from_iterator(i);
		const MGCurve* crv=edg->base_curve();
		if(crv){
			MGInterval rng=edg->range();
			area+=crv->curvilinear_integral
				(rng.low_point(), rng.high_point());
		}
	}
	return area;
}

//Copy loop data.
//This boundary data is cleared and loop's boundary is copied into this.
void MGLoop::copy_boundary(const MGBoundary& loop){
	assert(dynamic_cast<const MGLoop*>(&loop));
	const MGLoop* lpp=dynamic_cast<const MGLoop*>(&loop);
	MGBoundary::copy_boundary(loop);
	m_area=lpp->m_area;	m_kind=lpp->m_kind;
}

//Copy boundary data, but does not copy the binders.
void MGLoop::copy_boundary_without_binders(const MGBoundary& loop){
	assert(dynamic_cast<const MGLoop*>(&loop));
	const MGLoop* lpp=dynamic_cast<const MGLoop*>(&loop);
	MGBoundary::copy_boundary_without_binders(loop);
	m_area=lpp->m_area;	m_kind=lpp->m_kind;
}

//Obtain vector of curves(TrimmedCurve) of the loop.
//The curves are of parameter space expression.
MGPvector<MGCurve> MGLoop::curves()const{
	MGComplex::const_cellItr i=pcell_begin(),ie=pcell_end();
	size_t n=number_of_pcells(), j=0;
	MGPvector<MGCurve> crvs(n);
	for(; i!=ie; i++,j++){
		const MGEdge* edg=edge_from_iterator(i);
		crvs.reset(j,new MGTrimmedCurve(edg->trimmed_curve()));
	}
	return crvs;
}

//Obtain vector of curves(world coordinate expression) of the loop.
//Output curves are MGTrimmedCurve of edge's binders.
//When some of the edges do not have binders, they will be created.
MGPvector<MGCurve> MGLoop::curves_world()const{
	assert(surface());	//Surface expression must exist.

	MGPvector<MGCurve> crvs;
	MGComplex::const_cellItr i=pcell_begin(),ie=pcell_end();
	for(; i!=ie; i++){
		const MGEdge* eip=edge_from_iterator(i);
		MGEdge* bedg=eip->make_binder_with_curve();
		MGTrimmedCurve trim(bedg->trimmed_curve());
		MGCurve* crv=trim.clone();
		//if(!eip->equal_direction_to_binder())
		//	crv->negate();
		crvs.push_back(crv);
	}
	return crvs;
}

//Return i-th edge pointer.
MGEdge* MGLoop::edge(size_t i){
	m_kind=UNDEFINED;
	MGCellNB* cell=pcelli(i);
	return static_cast<MGEdge*>(cell);
}
const MGEdge* MGLoop::edge(size_t i) const{
	const MGCellNB* cell=pcelli(i);
	return static_cast<const MGEdge*>(cell);
}

//Get edge pointer from its iterator in MGComplex of MGBoudarynD.
MGEdge* edge_from_iterator(MGComplex::pcellItr i){
	return static_cast<MGEdge*>(*i);
}
const MGEdge* edge_from_iterator(MGComplex::const_pcellItr i){
	return static_cast<const MGEdge*>(*i);
}

//Get edge number in this loop.
//If e is not a member of this loop, 0 will  be returned.
size_t MGLoop::edge_num(const MGEdge* e)const{
	const_pcellItr cell=pcell_begin();
	const_pcellItr celle=pcell_end();
	size_t i=0;
	for(; cell!=celle; cell++, i++)
		if((*cell)==dynamic_cast<const MGCellNB*>(e)) return i;
	return 0;
}

//Return end point of this loop as MGLEPoint.
//loop must include at least one edge, or this output is undefined.
MGLEPoint MGLoop::end_LPoint() const{
	double t=0.;
	const_pcellItr cb=pcell_begin(), ce=pcell_end();
	if(cb!=ce){ce--; t=edge_from_iterator(ce)->param_e();}
	return MGLEPoint(ce,t);
}

//Return end point of this loop as MGPosition.
//The point is of parameter space.
MGPosition MGLoop::end_point() const{
	MGPosition P;
	if(edge_exist()){
		const MGEdge* e=last_edge();
		double t=e->param_e();
		P=e->base_curve()->eval(t);
	}
	return P;
}

//Get error of this loop.
//This is obtained from the box of the loop and relative zero.
double MGLoop::error()const{
	const MGCellNB* parent=star();
	double err;
	if(parent) err=parent->parameter_error();
	else{
		double err1,err2; MGBox prange=box();
		err1=prange(0).relative_error();
		err2=prange(1).relative_error();
		err=sqrt(err1*err1+err2*err2);
	}
	return err;
}

//Evaluation of the loop at the point t.
//When nderi=0, get a parameter (u,v) of the surface at the boundary point.
MGVector eval(const MGLEPoint& t, size_t nderi){
	return t.eval(nderi);
}

//Evaluation of the loop at the point t.
MGVector MGLoop::eval(const MGLPoint& t, size_t nderi)const{
	MGVector data;
	const MGCurve* crv=edge(t.edge_num())->base_curve();
	if(crv) data=crv->eval(t.param(), nderi);
	return data;
}

//Evaluation of the loop at i-th edge's parameter t.
MGVector MGLoop::eval(size_t i, double t, size_t nderi)const{
	MGVector data;
	const MGCurve* crv=edge(i)->base_curve();
	if(crv) data=crv->eval(t, nderi);
	return data;
}

//Return pointer of the face. If the loop is not a boundary of any face,
//null will be returned.
const MGFace* MGLoop::face() const{
	const MGCellNB* cell=star();
	return dynamic_cast<const MGFace*>(cell);
}
MGFace* MGLoop::face(){
	MGCellNB* cell=star();
	return dynamic_cast<MGFace*>(cell);
}

//Return edge pointer of the first edge.
const MGEdge* MGLoop::first_edge() const{
	const MGCellNB* cell=first_pcell();
	return static_cast<const MGEdge*>(cell);
}
MGEdge* MGLoop::first_edge(){
	m_kind=UNDEFINED;
	MGCellNB* cell=first_pcell();
	return static_cast<MGEdge*>(cell);
}

//Compute loop kind.
void MGLoop::compute_kind(const MGFSurface* srf)const{
	if(edge_exist()){
		const MGEdge*  pre=first_edge()->pre_edge();
		if(pre==last_edge()){//When closed.
			m_area=compute_area();
			if(m_area<0.)
				m_kind=INNER_LOOP;	//Inner boundary loop.
			else
				m_kind=OUTER_LOOP;	//Outer boundary loop.
		}else{
			if(get_perimeter_num(m_perim_num_s,m_perim_num_e,srf))
				m_kind=PERIMITER_LOOP;	//Perimeter boundary loop.
			else
				m_kind=INACTIVE;	//Inactive loop.
		}
	}else
		m_kind=UNDEFINED;
}

//Get the loop id of this loop in the star face baoundary.
//Let face=face(), then face->loop(get_loop_id_in_face())=this;
//When this does not have star face, or this is not a boundary of a face,
//-1 will be returned.
int MGLoop::get_loop_id_in_face()const{
	const MGFace* f=face();
	if(!f)
		return -1;

	size_t n=f->number_of_loops();
	for(size_t i=0; i<n; i++){
		if(f->loop(i)==this)
			return int(i);
	}
	return -1;
}

//Get loop kind.
void MGLoop::get_kind(const MGFSurface* srf) const{
	if(m_kind>=0)
		return;
	compute_kind(srf);
}

//Get perimeter number where start or end point lies on.
//Function's return value is:
// true if both end on perimeter,
// false if any of start or end point not on perimeter.
//pid_s, _e are valid only when true is returned.
bool MGLoop::get_perimeter_num(size_t& pid_s, size_t& pid_e,const MGFSurface* surf) const{
	pid_s=pid_e=0;
	const MGSurface* srf=0;
	if(surf)
		srf=surf->get_surface_pointer();
	if(!srf)
		srf=surface();
	if(!srf)
		return false;

	MGPosition p=start_point();
	double u=p.ref(0), v=p.ref(1);
	double rczero=MGTolerance::rc_zero();
	MGTolerance::set_rc_zero(rczero*50.);
	bool result;
	if(srf->on_a_perimeter(u,v,pid_s)){
		p=end_point();
		u=p.ref(0), v=p.ref(1);
		result=srf->on_a_perimeter(u,v,pid_e)!=0;
	}else result=false;
	MGTolerance::set_rc_zero(rczero);
	return result;
}

//Test if parameter value (u,v) is inside this loop or not.
//inside means inside face, that is,
//if the loop is inner, outside inner loops and
//if the loop is outer boundary loop, inside the outer boundary loop.
//This can be used for perimeter boundary loops.
//Returned is:
//  0:outside(not on the loop)
//  1:unknown
//  2:inside(not on the loop)
// otherwise:on the loop(size_t(MGEdge* of parameter edge))will be returned.
size_t MGLoop::inside(double u, double v) const{
	return inside(MGPosition(u,v));
}
size_t MGLoop::inside(const MGPosition& uv) const{
	double err=error(); err*=.2;
	size_t out;
	if(box().includes(uv))
		out=inside_test(err,uv);
	else{
		if(is_inner_boundary())
			out=2;
		else
			out=0;
	}
	return out;
}

//Compute intersections of this loop(parameter rep of the face)
//and param_curve.
MGLCisect_vector MGLoop::isect(const MGCurve& param_curve)const{
	MGLCisect_vector lcis(*this);
	if(!has_common(param_curve)) return lcis;

	double err=error();
	double errsave=MGTolerance::set_wc_zero(err);

	MGComplex::const_pcellItr
//		i=complex()->pcell_begin(),	ie=complex()->pcell_end();
		i=pcell_begin(),	ie=pcell_end();
	//Compute intersection points.
	for(; i!=ie ; i++){
		const MGEdge* e=edge_from_iterator(i);
		MGTrimmedCurve curve(e->trimmed_curve());
		MGCCisect_list list=curve.isect(param_curve);
		MGCCisect_list::CCiterator itr=list.begin(), itrend=list.end();
		for(;itr!=itrend; itr++)
			lcis.append(
				MGLCisect(
					MGLEPoint(i,(*itr).param1()),
					(*itr).param2(),
					(*itr).point()
					)
			);
	}

	MGTolerance::set_wc_zero(errsave);
	return lcis;
}

//Compute intersections of two loops.
MGLLisect_vector MGLoop::isect(const MGLoop& loop2) const{
	MGLLisect_vector vec(*this);
	if(!has_common(loop2))
		return vec;

	MGPosition P,Q; double len;
	double err=error();
	//1. Check end points of this loop.
	Q=start_point();
	MGLEPoint t=loop2.closest(Q,len);
	if(len<=err)
		vec.append(Q,start_LPoint(),t);
	Q=end_point();
	t=loop2.closest(Q,len);
	if(len<=err)
		vec.append(Q,end_LPoint(),t);
	//2. Check end points of loop2.
	Q=loop2.start_point();
	t=closest(Q,len);
	if(len<=err)
		vec.append(Q,t,loop2.start_LPoint());
	Q=loop2.end_point();
	t=closest(Q,len);
	if(len<=err)
		vec.append(Q,t,loop2.end_LPoint());

	double errsave=MGTolerance::set_wc_zero(err);
	//3. Compute normal intersection point.
	MGComplex::const_pcellItr i=loop2.pcell_begin(), ie=loop2.pcell_end();
	for(; i!=ie; i++){
		const MGEdge* edg=edge_from_iterator(i);
		MGLCisect_vector lcvec=isect(edg->trimmed_curve());
		size_t m=lcvec.entries();
		for(size_t j=0; j<m; j++){
			MGLCisect lc=lcvec.m_lcisects[j];
			vec.append(lc.uv(),lc.lp(),MGLEPoint(i,lc.t()));
		}
	}
	MGTolerance::set_wc_zero(errsave);
	return vec;
}

//Compute intersection points of 1D sub curves of the original loop.
//Parameter values of intersection points(MGLEPoint's) will be returned.
std::vector<MGLEPoint> MGLoop::isect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
) const{
	std::vector<MGLEPoint> vec;
	const MGBox& lpbx=box();
	if(lpbx[coordinate]<<f) return vec;

	double err=error();
	double tol=err*3.;
	double errsave=MGTolerance::set_wc_zero(err);
	MGPosition_list Plist;
		//This is used to avoid the answers of multiple points at a connection.
	MGPosition_list::iterator itr;
	MGComplex::const_pcellItr ei=pcell_begin(), ee=pcell_end();
	for(; ei!=ee; ei++){
		const MGEdge* edg=edge_from_iterator(ei);
		MGTrimmedCurve crv=edg->trimmed_curve();
		MGCParam_list tlist=crv.isect_1D(f, coordinate);
		MGCParam_list::Citerator ti=tlist.begin(), te=tlist.end();
		for(;ti!=te; ti++){
			MGPosition Q=crv.eval(*ti);
			if(!Plist.in(MGBox(Q,tol),itr)){
				Plist.push_back(Q);
				vec.push_back(MGLEPoint(ei, *ti));
			}
		}
	}
	MGTolerance::set_wc_zero(errsave);
	return vec;
}

//Compute intersection points of 1D sub curves of the original loop.
//Parameter values of intersection points(MGLEPoint's) will be returned.
//This is for tessellation and intersection with perimeter boudary edges will
//be excluded.
std::vector<MGLEPoint> MGLoop::isect_1D_tess(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
) const{
	std::vector<MGLEPoint> vec;
	const MGBox& lpbx=box();
	if(lpbx[coordinate]<<f) return vec;

	const MGSurface& surf=*(face()->surface());
	double u0=surf.param_s_u(), u1=surf.param_e_u()
			,v0=surf.param_s_v(), v1=surf.param_e_v();
	double err=error();
	double tol=err*3.;
	double errsave=MGTolerance::set_wc_zero(err);
	MGPosition_list Plist;
		//This is used to avoid the answers of multiple points at a connection.
	MGPosition_list::iterator itr;
	MGComplex::const_pcellItr ei=pcell_begin(), ee=pcell_end();
	for(; ei!=ee; ei++){
		const MGEdge* edg=edge_from_iterator(ei);
		const MGCurve* crv=edg->base_curve();
		const MGStraight* sl=dynamic_cast<const MGStraight*>(crv);
		if(sl){
			MGPosition sluvs=edg->start_point(), sluve=edg->end_point();
			if(sluvs[0]==sluve[0]) if(sluvs[0]==u0 || sluvs[0]==u1) continue;
			if(sluvs[1]==sluve[1]) if(sluvs[1]==v0 || sluvs[1]==v1) continue;
		}else{
			const MGLBRep* lb=dynamic_cast<const MGLBRep*>(crv);
			if(lb){
				if(lb->order()==2 && lb->bdim()==2){
					MGPosition sluvs=edg->start_point(), sluve=edg->end_point();
					if(sluvs[0]==sluve[0]) if(sluvs[0]==u0 || sluvs[0]==u1) continue;
					if(sluvs[1]==sluve[1]) if(sluvs[1]==v0 || sluvs[1]==v1) continue;
				}
			}
		}
		const MGInterval& frng=(edg->box())[coordinate];
		if(f<frng[0] || frng[1]<f)
			continue;
		MGTrimmedCurve tcrv=edg->trimmed_curve();
		MGCParam_list tlist=tcrv.isect_1D(f, coordinate);
		MGCParam_list::Citerator ti=tlist.begin(), te=tlist.end();
		for(;ti!=te; ti++){
			MGPosition Q=tcrv.eval(*ti);
			if(!Plist.in(MGBox(Q,tol),itr)){
				Plist.push_back(Q);
				vec.push_back(MGLEPoint(ei, *ti));
			}
		}
	}
	MGTolerance::set_wc_zero(errsave);
	return vec;
}

//Compute intersections of this loop's binder edges with f.
//If some parameter edges did not have a binder, isect_binder generates them.
MGLSPoint_vector MGLoop::isect_binder(
	const MGFSurface& f
)const{
	assert(face());//This loop must be a boundary of a face.

	double error=MGTolerance::wc_zero();
	double error2=MGTolerance::line_zero()*.5;
	MGLSPoint_vector leps;
	const_pcellItr i=pcell_begin(), ie=pcell_end();//To traverse all the edges.
	for(; i!=ie; i++){
		const MGEdge& pedge=*(edge_from_iterator(i));
		//std::cout<<pedge<<std::endl;////**********
		MGEdge& bedge=*(pedge.make_binder_with_curve());
		MGTrimmedCurve crv=bedge.trimmed_curve();
		MGCSisect_list csilist;
			MGTolerance::set_wc_zero(error2);
			csilist=f.isect(crv);
			MGTolerance::set_wc_zero(error);
		size_t m=csilist.entries(); 
		for(size_t j=0; j<m; j++){
			MGCSisect cs=csilist.removeFirst();
			MGPosition uv=cs.param_surface();
			double tb=cs.param_curve();		//tb is binder edge's parameter.
			leps.append(MGLSPoint(&pedge,tb,uv[0], uv[1]));
			//std::cout<<" tb="<<tb<<", "<<bedge.eval(tb)<<std::endl;
		}
	}
	return leps;
}

//Compute intersections of this loop(parameter rep of the face)
//and param_curve.
//isect_with_endpoints() includes endpoints of both if they are close enough
//(within tolerance).
MGLCisect_vector MGLoop::isect_with_endpoints(
	const MGCurve& param_curve)const{

	MGLCisect_vector lcis(*this);
	size_t npc=number_of_pcells();
	if(npc==0)
		return lcis;

	double err=error()*8.;
	double error_sqr=err*err;
	double errsave=MGTolerance::set_wc_zero(err);

	MGPosition P,Q; double t,len; MGVector dif; MGLEPoint p1;
	//1. Check end points of the param_curve.
	MGInterval range=param_curve.param_range();
	if(range.finite_below()){
		t=range.low_point(); P=param_curve.eval(t);
		p1=closest(P,len);
		if(len<=err)
			lcis.append(p1,t,P);
	}
	if(range.finite_above()){
		t=range.high_point(); P=param_curve.eval(t);
		p1=closest(P,len);
		if(len<=err)
			lcis.append(p1,t,P);
	}

	//2. Check end points of this loop.
	Q=start_point();
	param_curve.on(Q,t); P=param_curve.eval(t);
	dif=(P-Q);
	if(dif%dif<=error_sqr)
		lcis.append(start_LPoint(),t,P);

	Q=end_point();
	param_curve.on(Q,t); P=param_curve.eval(t);
	dif=(P-Q);
	if(dif%dif<=error_sqr)
		lcis.append(end_LPoint(),t,P);

	//3. Compute normal intersection point.
	MGLCisect_vector lcis2=isect(param_curve);
	lcis.append(lcis2);

	MGTolerance::set_wc_zero(errsave);
	return lcis;
}

//Test if this loop is inactive or not.
bool MGLoop::is_inactive(const MGFSurface* srf) const{
	get_kind(srf);
	return m_kind==INACTIVE;
}

//Test if this loop is inner boundary.
//Inner boundary is:
//  (1) closed loop. (2) the direction is clockwise.
bool MGLoop::is_inner_boundary(const MGFSurface* srf) const{
	get_kind(srf);
	return m_kind==INNER_LOOP;
}

//Test if this loop is outer boundary.
//outer boundary is:
//  (1) closed loop. (2) the direction is anti-clockwise.
bool MGLoop::is_outer_boundary(const MGFSurface* srf) const{
	get_kind(srf);
	return m_kind==OUTER_LOOP;
}

//Test if this loop is perimeter boundary.
//Perimeter boundary is:
//  both ends are on the surface perimeter.
bool MGLoop::is_perimeter_boundary(const MGFSurface* srf) const{
	get_kind(srf);
	return m_kind==PERIMITER_LOOP;
}

//Join two loops.
//start indicates which end of this loop loop2 should be connected to.
//start=true: connect start of this and loop2's end.
//start=false: connect end of this and loop2's start.
void MGLoop::join(bool start, const MGLoop& loop2){
	size_t npc=loop2.number_of_pcells();
	if(start){
		for(size_t i=1; i<=npc; i++) prepend(new MGEdge(*(loop2.edge(npc-i))));
	}else{
		const_pcellItr i=loop2.pcell_begin(), ie=loop2.pcell_end();
		for(; i!=ie; i++) append(new MGEdge(*(edge_from_iterator(i))));
	}
}
void MGLoop::join(bool start, MGLoop* loop2){
	if(!loop2) return;
	if(start){
		size_t npc=loop2->number_of_pcells();
		for(size_t i=1; i<=npc; i++){
			MGEdge* e=loop2->edge(npc-i);
			e->free_from_parent();
			prepend(e);
		}
	}else{
		pcellItr i=loop2->pcell_begin(), ie=loop2->pcell_end();
		while(i!=ie){
			MGEdge* e=edge_from_iterator(i++);
			e->free_from_parent();
			append(e);
		}
	}
	delete loop2;
}
void MGLoop::join(bool start, std::auto_ptr<MGLoop>& loop2){
	join(start,loop2.release());
}

//Return edge pointer of the last edge.
const MGEdge* MGLoop::last_edge() const{
	const MGCellNB* cell=last_pcell();
	return static_cast<const MGEdge*>(cell);
}
MGEdge* MGLoop::last_edge(){
	MGCellNB* cell=last_pcell();
	m_kind=UNDEFINED;
	return static_cast<MGEdge*>(cell);
}

//Reverse the direction of the boundary.
//(Coordinate transformation is not performed.)
void MGLoop::negate(){
	MGBoundary::negate();
	m_kind=UNDEFINED;
}

//Negate the boundary according to the parent cell negation.
//That is,
//1. Transform the coordinates of the bondary cell.
//(This transfromation depends on how the parent cell is transformed
//when negate() is invoked. So, the member cells of this boundary
//are transformed by negate_transoform of the parent cell.)
//2. Reverse the direction of the parameter cells(negate each cell).
//3. Reverse the ordering of the parameter cells.
//4. Negate the binders.
void MGLoop::negate_as_boundary(
	const MGCellNB* parent
){
	MGBoundary::negate_as_boundary(parent);
	if(m_kind==PERIMITER_LOOP)
		get_perimeter_num(m_perim_num_s,m_perim_num_e);
}

//Test if end point of the loop is on perimeter of the surface.
//Returned is true when on perimeter, false if not.
//When on perimeter, perimeter number of the surface is returned in pid_e.
bool MGLoop::on_perimeter_end(size_t& pid_e,const MGFSurface* surf)const{
	const MGSurface* srf=surf->get_surface_pointer();
	if(!srf)
		srf=surface();
	if(!srf) return false;

	MGPosition p=end_point();
	double u=p.ref(0), v=p.ref(1);
	return (srf->on_a_perimeter(u,v,pid_e))!=0;
}

//Test if start point of the loop is on perimeter of the surface.
//Returned is true when on perimeter, false if not.
//When on perimeter, perimeter number of the surface is returned in pid_s.
bool MGLoop::on_perimeter_start(size_t& pid_s,const MGFSurface* surf)const{
	const MGSurface* srf=0;
	if(surf)
		srf=surf->get_surface_pointer();
	if(!srf)
		srf=surface();
	if(!srf)
		return false;

	MGPosition p=start_point();
	double u=p.ref(0), v=p.ref(1);
	return (srf->on_a_perimeter(u,v,pid_s))!=0;
}

//Test if all the edges included are on a surface perimeter.
bool MGLoop::on_surface_perimeter(const MGFace& face)const{
	const_pcellItr i=pcell_begin(), ie=pcell_end();//To traverse all the edges.
	for(; i!=ie; i++){
		const MGEdge& pedge=*(edge_from_iterator(i));
		if(pedge.on_surface_perimeter(face)) continue;
		else return false;
	}
	return true;
}

// Output virtual function.
std::ostream& MGLoop::out(std::ostream& ostrm) const{
	ostrm<<"<<Loop="<<this;
	if(is_outer_boundary())
		ostrm<<", outer_boundary=";
	else if(is_inner_boundary())
		ostrm<<", inner_boundary=";
	else if(is_perimeter_boundary())
		ostrm<<", perimeter_boundary=";
	else if(is_network())
		ostrm<<", network=";
	else if(!active())
		ostrm<<", not active=";

	size_t nedge=number_of_edges();
	ostrm<<std::endl<<"Binder Edges::"<<std::endl;
	for(size_t i=0; i<nedge; i++){
		const MGEdge* pedge=edge(i); MGEdge* bedge=pedge->binder_edge();
		ostrm<<"Edge"<<i<<"=";
		if(bedge)
			ostrm<<(*bedge);
	}
	MGBoundary::out(ostrm);
	ostrm<<"=Loop>>"<<std::endl;
	return ostrm;
}

//Prepend edge to the start of the loop.
void MGLoop::prepend(MGEdge* e){
	if(edge_exist()) first_edge()->join(true,e);
	else append_pcell(e);
	m_kind=UNDEFINED;
}

//Return start point of this loop as MGLEPoint.
MGLEPoint MGLoop::start_LPoint() const{
	double t=0.;
	const_pcellItr cb=pcell_begin(), ce=pcell_end();
	if(cb!=ce){t=edge_from_iterator(cb)->param_s();}
	return MGLEPoint(cb,t);
}

//Return start point of this loop as MGPosition.
//The point is of parameter space.
MGPosition MGLoop::start_point() const{
	MGPosition P;
	if(edge_exist()){
		const MGEdge* e=first_edge();
		P=e->eval(e->param_s());
	}
	return P;
}

//Obtain parent surface pointer.
const MGSurface* MGLoop::surface()const{
	const MGSurface* surf=0;
	const MGCellNB* cell=star();
	if(cell){
		const MGGeometry* geo=cell->extent();
		if(geo) surf=dynamic_cast<const MGSurface*>(geo);
	}
	return surf;
}

//Test if (u,v) is in, on, or out of  curves by obtaining intersetion points
//of sl with the curves. The loop has a direction and if this is
//inner, inside means outside of the closed curves.
//Returned is:
//  0:outside(not on the loop)
//  1:unknown
//  2:inside(not on the loop)
// otherwise:on the loop(size_t(MGEdge* of parameter edge))will be returned.
//When 1 is returned, try again by changing the sample straight line sl.
size_t MGLoop::in_on_curves(
	const MGPosition& uv,	//Test point
	const double error[3],
		//error[0]: Tolerance allowed. Used for the test of on curve.
		//error[1]: loop's u span length
		//error[2]: loop's v span length
	const MGPvector<MGCurve>& crvs,//Vector of curves that make loop.
		//This has the direction as a face boundary, and the in-ness is
		//determined according to the direcion.
	const MGStraight& sl)	//Sample straight line to test the inclusion.
const{
	double errsave=MGTolerance::wc_zero();	//Save error.
	const double& err=error[0];

	size_t n=crvs.size();
	double len,lenmin=0.;
	const MGVector& dir=sl.sl_direction();
	const MGCurve *crvSave=0; const MGCurve *crv;
	MGCCisect isec, isecSave;
	MGVector difmin;
	//Compute nearest intersection point from uv to crvs.
	//If no intersection found, crvSave==0 is kept.
	MGTolerance::set_wc_zero(err*.5);	//Make error smaller.
	size_t isave;
	for(size_t i=0; i<n; i++){
		crv=crvs[i];
		MGCCisect_list list=crv->isect(sl);
		size_t n_is=list.entries();
		for(size_t j=0; j<n_is; j++){
			isec=list.removeFirst();
			MGVector dif(isec.point(), uv); len=dif%dif;
			if(!crvSave || len<lenmin){
				isave=i;
				isecSave=isec; crvSave=crv;
				lenmin=len; difmin=dif;
			}
		}
	}
	MGTolerance::set_wc_zero(err);	//Set specified error.

	if(!crvSave){
		//When no intersetion found.
		MGTolerance::set_wc_zero(errsave);	//Restore original error.
		if(is_inner_boundary())
			return 2;
		else
			return 0;
	}

	size_t out;
	double rcz=MGTolerance::rc_zero();
	double uerr=error[1]*rcz, verr=error[2]*rcz;
	if(fabs(difmin[0])<=uerr && fabs(difmin[1])<=verr){
		out=size_t(edge(isave));	//When intersection on a curve.
		MGTolerance::set_wc_zero(errsave);	//Restore original error.
		return out;
	}

	const MGFace* f=face();
	if(f){
		MGVector V0=f->eval(uv), V1=f->eval(isecSave.point());
		V0-=V1;
		if(V0%V0<=errsave*errsave){
			out=size_t(edge(isave));	//When intersection on a curve.
			MGTolerance::set_wc_zero(errsave);	//Restore original error.
			return out;
		}
	}

	MGVector dif(isecSave.point(), uv);
	MGVector tangen=crvSave->eval(isecSave.param1(),1);
	MGVector nrmal=dif*tangen;
	//cout<<"nrmal="<<nrmal<<" dif="<<dif<<" tangen="
	//<<tangen<<isecSave.point()<<endl;
	double nrmalz=nrmal.ref(2);
	if(fabs(nrmalz)<err*5.){
		//When vetor dif is parallel to tangen.
		double t; crvSave->on(uv,t);
		dif=MGVector(crvSave->eval(t),uv);
		tangen=crvSave->eval(t,1);
		nrmal=dif*tangen; nrmalz=nrmal.ref(2);
		if(nrmalz>0.) out=2; else out=0;
	}else{
		if(nrmalz>0.) out=2; else out=0;
	}

	MGTolerance::set_wc_zero(errsave);	//Restore original error.
	return out;
}

const static MGStraight loopTestSL1(MGSTRAIGHT_UNLIMIT,MGVector( .866,.5),MGPosition(0.,0.));
const static MGStraight loopTestSL2(MGSTRAIGHT_UNLIMIT,MGVector(-.866,.5),MGPosition(0.,0.));
const static MGStraight loopTestSL3(MGSTRAIGHT_UNLIMIT,MGVector(.5, .866),MGPosition(0.,0.));
const static MGStraight loopTestSL4(MGSTRAIGHT_UNLIMIT,MGVector(-.5,.866),MGPosition(0.,0.));

//Test if (u,v) is in, on, or out of this loop.
//inside_test() does not use any face data of this loop. That is, this loop
//does not have to have parent face. The loop has a direction and if this is
//inner, inside means outside of the closed curves.
//This loop does not have to make a closed loop.
//returned is:
//  0:outside(not on the loop)
//  1:unknown
//  2:inside(not on the loop)
// otherwise:on the loop(size_t(MGEdge* of parameter edge))will be returned.
size_t MGLoop::inside_test(
	double err,				//Error allowed to get intersection with loop.
	const MGPosition& uv	//Point to test if inside of the loop.
)const{
	const MGBox& lpbox=box();
	double error[3]={err, lpbox[0].length().value(), lpbox[1].length().value()};

	//Obtain curves of the loop.
	MGPvector<MGCurve> crvs = curves();

	MGStraight SL1=loopTestSL1+uv;
	size_t in1=in_on_curves(uv,error,crvs,SL1);
	if(in1>2)
		return in1;

	MGStraight SL2=loopTestSL2+uv;
	size_t in2=in_on_curves(uv,error,crvs,SL2);
	if(in2>2)
		return in2;

	if(in1==in2 && in1!=1)
		return in1;

	MGStraight SL3=loopTestSL3+uv;
	size_t in3=in_on_curves(uv,error,crvs,SL3);
	if(in3>2) return in3;
	if(in1==in3 && in1!=1)
		return in1;
	if(in2==in3 && in2!=1)
		return in2;

	//Now in1=1, in2=1,  or in1!=in2!=in3;
	if(in1==1){
		if(in2==1 && in3!=1)
			return in3;
		if(in2!=1 && in3==1)
			return in2;
	}else
		if(in2==1 && in3==1)
			return in1;

	//Final try.
	MGStraight SL4=loopTestSL4+uv;
	size_t in4=in_on_curves(uv,error,crvs,SL4);
	if(in4>2)
		return in4;
	if(in1==1 && in2==1 && in3==1)
		return in4;//in this case, 1 may be returned.
	if(in4==1){
		if(in1!=1)
			return in1;
		if(in2!=1)
			return in2;
		return in3;//in this case, maybe 1 is returned.
	}
	return in4;//in this case, 1 may be returned.
}

//Compute the mid point of this loop.
//Mid point is the point of the mid point of m-th edge,
//where m=number_edges()/2.
MGPosition MGLoop::mid_point()const{
	size_t m=number_of_edges();
	if(!m)
		return MGPosition();

	m/=2;
	return edge(m)->mid_point();
}

size_t count_perimeter_boundaries(
  const std::vector<const MGLoop*>& loops
){
	size_t n=loops.size();
	size_t i=0;
	while(i<n && loops[i]->is_perimeter_boundary()) i++;
	return i;
}

//Test if (u,v) is inside the outer boundary(of std::vector<MGLoop*>& loops).
//Inside the outer boundary means that inside outer_boudary_param() or not.
//This must not be used for faces that do not have perimeter or outer boundary loop.
//Function's return value is:
//  0:outside the outer boundary(not on a loop)
//  1:unknown
//  2:inside the outer boundary(not on a loop)
// otherwise:on the outer boundary loop 
size_t inside_outer_loop(
	const std::vector<const MGLoop*>& loops,
	const MGSurface& surf,
	const MGPosition& uv
){
	double err=surf.parameter_error();
	size_t i, n=count_perimeter_boundaries(loops);
	if(n){//When there exist perimeter boundaries.
		for(size_t i=0; i<n; i++){
			const MGLoop& lp=*(loops[i]);
			if(lp.box()<<uv)
				continue;
				//If outside the box of lp, test next loop.
			size_t in=lp.inside_test(err,uv);
			if(in==0)
				return 0;
			if(in==1)
				return 2;
			if(in>2)
				return in;
		}
		return 2;
	}

	//When outer boundary.
	const MGLoop& olp=*(loops[0]);
	const MGBox& olpBox=olp.box();
	if(olpBox<<uv)
		return 0;

	const MGInterval& olpu=olpBox[0];
	const MGInterval& olpv=olpBox[1];
	if(olpu>surf.param_s_u() && olpu<surf.param_e_u()
		&& olpv>surf.param_s_v() && olpv<surf.param_e_v())
		return olp.inside_test(err,uv); 

	size_t ne=olp.number_of_edges();
	std::vector<int> pedge(ne);
	MGComplex::const_pcellItr itrstart=olp.pcell_begin(), itr,itrsave;
	size_t isave=ne;
	for(i=0, itr=itrstart; i<ne; i++,itr++){
		pedge[i]=(*(edge_from_iterator(itr))).surface_perimeter(surf);
		if(pedge[i]>=0){
			itrsave=itr;
			isave=i;
		}
	}
	if(isave>=ne)
		return olp.inside_test(err,uv);

	size_t perim_num;
	double puv[2]={uv[0], uv[1]};
	if(surf.on_a_perimeter(puv[0], puv[1], perim_num)){
		int aftperim=perim_num+1, preperim=perim_num+3;
		aftperim%=4; preperim%=4;
		if(!surf.on_the_perimeter(aftperim,puv[0],uv[1])) aftperim=-1;
		if(!surf.on_the_perimeter(preperim,puv[0],uv[1])) preperim=-1;
		double rczero=MGTolerance::rc_zero()*2.;
		for(i=0; i<ne; i++){
			int edgeperi=pedge[i];
			if(edgeperi<0)
				continue;
			if(edgeperi!=perim_num && edgeperi!=aftperim && edgeperi!=preperim)
				continue;
			const MGEdge& ei=*(olp.edge(i));
			size_t ii=edgeperi%2;
			double x0=ei.start_point()[ii], x1=ei.end_point()[ii];
			if(x0>x1){
				double save=x0; x0=x1; x1=save;
			}
			double error=(x1-x0)*rczero;
			if(x0-error<=puv[ii] && puv[ii]<=x1+error)
				return size_t(olp.edge(i));
		}
		return 0;
	}else{
		MGBox lpbx;
		for(i=0, itr=itrsave; i<ne; i++,isave++,itr++){
			if(isave==ne){isave=0; itr=itrstart;}
			const MGEdge& ei=*(edge_from_iterator(itr));
			if(pedge[isave]>=0){
				if(lpbx.is_null())
					continue;
				else if(lpbx>>uv)
					return olp.inside_test(err,uv);
				lpbx.set_null();
			}else
				lpbx|=ei.box();
		}

		if(lpbx>>uv)
			return	olp.inside_test(err,uv);
		return 2;
	}
}
