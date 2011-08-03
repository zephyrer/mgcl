/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/Pvector.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/LEPoint.h"
#include "topo/LPoint.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

//
//Implements MGFace Class.
//MGFace is an instance of MGCell.

/////// Constructor ///////

//Fundamental constructor.
//Construct a face from geometry of manifold dimension 2
//and the boundaries.
//The constructor takes the ownership of geo and MGBoundary in boundaries.
//boundaries must be loops.
MGFace::MGFace(
	MGSurface* geo,
	std::vector<MGBoundary*>& boundaries,
	MGCell* binder
):MGCell(geo, boundaries, binder){
	assert(geo->manifold_dimension()==2);
}

//Conversion constructor from MGFSurface to MGFace.
MGFace::MGFace(const MGFSurface& surf)
:MGCell(static_cast<MGGeometry*>(surf.get_surface_pointer()->clone())),m_box_param(surf.param_range())
{
	const MGFace* f=dynamic_cast<const MGFace*>(&surf);
	if(!f)
		return;

	copy_all_boundaries(*f);
	if(f->is_pcell()){
		const MGFace* bindr=f->binder_face();
		if(bindr){
			MGFace* fb=new MGFace(*bindr);	//Copy binder.
			set_binder(*fb);
		}
	}
}

//Construct a face by copying boundaries(only parameter rep of the boundary)
//from argument boundaries.
MGFace::MGFace(
	const MGSurface& surf,
	const std::vector<MGBoundary*>& boundaries
):MGCell(surf){
	std::vector<MGBoundary*>::const_iterator i=boundaries.begin(),
		iend=boundaries.end();
	for(;i!=iend; i++){
		const MGLoop* loopi=dynamic_cast<const MGLoop*>(*i);
		if(loopi) append_boundary(loopi->clone_without_binders());
	}
}

//This form is to input a newed surface. The constructor takes the ownership
//of the surf.
MGFace::MGFace(
	MGSurface* surf,
	const std::vector<MGBoundary*>& boundaries
):MGCell(surf){
	std::vector<MGBoundary*>::const_iterator i=boundaries.begin(),
		iend=boundaries.end();
	for(;i!=iend; i++){
		const MGLoop* loopi=dynamic_cast<const MGLoop*>(*i);
		if(loopi) append_boundary(loopi->clone_without_binders());
	}
}

//Copy constructor.
MGFace::MGFace(const MGFace& face, bool copy_boundary, bool no_binder)
:MGCell(face),m_box_param(face.m_box_param){
	if(copy_boundary) copy_all_boundaries(face);
	if(no_binder)
		return;
	
	if(is_pcell()){
		const MGFace* bindr=face.binder_face();
		if(bindr){
			MGFace* fb=new MGFace(*bindr);	//Copy binder.
			set_binder(*fb);
		}
	}
}

/////// operator overload///////

// Faceに平行移動を行ないオブジェクトを生成する。
//Translation of the Face
MGFace MGFace::operator+ (const MGVector& v) const{
	MGFace face(*this);
	face+=v;
	return face;
}

// Faceに逆方向の平行移動を行ないオブジェクトを生成する。
//Translation of the Face
MGFace MGFace::operator- (const MGVector& v) const{
	MGFace face(*this);
	face-=v;
	return face;
}

//Faceのスケーリングを行い，Faceを作成する。
//Scaling of the Face by a double.
MGFace MGFace::operator* (double s) const{
	MGFace face(*this);
	face*=s;
	return face;
}

//Faceのスケーリングを行い，Faceを作成する。
//Scaling of the Face by a double.
MGFace operator* (double s, const MGFace& face){
	MGFace tface(face);
	tface*=s;
	return tface;
}

// 与えられた変換でFaceの変換を行い，Faceを作成する。
//Transformation of the Face by a matrix.
MGFace MGFace::operator* (const MGMatrix& mat) const{
	MGFace face(*this);
	face*=mat;
	return face;
}

// 与えられた変換によってトランスフォームをおこないFaceを生成する。
//Transformation of the Face by a MGTransf.
MGFace MGFace::operator* (const MGTransf& tr) const{
	MGFace face(*this);
	face*=tr;
	return face;
}

//Assignment.
//When the leaf object of this and cell2 are not equal, this assignment
//does nothing.
MGFace& MGFace::operator=(const MGFace& face){
	if(this==&face)
		return *this;

	set_cell(face);
	m_box_param=face.m_box_param;
	return *this;
}
MGFace& MGFace::operator=(const MGGel& gel2){
	const MGFace* gel2_is_this=dynamic_cast<const MGFace*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}
bool MGFace::operator<(const MGGel& gel2)const{
	const MGFace* gel2_is_this=dynamic_cast<const MGFace*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

std::ostream& MGFace::out(std::ostream& ostrm) const{
	ostrm<<"<<Face="<<this;
	ostrm<<",m_box_param="<<m_box_param;
	MGCell::out(ostrm);
	ostrm<<"=Face>>"<<endl;
	return ostrm;
}

///////Member Function///////

//Generate arrow data of the tangent along u and v and the normal
//at the parameter value (u,v) of the surface.
//data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
//data[2], [3] are two bottoms of u-tangent arrowhead.
//data[0], [4], [5], [6] are the points of v-tangent arrow.
//data[0], [7], [8], [9] are the points of v-tangent arrow.
void MGFace::arrow(double u,double v, MGPosition data[10])const{
	const MGSurface* f=surface();
	f->arrow(box(),u,v,data);
}

//Add a new loop to this face as aboundary.
//When the old inner loops that are outside the nloop will be removed from this.
//nloop can be inner or outer.
void MGFace::add_boundary(MGLoop* nloop){
	if(nloop->is_outer_boundary()){
		if(hasOuterBoundaryLoop())
			erase_boundary(size_t(0));
	}

	int nlps=number_of_loops();
	for(int i=nlps-1; i>=0; i--){
		if(!nloop->inside(loop(i)->start_point()))
			erase_boundary(i);
	}

	m_box_param.set_null();
	if(nloop->is_outer_boundary())
		MGCell::prepend_boundary(nloop);
	else
		MGCell::append_boundary(nloop);

	shrink_base_surface_to_knot();
}

//Append new one boundary to boundary vectors.
//Returned is the number of boudaries after appending.
//bound must be a newed MGLoop, and the ownership is transfered to this.
//*** append_boundary does not check validity with other loops
//(e.x. already existed loops will be outside the new boudanry bound).
//If the validity check is necessary, use add_boudanry().
size_t MGFace::append_boundary(MGBoundary* bound){
	m_box_param.set_null();
	return MGCell::append_boundary(bound);
}

//Obtain binder face pointer.
//Null when this does not have binder.
MGFace* MGFace::binder_face() const{
	MGCellNB* bcell=binder();
	return static_cast<MGFace*>(bcell);
}

//Get the clone of this as a MGFace.
//If this is MGSurface, it is converted to MGFace.
MGFace* MGSurface::clone_as_face()const{
	return new MGFace(copy_surface());
}

//Obtain all the boundary curves(world coordinates representation)
//of the face.
//That is, all of the outer boundaries and all of the inner boundaries.
MGPvector<MGCurve> MGFace::face_boundaries()const{
	MGPvector<MGCurve> crvs=outer_boundary();
	size_t n=number_of_inner_boundaries(),i;
	for(i=0; i<n; i++){
		MGPvector<MGCurve> crvsi=inner_boundary(i);
		crvs.push_back(crvsi);
	}
	return crvs;
}

//Return box of the parameter space of the face.
//After trimmed one.
const MGBox& MGFace::box_param() const{
	if(m_box_param.is_null())
		compute_box_param();
	return m_box_param;
}

//Make a clone of the cell.
MGFace* MGFace::clone() const{
	return new MGFace(*this, true,true);
}

//Generate a new MGCellBase pointer by newing the original MGCellBase.
//This is a proprietry routine of MGComplex copy.
//Copy all boundary data, (but does not copy own binder cell relation)
//and register boundary binder association of new and old into cmap.
MGFace* MGFace::clone(MGCellMap& cmap)const{
	MGFace* cell=clone_without_boundaries();
	cell->copy_all_boundaries(*this,cmap);
	cell->copy_box(*this);
	cell->copy_perror(*this);
	return cell;
}

//Make a clone without boundaries of this face.
//clone_without_boundaries() does not copy the binder cell relation.
MGFace* MGFace::clone_without_boundaries() const{
	return new MGFace(*this, false,true);
}

//Make a clone of this(this is a binder), and set binder and parameter cell
//relation between the new binder and  the parameter cell f.
MGFace* MGFace::clone_binder(const MGCellBase& f) const{
	MGFace* bface=new MGFace(*this);
	f.set_binder(*bface);
	return bface;
}

//Compute closest point from a point.
//Returned is the parameter value of the face that is closest to point.
MGPosition MGFace::closest(const MGPosition& point) const{
	const MGSurface* srf=surface();
	MGPosition P,P1; double dist,dist1; MGPosition uv;
	MGVector dif;
	//(P,uv,dist) will be the colsest pair of the face to point.
	//uv:parameter value of P, dist:distance between P and point.

	//First, compute closest candidates on the surface.
	int n=perp_one(point,uv);

	MGPvector<MGCurve> prmtr=outer_boundary();
	if(n==0){P=prmtr[0]->start_point(); uv=MGPosition();}
	else P=eval(uv);
	dif=(P-point); dist=dif%dif;

	//Second, compute closest candidates on the boundaries.
	//	2.1 Outer boundary.
	size_t pnum,i;
	pnum=prmtr.size();
	for(i=0; i<pnum; i++){
		 double t=prmtr[i]->closest(point);
		 P1=prmtr[i]->eval(t); dif=P1-point; dist1=dif%dif;
		 if(dist1<dist){dist=dist1; P=P1; uv=MGPosition();}
	}

	//	2.2 Inner boundaries.
	size_t m=number_of_inner_boundaries();
	for(size_t j=0; j<m; j++){
		MGPvector<MGCurve> prmtrj=inner_boundary(j);
		pnum=prmtrj.size();
		for(i=0; i<pnum; i++){
			 double t=prmtrj[i]->closest(point);
			 P1=prmtrj[i]->eval(t); dif=P1-point; dist1=dif%dif;
			 if(dist1<dist){dist=dist1; P=P1; uv=MGPosition();}
		}
	}

	if(uv.sdim()==0) uv=param(P);
	return uv;
}

//Compute closest point from a line to the boundary of the MGFSurface.
//Returned is the parameter value of the FSurface that is closest to point.
MGPosition MGFace::closest_on_boundary(const MGStraight& sl) const{
	MGPosition P,P1; double dist,dist1; MGPosition uv;
	MGVector dif;
	//(P,uv,dist) will be the colsest pair of the face to point.
	//uv:parameter value of P, dist:distance between P and point.

	const MGSurface& srf=*surface();
	uv=MGPosition(srf.param_s_u(), srf.param_s_v());
	P=srf.eval(uv);
	dif=(P-sl.eval(sl.closest(P))); dist=dif%dif;
	//Above are the candidate of the closest.

	//Second, compute closest candidates on the boundaries.
	MGPvector<MGCurve> bcrvs=face_boundaries();
	size_t pnum,i;
	pnum=bcrvs.size();
	for(i=0; i<pnum; i++){
		 MGPosition tt=bcrvs[i]->closest(sl);
		 P1=bcrvs[i]->eval(tt[0]);
		 dif=P1-sl.eval(tt[1]); dist1=dif%dif;
		 if(dist1<dist){dist=dist1; P=P1; uv=MGPosition();}
	}
	if(uv.sdim()==0) uv=param(P);
	return uv;
}

//compute box of the cell in m_box.
//Currently this does not compute correct box, compute m_extent box.
void MGFace::compute_box() const{
	const MGSurface* srf=surface();
	if(srf)
		m_box=srf->box_limitted(box_param());
	else
		set_box_as_null();
}

//Compute parameter range box.
void MGFace::compute_box_param() const{
	//Compute parameter range box.
	m_box_param=surface()->param_range();
	if(no_outer_boundaries())
		return;
	MGBox temp;
	if(hasOuterBoundaryLoop()){
		const MGLoop* lp=loop(size_t(0));
		temp=lp->box();
	}else{
		MGPvector<MGCurve> crvs=outer_boundary_param();
		vector<MGCurve*>::const_iterator i=crvs.begin(), iend=crvs.end();
		if(i==iend){
			return;
		}else{
			for(; i!=iend; i++) temp|=(*i)->box();
		}
	}
	m_box_param&=temp;
}

//set box as null(to set the box as initial)
void MGFace::set_box_as_null()const{
	MGCell::set_box_as_null();
	m_box_param.set_null();
}

//Test if directions of parameter curve and world curve of the face boundary
//is equal or not. This function can be used to test the pair of
//the output of outer_boundary() and outer_boundary_param(), or the pair of
//inner_boundary() and inner_boundary_param().
//Return is:
//true if equal direction, false if opposite direction.
bool MGFace::equal_direction(
	const MGPvector<MGCurve>& wcurves,
		//output of outer_boundary() or inner_boundary().
	const MGPvector<MGCurve>& pcurves,
		//output of outer_boundary_param() or inner_boundary_param().
	size_t i
		//id of the curve in wcurves and pcurves to test the direction.
		) const{
	assert(i<wcurves.size() && wcurves.size()==pcurves.size());

	MGSurfCurve scrv(*(surface()),*(pcurves[i]));
	double t1=(scrv.param_s()+scrv.param_e())*0.5;
	double t2=(wcurves[i]->param_s()+wcurves[i]->param_e())*0.5;
	return (scrv.direction(t1))%(wcurves[i]->direction(t2))>0.;
}

//Erase i-th loop.
//void MGFace::erase_boundary(size_t i){
//	bool is_in=(*(loop(i))).is_inner_boundary();
//	MGCell::erase_boundary(i);
//}

//Evaluate.
//Input parameter value is not checked if it is in_range() or not.
//Even if it is not in_range(), surface evaluation will be executed.
MGVector MGFace::eval(double u, double v,	//Face parameter value(u,v)
			  size_t ndu, size_t ndv) const//Order of derivative.
{
	const MGSurface* srf=surface();
	return srf->eval(u,v,ndu,ndv);
}
MGVector MGFace::eval(const MGPosition& uv,	//Face parameter value(u,v)
			  size_t ndu, size_t ndv) const//Order of derivative.
{
	const MGSurface* srf=surface();
	return srf->eval(uv,ndu,ndv);
}

//Get inner_aboundary loops included in the input box.
std::vector<const MGLoop*> MGFace::get_inner_boundary_loops(
	const MGBox& uvbox
)const{
	std::vector<const MGLoop*> lps;
	size_t is;
	size_t n=number_of_inner_boundaries(is);
	if(n==0) return lps;
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	for(size_t i=0; i<n; i++){
		const MGLoop* lpi=loop(i+is);
		const MGBox& bx=lpi->box();
		const MGInterval& urange=bx[0];
		if(urange<u0) continue;
		if(urange>u1) continue;
		const MGInterval& vrange=bx[1];
		if(vrange<v0) continue;
		if(vrange>v1) continue;
		lps.push_back(lpi);
	}
	return lps;
}

//Test if this face has boundary loops or not in the specified box.
//If this has one, return true.
bool MGFace::hasLoop(const MGBox& uvbox) const{
	size_t n=number_of_loops();
	if(n==0) return false;
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	for(size_t i=0; i<n; i++){
		const MGBox& bx=loop(i)->box();
		const MGInterval& urange=bx[0];
		const MGInterval& vrange=bx[1];
		if(urange>=u0 && urange<=u1 && vrange>=v0 && vrange<=v1) return true;
	}
	return false;
}

//Test if this face has the outer boundary loop instead of perimeter boundary
//loops. If this has the outer boundary loop and has not perimeter boundary
//loops, return true.
bool MGFace::hasOuterBoundaryLoop() const{
	size_t n=number_of_boundaries();
	if(!n) return false;
	if(loop(size_t(0))->is_outer_boundary()) return true;
	return false;
}

//Test if this face has perimeter boundary loops or not.
//If this has one, return true.
bool MGFace::hasPerimeterBoundaryLoop() const{
	size_t n=number_of_boundaries();
	if(!n) return false;
	if(loop(size_t(0))->is_perimeter_boundary()) return true;
	return false;
}

//Obtain i-th inner_boundary curves(world coordinates representation)
//of the face. Let the output of inner_boundary(i) be wcurves and
//of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
//to pcurves[j] one to one. Number of inner_boundary can be obtained
//by the function number_of_inner_boundaries().
MGPvector<MGCurve> MGFace::inner_boundary(size_t i)const{
	size_t start;
	size_t nb=number_of_inner_boundaries(start);
	MGPvector<MGCurve> crvs;
	if(nb && (i<nb)){
		size_t j=start+i;
		const MGLoop* lp=loop(j);
		MGPvector<MGCurve> crvsj=lp->curves_world();
		crvs.push_back(crvsj);
	}
	return crvs;
}

//Obtain i-th inner_boundary curves(parameter space representation)
//of the face. Let the output of inner_boundary(i) be wcurves and
//of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
//to pcurves[j] one to one. Number of inner_boundary can be obtained
//by the function number_of_inner_boundaries().
MGPvector<MGCurve> MGFace::inner_boundary_param(size_t i)const{
	size_t start;
	MGPvector<MGCurve> crvs;
	size_t nb=number_of_inner_boundaries(start);
	if(nb && (i<nb)){
		size_t j=start+i;
		const MGLoop* lp=loop(j);
		MGComplex::const_pcellItr i=lp->pcell_begin(),ie=lp->pcell_end();
		for(; i!=ie; i++){
			const MGEdge* edg=edge_from_iterator(i);
			crvs.push_back(new MGTrimmedCurve(edg->trimmed_curve()));
		}		
	}
	return crvs;
}

//Test if (u,v) is inside inner boundary. inside means not on the
//boundary and not included inside the face.
//If true is returned, the id of m_boundaies is returned.
//Function's return value is:
//  0:outside of all the inner loops(not on the loop)
//  1:unknown
//  2:inside an inner loop(not on the loop), and the loop id is returned in id.
// otherwise:on the loop(size_t(MGEdge* of parameter edge))will be returned, 
//			 and the loop id is returned in id.
size_t MGFace::inside_inner_boundary(
	const MGPosition& uv, size_t& id
)const{
	size_t start;
	size_t n=number_of_inner_boundaries(start);
	for(size_t i=0; i<n; i++){
		id=start+i;
		size_t in_the_loop=loop(id)->inside(uv);
		if(!in_the_loop)
			return 2;
		if(in_the_loop>2)
			return in_the_loop;
	}
	return 0;
}

//Test if (u,v) is inside the outer boundary.
//Inside the outer boundary means that inside outer_boudary_param() or not.
//***Caution***
//(1)This must not be used for faces that do not have perimeter or outer boundary
//loop.
//(2)inside_outer_boundary does not check about the inner loops.
//
//Function's return value is:
//  0:outside the outer boundary(not on a loop)
//  1:unknown
//  2:inside the outer boundary(not on a loop)
// otherwise:on the outer boundary loop(including perimeter loops)
size_t MGFace::inside_outer_boundary(
	const MGPosition& uv
)const{
	std::vector<const MGLoop*> loops;
	extract_loops(loops);
	return inside_outer_loop(loops,*(surface()),uv);
}

//Test if parameter value (u,v) is in the range of the face parameter.
bool MGFace::in_range(double u, double v)const{
	return in_range(MGPosition(u,v));
}
bool MGFace::in_range(const MGPosition& uv)const{
	if(box_param()<<uv) return false;
	else if(!number_of_boundaries())
		return true;
	else{
		size_t dummy;
		size_t in=inside_inner_boundary(uv,dummy);
		if(in>2)
			return true;
		else if(in==2)
			return false;
		else if(no_outer_boundaries())
			return true;
		else
			return inside_outer_boundary(uv)>=2;
	}
}

//Test if (u,v) is inside the face.
//Function's return value is:
//  0:outside the face.
//  1:unknown.
//  2:inside the face, not on a boundary.
//  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
//  4:(u,v) is on the outer boundary(including perimeter loops).
//  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int MGFace::in_range_with_on(
	const MGPosition& uv
)const{
	if(box_param()<<uv)
		return 0;

	size_t loop_id;
	size_t in=inside_inner_boundary(uv,loop_id);
	if(in>2)
		return -int(loop_id);
	else if(in==2)
		return 0;

	if(no_outer_boundaries()){
		const MGSurface& srf=*(surface());
		return srf.in_range_with_on(uv);
	}

	in=inside_outer_boundary(uv);
	if(in>2)
		return 4;
	else
		return in;
}

//Obtain i-th boundary loop of the face.
MGLoop* MGFace::loop(size_t i){
	MGBoundary* bnd=boundary(i);
	return static_cast<MGLoop*>(bnd);
}
const MGLoop* MGFace::loop(size_t i)const{
	const MGBoundary* bnd=boundary(i);
	return static_cast<const MGLoop*>(bnd);
}
MGLoop* MGFace::loop(iterator i){
	return static_cast<MGLoop*>(*i);
}
const MGLoop* MGFace::loop(const_iterator i)const{
	return static_cast<const MGLoop*>(*i);
}

//Make a binder cell of this parameter cell.
//Returned is the binder pointer generated by new.
//The binder has no geometry, only has binder and parameter cell relationship.
MGCellNB* MGFace::make_binder() const{
	MGCellNB* cell=binder();
	if(cell) return cell;
	MGFace* bface=new MGFace();
	set_binder(*bface);
	return bface;
}

//Make this cell's binder cell's extent expression.
//Returned is a MGGeometry pointer generated by new.
//When this cell does not have star cell, null pointer will be returned.
//make_binder_extent() only makes the expression, and does nothing to
//the topology structure.
MGGeometry* MGFace::make_binder_extent() const{
//Currently this does nothing since class Solid is not supported.
	return 0;
}

//Make sure that this has an extent expression.
//When this did not have an extent, make the extent from the partner
//member's parameter expression and the star cell.
//This must be a binder cell that has partner members that are
//boundaries. When this is not the case or this had an extent already,
//it does nothing.
void MGFace::make_extent() const{
	if(extent()) return;
	if(!is_bcell()) return;

//	const MGCellBase* partner=m_partners[0];
//	const MGFace* pface=dynamic_cast<const MGFace*>(partner);
	return;
//Currently this does nothing since class Solid is not supported.
/*	assert(pface->solid());
	assert(pface->solid()->xxxxxxx());

	MGCurve* bcrv;
	//If binder edge did not have curve representation.
	bcrv=new MGSurfCurve(*srf,pedge->trimmed_curve());
	MGEdge* thisE=const_cast<MGEdge*>(this);
	thisE->set_extent(bcrv);
	*/
}

//Obtain the i-th member partner face.
const MGFace* MGFace::member_partner_face(size_t i)const{
	return static_cast<const MGFace*>(member_partner(i));
}

//Negate the face.
void MGFace::negate(){
	MGCellNB::negate();
	if(!m_box_param.is_null())
		m_box_param=MGBox(2,m_box_param,1,0);
		//Exchange parameter range of(u,v).
}

//Test if no outer boundary except the surface perimeters.
//That is, test if the following two conditions are satisfied:
//         1. no perimeter boundaries.
//         2. no outer boundary.
bool MGFace::no_outer_boundaries()const{
	size_t n=number_of_boundaries();
	if(!n) return true;
	const MGLoop* lp_first=loop(size_t(0));
	if(lp_first->is_perimeter_boundary()) return false;
	if(lp_first->is_outer_boundary()) return false;
	return true;
}

//Get number of inner boundaries.
//Returned i is the id of the first inner boundary loop if inner boundaries
//exist.
size_t MGFace::number_of_inner_boundaries(size_t& i)const{
	size_t n=number_of_boundaries();
	i=0;
	while(i<n && loop(i)->is_perimeter_boundary()) i++;
	while(i<n && loop(i)->is_outer_boundary()) i++;
	size_t j=i;
	while(j<n && loop(j)->is_inner_boundary()) j++;
	return j-i;
}

//Compute number of active loops.
size_t MGFace::number_of_loops()const{
	size_t n=number_of_boundaries();
	size_t i=0;
	while(i<n && loop(i)->is_perimeter_boundary()) i++;
	while(i<n && loop(i)->is_outer_boundary()) i++;
	while(i<n && loop(i)->is_inner_boundary()) i++;
	return i;
}

//Get number of perimeter boundary loop.
size_t  MGFace::number_of_perimeter_boundaries()const{
	size_t n=number_of_boundaries();
	size_t i=0;
	while(i<n && loop(i)->is_perimeter_boundary()) i++;
	return i;
}

//Test if a point P is on the face.
//Returned is true if the point P is on the face.
//false if P was not on the face.
bool MGFace::on(const MGPosition& P,
		MGPosition& uv	//Parameter value of the face is returrned.
						//Even if P is not on the face, nearest point
						//parameter value will be returned.
		) const{
	uv=closest(P);
	return P==eval(uv);
}

//Test if input (u,v) is parameter value on a perimeter of the base surface.
//If u or v is on a perimeter, they will be updated to the perimeter value.
bool MGFace::on_a_perimeter(
	double& u, double& v,		//Surface parameter (u,v)
	size_t& perim_num	//if function returns true,
						//the perimete number is output.
)const{
	return surface()->on_a_perimeter(u,v,perim_num);
}

//Obtain parameter value of the face whose world coordinates are P.
MGPosition MGFace::param(const MGPosition& P)const{
	const MGSurface* srf=surface();
	MGPosition uv=srf->param(P);
	return range(uv);
}

// Return ending parameter value.
double MGFace::param_e_u()const{
	const MGBox& uvbox=box_param();
	return uvbox[0].high_point();
}
double MGFace::param_e_v()const{
	const MGBox& uvbox=box_param();
	return uvbox[1].high_point();
}

// Return starting parameter value of the base surface.
double MGFace::param_s_u()const{
	const MGBox& uvbox=box_param();
	return uvbox[0].low_point();
}
double MGFace::param_s_v()const{
	const MGBox& uvbox=box_param();
	return uvbox[1].low_point();
}

//Return the foot of the perpendicular straight line from P.
//Computation is done from the guess parameter value.
//Function's return value is whether point is obtained(true) or not(false).
bool MGFace::perp_guess(
	const MGPosition& P,		//Point
	const MGPosition& uvguess,	// guess parameter value of the shell
	MGPosition& uv				// Parameter value will be returned.
) const{
	const MGBox& pbox=box_param();	//Parameter range of the face.
	MGPosition uv0(pbox[0].low_point(), pbox[1].low_point());
	MGPosition uv1(pbox[0].high_point(), pbox[1].high_point());

	bool obtained=false;
	if(surface()->perp_guess(uv0,uv1,P,uvguess,uv)) obtained=in_range(uv);
	return obtained;
}

//Compute perpendicular points of a curve and the face, given
//guess starting paramter values.
//Function's return value is:
//   perp_guess=true if perpendicular points obtained,
//   perp_guess=false if perpendicular points not obtained,
bool MGFace::perp_guess(
	const MGCurve& curve,	//curve.
	const MGPosition& uvguess,	//Guess parameter value of the face.
	double tguess,			//Guess parameter value of the curve.
	MGPosition& uv,			//perpendicular point's parameter values of the shell
	double& t				//will be output.
) const{
	const MGBox& pbox=box_param();	//Parameter range of the face.
	MGPosition uv0(pbox[0].low_point(), pbox[1].low_point());
	MGPosition uv1(pbox[0].high_point(), pbox[1].high_point());

	bool obtained=false;
	MGPosition tuvguess(3), tuv;
	tuvguess(0)=tguess;tuvguess(1)=uvguess[0];tuvguess(2)=uvguess[1];
	if(surface()->perp_guess(uv0,uv1,curve,1.,-1.,tuvguess,tuv)){
		t=tuv[0]; uv=MGPosition(2,tuv,0,1);
		obtained=in_range(uv);
	}
	return obtained;
}

//指定点から最も近い、垂線の足とパラメータ値を返す。
//Return the foot of the perpendicular straight line from p that is 
//nearest to point p.
// Function's return value is whether point is obtained(1) or not(0)
int MGFace::perp_point (
	const MGPosition& p,		// 指定点(point)
	MGPosition& uv,		//Parameter value of the surface will be returned.
	const MGPosition* uvguess	// guess parameter value of surface
	) const
{
	MGPosition uv0(1.,1.), uv1(0.,0.);
	if(uvguess) return perp_guess(p,*uvguess,uv);
	else        return perp_one(p,uv);
}

//Compute perpendicular points on the face from a point P((x,y,z)).
//MGPosition uv in the MGPosition_list is:
//uv(0): u parameter, and uv(1): v parameter of the face.
//Generally number of uv are more than one.
MGPosition_list MGFace::perps(const MGPosition& P) const
{
	MGPosition_list uvs;
	const MGSurface* srf=surface(); if(!srf) return uvs;
	uvs=srf->perps(P);
	remove_outside_param(uvs);
	return uvs;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGFace::pick_closest(const MGStraight& sl)const{
	const MGSurface* srf=surface();
	if(no_outer_boundaries()) return srf->pick_closest(sl);

	MGCSisect_list ises=srf->isect(sl);
	size_t n=ises.size();
	double t; MGPosition uv;
	MGCSisect_list::CSiterator is=ises.begin(), ie=ises.end();
	for(; is!=ie; is++){
		const MGPosition& uv2=(*is).param_surface();
		if(in_range(uv2)){
			double t2=(*is).param_curve();
			if(uv.is_null()){
				t=t2; uv=uv2;
			}else{
				if(t2>t){
					t=t2; uv=uv2;
				}
			}
		}
	}
	if(!uv.is_null()) return uv;

	//Second, compute closest candidates on the boundaries.
	//	2.1 Outer boundary.
	MGPvector<MGCurve> prmtr=outer_boundary();
	MGPosition P=prmtr[0]->start_point(), P1;
	double tsl=sl.closest(P);
	MGVector dif=(P-sl.eval(tsl));
	double dist=dif%dif;

	size_t pnum,i;
	pnum=prmtr.size();
	for(i=0; i<pnum; i++){
		 MGPosition tt=prmtr[i]->closest(sl);
		 P1=prmtr[i]->eval(tt[0]);
		 dif=P1-sl.eval(tt[1]);
		 double dist1=dif%dif;
		 if(dist1<dist){dist=dist1; P=P1;}
	}

	//	2.2 Inner boundaries.
	size_t m=number_of_inner_boundaries();
	for(size_t j=0; j<m; j++){
		MGPvector<MGCurve> prmtrj=inner_boundary(j);
		pnum=prmtrj.size();
		for(i=0; i<pnum; i++){
			 MGPosition tt=prmtrj[i]->closest(sl);
			 P1=prmtrj[i]->eval(tt[0]);
			 dif=P1-sl.eval(tt[1]);
			 double dist1=dif%dif;
			if(dist1<dist){dist=dist1; P=P1;}
		}
	}
	return param(P);;
}

//Dedicated function of range.
//Will check if point (u,v) is inside inner boundary or not.
//If inside an inner boundary, obtain the closest point to the boundary.
MGPosition MGFace::range_check_inner_boundary(
	const MGPosition& uv
)const{
	size_t loop_id;
	double len;
	size_t in=inside_inner_boundary(uv,loop_id);
	if(in>=2){
		const MGLoop* cloop=loop(loop_id);
		MGLEPoint lp=cloop->closest(uv,len);
		return cloop->eval(lp);
	}else{
		return uv;
	}
}

//Round the input parameter (u,v) of the face to the nearest point of
//the face parameter range.
//**********This algorithm should be improved later.**********
//(Especially about using MGClosest_to_curves)
MGPosition MGFace::range(const MGPosition& uv) const{
	const MGSurface* srf=surface();
	if(!srf) return uv;

	double errsave=MGTolerance::wc_zero();
	MGTolerance::set_wc_zero(parameter_error());
	MGPosition uv_new;
	if(!srf->in_range(uv)){
		if(no_outer_boundaries()) uv_new=srf->range(uv);
		else uv_new=MGClosest_to_curves(uv,outer_boundary_param());
		goto end_process;
	}

	//Now uv is in the parameter range of surface.
	if(!number_of_boundaries()){
		uv_new=uv;	//When no boundaries.
		goto end_process;
	}
	if(no_outer_boundaries()){
		uv_new=range_check_inner_boundary(uv);
		goto end_process;
	}
	//When there exists perimeter boundary.
	if(inside_outer_boundary(uv))
		uv_new=range_check_inner_boundary(uv);
	else
		uv_new=MGClosest_to_curves(uv,outer_boundary_param());

end_process:
	MGTolerance::set_wc_zero(errsave);
	return uv_new;
}

//Test if this face has an inactive loop.
//If this has one, return true.
bool MGFace::hasInactiveLoop()const{
	size_t n=number_of_boundaries();
	for(int i=n-1; i>=0; i--){
		if(loop(i)->is_inactive(this))
			return true;
	}
	return false;
}

//Remove inactive loops from this face.
void MGFace::remove_inactive_loops(){
	size_t n=number_of_boundaries();
	for(int i=n-1; i>=0; i--){
		if(loop(i)->is_inactive(this))
			erase_boundary(i);
	}
}

//Remove parameter uv from uvs that is outside face parameter range.
void MGFace::remove_outside_param(MGPosition_list& uvs)const{
	MGPosition_list::iterator i=uvs.begin(), iend=uvs.end(), i1;
	while(i!=iend){
		i1=i; i1++;
		if(!in_range(*i)) uvs.removeAt(i);
		i=i1;
	}
}

//Get surface pointer.
MGSurface* MGFace::surface(){
	MGGeometry* cell=extent();
	return dynamic_cast<MGSurface*>(cell);
}
const MGSurface* MGFace::surface() const{
	const MGGeometry* cell=extent();
	return dynamic_cast<const MGSurface*>(cell);
}

//Obtain the closest point from point uv to vector of curves.
//MGClosest_to_curves does not change wc_zero, and so calling program of
//MGClosest_to_curves should change it if necessary.
MGPosition MGClosest_to_curves(
	const MGPosition& uv,				//Point.
	const MGPvector<MGCurve>& curves	//vector of curves.
){
	size_t n=curves.size(); if(n==0) return uv;
	double t=curves[0]->closest(uv);
	MGPosition uv_min=curves[0]->eval(t), uv1;
	MGVector dif=uv-uv_min;
	double lenmin=dif%dif, len1;
	for(size_t i=1; i<n; i++){
		t=curves[i]->closest(uv);
		uv1=curves[i]->eval(t);
		dif=uv-uv1; len1=dif%dif;
		if(len1<lenmin){lenmin=len1; uv_min=uv1;}
	}
	return uv_min;
}
