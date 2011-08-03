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
#include "mg/Plane.h"
#include "mg/FSurface.h"
#include "mg/SBRep.h"
#include "mg/isects.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/HHisect_vector.h"
#include "topo/FOuterCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//
//Implements MGFace Class.
//MGFace is an instance of MGCell.

///////Member Function///////

//Compute intersection points of an inner parameter line of this face and face2.
//The intersection point is used to compute surface to surface intersection lines.
//Function's return value is at most one intersection point in uvuv_list.
//One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
//this surface and (u2,v2) is a parameter of surf.
MGPosition_list MGFace::intersectInner(
	const MGFace& face2		//The second surface.
) const{
	MGPosition_list uvuv_list;
	const MGBox& pbox1=box_param(); const MGBox& pbox2=face2.box_param();

	double u10=pbox1[0].low_point(), u11=pbox1[0].high_point();
	double v10=pbox1[1].low_point(), v11=pbox1[1].high_point();
	const MGKnotVector& t1u=knot_vector_u();
	const MGKnotVector& t1v=knot_vector_v();
	size_t nspan1u=t1u.locate(u11,1)+1-t1u.locate(u10);
	size_t nspan1v=t1v.locate(v11,1)+1-t1v.locate(v10);

	double u20=pbox2[0].low_point(), u21=pbox2[0].high_point();
	double v20=pbox2[1].low_point(), v21=pbox2[1].high_point();
	const MGKnotVector& t2u=face2.knot_vector_u();
	const MGKnotVector& t2v=face2.knot_vector_v();
	size_t nspan2u=t2u.locate(u21,1)+1-t2u.locate(u20);
	size_t nspan2v=t2v.locate(v21,1)+1-t2v.locate(v20);

	int maximum;
	if(nspan1u<nspan1v){
		if(nspan1v<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1v<nspan2v) maximum=3; else maximum=1;
		}
	}else{
		if(nspan1u<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1u<nspan2v) maximum=3; else maximum=0;
		}
	}
	const MGFace *f1, *f2;
	int isU=1; if(maximum%2) isU=0;
	double t0,t1;
	int nspan;
	if(maximum<=1){
		f1=this; f2=&face2;
		if(isU){ t0=u10; t1=u11; nspan=nspan1u;}
		else{ t0=v10; t1=v11; nspan=nspan1v;}
	}else{
		f2=this; f1=&face2;
		if(isU){ t0=u20; t1=u21; nspan=nspan2u;}
		else{ t0=v20; t1=v21; nspan=nspan2v;}
	}
	double tau=(t0+t1)*.5;
	double delta=(t1-t0)/double(nspan);
	int sign=-1;
	MGCSisect_list csiList;
	for(int i=0; i<nspan; i++){
		tau+=double(i*sign)*delta;
		MGPvector<MGCurve> crvs;
		if(tau>t0 && tau<t1)
			crvs=f1->parameter_curves(isU,tau);
		for(size_t j=0; j<crvs.size(); j++){
			csiList=crvs[j]->isect(*f2);
			if(csiList.size()) break;
		}
		if(csiList.size()) break;
		sign*=-1;
	}
	if(!csiList.size()) return uvuv_list;

	double u1,v1;
	u1=tau; v1=csiList.front().param_curve();
	if(!isU){ u1=v1; v1=tau;}
	const MGPosition& uv2=csiList.front().param_surface();
	MGPosition uvuv(4);
	if(f1==this){
		uvuv(0)=u1; uvuv(1)=v1; uvuv.store_at(2,uv2,0,2);
	}else{
		uvuv(2)=u1; uvuv(3)=v1; uvuv.store_at(0,uv2,0,2);
	}
	uvuv_list.append(uvuv);
	return uvuv_list;
}

//Compute intersection points of an inner parameter line of this face and sf2.
//The intersection point is used to compute surface to surface intersection lines.
//Function's return value is at most one intersection point in uvuv_list.
//One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
//this surface and (u2,v2) is a parameter of surf.
MGPosition_list MGFace::intersectInner(
	const MGSurface& sf2		//The second surface.
) const{
	MGPosition_list uvuv_list;
	const MGBox& pbox1=box_param();

	double u10=pbox1[0].low_point(), u11=pbox1[0].high_point();
	double v10=pbox1[1].low_point(), v11=pbox1[1].high_point();
	const MGKnotVector& t1u=knot_vector_u();
	const MGKnotVector& t1v=knot_vector_v();
	size_t nspan1u=t1u.locate(u11,1)+1-t1u.locate(u10);
	size_t nspan1v=t1v.locate(v11,1)+1-t1v.locate(v10);

	double u20=sf2.param_s_u(), u21=sf2.param_e_u();
	double v20=sf2.param_s_v(), v21=sf2.param_e_v();
	const MGKnotVector& t2u=sf2.knot_vector_u();
	const MGKnotVector& t2v=sf2.knot_vector_v();
	size_t nspan2u=t2u.locate(u21,1)+1-t2u.locate(u20);
	size_t nspan2v=t2v.locate(v21,1)+1-t2v.locate(v20);

	int maximum;
	if(nspan1u<nspan1v){
		if(nspan1v<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1v<nspan2v) maximum=3; else maximum=1;
		}
	}else{
		if(nspan1u<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1u<nspan2v) maximum=3; else maximum=0;
		}
	}
	int isU=1; if(maximum%2) isU=0;
	double t0,t1,tau;
	int nspan;
	int sign=-1;
	MGCSisect_list csiList;
	if(maximum<=1){
		if(isU){ t0=u10; t1=u11; nspan=nspan1u;}
		else{ t0=v10; t1=v11; nspan=nspan1v;}
		tau=(t0+t1)*.5;
		double delta=(t1-t0)/double(nspan);
		for(int i=0; i<nspan; i++){
			tau+=double(i*sign)*delta;
			MGPvector<MGCurve> crvs;
			if(tau>t0 && tau<t1) crvs=parameter_curves(isU,tau);
			for(size_t j=0; j<crvs.size(); j++){
				csiList=sf2.isect(*(crvs[j]));
				if(csiList.size()) break;
			}
			if(csiList.size()) break;
			sign*=-1;
		}
	}else{
		if(isU){ t0=u20; t1=u21; nspan=nspan2u;}
		else{ t0=v20; t1=v21; nspan=nspan2v;}
		tau=(t0+t1)*.5;
		double delta=(t1-t0)/double(nspan);
		for(int i=0; i<nspan; i++){
			tau+=double(i*sign)*delta;
			if(tau>t0 && tau<t1){
				MGCurve* crv=sf2.parameter_curve(isU,tau);
				csiList=isect(*crv); delete crv;
				if(csiList.size()) break;
			}
			if(csiList.size()) break;
			sign*=-1;
		}
	}
	if(!csiList.size()) return uvuv_list;

	double u1,v1;
	u1=tau; v1=csiList.front().param_curve();
	if(!isU){ u1=v1; v1=tau;}
	const MGPosition& uv2=csiList.front().param_surface();
	MGPosition uvuv(4);
	if(maximum<=1){
		uvuv(0)=u1; uvuv(1)=v1; uvuv.store_at(2,uv2,0,2);
	}else{
		uvuv(2)=u1; uvuv(3)=v1; uvuv.store_at(0,uv2,0,2);
	}
	uvuv_list.append(uvuv);
	return uvuv_list;
}

//Curve to face intersection.
MGCSisect_list MGFace::isect(const MGCurve& curv)const{
	return curv.isect(*this);
}

MGSSisect_list MGFace::isect(const MGFSurface& surf2)const{
	return intersect(surf2);
}
MGSSisect_list MGFace::isect(const MGFace& surf2)const{
	return intersect(surf2);
}
MGSSisect_list MGFace::isect(const MGSurface& surf2)const{
	return intersect(surf2);
}

//Face to Face intersection.
MGSSisect_list MGFace::intersect(const MGFSurface& face2) const{
	if(!face2.has_commonFS(*this))
		return MGSSisect_list(this,&face2);

	MGPosition_list uvuv_list;

	//Compute intersection points of this face's boundary and face2.
	//Will be stored in uvuv_list. The one member is uvuv(7):
	//0,1:(u,v) of this, 2,3:(u,v) of face2, 4-6:direction of isect.
	intersect12Boundary(face2,uvuv_list);
	//cout<<"MGFace::isect(MGSurface), uvuv_list="<<uvuv_list<<endl;
	return isectEP(uvuv_list,face2);
}

//Obtain parameter u(kcod=1), or v(kcod=0) of the intersection point of
//v=x(const) line(kcod=1), or u=x(const) line (kcod=0) with
//all the face boundaries.
std::vector<double> MGFace::isect1D_with_boundaries(
	double x,							//coordinate value of kcod.
	size_t kcod)						//Coordinate kind, =0:u, =1: v.
const{
	MGPosition_list plist;
	double error=parameter_error();
	double error_save=MGTolerance::wc_zero();
	MGTolerance::set_wc_zero(error);

	//1. outer boundaries.
	MGPvector<MGCurve> curves=outer_boundary_param();
	//cout<<curves<<endl;
	size_t n=curves.size(),i;
	for(i=0; i<n; i++){
		MGCParam_list list=curves[i]->isect_1D(x,kcod);
		MGCParam_list::Citerator ci=list.begin(), ciend=list.end();
		while(ci!=ciend){
			MGPosition uv=curves[i]->eval(*ci++);
			plist.append(*this, uv);
		}
	}

	//2. inner boundaries.
	size_t n_inner=number_of_inner_boundaries();
	for(size_t j=0; j<n_inner; j++){
		MGPvector<MGCurve> curvesj=inner_boundary_param(j);
		n=curvesj.size();
		for(i=0; i<n; i++){
			MGCParam_list list=curvesj[i]->isect_1D(x,kcod);
			MGCParam_list::Citerator ci=list.begin(), ciend=list.end();
			while(ci!=ciend){
				MGPosition uv=curvesj[i]->eval(*ci++);
				plist.append(*this, uv);
			}
		}
	}

	size_t other_cod=(kcod+1)%2;
	size_t np=plist.size();
	std::vector<double> ivec(np);
	MGPosition_list::iterator pi=plist.begin();
	for(i=0; i<np; i++)	ivec[i]=(*pi++).ref(other_cod);
	std::sort(ivec.begin(),ivec.end());

	MGTolerance::set_wc_zero(error_save);
	return ivec;
}

//Compute intersection points of this face's boundary(outer and inners) with
//face2. If intersection points are found and the boundary is a loop,
//the point's edge pointer(of this) will be stored in a member uvuv of uvuvs.
//uvuv[7] is the edge pointer. If the boundary is not a loop(that is, a perimeter of
//Surfaces), uvuv.sdim()==7 and an edge pointer is not returned.
//When uvuv.sdim()==8, the edge pointer of uvuv[7] is accessed through union mgEdgeP.
//uvuvs[i] is i-th intersection points.
size_t MGFace::isect_boundary(
	const MGFSurface& face2,
	MGPosition_list& uvuvs,
	//id1 and id2 are the ids of uvuv where this face's and f2's parameters
	//are to be stored in a member of uvuvs.
	//This face's (u,v) is stored in uvuv(id1) and (id1+1).
	//f2's (u,v) is stored in uvuv(id2) and (id2+1).
	//id2=0 if id1=2, and id2=2 if id1=0.
	size_t id1
)const{
	size_t inum=0;

	//Compute intersection points of the outer boundary of this face and face2.
	//Will be stored in uvuvs. The one member is uvuv(8):
	//0,1:(u,v) of face, 2,3:(u,v) of surf, 4-6:direction of isect.
	//7:for edge pointer at the intersection point.
	inum+=isect_outcurves(face2,uvuvs,id1);

	//Compute intersection points of inner boundaries of face fi and surf.
	size_t ibStrt;
	size_t nin=number_of_inner_boundaries(ibStrt);
	for(size_t k=0; k<nin; k++)
		inum+=isect_lpcurves(face2,ibStrt+k,uvuvs,id1);
	return inum;
}

//Compute intersection points between the boundary of iid-th inner boundary
//of this face and face2 to compute intersections of face with face2.
//Function's return value is the number of ip's obtained before appending
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t MGFace::isect_incurves(
	const MGFSurface& face2,
	size_t iid,				//Loop id of this face.
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 8,
		//(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	size_t start_id;
	size_t num=number_of_inner_boundaries(start_id);
	return isect_lpcurves(face2,start_id+iid,uvuv_list,id1);
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGFace::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,	//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	size_t incr			//Incremental valuse of B-coef's id.
)const{
	const MGSurface* surf=surface();
	return surf->isect_incr_pline(uv,kdt,du,dv,u,v,incr);
}

//Compute intersection points between loop lp of this face and face2
//to compute intersections of face with face2.
//Function's return value is the number of ip's obtained before appending
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t MGFace::isect_lpcurves(
	const MGFSurface& face2,		//srf!=null, face2=null.
	const MGLoop& lp,				//Loop id of this face.
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 8,
		//(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	size_t id2=2; if(id1==2) id2=0;
	double zero_angle=MGTolerance::angle_zero();

	mgEdgeP edgeP;//Union to access MGEdge* in MGPostion.
	MGLSPoint_vector lsps=lp.isect_binder(face2);
	size_t m=lsps.entries();
	for(size_t j=0; j<m; j++){
		const MGLSPoint& lsp=lsps.m_lspoints[j];
		const MGPosition& suv=lsp.surface_param();
		MGPosition uvuv(8,suv,id2,0);
		double tb=lsp.binder_param();	//tb is binder edge's curve parameter value.
		const MGEdge* pedge=lsp.parameter_edge();
		double tp=pedge->param_pcell(tb);//tp is parameter edge's parameter.
		MGPosition uv=pedge->eval(tp);
		uvuv.store_at(id1,uv,0,2);
		MGVector N2=face2.normal(suv);
		MGVector N=normal(uv), V=pedge->eval_star(tp,1);
		double vn2angle=V.cangle(N2); if(vn2angle<.0) vn2angle*=-1.;
		if(vn2angle<=zero_angle){
			//vn2angle<=zero_angle means V and srf(or face2) are almost parallel.
			//double tps=pedge->param_s(), tpe=pedge->param_e();
			//if(tp>(tps+tpe)*.5) V*=-1.;
			//uvuv.store_at(4,V,0,3);
			uvuv.store_at(4,MGDefault::zero_vector(),0,3);
		}else uvuv.store_at(4,N2*(N*V)*N2,0,3);
			//cout<<"isect_outcurves:N="<<N<<",V="<<V<<" "<<uvuv<<endl;//////
			//cout<<eval(uvuv[id1],uvuv[id1+1])<<endl;
		edgeP.pointer=pedge; uvuv(7)=edgeP.doubleV;
		//cout<<"isect_lpcurves:N="<<N<<"V="<<V<<" "<<uvuv<<endl;//////
		uvuv_list.append(*this,face2,uvuv);
	}
	return m;
}

//Compute intersection points between loop lpid of this face and face2
//to compute intersections of face with face2.
//Function's return value is the number of ip's obtained before appending
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t MGFace::isect_lpcurves(
	const MGFSurface& face2,		//srf!=null, face2=null.
	size_t lpid,				//Loop id of this face.
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 8,
		//(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	size_t id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	const MGLoop& lp=*(loop(lpid));
	return isect_lpcurves(face2,lp,uvuv_list,id1);
}

//Compute intersection points of outer boundary curves of this face 
//with face2 to compute intersections.
//Function's return value is the number of ip's obtained(appended)
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t MGFace::isect_outcurves(
	const MGFSurface& face2,
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 7,
		//and the last three elements are the ip direction vector.
	size_t id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	size_t id2=2; if(id1==2) id2=0;
	const MGSurface* srf1=surface();
	MGPvector<MGCurve> perim(srf1->perimeter_num());
		//To save whole perimeter curve construction
		//at multiple processes of same perimeters.

	double zero_angle=MGTolerance::angle_zero();
	std::vector<MGFOuterCurve> cid=outer_curve();
	
	//Compute intersection points.
	size_t inum=0;
	std::vector<MGFOuterCurve>::iterator i=cid.begin(), ie=cid.end();
	for(; i!=ie; i++){
		if(i->is_loop()){
			//When loop.
			const MGLoop& lpi=*(i->loop());
			inum+=isect_lpcurves(face2,lpi,uvuv_list,id1);
		}else{
			//When perimeter.
			size_t pid=i->perimeter_id();
			if(perim[pid]==0) perim.reset(pid,srf1->perimeter_curve(pid));
			double t0,t1; i->range(t0,t1);
			//if(t0>t1){ double save=t0; t0=t1; t1=save;}
			MGTrimmedCurve crvi(*(perim[pid]), t0,t1);
			double error=MGTolerance::wc_zero();
			MGTolerance::set_wc_zero(error*.1);
			MGCSisect_list csilist=face2.isect(crvi);
			MGTolerance::set_wc_zero(error);
			size_t m=csilist.entries(); inum+=m;
			for(size_t j=0; j<m; j++){
				MGCSisect cs=csilist.removeFirst();
				MGPosition uvuv(4,cs.param_surface(),id2,0);
				double t=cs.param_curve();	//Perimeter parameter.
				MGPosition uv=srf1->perimeter_uv(pid,t);
				uvuv.store_at(id1,uv,0,2);
				uvuv_list.prepend(*this,face2,uvuv);
			}
		}
	}
	return inum;
}

//Compute intersection lines, given end points of the i.l.
MGSSisect_list MGFace::isectEP(
	MGPosition_list& uvuv_list,	//End points list of the intersection.
		//On return, uvuv_list.size() will be 0.
	const MGFSurface& fsrf2		//The second surface.
) const{
	MGSSisect_list lst(this,&fsrf2);
	if(uvuv_list.size()==0) return lst;

	MGSSisect ssi;
	MGPosition_list::iterator uvuv_id;
	int obtained;
	const MGPlane* pl1=dynamic_cast<const MGPlane*>(surface());
	const MGPlane* pl2=dynamic_cast<const MGPlane*>(fsrf2.get_surface_pointer());
	if(pl1 && pl2){//When both are planes.
		while(uvuv_list.entries()){
			MGPosition uvuvS=uvuv_list.removeFirst();
			if(obtained=pl1->isect_startHPL(uvuvS,uvuv_list,*pl2,ssi,uvuv_id)){
				lst.append(ssi);
				if(obtained==3) uvuv_list.removeAt(uvuv_id);
			}
		}
	}else if(pl2){//When this is not a plane and srf2 is a plane.
		//Compute intersection lines using isect_with_plane.
		lst=isect_with_plane(uvuv_list,*pl2,fsrf2);
	}else if(pl1){//When this is a plane and srf2 is not a plane.
		MGPosition_list uvuvlist2;
		MGPosition_list::iterator i=uvuv_list.begin(), ie=uvuv_list.end();
		for(; i!=ie; i++){
			MGPosition uvuv2(*i);
			uvuv2(0)=(*i)[2]; uvuv2(1)=(*i)[3];
			uvuv2(2)=(*i)[0]; uvuv2(3)=(*i)[1];
			uvuvlist2.append(uvuv2);
		}
		lst=fsrf2.isect_with_plane(uvuvlist2,*pl1,*this);
		lst.replace12();
		uvuv_list.clear();
	}else{//When both are not planes.
	//Compute intersection line using isect_with_surf.
		lst=isect_with_surf(uvuv_list,fsrf2);
	}
	return lst;
}

//Compute the intersections of two objects.
MGisects MGCurve::intersection(const MGFace& obj2)const{
	MGCSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}

MGisects MGSurface::intersection(const MGFace& obj2)const{
	MGSSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}

//Compute the intersections of two objects.
MGisects MGFace::intersection(const MGObject& obj2)const{
	MGisects isects=obj2.intersection(*this);
	isects.exchange12();
	return isects;
}
MGisects MGFace::intersection(const MGCurve& obj2)const{
	MGCSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGFace::intersection(const MGFSurface& obj2)const{
	MGSSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGFace::intersection(const MGSurface& obj2)const{
	MGSSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGFace::intersection(const MGFace& obj2)const{
	MGSSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGFace::intersection(const MGShell& obj2)const{
	MGHHisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}
