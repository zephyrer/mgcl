/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <iterator>
#include "mg/Box.h"
#include "mg/Position_list.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect.h"
#include "mg/isects.h"
#include "mg/Tolerance.h"
#include "topo/LEPoint.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "topo/HHisect_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGShell Class.

///////Member Function///////

//Compute the closest point from a point to this shell.
MGFPoint MGShell::closest(const MGPosition& point) const{
	size_t n=number_of_pcells(); if(!n) return MGFPoint();
	double* dist=new double[n], min, disti;
	size_t j;

	const_pcellItr i=pcell_begin(), ie=pcell_end();
	const MGFace* f=0;
	for(j=0; j<n; i++, j++){
		const MGFace* fi=dynamic_cast<const MGFace*>(*i);
		if(!fi)
			continue;
		MGBox boxi=fi->box();dist[j]=disti=boxi.distance(point);
		if(!f){//For the first fi.
			f=fi; min=disti;
		}else{
			if(disti<min){
				f=fi; min=disti;
			}
		}
	}

	MGPosition uv=f->closest(point);
	min=(f->eval(uv)-point).len();
	const MGFace* fmin=f;
	if(!MGAZero(min)){ 

	for(i=pcell_begin(), j=0; j<n; i++, j++){
		const MGFace* fi=dynamic_cast<const MGFace*>(*i);
		if(!fi || fi==f)
			continue;
		if(min<dist[j])
			continue;
		MGPosition uvi=fi->closest(point);
		disti=(fi->eval(uvi)-point).len();
		if(disti<min){
			min=disti; fmin=fi; uv=uvi;
			if(MGAZero(min))
				break;
		}
	}

	}

	delete[] dist;
	return MGFPoint(*fmin, uv);
}

//Intersection of a shell and a curve.
MGCFisect_vector MGCurve::isect(const MGShell& shl) const{
	return shl.isect(*this);
}

//Intersection of a shell and a curve.
MGCFisect_vector MGShell::isect(const MGCurve& curve) const{
	MGCFisect_vector is(&curve);
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=dynamic_cast<const MGFace*>(*i);
		if(!fi)
			continue;
		MGCSisect_list list=fi->isect(curve);
		while(!list.empty())
			is.push_back(MGCFisect(list.removeFirst(),*fi));
	}
	return is;
}

//Intersection of a shell and a surface.
MGHHisect_vector MGShell::isect(const MGSurface& surf) const{
	MGHHisect_vector vec;
	size_t nf=number_of_pcells();//nf=num of faces included.
	if(!nf) return vec;
	if(!has_common(surf))
		return vec;

	MGHHisect_vector lines2;
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i);
		if(!fi)
			continue;
		MGSSisect_list ilst=fi->isect(surf);
		int ncrvs=ilst.size();
		for(int j=0; j<ncrvs; j++){
			lines2.push_back(MGHHisect(fi,0,ilst.front()));
			ilst.pop_front();
		}
	}

	size_t nlines=lines2.size();
	int nlines2=0;
	for(size_t k=0; k<nlines; k++){
		MGHHisect is(lines2[k]);
		if(is.is_null()) continue;
		is.build_one(lines2);
		vec.push_back(is);nlines2++;
	}
	return vec;
}

MGHHisect_vector MGFace::isect(const MGShell& shell) const{
	MGHHisect_vector hhi=shell.isect(*this);
	hhi.replace12();
	return hhi;
}
MGHHisect_vector MGSurface::isect(const MGShell& shell) const{
	MGHHisect_vector hhi=shell.isect(*this);
	hhi.replace12();
	return hhi;
}

//Intersection of a shell and a face.
MGHHisect_vector MGShell::isect(const MGFace& face2) const{
	MGHHisect_vector vec;
	size_t nf=number_of_faces();//nf=num of faces included.
	if(!nf) return vec;
	if(!has_common(face2)) return vec;

	MGHHisect_vector lines2;
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		MGSSisect_list ilst=fi->isect(face2);
		int ncrvs=ilst.size();
		for(int j=0; j<ncrvs; j++){
			lines2.push_back(MGHHisect(fi,&face2,ilst.front()));
			ilst.pop_front();
		}
	}

	size_t nlines=lines2.size();
	int nlines2=0;
	for(size_t k=0; k<nlines; k++){
		MGHHisect is(lines2[k]);
		if(is.is_null()) continue;
		is.build_one(lines2);
		vec.push_back(is);nlines2++;
	}
	return vec;

/***********************************OLD ALGORITHM***********
	const MGSurface& surf2=*(face2.surface());
	MGPosition_list* uvuv_list=new  MGPosition_list[nf];
	//In uvuv_list[j] is the list of intersection points of j-th face's
	//with the surface.
	//Let uvuv=uvuv_list[j], then uvuv[0-1]:(u,v) of this j-th face,
	//uvuv[2-3]:(u,v) of surf, uvuv[4-6]: directio vector of the ip., 
	//uvuv[7]:edge pointer of j-th face(is accessed through mgEdgeP).

	size_t inS;
	size_t nin2=face2.number_of_inner_boundaries(inS);
	size_t j;
	for(j=0; j<nf; j++){
		const MGFace* fj=face(j); if(!fj) continue;
		if(!face2.has_common(*fj)) continue;

		//Compute intersection points of fj's boundary(outer and inners) 
		//and face2. If intersection points are found,
		//the point's edge pointer(of fj) will be stored in uvuv_list[j].
		//i-th member of uvuv_list[j] is i-th intersection points
		//and the edge pointer of fj's boundary.
		fj->isect_boundary(&face2, uvuv_list[j]);
			//A member of uvuv_list is of sdim 8.

		//Compute intersection points of face2 boundary and this fj face.
		face2.isect_outcurves(0,fj,uvuv_list[j],2);
			//A member of uvuv_list is of sdim 7.
		for(size_t k2=0; k2<nin2; k2++){
			face2.isect_lpcurves(0,fj,inS+k2,uvuv_list[j],2);
		}

		//Compute a loop intersection line of fj.
		if(uvuv_list[j].entries()) continue;
		MGPosition_list listuvuv=fj->intersectInner(face2);
		if(!listuvuv.size()) continue;
		MGSSisect_list ssiList=fj->isectEP(listuvuv,surf2,&face2);
		if(!ssiList.size()) continue;
		MGSSisect iis=ssiList.removeFirst();
		MGHHisect hhi(fj,&face2,iis);
		vec.push_back(hhi);
	}

	//Compute intersection line using isect_startH.
	int obtained;
	MGPosition_list::iterator uvuv_id;
	MGSSisect ssi;

	for(j=0; j<nf; j++){//Loop for all faces of this shell.

	while(uvuv_list[j].entries()){//Loop for all ips of boundaries.
		uvuv_id=uvuv_list[j].begin();
		MGPosition uvuvSS=*uvuv_id;//Save the starting point of j-th face.

		MGPosition uvuvS, uvuvE;
		MGHHisect hhi;
		//uvuv_id, obtained, and jtemp are loop variables of while below.
	//1. Get the ip to the end direction.
		size_t jtemp=j;
		do{

		const MGFace* fj=face(jtemp);
		const MGSurface* fjsrf=fj->surface();
		uvuvS=*uvuv_id; uvuv_list[jtemp].erase(uvuv_id);
		obtained=fjsrf->isect_startH(uvuvS,uvuv_list[jtemp],surf2,ssi,uvuv_id);
		if(obtained){
			hhi.connect_line_to_end(ssi.release_line(),
				fj,ssi.release_param1(), &face2,ssi.release_param2());
			if(obtained==3){
					//Trace to the end direction from its partner.
				uvuvE=*uvuv_id; uvuv_list[jtemp].erase(uvuv_id);
				if(!isect_partner(uvuvE,uvuv_list,jtemp,uvuv_id)) break;
					//To halt the while loop if the end point did not have a partner.
			}
		}

		}while(obtained==3);

		uvuvE=uvuvSS;
	//2. Get the ip to the start direction.
		do{
			if(!isect_partner(uvuvE,uvuv_list,jtemp,uvuv_id)) break;
				//To halt the while loop if the end point did not have a partner.
			const MGFace* fj=face(jtemp);
			const MGSurface* fjsrf=fj->surface();
			uvuvS=*uvuv_id; uvuv_list[jtemp].erase(uvuv_id);
			obtained=fjsrf->
				isect_startH(uvuvS,uvuv_list[jtemp],surf2,ssi,uvuv_id);
			if(obtained){
				hhi.connect_line_to_start(ssi.release_line(),
					fj,ssi.release_param1(), &face2,ssi.release_param2());
				if(obtained==3){//Trace to the start direction from its partner.
					uvuvE=*uvuv_id; uvuv_list[jtemp].erase(uvuv_id);
				}
			}
		}while(obtained==3);

		vec.push_back(hhi);
	}

	}

	delete[] uvuv_list;
	return vec;
********************************************/
}

//Intersection of two shells.
MGHHisect_vector MGShell::isect(const MGShell& shell2)const{
	MGHHisect_vector vec;
	size_t nf1=number_of_faces(), nf2=shell2.number_of_faces();
		//nf1,2=num of faces included in shell1,2.
	if(!nf1 || !nf2) return vec;
	if(!has_common(shell2)) return vec;

	MGHHisect_vector lines2;
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		const_pcellItr i2=shell2.pcell_begin(), ie2=shell2.pcell_end();
		for(; i2!=ie2; i2++){
			const MGFace* fi2=shell2.face(i2); if(!fi2) continue;
			MGSSisect_list ilst=fi->isect(*fi2);
			int ncrvs=ilst.size();
			for(int j=0; j<ncrvs; j++){
				lines2.push_back(MGHHisect(fi,fi2,ilst.front()));
				ilst.pop_front();
			}
		}
	}

	size_t nlines=lines2.size();
	int nlines2=0;
	for(size_t k=0; k<nlines; k++){
		MGHHisect is(lines2[k]);
		if(is.is_null()) continue;
		is.build_one(lines2);
		vec.push_back(is);nlines2++;
	}
	return vec;
}

//Intersection of a shell and a face.
//This shell's face is face1 in HHisect and face2 is face.
MGHHisect_vector MGShell::isect(const MGFSurface& face) const{
	MGHHisect_vector isects=face.isect(*this);
	isects.replace12();
	return isects;
}

//Get the partner point uvuv_id of the uvuv_list[j].
//If found, return true, if not, return false.
//Shell*(face or surface) version.
//The other one must not be shell, a surface or a face.
bool MGShell::isect_partner(
	const MGPosition& uvuv,//The point to get the partner.
	MGPosition_list* uvuv_list,
	size_t& j,		//Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr
		//Partner point's iterator of uvuv_list[j] will be returned.
)const{
	if(uvuv.sdim()<8)
		return false;

	mgEdgeP edgP;	//Work area to access to MGEdge* in uvuv_list.
	const MGEdge* e; edgP.doubleV=uvuv[7]; e=edgP.pointer;
	const MGEdge* pe=e->first_partner();
	if(!pe) return false;//If e had no partnres.

	const MGCellNB* pf=pe->star();
	const_pcellItr fitr=std::find(pcell_begin(), pcell_end(), pf);
	if(fitr==pcell_end())
		return false;

	j=std::distance(pcell_begin(), fitr);
		//j is the id of uvuv_list of face pf.		

	MGPosition P=e->face()->eval(uvuv[0],uvuv[1]);

	//Find uvuv_list[j]'s uvuv that has the same edge and nearest to P.
	MGPosition_list::iterator itr=uvuv_list[j].begin(), itre=uvuv_list[j].end();
	double dlen=-1.;
	for(; itr!=itre; itr++){
		if(itr->sdim()<8) break;
			//We can break since elements of sdim()<8 are set last in the array.
		const MGEdge* et; edgP.doubleV=itr->ref(7); et=edgP.pointer;
		if(et==pe){
			MGPosition Q=et->face()->eval((*itr)[0],(*itr)[1]);
			if(dlen<0.){
				dlen=(P-Q).len(); uvuvItr=itr;
			}else{
				double dlen2=(P-Q).len();
				if(dlen>dlen2){
					dlen=dlen2; uvuvItr=itr;
				}
			}
		}
	}
	if(dlen<0.) return false;
	return true;
}

//Get the partner point uvuv_id of the uvuv_list[j].
//If found, return true, if not, return false.
//Shell*shell intersection version.
bool MGShell::isect_partner(
	const MGPosition& uvuv,//The point to get the partner.
	MGPosition_list* uvuv_list,
	size_t& i,
	size_t& j,		//Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr
		//Partner point's iterator of uvuv_list[i+nf1*j] will be returned.
)const{
	if(uvuv.sdim()<8)
		return false;

	mgEdgeP edgP;	//Work area to access to MGEdge* in uvuv_list.
	const MGEdge* e; edgP.doubleV=uvuv[7]; e=edgP.pointer;
	const MGEdge* pe=e->first_partner();
	if(!pe) return false;//If e had no partnres.
	const MGCellNB* pf=pe->star();
	const MGComplex* shl=pf->parent_complex();
		//Shell of pf and e->star() is the same.
	const_pcellItr fitr=std::find(shl->pcell_begin(), shl->pcell_end(), pf);
	if(fitr==shl->pcell_end()) return false;

	size_t k, kp1;
		//Id of uvuv to evaluate the position(this shell or the other shell).
	if(shl==this){
		i=std::distance(pcell_begin(), fitr);
		k=0; 
	}else{
		j=std::distance(shl->pcell_begin(), fitr);
		k=2; 
	}
	kp1=k+1;
	MGPosition P=e->face()->eval(uvuv[k],uvuv[kp1]);
	
	//i+j*number_faces() is the id of uvuv_list of face pf.
	size_t id=i+j*number_of_faces();
	//Find uvuv_list[j]'s uvuv that has the same edge and nearest to P.
	MGPosition_list::iterator itr=uvuv_list[id].begin(), itre=uvuv_list[id].end();
	double dlen=-1.;
	for(; itr!=itre; itr++){
		const MGEdge* et; edgP.doubleV=itr->ref(7); et=edgP.pointer;
		if(et==pe){
			MGPosition Q=et->face()->eval((*itr)[k],(*itr)[kp1]);
			if(dlen<0.){
				dlen=(P-Q).len(); uvuvItr=itr;
			}else{
				double dlen2=(P-Q).len();
				if(dlen>dlen2){
					dlen=dlen2; uvuvItr=itr;
				}
			}
		}
	}
	if(dlen<0.) return false;
	return true;
}

//Test if a point is on the shell or not.
//Function's return value is true if the point is on the shell, and false if not.
//The point parameter of the shell is returned in fp if true is returned.
//If false is returned, the closest point of the shell will be returned in fp.
bool MGShell::on(const MGPosition& point,
	MGFPoint& fp				//Shell's point parameter value.
)const{
	fp=closest(point);
	return MGAZero((fp.eval()-point).len())!=0;
}

//Obtain perpendicular points of a shell from a point.
std::vector<MGFPoint> MGShell::perps(const MGPosition& point) const{
	double error=MGTolerance::wc_zero(); error*=4.;
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	MGPosition_list Ps;
		//This is used for the same points to be not added in the answer. 
	std::vector<MGFPoint> fps;//return value.
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		MGPosition_list uvs=fi->perps(point);
		MGPosition_list::iterator j=uvs.begin(), je=uvs.end(), jwk;
		for(;j!=je; j++){
			MGPosition P=fi->eval(*j);
			MGBox box(P,error);
			if(!Ps.in(box,jwk)){
				Ps.push_back(P);
				fps.push_back(MGFPoint(*fi,*j));
			}
		}
	}
	return fps;
}

//Obtain the projected curves of a curve onto the shell.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the shell if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves(iline() of the MGHHisect members of lines),
//and the other is (u,v) curves of the parameter space of the surfaces
//(uvline1() of the MGHHisect members of lines ).
//*** uvline2() of the MGHHisect members of lines is a deque of length zero
//(not used).
//Function's return value is the number of curves obtained.
//When <0 is returned, some internal error occured.
int MGShell::project(
	const MGCurve& crv,		//given world coordinate curve to project.
	MGHHisect_vector& lines,
			//World coordinates (x,y,z) lines and (u,v) lines of the projected
			//curves will be returned in lines.
	const MGVector& vec	//projection direction.
			//if vec = NULL, then projection that is normal to the shell.
)const{
	MGHHisect_vector lines2;
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		MGPvector<MGCurve> crv_uvs;		//uv projection curve.
		MGPvector<MGCurve> crvs;		//projection curve.
		int ncrvs=fi->project(crv,crv_uvs,crvs,vec);
		if(ncrvs<0)
			continue;//since error detected.
		for(int j=0; j<ncrvs; j++)
			lines2.push_back(MGHHisect(crvs.release(j),fi,crv_uvs.release(j)));
	}
	int nlines=lines2.size();
	if(!nlines) return 0;

	int k, nlines2=0;

	//Search hhisect closest to the start point of crv.
	MGPosition PS=crv.start_point();
	double dist=-1.;
	int nmindif=0;
	for(k=0; k<nlines; k++){
		MGVector diff=lines2[k].iline().start_point()-PS;
		double dist2=diff%diff;
		if(dist<0. || dist2<dist){
			nmindif=k;
			dist=dist2;
		}
	}

	MGHHisect is(lines2[nmindif]);
	for(k=0; k<nlines; k++){
		if(k)
			is=lines2[k];
		if(is.is_null()) continue;
		is.build_one(lines2);
		lines.push_back(is);nlines2++;
	}
	return nlines2;
}

MGisects MGCurve::intersection(const MGShell& obj2)const{
	MGCFisect_vector isects2=obj2.isect(*this);
	MGisects isects(isects2);
	return isects;
}

MGisects MGSurface::intersection(const MGShell& obj2)const{
	MGHHisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}

//Compute the intersections of two objects.
MGisects MGShell::intersection(const MGObject& obj2)const{
	MGisects isects=obj2.intersection(*this);
	isects.exchange12();
	return isects;
}
MGisects MGShell::intersection(const MGCurve& obj2)const{
	MGCFisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGShell::intersection(const MGFSurface& obj2)const{
	MGHHisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGShell::intersection(const MGSurface& obj2)const{
	MGHHisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGShell::intersection(const MGFace& obj2)const{
	MGHHisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}
MGisects MGShell::intersection(const MGShell& obj2)const{
	MGHHisect_vector isects2=isect(obj2);
	return MGisects(isects2);
}
