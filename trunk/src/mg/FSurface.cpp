/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"
#include "mg/CParam_list.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"

using namespace std;

//Comparison operator.
bool MGFSurface::operator< (const MGFSurface& f2) const{
	const MGFace* face1=dynamic_cast<const MGFace*>(this);
	const MGFace* face2=dynamic_cast<const MGFace*>(&f2);
	if(face1 && face2)
		return (*face1)<(*face2);

	if((!face1) && (!face2))
		return *(get_surface_pointer())< *(f2.get_surface_pointer());

	return false;
}

//Get the box of the object.
const MGBox& MGFSurface::get_box() const{
	const MGObject* obj=object_pointer();
	return obj->box();
}

/**
 *  @brief eval_discrete_deviationの下請け関数
 *  @param face1 一方の面データ
 *  @param face2 もう一方の面データ
 *  @param wcrv1 face1側エッジのワールドカーブ(must be trimmed into evaluatin range).
 *  @param wcrv2 face2側エッジのワールドカーブ(must be trimmed into evaluatin range).
 *  @param pcrv1 face1側エッジの面上パラメータカーブ(must be the same range as wcrv1).
 *  @param pcrv2 face2側エッジの面上パラメータカーブ(must be the same range as wcrv2).
 *
 *  wcrv1とwcrv2はほぼ共線となっている。
 */
void deviation(
	const MGFSurface& face1,
	const MGFSurface& face2,
	const MGCurve&    wcrv1,
	const MGCurve&    wcrv2,
	const MGCurve&    pcrv1,
	const MGCurve&    pcrv2,
	int npoint,		//num of discrete points.
	std::vector<MGPosition>& uvuvs
){
	MGPosition wcrv1s = wcrv1.start_point(), wcrv1e = wcrv1.end_point();
	MGPosition uv1s=pcrv1.start_point(), uv1e=pcrv1.end_point();
	bool w1_equal_to_p1=face1.get_surface_pointer()->equal_direction(pcrv1,wcrv1)==1;

	MGPosition wcrv2s = wcrv2.start_point(), wcrv2e = wcrv2.end_point();
	MGPosition uv2s=pcrv2.start_point(), uv2e=pcrv2.end_point();
	bool w2_equal_to_p2=face2.get_surface_pointer()->equal_direction(pcrv2,wcrv2)==1;
	//cout << "端点:\n";
	//cout << "f1S:uv="<<uv1s<<",world="<<wcrv1s<<",eval="<< face1.eval(uv1s)<<",normal="<<face1.unit_normal(uv1s)<<endl;
	//cout << "f1E:uv="<<uv1e<<",world="<<wcrv1e<<",eval="<< face1.eval(uv1e)<<",normal="<<face1.unit_normal(uv1e)<<endl;
	//cout << "f2S:uv="<<uv2s<<",world="<<wcrv2s<<",eval="<< face2.eval(uv2s)<<",normal="<<face2.unit_normal(uv2s)<<endl;
	//cout << "f2E:uv="<<uv2e<<",world="<<wcrv2e<<",eval="<< face2.eval(uv2e)<<",normal="<<face2.unit_normal(uv2e)<<endl;
	//cout<<" p,w crv 1="<<pcrv1<<wcrv1<<endl;
	//cout<<" p,w crv 2="<<pcrv2<<wcrv2<<endl;

	MGPosition uvuv(4);
	if(w1_equal_to_p1)
		uvuv.store_at(0,uv1s);
	else
		uvuv.store_at(0,uv1e);
	const double lenss = (wcrv1s-wcrv2s).len();
	const double lense = (wcrv1s-wcrv2e).len();
	if(lenss<lense){
		if(w2_equal_to_p2)
			uvuv.store_at(2,uv2s);
		else
			uvuv.store_at(2,uv2e);
	}else{
		if(w2_equal_to_p2)
			uvuv.store_at(2,uv2e);
		else
			uvuv.store_at(2,uv2s);
	}
	uvuvs.push_back(uvuv);

	double pcrv1prms = pcrv1.param_s(), pcrv1prme = pcrv1.param_e();
	double pcrv2prms = pcrv2.param_s(), pcrv2prme = pcrv2.param_e();
	
	double dwt1=wcrv1.param_span()/npoint, dwt2, 
			dpt1=(pcrv1prme-pcrv1prms)/npoint, dpt2, 
			wt1= wcrv1.param_s(), wt2, 
			pt1 = pcrv1prms, pt2;
	if(!w1_equal_to_p1){
		dpt1*=-1.;
		pt1=pcrv1prme;
	}

	if(lenss < lense){
		wt2 = wcrv2.param_s();
		pt2 =w2_equal_to_p2?pcrv2prms:pcrv2prme;
		dwt2 = wcrv2.param_span()/npoint;
		dpt2 = (pcrv2prme - pcrv2prms)/npoint;
		if(!w2_equal_to_p2)
			dpt2*=-1.;
	}else{
		wt2 = wcrv2.param_e();
		pt2 =w2_equal_to_p2?pcrv2prme:pcrv2prms;
		dwt2 = - wcrv2.param_span() / npoint;
		dpt2 = (pcrv2prms - pcrv2prme) / npoint;
		if(!w2_equal_to_p2)
			dpt2*=-1.;
	}

	// pos1, pos2の、それぞれの最初と最後の要素以外のすべての要素を計算する
	const MGSurface& srf1=*(face1.get_surface_pointer());
	const MGSurface& srf2=*(face2.get_surface_pointer());
	for(int i = 1; i < npoint; i++){
		wt1 += dwt1; wt2 += dwt2;
		pt1 += dpt1; pt2 += dpt2;
		MGPosition P1=wcrv1.eval(wt1);
		if(!wcrv2.perp_guess(0.,-1., P1, wt2, wt2)){
			wt2 = wcrv2.closest(P1);
		}

		pt1=srf1.param_of_pcurve(wt1,wcrv1,pcrv1,&pt1);
		pt2=srf2.param_of_pcurve(wt2,wcrv2,pcrv2,&pt2);
		uvuv.store_at(0,pcrv1.eval(pt1),0,2);
		uvuv.store_at(2,pcrv2.eval(pt2),0,2);
		uvuvs.push_back(uvuv);
	}
	// pos1, pos2それぞれの最後の点を決定する
	if(w1_equal_to_p1)
		uvuv.store_at(0,uv1e,0,2);
	else
		uvuv.store_at(0,uv1s,0,2);
	if(lenss<lense){
		if(w2_equal_to_p2)
			uvuv.store_at(2,uv2e);
		else
			uvuv.store_at(2,uv2s);
	}else{
		if(w2_equal_to_p2)
			uvuv.store_at(2,uv2s);
		else
			uvuv.store_at(2,uv2e);
	}
	uvuvs.push_back(uvuv);
}

//Evaluate deviations of two faces(this and face2) at npoint
//discrete points.
//(1)Search the common edges which have the distance within tolerance.
//(2)Compute the nearest points from npoint discrete points of this to face2.
//Let uvuvi=uvuvs[i], then
//uvuvi[0], [1] are this face's parameter value(u1,v1), and uvuvi[2], [3] are
//parameter value(u2,v2) of face2 which is the nearest point from the point (u1, v1).
void MGFSurface::eval_discrete_deviation(
	const MGFSurface& face2,
	std::vector<MGPosition>& uvuvs,
	int npoint,		//indicates how many discrete points be obtained.
	double tolerance//tolerance to get two edge to compute deviation.
)const{

	const MGFSurface& face1=*this;
	const MGSurface& srf1=*(face1.get_surface_pointer());
	const MGSurface& srf2=*(face2.get_surface_pointer());

	// 1. faceの境界線を取得する。
	MGPvector<MGCurve> outer1 = face1.outer_boundary();
	MGPvector<MGCurve> outer2 = face2.outer_boundary();
	MGPvector<MGCurve> outerp1 = face1.outer_boundary_param();
	MGPvector<MGCurve> outerp2 = face2.outer_boundary_param();
	/*for(size_t iii=0; iii<4; iii++){
		std::cout<<iii<<" 1W="<<*(outer1[iii])<<std::endl;
		std::cout<<iii<<" 1P="<<*(outerp1[iii])<<std::endl;
		std::cout<<iii<<" 2W="<<*(outer2[iii])<<std::endl;
		std::cout<<iii<<" 2P="<<*(outerp2[iii])<<std::endl;
	}*/

	// 2. エッジ同士の離れをチェックし、近いようであれば離れ計算の対象にする。
	// ループスタート
	MGTolerance& tol = MGTolerance::instance();
	for(size_t i=0; i<outer1.size(); i++){
		// face1側のエッジのworld curve
		const MGCurve& wcrv1 = *outer1[i];
		const MGCurve& pcrv1 = *outerp1[i];
		MGPosition crv1s = wcrv1.start_point();
		MGPosition crv1e = wcrv1.end_point();
		for(size_t j=0; j<outer2.size(); j++){
			// face2側のエッジのworld curve
			const MGCurve& wcrv2 = *outer2[j];
			const MGCurve& pcrv2 = *outerp2[j];

			std::vector<double> spans;
			// 一時的にトレランスを緩くしてからMGCurve::common関数を呼ぶ。
			double savedWorld = tol.set_wc_zero(tolerance);
			double savedLine = tol.set_line_zero(tolerance);
			int cmnret = wcrv1.common(wcrv2, spans);
			tol.set_line_zero(savedLine);
			tol.set_wc_zero(savedWorld);

			// 共線でなければ次のエッジの組み合わせへジャンプ。
			if(cmnret <= 0)
				continue;

			int k, nspan=spans.size()/4;
			double totalLen=fabs(spans[1]-spans[0]);
			for(k=1; k<nspan; k++)
				totalLen+=fabs(spans[4*k+1]-spans[4*k]);

			for(k=0; k<nspan; k++){
				int fk=4*k;
				// face1側エッジの共線部分曲線を取り出す
				double sp0 = spans[fk], sp1 = spans[fk+1];
				MGTrimmedCurve crv1trmd(wcrv1,sp0,sp1);
				double spara1=srf1.param_of_pcurve(sp0,wcrv1,pcrv1);
				double epara1=srf1.param_of_pcurve(sp1,wcrv1,pcrv1);
				MGTrimmedCurve pcrv1trmd(pcrv1,spara1,epara1);

				// face2側エッジの共線部分曲線を取り出す
				double sp2 = spans[fk+2], sp3 = spans[fk+3];
				MGTrimmedCurve crv2trmd(wcrv2,sp2,sp3);
				double spara2=srf2.param_of_pcurve(sp2,wcrv2,pcrv2);
				double epara2=srf2.param_of_pcurve(sp3,wcrv2,pcrv2);
				MGTrimmedCurve pcrv2trmd(pcrv2,spara2,epara2);

				// 補助関数に丸投げ
				int npointk=int(npoint*fabs(sp1-sp0)/totalLen)+1;
				deviation(face1,face2,crv1trmd,crv2trmd,pcrv1trmd,pcrv2trmd,npointk,uvuvs);
			}
		}
	}
}

//Obtain all the boundaries(i.e., outer boundary and all the inner boundaries)
MGPvector<MGCurve> MGFSurface::get_all_boundaries(void)const{
	MGPvector<MGCurve> crvs=outer_boundary();
	size_t n=number_of_inner_boundaries();
	for(size_t i=0; i<n; i++){
		crvs.push_back(inner_boundary(i));
	}
	return crvs;
}

//Obtain parameter space error.
double MGFSurface::param_error() const{
	const MGSurface& f=*(get_surface_pointer());
	return f.param_error();
}
double MGFSurface::param_error_u() const{
	const MGSurface& f=*(get_surface_pointer());
	return f.param_error_u();
}
double MGFSurface::param_error_v() const{
	const MGSurface& f=*(get_surface_pointer());
	return f.param_error_v();
}

//指定点から最も近い、垂線の足とパラメータ値を返す。
//Return the foot of the perpendicular straight line from p that is 
//nearest to point P.
//Function's return value is whether point is obtained(>0) or not(0)
int MGFSurface::perp_one(
	const MGPosition& P, // 指定点(point)
	MGPosition& uv 		//Parameter value of the surface will be returned.
)const{
	MGPosition_list list=perps(P);
	size_t n=list.entries();
	if(n){	//Compute the nearest point.
		MGPosition_list::iterator i=list.begin();
		uv=*i++;
		MGPosition uv2;
		double dist=(eval(uv)-P).len(), dist2;
		for(;i!=list.end();i++){
			uv2=(*i);
			dist2=(eval(uv2)-P).len();
			if(dist2<dist){uv=uv2; dist=dist2;}
		}
	}
	return n;
}

//Obtain main parameter lines of the FSurface without boundaries.
//inner_skeleton includes only inner parameter lines without boundaries.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
MGPvector<MGCurve> MGFSurface::inner_skeleton(int density) const{
	MGPvector<MGCurve> crv_list;
	if(density>0){
		if(density>10)
			density=10;
		MGBox prange = param_range();
		const MGInterval& uspan=prange[0];
		const MGInterval& vspan=prange[1];
		double u0=uspan[0].value(), u1=uspan[1].value();
		double ulen=u1-u0;
		double v0=vspan[0].value(), v1=vspan[1].value();
		double vlen=v1-v0;
		double divider=double(density+1);
		for(int i=1; i<=density; i++){
			double ith=double(i)/divider;
			double u=u0+ulen*ith;
			crv_list.push_back(parameter_curves(1,u));
			double v=v0+vlen*ith;
			crv_list.push_back(parameter_curves(0,v));
		}
	}
	return crv_list;
}

//Obtain boundary and main parameter lines of the FSurface.
//skeleton includes boundary() and inner parameter lines.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
MGPvector<MGCurve> MGFSurface::skeleton(int density) const{
	MGPvector<MGCurve> crv_list=get_all_boundaries();
	if(density<0)
		density=0;

	crv_list.push_back(inner_skeleton(density));
	return crv_list;
}

//Obtain all the parameter curves at knots of u and v knot vector.
MGPvector<MGCurve> MGFSurface::skeleton_at_knots()const{
	MGPvector<MGCurve> crvs=get_all_boundaries();
	const MGSurface* srf=get_surface_pointer();
	const MGBox& pbox=box_param2();
	double u0=pbox[0].low_point(), u1=pbox[0].high_point();
	double v0=pbox[1].low_point(), v1=pbox[1].high_point();

	//all the parameter lines are necessary.
	const MGKnotVector& tu=srf->knot_vector_u();
	if(tu!=mgNULL_KNOT_VECTOR){
		size_t num1=tu.bdim()-1;
		size_t ku=tu.order();
		double uold=tu[ku]-1.;
		for(size_t i=ku; i<=num1; i++){
			double u=tu[i];
			if(u==uold)
				continue;
			if(u0<u && u<u1){
				MGPvector<MGCurve> pcrvsi=parameter_curves(1,u);
				crvs.push_back(pcrvsi);
				uold=u;
			}
		}
	}
	
	const MGKnotVector& tv=srf->knot_vector_v();
	if(tv!=mgNULL_KNOT_VECTOR){
		size_t nvm1=tv.bdim()-1;
		size_t kv=tv.order();
		double vold=tv[kv]-1.;
		for(size_t i=kv; i<=nvm1; i++){
			double v=tv[i];
			if(v==vold)
				continue;
			if(v0<v && v<v1){
				MGPvector<MGCurve> pcrvsi=parameter_curves(0,v);
				crvs.push_back(pcrvsi);
				vold=v;
			}
		}
	}
	return crvs;
}

void trimProject(
	const std::vector<const MGCurve*>& trimmers,//Trimmer curves
	const MGVector&         dir, // dir == mgNULL_VEC means pull-projection
	const MGFSurface&       surf,
	MGPvector<MGCurve>&     result_world,
	MGPvector<MGCurve>&     result_param
){
	std::vector<const MGCurve*>::const_iterator i=trimmers.begin(), iend=trimmers.end();
	for(; i != iend; ++i){
		surf.project(**i, result_param, result_world, dir);
			// dir == mgNULL_VEC means pull-projection
	}
}

//Build networks of surf, given parameter curves vector.
void build_networks(
	const MGFSurface& surf,		//The objective surface
	const MGPvector<MGCurve>& pcurves,//(u,v) 2D parameter curves of surf.
	MGPvector<MGLoop>& networks	//Built networks
){
	double err=surf.param_error()*8.;
	int n=pcurves.size();
	std::vector<const MGCurve*> pcurves2(n);
	for(int i=0; i<n; i++)
		pcurves2[i]=pcurves[i];

	while(true){
		int npcurves=pcurves2.size();
		if(!npcurves)
			break;;

		//Search non processed curve.
		const MGCurve* curve=0;
		while(curve==0 && npcurves){
			npcurves--;
			curve=pcurves2[npcurves];pcurves2.pop_back();
		}
		if(!curve)
			break;//If non processed curves not found.

		
		int j, nnet=networks.size();
		for(j=nnet-1; j>=0; j--){
			MGLoop* netj=networks[j];
			if(netj->merge_network(*curve)){
				break;
			}
		}
		if(j>=0)
			continue;

		std::auto_ptr<MGLoop> networki(new MGLoop(new MGEdge(curve->clone())));
		//std::cout<<(*curve)<<std::endl;
		double errsave=MGTolerance::set_wc_zero(err);
		if(curve->start_point()==curve->end_point())
			networki->make_close();
		MGTolerance::set_wc_zero(errsave);
		npcurves=pcurves2.size();
		for(int i=npcurves-1; i>=0; i--){
			if(pcurves2[i]==0)
				continue;
			bool mergedi=networki->merge_network(*(pcurves2[i]));
			if(mergedi){
				pcurves2[i]=0;
			}
		}

		//std::cout<<*networki<<std::endl;
		networks.push_back(networki.release());
	}

	int nnet=networks.size();
	for(int j=nnet-1; j>=0; j--){
		MGLoop* netj=networks[j];
		netj->remove_pendent_edge(surf);
		if(netj->number_of_edges()==0)
			networks.removeAt(j);
	}
}

//Trim this fsurface with trimmers. trimmers are 3D curves and will be projected
//onto this surface tword the direction dir. If dir is null vector, surface normal
//prjection will be performed. Trimming is so performed that the smallest region enclosed
//by trimmers that includes the surface point uv will be removed. 
void MGFSurface::trim(
	const std::vector<const MGCurve*>& trimmers,//Trimmer curves
	const MGVector&  dir,	//trimmers projection direction.
	const MGPosition& uv,	//surface parameter (u,v) that indicates the region to remove.
							//The smallest region that inclued uv will be removed.
	MGPvector<MGFace>& faces//Result trimmed face(s) will be appended.
			//If no trimming was performed, no faces will be appended.
)const{
	// 1 - get surf-curve by project
	MGPvector<MGCurve> wcurves, pcurves;
	trimProject(trimmers,dir,*this,wcurves,pcurves);
	//for(size_t iii=0; iii<trimmers.size(); iii++)std::cout<<*(trimmers[iii])<<std::endl;//////////////////::::::*******
	//for(size_t iii=0; iii<pcurves.size(); iii++)std::cout<<*(pcurves[iii])<<std::endl;//////////////////::::::*******

	MGPvector<MGLoop> networks;	//Built networks
	build_networks(*this,pcurves,networks);
	//std::cout<<networks<<std::endl;

	MGFace face(*this);
	face.remove_inactive_loops();
	face.make_outer_boundary();
	face.trim(networks,uv,faces);
}

//Extract a sub surface with trimmers. trimmers are 3D curves and will be projected
//onto this surface tword the direction dir. If dir is null vector, surface normal
//prjection will be performed. Extraction is so performed that the smallest region
//enclosed by trimmers that includes the surface point uv will be extracted. 
void MGFSurface::extract(
	const std::vector<const MGCurve*>& trimmers,	//Trimmer curves
	const MGVector&  dir,	//trimmers projection direction.
	const MGPosition& uv,	//surface parameter (u,v) that indicates the region to extract.
							//The smallest region that inclued uv will be extracted.
	std::auto_ptr<MGFace>& eface//Result extracted face will be output.
)const{
	// 1 - get surf-curve by project
	MGPvector<MGCurve> wcurves, pcurves;
	trimProject(trimmers,dir,*this,wcurves,pcurves);
	//for(size_t iii=0; iii<trimmers.size(); iii++)std::cout<<*(trimmers[iii])<<std::endl;//////////////////::::::*******
	//for(size_t iii=0; iii<pcurves.size(); iii++)std::cout<<*(pcurves[iii])<<std::endl;//////////////////::::::*******

	MGPvector<MGLoop> networks;	//Built networks
	build_networks(*this,pcurves,networks);

	MGFace face(*this);
	face.remove_inactive_loops();
	face.make_outer_boundary();
	face.extract_sub_face(networks,uv,eface);
}

///split this fsurface with splitters. splitters are 2D (u,v) surfaces's parameter curves.
void MGFSurface::split(
	const std::vector<const MGCurve*>& splitters,	//splitter world curves.
	const MGVector&  dir,	//splitter projection direction.
							//If dir.is_null(), normal projection will be performed.
	MGPvector<MGFace>& faces//Result splitted face(s) will be appended.
			//If no splitting was performed, no faces will be appended.
)const{
	// 1 - get surf-curve by project
	MGPvector<MGCurve> wcurves, pcurves;
	trimProject(splitters,dir,*this,wcurves,pcurves);
	split(pcurves,faces);
}

///split this fsurface with splitters. splitters are 2D (u,v) surfaces's parameter curves.
void MGFSurface::split(
	const MGPvector<MGCurve>& splitters,//splitter (u,v) curves.
	MGPvector<MGFace>& faces//Result splitted face(s) will be appended.
			//If no splitting was performed, no faces will be appended.
)const{
	MGPvector<MGLoop> networks;	//Built networks
	build_networks(*this,splitters,networks);

	MGFace face(*this);
	face.remove_inactive_loops();
	face.make_outer_boundary();
	face.split(networks,faces);
}
