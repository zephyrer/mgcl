/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/SurfCurve.h"
#include "mg/FSurface.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mg/Group.h"
#include "Tl/TLDataVector.h"
#include "Tl/TLInputParam.h"

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//mgTLDataVector is a vector of mgTLData which is a tessellated data of one face or surface.
//If two (or more than two) faces share a binder edge, the faces are guaranteed to
//be tessallated from the same one straight line.

typedef std::pair<MGEdge*,MGCurve*> edge_curve;
typedef std::vector<edge_curve> edge_curves;

////////// constructor /////////////

//construct an mgTLDataVector from an object.
mgTLDataVector::mgTLDataVector(
	const MGObject& obj,	//obj must be a MGSurface, MGFace, or MGShell
	const mgTLInputParam& param//parameter for the tessellation.
){
	push_back(obj, param);
}

mgTLDataVector::mgTLDataVector(
	const MGObject& obj,	//obj must be a MGSurface, MGFace, or MGShell
	double crvTol,			//バウンダリのトレランス
	double surfTol,			//平面とみなすトレランス
	double max_ratio,		//最大アスペクト比
	MGCL::fan_kind fk,
		//fk=SINGLE_TRIANGLE:   1 triangle/FAN
		//fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri,		//Specify minimum number of triangles.
	double max_edge_len
){
	push_back(obj,crvTol,surfTol,max_ratio,fk,minimum_tri,max_edge_len);
}

//construct an mgTLDataVector from an object.
mgTLDataVector::mgTLDataVector(
	const MGGroup& group,	//From group, a MGSurface, MGFace, or MGShell will be extracted.
	double crvTol,			//バウンダリのトレランス
	double surfTol,			//平面とみなすトレランス
	double max_ratio,	//最大アスペクト比
	MGCL::fan_kind fk,
		//fk=SINGLE_TRIANGLE:   1 triangle/FAN
		//fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri,	//Specify minimum number of triangles.
	double max_edge_len	//when max_edge_len<=0, this means no limits on an edge length.
){
	push_back(group,crvTol,surfTol,max_ratio,fk,minimum_tri,max_edge_len);
}

//compute all the texture coordinates of this m_tlpoints.
//compute_texture_coordinates will find f in this vector,
//then compute texture coordinates of f.
//If f is a constituent face of a shell, all the textures of the faces of the shell
//will be  computed.
void mgTLDataVector::compute_texture_coordinates(
	const MGFSurface& f,	//face or surface of uv.
	const MGPosition& uv,	//surface parameter of the point st.
	const MGPosition& st,	//texture coordinate of uv.
	const MGUnit_vector& texXaxis,//world coordinate x axis vector.
 	double distortion_ratio//When >1. distorted points(points whose ratio are more
							//than distortion_ratio or less than 1/distortion_ratio will be
							//stored in m_tex_distorted.
){
	const MGSurface* surf=f.get_surface_pointer();
	//1. search the mgTLData where f belongs to.
	mgTLDataVector::iterator datai=begin(), dataie=end();
	for(; datai!=dataie; datai++){
		if(surf==&((*datai).surface()))
			break;
	}
	if(datai==dataie)
		return;//when not found, return.

	//2. set the texture coordinates of the found mgTLData.
	mgTLData& tldata=*datai;
	tldata.compute_texture_coordinates(uv,st,texXaxis,distortion_ratio);

	compute_shell_texture_coordinates(distortion_ratio);
}

void project_fp(
	const std::vector<MGFPoints>& fpointsVec,
				//vector of vector of MGFPoints for high priority line.
	MGHHisect_vector& hhisVec//Projected lines will be output.
){
	size_t nvec=fpointsVec.size();
	if(!nvec)
		return;

	const MGFPoints& fpoints0=fpointsVec[0];
	if(fpoints0.empty())
		return;

	const MGFPoint& fp00=fpoints0[0];
	const MGFSurface& f=fp00.fsurface();
	const MGFace* face=f.get_face_pointer();
	const MGShell* shl=0;
	if(face)
		shl=face->shell();

	//1. Convert fpointsVect to HHisect_vector by projecting the polylines onto the shell.
	if(shl){
		for(size_t i=0; i<nvec; i++){
			const MGFPoints& fpointsi=fpointsVec[i];

			MGHHisect hhis;
			size_t nfps=fpointsi.size();
			for(size_t j=1; j<nfps; j++){
				const MGFPoint& fpjm1=fpointsi[j-1];
				const MGFPoint& fpj=fpointsi[j];
				MGStraight sl(fpj.eval(),fpjm1.eval());
				MGHHisect_vector hhisj;
				shl->project(sl,hhisj);
				//std::cout<<std::endl<<" ***Before buidl one=="<<hhis<<std::endl;
				//std::cout<<hhisj<<std::endl;
				hhis.build_one(hhisj);
			}
			if(hhis.num_of_uvline())
				hhisVec.push_back(hhis);
		}
	}else{
		size_t nvec=fpointsVec.size();
		for(size_t i=0; i<nvec; i++){
			const MGFPoints& fpointsi=fpointsVec[i];

			MGHHisect hhis;
			size_t nfps=fpointsi.size();
			for(size_t j=1; j<nfps; j++){
				const MGFPoint& fpjm1=fpointsi[j-1];
				const MGFPoint& fpj=fpointsi[j];
				MGStraight sl(fpj.eval(),fpjm1.eval());
				MGPvector<MGCurve> uvcrvs, wcrvs;
				f.project(sl,uvcrvs,wcrvs);
				size_t ncrvs=uvcrvs.size();
				for(size_t k=0; k<ncrvs; k++){
					MGHHisect is(wcrvs.release(k),&f,uvcrvs.release(k));
					hhis.build_one(is);
				}
			}
			if(hhis.num_of_uvline())
				hhisVec.push_back(hhis);
		}
	}
}

#define MAX_PROPAGATION 4
//compute all the texture coordinates of this m_tlpoints.
//compute_texture_coordinates will find f in this vector,
//then compute texture coordinates of f.
//If f is a constituent face of a shell, all the textures of the faces of the shell
//will be  computed.
void mgTLDataVector::compute_texture_coordinates(
	const std::vector<MGFPoints>& fpointsVec,
				//vector of vector of MGFPoints for high priority line.
	const MGFSurface& f,	//face or surface of uv.
	const MGPosition& uv,	//surface parameter of the point st.
	const MGPosition& st,	//texture coordinate of uv.
	const MGUnit_vector& texXaxis,//world coordinate x axis vector.
 	double distortion_ratio//When >1. distorted points(points whose ratio are more
							//than distortion_ratio or less than 1/distortion_ratio will be
							//stored in m_tex_distorted.
){
	mgTLData* tld=find_tldata(&f);
	if(!tld)
		return;

	//1. convert MGFPoints to MGHHisect_vector.
	MGHHisect_vector hhisVec;
	project_fp(fpointsVec,hhisVec);
	//std::cout<<hhisVec<<std::endl;

	//2. set the texture coordinates of the found mgTLData partially.
	int hhisVsize=hhisVec.size();
	if(!hhisVsize){
		compute_texture_coordinates(f,uv,st,texXaxis,distortion_ratio);

		//6. set the texture coordinates of all the non-textured mgTLData
		//connected to texturedTLDs by shell structure.
		compute_shell_texture_coordinates(distortion_ratio);
		return;
	}

	//3. set the texture coordinates of rects connected by hhisVec.
	std::vector<mgTLData*> texturedTLDs;
	MGPosition uv2=uv;//initial surface parameter of the point st.
	MGPosition st2=st;//initial texture coordinate of uv.
	MGUnit_vector texXaxis2=texXaxis;//initial world coordinate x axis vector.
	for(int i=0; i<hhisVsize; i++){
		MGHHisect& hhisi=hhisVec[i];//std::cout<<hhisi<<std::endl;////////*******
		if(i){
			const MGFPline& fpl0=hhisi.uvline1(0);
			const MGFSurface* facei=fpl0.face();
			uv2=fpl0.uvline().start_point();
			uv2=facei->range(uv2);
			tld=find_tldata(facei);//The next face's tld.
			if(!tld)
				continue;
			if(!tld->get_texXaxis(uv2,st2,texXaxis2))
				continue;
		}
		int next=0;
		int hhisiSize=hhisi.num_of_uvline();
		while(true){
			int nextOld=next;
			if(tld->compute_texture_coordinates_partially(
				next,hhisi,uv2,st2,texXaxis2,distortion_ratio))
				texturedTLDs.push_back(tld);
			if(next>=hhisiSize)
				break;
			if(nextOld==next)
				break;

			const MGFPline& fpl_before=hhisi.uvline1(next-1);
			const MGFSurface* faceB=fpl_before.face();
			tld=find_tldata(faceB);//The next face's tld.
			if(!tld)
				break;
			MGPosition uv3=fpl_before.uvline().end_point();
			uv3=faceB->range(uv3);
			if(!tld->get_texXaxis(uv3,st2,texXaxis2))
				break;
			const MGFPline& fpl_next=hhisi.uvline1(next);
			const MGFSurface* faceN=fpl_next.face();
			uv2=fpl_next.uvline().start_point();
			uv2=faceN->range(uv2);
			tld=find_tldata(fpl_next.face());//The next face's tld.
			if(!tld)
				break;
		}
	}

	//4. Propagation to MAX_PROPAGATION.
	int level;
	for(level=3; level<=MAX_PROPAGATION; level++){
		std::vector<mgTLData*>::iterator itld=texturedTLDs.begin(),
											itlde=texturedTLDs.end();
		for(;itld!=itlde; itld++){
			(**itld).compute_level_texture(level);
		}
	}

	//5. texture all the mgTLData in texturedTLDs by raising level number.
	std::vector<mgTLData*>::iterator itld=texturedTLDs.begin(),
											itlde=texturedTLDs.end();
	for(;itld!=itlde; itld++){
		level=MAX_PROPAGATION+1;
		mgTLData& tldi=**itld;
		while(tldi.compute_level_texture(level))
			level++;
		tldi.set_textured(mgTLData::TOTALLY_TEXTURED);
	}

	//6. set the texture coordinates of all the non-textured mgTLData
	//connected to texturedTLDs by shell structure.
	compute_shell_texture_coordinates(distortion_ratio);
}

//This DataVector is a vector of members of a shell, and at least one of
//mgTLData(that is, at least one of the member face) is textured.
//Then compute_shell_texture_coordinates textures all of non-textured mgTLData.
void mgTLDataVector::compute_shell_texture_coordinates(
	double distortion_ratio
){
	// set the texture coordinates of all the non-textured mgTLData
	//connected to texturedTLDs by shell structure.
	int ntextured;
	do{
		ntextured=0;
		iterator datai=begin(), dataie=end();
		for(; datai!=dataie; datai++){
			if(!datai->textured()){
				if(datai->find_textured_and_compute(*this,distortion_ratio))
					ntextured++;
			}
		}
	}while(ntextured>0);
}

std::ostream& operator<<(std::ostream& out, const mgTLDataVector& datas){
	size_t n=datas.size();
	out<<"mgTLDataVector="<<(&datas)<<"::"<<std::endl;
	for(size_t i=0; i<n; i++) out<<datas[i];
	return out;
}

//Find mgTLData that is for the input face in this vector.
//If found the pointer will be returned,
//else null will be returned.
mgTLData* mgTLDataVector::find_tldata(const MGFSurface* face){
	const MGSurface* surf=face->get_surface_pointer();
	iterator datai=begin(), dataie=end();
	for(; datai!=dataie; datai++){
		if(surf==&((*datai).surface()))
			break;
	}
	if(datai==dataie)
		return 0;//when not found, return.

	mgTLData& tldi=*datai;
	return &tldi;
}

//Compute minimum and maximum curvature of MGCL::SURFACE_CURVATURE_KIND kind.
void mgTLDataVector::get_min_max_curvatures(
	MGCL::SURFACE_CURVATURE_KIND kind,
	double& minimum,
	double& maximum
)const{
	minimum=maximum=0.;
	const_iterator i=begin(), iend=end();
	if(i==iend)
		return;

	(*i++).get_min_max_curvatures(kind,minimum,maximum);
	for(;i!=iend;i++){
		double mini,maxi;
		(*i).get_min_max_curvatures(kind,mini,maxi);
		if(mini<minimum)
			minimum=mini;
		else if(maxi>maximum)
			maximum=maxi;
	}
}

//re-construct mgTLDataVector by the input parameters.
void mgTLDataVector::push_back(
	const MGObject& obj,//obj must be a MGSurface, MGFace, or MGShell
	double crvTol,		//バウンダリのトレランス
	double surfTol,		//平面とみなすトレランス
	double max_ratio,	//最大アスペクト比
	MGCL::fan_kind fk,
		//fk=SINGLE_TRIANGLE:   1 triangle/FAN
		//fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri,	//Specify minimum number of triangles.
	double max_edge_len	//when max_edge_len<=0, this means no limits on an edge length.
){
	mgTLInputParam param(crvTol,surfTol,max_ratio,fk,minimum_tri,max_edge_len);
	push_back(obj,param);
}

void mgTLDataVector::push_back(
	const MGObject& obj,	//obj must be a MGSurface, MGFace, or MGShell
	const mgTLInputParam& param//parameter for the tessellation.
){
	if(obj.manifold_dimension()!=2)
		return;//If not MGSurface, MGFace, or MGShell.

	const MGSurface* sf=obj.surf();
	const MGFace* f=0;
	if(!sf){
		f=obj.face();
		if(f)
			sf=f->surface();
	}
	if(sf){//if MFSurface or MGFace.
		if(f)
			m_datas.push_back(mgTLData(*f,param));
		else
			m_datas.push_back(mgTLData(*sf,param));
		return;
	}

	const MGShell* shell=dynamic_cast<const MGShell*>(&obj);
	if(!shell) return;
/*
	MGShell* pshell=const_cast<MGShell*>(shell);

	double crvTol=param.crvTol();
	//1. change the original edge curves to straight lines.
	double dPreLineTol = MGTolerance::set_line_zero(crvTol);
	MGComplex::bcellItr i=pshell->bcell_begin(), ie=pshell->bcell_end();
	edge_curves ecs;//original edge-curve will be saved in this ecs.
	//size_t ii=0;
	for(; i!=ie; i++){
		MGEdge* bedge=dynamic_cast<MGEdge*>(*i);
		//cout<<ii++<<":"<<*(bedge);
		if(!bedge) continue;

		MGCurve* bcrv=bedge->base_curve();
		ecs.push_back(edge_curve(bedge,bcrv));//save the original edge curve.

		const MGEdge* pedge=bedge->member_partner_edge(0);
		MGTrimmedCurve pcrv=pedge->trimmed_curve();
		double pts=pcrv.param_s(), pte=pcrv.param_e();
		MGPosition P0=pedge->eval_star(pts), P1=pedge->eval_star(pte);
		if(P0==P1){
			bedge->m_extent=new MGStraight(P1,P0);
		}else if(bcrv){
			double ts=bedge->param_s(), te=bedge->param_e();
			MGLBRep* lb=new MGLBRep(bedge->trimmed_curve(), 2);
			lb->change_range(ts,te);
			bedge->m_extent=lb;
		}else{
			//If binder edge did not have curve representation.
			assert(pedge->face());assert(pedge->face()->surface());
			const MGSurface* srf=pedge->face()->surface();
			bedge->m_extent=new MGLBRep(MGSurfCurve(*srf,pedge->trimmed_curve()),2);
		}
	}
	MGTolerance::set_line_zero(dPreLineTol);		
*/
	//2. perform tessellation to each face.
	//size_t k = m_datas.size();
	//m_datas.reserve(shell->number_of_faces());
	MGComplex::const_pcellItr j=shell->pcell_begin(), je=shell->pcell_end();
	for(; j!=je; j++){
		const MGFace* ff=static_cast<const MGFace*>(*j);
		m_datas.push_back(mgTLData(*ff,param));
	}
/*
	//3. restore the original edge curves.
	size_t necs=ecs.size();
	for(k=0; k<necs; k++){
		MGEdge* bedge=ecs[k].first;
		delete bedge->m_extent;
		bedge->m_extent=ecs[k].second;
	}
*/
}

//construct mgTLData by the input parameters and push backt the data.
void mgTLDataVector::push_back(
	const MGGroup& group,	//From group, a MGSurface, MGFace, or MGShell will be extracted.
	const mgTLInputParam& param//parameter for the tessellation.
){
	MGGroup::const_iterator i=group.begin(), ie=group.end();
	for(; i!=ie; i++){
		const MGObject* obj=dynamic_cast<const MGObject*>(*i);
		if(obj) push_back(*obj,param);
		else{
			const MGGroup* grpi=dynamic_cast<const MGGroup*>(*i);
			if(grpi) push_back(*grpi,param);
		}
	}
}

//construct mgTLData by the input parameters and push back the data.
void mgTLDataVector::push_back(
	const MGGroup& group,//From group, a MGSurface, MGFace, or MGShell will be extracted.
	double crvTol,		//バウンダリのトレランス
	double surfTol,		//平面とみなすトレランス
	double max_ratio,	//最大アスペクト比
	MGCL::fan_kind fk,
		//fk=SINGLE_TRIANGLE:   1 triangle/FAN
		//fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri,	//Specify minimum number of triangles.
	double max_edge_len	//when max_edge_len<=0, this means no limits on an edge length.
){
	MGGroup::const_iterator i=group.begin(), ie=group.end();
	for(; i!=ie; i++){
		const MGObject* obj=dynamic_cast<const MGObject*>(*i);
		if(obj) push_back(*obj,crvTol,surfTol,max_ratio,fk,minimum_tri,max_edge_len);
		else{
			const MGGroup* grpi=dynamic_cast<const MGGroup*>(*i);
			if(grpi) push_back(*grpi,crvTol,surfTol,max_ratio,fk,minimum_tri,max_edge_len);
		}
	}
}

//Compute maximum  width and height out of all the tex_planes.
void mgTLDataVector::tex_plane_maximum(double& width, double& height)const{
	width=height=0.;
	mgTLDataVector::const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		double wi, hi;
		i->tex_plane_maximum(wi,hi);
		if(wi>width)
			width=wi;
		if(hi>height)
			height=hi;
	}
}
