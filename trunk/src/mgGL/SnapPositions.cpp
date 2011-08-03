/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Plane.h"
#include "mg/Curve.h"
#include "mg/Point.h"
#include "mg/Group.h"
#include "mg/PickObjects.h"
#include "topo/BVertex.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGL/SnapPositions.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGSnapPositions Class.
//MGSnapPositions is a class to store array(vector) of MGPosition's,
//used for MGLocateTool specifically.
//

//append this positions in m_positions into points.
void MGSnapPositions::append_position_data(std::vector<MGPosition>& points)const{
	points.insert(points.end(), m_positions.begin(), m_positions.end());
}

//Extract position data.
void MGSnapPositions::extract(
	const MGPoint& point	//the curve to extract.
){
	switch(m_snap_kind){
	case endpos:
		m_obj_nums.push_back(obj_num(&point,1));
		m_positions.push_back(point.position());
		return;
	default:
		;
	}
}

//Extract position data.
void MGSnapPositions::extract(
	const MGCurve& crv	//the curve to extract.
){
	switch(m_snap_kind){
	case endpos:
		m_obj_nums.push_back(obj_num(&crv,2));
		m_positions.push_back(crv.start_point());
		m_positions.push_back(crv.end_point());
		return;
	case knotpos:
		{
		const MGKnotVector& t=crv.knot_vector();
		size_t k=t.order(), n=t.bdim();
		int numpoint=0;
		for(size_t i=k; i<n; i++){
			if(t[i]==t[i-1]) continue;
			m_positions.push_back(crv.eval(t[i]));
			numpoint++;
		}
		m_obj_nums.push_back(obj_num(&crv,numpoint));
		return;
		}
	case centerpos:
		m_obj_nums.push_back(obj_num(&crv,1));
		m_positions.push_back(crv.center());
		return;
	default:
		;
	}
}

//Extract position data.
void MGSnapPositions::extract(
	const MGSurface& srf//the surface to extract.
){
	double u0=srf.param_s_u(), u1=srf.param_e_u();
	double v0=srf.param_s_v(), v1=srf.param_e_v();
	switch(m_snap_kind){
	case endpos:
		{
		const MGPlane* pl=dynamic_cast<const MGPlane*>(&srf);
		if(pl) return;
		m_obj_nums.push_back(obj_num(&srf,4));
		m_positions.push_back(srf.eval(u0,v0));
		m_positions.push_back(srf.eval(u0,v1));
		m_positions.push_back(srf.eval(u1,v0));
		m_positions.push_back(srf.eval(u1,v1));
		return;
		}

	case knotpos:
		{
		const MGKnotVector& tu=srf.knot_vector_u();
		size_t ku=tu.order(), nu=tu.bdim();
		const MGKnotVector& tv=srf.knot_vector_v();
		size_t kv=tv.order(), nv=tv.bdim();
		int numpoint=0;
		for(size_t i=ku; i<nu; i++){
			if(tu[i]==tu[i-1]) continue;
			m_positions.push_back(srf.eval(tu[i],v0));
			m_positions.push_back(srf.eval(tu[i],v1));
			numpoint+=2;
		}
		for(size_t j=kv; j<nv; j++){
			if(tv[j]==tv[j-1]) continue;
			m_positions.push_back(srf.eval(u0,tv[j]));
			m_positions.push_back(srf.eval(u1,tv[j]));
			numpoint+=2;
		}
		m_obj_nums.push_back(obj_num(&srf,numpoint));
		return;
		}

	case centerpos:
		m_obj_nums.push_back(obj_num(&srf,1));
		m_positions.push_back(srf.center());
		return;

	default:
		;
	}
}

//Extract vertex of a topology.
void MGSnapPositions::extract(
	const MGFace& face	//the curve to extract.
){
	switch(m_snap_kind){
	case vertexpos:
		{
		size_t num_vertex=0;
		size_t nlp=face.number_of_boundaries();
		for(size_t i=0; i<nlp; i++){
			const MGLoop& lp=*(face.loop(i));
			MGLoop::const_bcellItr j=lp.bcell_begin(), je=lp.bcell_end();
			for(; j!=je; j++){
				const MGBVertex* v=dynamic_cast<const MGBVertex*>(*j);
				if(v){
					MGPosition uv=v->position();
					m_positions.push_back(face.eval(uv));
					num_vertex++;
				}
			}
		}
		m_obj_nums.push_back(obj_num(&face,num_vertex));
		return;
		}
	default:
		;
	}
}

void MGSnapPositions::extract(
	const MGShell& shell	//the surface to extract.
){
	switch(m_snap_kind){
	case vertexpos:
		{
		MGShell::const_pcellItr i=shell.pcell_begin(), ie=shell.pcell_end();
		for(; i!=ie; i++){
			const MGFace* fi=shell.face(i);
			extract(*fi);
		}
		return;
		}
	default:
		;
	}
}

void MGSnapPositions::extract(
	const MGGel& gel	//the gel to extract.
){
	const MGCurve* crv=gel.curve();
	if(crv){ extract(*crv); return;}

	const MGFace* f=gel.face();
	if(f){ extract(*f); return;}

	const MGSurface* srf=gel.surf();
	if(srf){ extract(*srf); return;}

	const MGPoint* P=gel.point();
	if(P){	extract(*P); return;}

	const MGShell* shl=gel.shell();
	if(shl){  extract(*shl); return;}

	const MGGroup* group=gel.group();
	if(group){ extract(group->m_gels.m_list); return;}
}

void MGSnapPositions::extract(
	const std::list<MGGel*>& gel_list	//the group to extract.
){
	std::list<MGGel*>::const_iterator i=gel_list.begin(), ie=gel_list.end();	
	for(; i!=ie; i++) extract(**i);
}

void MGSnapPositions::extract(
	const MGPickObjects& pobjs	//array of pick objects to extract.
){
	MGPickObjects::const_iterator i=pobjs.begin(), ie=pobjs.end();
	for(; i!=ie; i++) extract(*((**i).gel()));
}

void MGSnapPositions::generate_point_display_list()const{
	size_t nm=dlist_name();
	glNewList(nm, GL_COMPILE);
		glPointSize(1.);
		std::vector<obj_num>::const_iterator i=m_obj_nums.begin(), ie=m_obj_nums.end();
		MGSnapPositions::const_iterator j=m_positions.begin(), je=m_positions.end();
		for(; i!=ie; i++){
			int n=i->second;//number of points in this object.
			if(!n) continue;
			nm=i->first->dlist_name();//MGObject*
			glPushName(nm);
			for(int k=0; k<n; k++, j++){
				nm=size_t(j->dlist_name());//MGPosition*
				glPushName(nm);
				glBegin(GL_POINTS); glVertex3dv((*j).data()); glEnd();
				glPopName();
			}
			glPopName();
		}
	glEndList();
}

void MGSnapPositions::get_pick_data(
	const unsigned int* selectBuf, 
	MGPosition& point,	//point data will be output.
	const MGObject*& obj,//When function's return value is nearpos, end, or knot,
				//the point's parameter value of the object will be returned.
	double& t	//When center, only MGObject pointer will be rturned.
)const{
	obj=reinterpret_cast<const MGObject*>(selectBuf[3]);
	const MGPosition* pos=reinterpret_cast<const MGPosition*>(selectBuf[4]);
	point=*pos;
	const MGCurve* crv=obj->curve();
	if(!crv) return;

	switch(m_snap_kind){
	case endpos:
		{
		MGVector Ps=crv->start_point(), Pe=crv->end_point();
		MGVector dif1=Ps-point, dif2=Pe-point;
		if(dif1%dif1<= dif2%dif2) t=crv->param_s();
		else t=crv->param_e();
		return;
		}
	case knotpos:
		{
		const MGKnotVector& tv=crv->knot_vector();
		size_t k=tv.order(), n=tv.bdim();
		double error=MGTolerance::mach_zero();
		for(size_t i=k; i<n; i++){
			t=tv[i];
			MGVector dif=crv->eval(t)-point;
			if(dif%dif<error) return;
		}
		return;
		}
	default:
		;
	}
}
