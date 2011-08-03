/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Point.h"
#include "mg/Plane.h"
#include "mg/DrawFunc.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "Tl/TLInputParam.h"
#include "Tl/TLDataVector.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/GLDrawFunc.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Implements the drawWire functions of all the classes.

void MGComplex::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGCellNB& celli=**i;
		const MGObject* obj=celli.object();
		size_t dlname=obj->dlist_name();
		glPushName(dlname);
		celli.drawWire(span_length,line_density);	
		glPopName();
	}
}

//Draw 3D curve in the topology's star cell world coordinates.
//The object is converted to curve(s) and is drawn.
void MGComplex::drawWire_in_star(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	if(!star()) return;

	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGCellNB& celli=**i;
		//const MGObject* obj=celli.object();
		//glLoadName(obj->dlist_name());
		celli.drawWire_in_star(span_length,line_density);	
	}
}

//Draw 3D point(vertex) in world coordinates.
//The object is converted to point(s) and is drawn.
void MGComplex::draw3DVertex(
)const{
	const_bcellItr i=bcell_begin(), ie=bcell_end();
	for(; i!=ie; i++){
		const MGCellNB& celli=**i;
		const MGObject* obj=celli.object();
		glLoadName(obj->dlist_name());
		celli.draw3DVertex();	
	}
}

//Draw 3D point(vertex) in star cell's world coordinates.
//The object is converted to point(s) and is drawn.
void MGComplex::draw3DVertex_in_star(
)const{
	const_pcellItr i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGCellNB& celli=**i;
		//const MGObject* obj=celli.object();
		//glLoadName(obj->dlist_name());
		celli.draw3DVertex_in_star();
	}
}

//////////////////////////////////////////////

//Draw 3D curve in the topology's star cell world coordinates.
//The object is converted to curve(s) and is drawn.
void MGCellBase::drawWire_in_star(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	if(!star()) return;
	if(manifold_dimension()<=0) return;
	const MGCellNB* bndr=make_binder_with_extent();
	bndr->drawWire(span_length,line_density);
}

//Draw 3D curve in the topology's star cell world coordinates.
//The object is converted to curve(s) and is drawn.
void MGCellBase::draw3DVertex_in_star(
)const{
	if(!star()) return;
	if(manifold_dimension()>0) return;
	const MGCellNB* bndr=make_binder_with_extent();
	bndr->draw3DVertex_in_star();
}

//Draw the vertex in the parameter space of the star edge.
//We expect this drawWire will not be used by any chance.
void MGPVertex::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
//	double tau=t();
//	(*(MGDrawFunc::pfunc()))(tau,0.,0.);
}

//Draw 3D point(vertex) in world coordinates.
void MGPVertex::draw3DVertex(
)const{
	double tau=t();
	(*(MGDrawFunc::pfunc()))(tau,0.,0.);
}

void MGBVertex::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
//	const MGPosition& P=point()->position();
//	(*(MGDrawFunc::pfunc()))(P[0],P[1],P[2]);
}

void MGBVertex::draw3DVertex(
)const{
	const MGPoint* pnt=point();
	if(pnt){
		const MGPosition& P=pnt->position();
		(*(MGDrawFunc::pfunc()))(P[0],P[1],P[2]);
		return;
	}
	const MGCellBase* partner=m_partners[0];
	const MGEdge* edg=dynamic_cast<const MGEdge*>(partner->star());
	if(!edg) return;
	const MGPVertex* pv=dynamic_cast<const MGPVertex*>(partner);
	if(!pv) return;
	MGVector P=edg->eval(pv->t());
	(*(MGDrawFunc::pfunc()))(P[0],P[1],P[2]);
}

void MGEdge::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	const MGCurve* crv=base_curve();
	crv->drawSE(span_length, param_s(), param_e());
}

void MGEdge::draw3DVertex(
)const{
	MGDrawFunc::DrawPoint pf=*(MGDrawFunc::pfunc());
	if(active_start()){
		MGPosition Ps=start_point();
		pf(Ps[0],Ps[1],Ps[2]);
	}
	if(active_end()){
		MGPosition Pe=end_point();
		pf(Pe[0],Pe[1],Pe[2]);
	}
}

void MGFace::draw3DVertex(
)const{
	const MGBox& pbx=box_param();
	if(pbx.is_null()) return;

	//Draw boundaries' vertexes.
	size_t m=number_of_boundaries();
	for(size_t i=0; i<m; i++){
		const MGBoundary& loop=*(boundary(i));
		loop.draw3DVertex_in_star();
	}
}

//Shade the object in world coordinates.
void MGFace::shade(
	double span_length	//Line segment span length.
)const{
	mgTLData tess(*this,mgTLInputParam(*this,span_length));
	glShadeModel(GL_SMOOTH);// shading model
	mgGDL::shade(tess);
}

//Shade the object in world coordinates.
void MGShell::shade(
	double span_length	//Line segment span length.
)const{
	mgTLDataVector tess(*this,mgTLInputParam(*this,span_length));
	glShadeModel(GL_SMOOTH);// shading model

	mgTLDataVector::iterator i=tess.begin(), iend=tess.end();
	for(; i!=iend; i++)
		mgGDL::shade(*i);
}

//Draw an object in its parameter space.
//This is valid only for Surface, Face, Loop, Edge.
void MGDraw_in_parameter_space(
	const MGObject& obj,	//The object to draw.
	double span_length		//Line segment span length.
){
	MGDrawFunc::DrawPoint pf=*(MGDrawFunc::pfunc());
	MGDrawFunc::VertexF vf=*(MGDrawFunc::vfunc());
	const MGSurface* sf=dynamic_cast<const MGSurface*>(&obj);
	if(sf){
		const MGPlane* pl=dynamic_cast<const MGPlane*>(sf);
		if(pl) return;//Plane will not be drawn.
		MGBox uv=sf->param_range();
		double u0=uv[0].low_point(), u1=uv[0].high_point();
		double v0=uv[1].low_point(), v1=uv[1].high_point();
		(*(MGDrawFunc::bfunc()))();
		const MGDrawFunc::VertexF vertex=MGDrawFunc::vfunc();
		vf(u0, v0, 0.); vf(u1, v0, 0.);
		vf(u1, v1, 0.); vf(u0, v1, 0.);
		vf(u0, v0, 0.);
		(*(MGDrawFunc::efunc()))();
		pf(u0, v0, 0.);
		pf(u1, v0, 0.);
		pf(u1, v1, 0.);
		pf(u0, v1, 0.);
		return;
	}
	const MGFace* f=dynamic_cast<const MGFace*>(&obj);
	if(f){
		if(f->hasOuterBoundaryLoop()){
			const MGLoop& ol=*(f->loop(size_t(0)));
			ol.drawWire(span_length);
		}else{
			MGPvector<MGCurve> crvs=f->outer_boundary_param();
			size_t n=crvs.size();
			for(size_t i=0; i<n; i++) crvs[i]->drawWire(span_length);
		}
		size_t i,if0;
		size_t nib=f->number_of_inner_boundaries(if0);
		for(i=0; i<nib; i++,if0++) (f->loop(if0))->drawWire(span_length);

		size_t nlp=f->number_of_boundaries();
		for(; if0<nlp; if0++) (f->loop(if0))->drawWire(span_length);

		for(i=0; i<nlp; i++){
			const MGLoop& lp=*(f->loop(i));
			lp.draw3DVertex();
		}
		return;
	}
	const MGLoop* lp=dynamic_cast<const MGLoop*>(&obj);
	if(lp){
		lp->drawWire(span_length);
		lp->draw3DVertex();
		return;
	}
	const MGEdge* edg=dynamic_cast<const MGEdge*>(&obj);
	if(edg){
		edg->drawWire(span_length);
		edg->draw3DVertex();
		return;
	}
}

//Make a display list of only call lists of this shell's faces.
void MGShell::make_only_call_list(
	int mode	//=0 when wire, =1 when wire and no color, =2 when SHADING.
)const{
	bool no_color=false;//if true, color attribute will be neglected.
	if(mode==1)
		no_color=true;

	size_t glname=dlist_name()+mode;
	glNewList(glname, GL_COMPILE);
		size_t mask=get_draw_attrib_mask();
		if(mask){
			glPushAttrib(mask);
			draw_attribute(no_color);
		}
		size_t nfaces=number_of_faces();
		//Make display list of the faces.
		for(size_t i=0; i<nfaces; i++){
			const MGFace& facei=*(face(i));
			if(facei.no_display())
				continue;
			size_t name2=facei.dlist_name();
			name2+=mode;
			glCallList(name2);//call object.
		}
		if(mask)
			glPopAttrib();
	glEndList();
}

//Delete a display list of this gel.
void MGShell::delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be apppended.
)const{
	size_t nfaces=number_of_faces();
	for(size_t i=0; i<nfaces; i++){
		const MGFace& facei=*(face(i));
		if(facei.no_display())
			continue;
		facei.delete_display_list(glv,functions);
	}
	size_t nm=dlist_name();
	glDeleteLists(nm,3);
	glv.m_sysgllist.delete_lists_by_object_id(this,functions);
}

//Make a display list of this gel.
size_t MGShell::make_display_list(
	double span_length,//span length to approximate by polyline.
	int line_density//line density to draw surface in wire mode.
)const{
	size_t name=dlist_name();

	//make the wire mode display list to call list.
	make_only_call_list(0);

	//make the shading mode display list to call list.
	make_only_call_list(2);

	size_t nfaces=number_of_faces();
	//Make display list of the gel that has an object.
	for(size_t i=0; i<nfaces; i++){
		const MGFace& facei=*(face(i));
		if(facei.no_display())
			continue;
		facei.make_display_list(span_length,line_density);
	}
	return name;
}

///Make a display list without color of this gel.
///Return is the display list name.
size_t MGShell::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{
	size_t name=dlist_name()+1;

	//make the wire mode display list to call list.
	make_only_call_list(1);
	size_t nfaces=number_of_faces();
	//Make display list of the gel that has an object.
	for(size_t i=0; i<nfaces; i++){
		const MGFace& facei=*(face(i));
		if(facei.no_display())
			continue;
		facei.make_display_list_to_hilight(span_length,line_density);
	}
	return name;
}