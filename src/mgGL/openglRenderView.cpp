/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Group.h"
#include "mgGL/OpenGLRenderView.h"
#include "mgGL/GLDrawFunc.h"
#include "Tl/TLDataVector.h"

/////////////////////////////////////////////////////////////////////////////
// COpenGLRender ÉrÉÖÅ[

//////////////////////////////METHOD//////////////////

MGOpenGLRenderView::MGOpenGLRenderView(
	void* object,	//Object pointer who uses this MGOpenGLRenderView.
	const MGColor* Bclr,//Background color.
	HDC hdc,		//Device Context
    HGLRC hRC,	//Rendereing context
	Compute_eye_func eye_func,//eye position copmuting function.
					//this function is used wnly when the following function is invoked.
					//"void initialize_viewing_environment(const MGBox& box);"
	double fovy//viewing pyramid angle in degrees.
){
	set_object(object);
	if(Bclr)
		setBcolor(*Bclr);
	setDCRC(hdc,hRC);
	set_eye_func(eye_func);
	set_fovy(fovy);
	// initialize();
}

//object, dvvf, and mmf must be set valid data before use of MGOpenGLRenderView.

//Construct from MGglViewAttrib.
MGOpenGLRenderView::MGOpenGLRenderView(
	const mgTLInputParam& tessellate_param,
	const MGglViewAttrib& glatr,
	void* object,	//Object pointer who uses this MGOpenGLRenderView.
					//This pointer will be passed to dvvf as the 1st
					//parameter as (*dvvf)(object, cx, cy).
	double smooth,
	double pick_aperture,
	const MGColor* Bclr,//Background color.
	HDC hdc,	//Device Context
    HGLRC hRC,	//Rencereing context
	Compute_eye_func eye_func	//eye position copmuting function.
):MGOpenGLView(glatr){
	set_object(object);
	set_smooth(smooth);
	set_pick_aperture(pick_aperture);
	if(Bclr)
		setBcolor(*Bclr);
	setDCRC(hdc,hRC);
	set_eye_func(eye_func);
	make_RC_current();
	glClearDepth(0.);
	const float* Bcolr=m_Bcolor.color();
    glClearColor(Bcolr[0], Bcolr[1], Bcolr[2], Bcolr[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
//	glEnable(GL_LIGHT0);
//	glEnable(GL_LIGHTING);
}

///////////////////////////////////////////////////////////////////////

//Generate a display list of the object obj.
//GDLdraw will NOT invoke make_RC_current();
void MGOpenGLRenderView::GDLdraw(
	const MGObject& obj,//the object to draw.
	bool no_color	//if true, color attribute will be neglected.
){
	size_t nm=size_t(&obj);
	glNewList(nm, GL_COMPILE);
		glPushName(nm);
		if(obj.manifold_dimension()==2){
			mgTLDataVector tlds(obj,m_tessellate_param);
			std::for_each(tlds.begin(), tlds.end(), mgGDL::shade);
		}else
			obj.drawWire(span_length(),line_density());
		glPopName();
	glEndList();
}

//When inheritted class from this MGOpenGLRenderView needs to process GL attributes
//just after glNewList(nm, GL_COMPILE), define GDLexec_attributes() to
//write the necessary process. The default is line drawing attributes.
void MGOpenGLRenderView::GDLexec_attributes(
	const MGGroup& group,//the group to replace from.
	bool no_color		//if true, color attribute will be neglected.
){
	group.render_attribute();
}

//When inheritted class from this MGOpenGLRenderView needs to process GL attributes
//just after glNewList(nm, GL_COMPILE), define GDLexec_attributes() to
//write the necessary process. The default is line drawing attributes.
void MGOpenGLRenderView::GDLexec_attributes(
	const MGObject& obj,	//the object to process.
	bool no_color		//if true, color attribute will be neglected.
){
}

void MGOpenGLRenderView::initialize(){
	make_RC_current();
	glClearDepth(1.);
	const float* Bcolr=m_Bcolor.color();
    glClearColor(Bcolr[0], Bcolr[1], Bcolr[2], Bcolr[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
//	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
}