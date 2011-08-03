/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgOpenGLRenderView_HH_
#define _mgOpenGLRenderView_HH_

#include "mgGL/OpenGLView.h"
#include "Tl/TLInputParam.h"

///@cond

/// COpenGLWindow ビュー
class MGCLASS MGOpenGLRenderView:public MGOpenGLView{

public:

//////////////////////////////////METHOD/////////////////////

MGOpenGLRenderView(
	void* object=0,	///<Object pointer who uses this MGOpenGLRenderView.
		///<This pointer will be passed to dvvf as the 1st parameter as (*dvvf)(object).
	const MGColor* Bcolor=0,///<Background color.
	HDC hdc=0,		///<Device Context
    HGLRC hRC=0,	///<Rendereing context
	Compute_eye_func eye_func=0,///<eye position copmuting function.
					///<this function is used wnly when the following function is invoked.
					///<"void initialize_viewing_environment(const MGBox& box);"
	double fovy=45.///<viewing pyramid angle in degrees.
);
///object, dvvf, and mmf must be set valid data before use of MGOpenGLRenderView.

///Construct from MGglViewAttrib.
MGOpenGLRenderView(
	const mgTLInputParam& m_tessellate_param,
	const MGglViewAttrib& glatr,
	void* object=0,	///<Object pointer who uses this MGOpenGLRenderView.
					///<This pointer will be passed to dvvf as the 1st
					///<parameter as (*dvvf)(object, cx, cy).
	double smooth=.01,
	double pick_aperture=5.,
	const MGColor* Bcolor=0,///<Background color.
	HDC hdc=0,	///<Device Context
    HGLRC hRC=0,	///<Rencereing context
	Compute_eye_func eye_func=0	///<eye position copmuting function.
);

///~MGOpenGLRenderView();

/////////////////////////////////////////////////////////////////////////////////

///Generate a display list of the object obj.
///GDLdraw will NOT invoke make_RC_current();
void GDLdraw(
	const MGObject& obj,///<the object to draw.
	bool no_color=false	///<if true, color attribute will be neglected.
);

///When inheritted class from this MGOpenGLRenderView needs to process GL attributes
///just after glNewList(nm, GL_COMPILE), define GDLexec_attributes() to
///write the necessary process. The default is line drawing attributes.
void GDLexec_attributes(
	const MGGroup& group,///<the group to replace from.
	bool no_color=false		///<if true, color attribute will be neglected.
);

///When inheritted class from this MGOpenGLRenderView needs to process GL attributes
///just after glNewList(nm, GL_COMPILE), define GDLexec_attributes() to
///write the necessary process. The default is line drawing attributes.
void GDLexec_attributes(
	const MGObject& obj,	///<the object to process.
	bool no_color=false		///<if true, color attribute will be neglected.
);

void initialize();

void set_tessellate_param(const mgTLInputParam& tlparam){m_tessellate_param=tlparam;};
const mgTLInputParam& tessellate_param()const{return m_tessellate_param;};
mgTLInputParam& tessellate_param(){return m_tessellate_param;};

private:
/// アトリビュート

	mgTLInputParam m_tessellate_param;///<tessellation parameter.

};

///@endcond

#endif
