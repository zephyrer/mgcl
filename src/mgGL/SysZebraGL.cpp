/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// SysZebraGL.h: mgSysZebraGL クラスのインターフェイス
#include "MGCLStdAfx.h"
#include "mgGL/SysZebraGL.h"
#include "mgGL/GLDrawFunc.h"
#include "Tl/TLDataVector.h"

mgSysZebraGL::mgSysZebraGL(
	size_t function_code,//Function code.
	const MGGel* gel,	//This must be MGSurface, MGFace, or MGShell(mgAll_2Manifold).
	const mgTLInputParam& tlparam,//Tessellation parameter.
	unsigned char color[4],//zebra color
	float thickness,	//zebra thickness.
	bool is_vertical	//direction of the zera stripes.
):mgSysGL(function_code,gel),m_tlparam(tlparam),m_thickness(thickness),m_vertical(is_vertical){
	for(int i = 0; i < 4; ++i){
		m_color[i]= color[i];
	}
}

//Construct new object by copying to newed area.
//User must delete this copied object by "delete".
mgSysZebraGL* mgSysZebraGL::clone()const{
	return new mgSysZebraGL(*this);
}

//Draw this Sysgl.
//This draw is used to draw the pictures for Undo(, Redo) operations.
//When this draw is invoked, all of the functions of mgGDL can be
//used in that context.
void mgSysZebraGL::draw(MGOpenGLView& glv)const{
	glPushAttrib(GL_TEXTURE_BIT | GL_CURRENT_BIT | GL_POLYGON_BIT
				| GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT
				| GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);//Push Attrib//////

	glDisable(GL_BLEND);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_CULL_FACE);
	glDisable(GL_TEXTURE_2D);

	size_t tex_name=dlist_name();
	glBindTexture(GL_TEXTURE_1D, tex_name);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	const int border = TEXGEN_WIDTH/2;
	zebraData ubvImage[TEXGEN_WIDTH];
	for(int i = 0; i < TEXGEN_WIDTH; ++i){
		ubvImage[i].BYTE[0]= (GLubyte)((i <= border) ? 255 : m_color[0]);
		ubvImage[i].BYTE[1]= (GLubyte)((i <= border) ? 255 : m_color[1]);
		ubvImage[i].BYTE[2]= (GLubyte)((i <= border) ? 255 : m_color[2]);
		ubvImage[i].BYTE[3]= m_color[3];
	}

	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, TEXGEN_WIDTH, 0,
					GL_RGBA, GL_UNSIGNED_BYTE, ubvImage);
	glEnable(GL_TEXTURE_1D);
	glEnable(GL_TEXTURE_GEN_S);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glClearDepth(1.0);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	const float sc_mtrShininess = 128.0f;
	glMaterialf(GL_FRONT, GL_SHININESS, sc_mtrShininess);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	mgTLDataVector       tldata;   // tessellation 情報
	tldata.push_back(*(object_id()->object()),m_tlparam);
	std::for_each(tldata.begin(), tldata.end(), mgGDL::shade);

	glPopAttrib();//Pop Attrib////////
}

//Process necessary before transformation in MGOpenGLView::DrawScene.
void mgSysZebraGL::pre_transform_process()const{
	float fvPlane[4];
	if(m_vertical){
		fvPlane[0] = m_thickness;
		fvPlane[1] = 0.0f;
		fvPlane[2] = m_thickness;
		fvPlane[3] = 0.0f;
	}else{
		fvPlane[0] = 0.0f;
		fvPlane[1] = m_thickness;
		fvPlane[2] = m_thickness;
		fvPlane[3] = 0.0f;
	}
	glDisable(GL_COLOR_MATERIAL);
	glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR);
	glTexGenfv(GL_S, GL_EYE_PLANE, fvPlane);
}

// Output virtual function.
//Output to stream file:メンバデータを標準出力に出力する。
std::ostream& mgSysZebraGL::out(std::ostream& ostrm) const{
	ostrm<<"mgSysZebraGL::"<<this;
	mgSysGL::out(ostrm);
	ostrm<<",color=("<<m_color[0]<<","<<m_color[1]<<","<<m_color[2]<<","<<m_color[3]<<")";
	ostrm<<",thickness="<<m_thickness<<",vertical="<<m_vertical;
	return ostrm;
}
