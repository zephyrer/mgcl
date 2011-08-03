/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// mgSysGLList.cpp : mgSysGLList クラスのimplementation。
#include "MGCLStdAfx.h"
#include "mgGL/OpenGLView.h"
#include "mgGL/SysGL.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Copy constructor, replacing gel_old to gel_new.
mgSysGL::mgSysGL(
	mgSysGL& glold,
	const MGGel* gel_old,
	const MGGel* gel_new
):m_fucntion_id(glold.m_fucntion_id),m_gel(glold.m_gel){
	replace(gel_old,gel_new);	
}

//Construct new object by copying to newed area.
//User must delete this copied object by "delete".
mgSysGL* mgSysGL::clone()const{
	return new mgSysGL(*this);
}

bool mgSysGL::includes(const MGGel* gel)const{
	if(!gel || !m_gel)
		return false;
	return gel==m_gel;
}

//Make system display list in glv.
//This must be a newed object and the ownership will be transfered to
//glv(glv.m_sysgllist).
void mgSysGL::make_display_list(
	MGOpenGLView& glv
){
	size_t temporal_draw_id = glv.push_back_to_sysgl(this);
	// OpenGLコンテキストに依存する描画を行う
	::glNewList(temporal_draw_id, GL_COMPILE);
	draw(glv);
	::glEndList();
}

//replace gel_old to gel_new.
//If gel_old is not included in this, do nothing.
void mgSysGL::replace(
	const MGGel* gel_old, const MGGel* gel_new
){
	if(m_gel==gel_old)
		m_gel=const_cast<MGGel*>(gel_new);
}

// Output virtual function.
//Output to stream file:メンバデータを標準出力に出力する。
std::ostream& mgSysGL::out(std::ostream& ostrm) const{
//	ostrm.setf(ios::scientific,ios::floatfield);
//	ostrm.precision(10);
	ostrm<<"mgSysGL::"<<this;
	ostrm<<",m_fucntion_id="<<m_fucntion_id
		<<",m_gel="<<m_gel;
	return ostrm;
}

//////////// mgSysGL output ////////////
std::ostream& operator<<(std::ostream& outp, const mgSysGL& sysgl){
	sysgl.out(outp);
	return outp;
}