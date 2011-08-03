/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Fog.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGFog::MGFog(){;};
MGFog::MGFog(FOGMODE mode):MGGLAttrib(mode){;}
MGFog::MGFog(FOGMODE mode,float density,float start,float end, const MGColor& color)
:MGGLAttrib(mode), m_density(density),m_start(start),m_end(end){
	setColor(color);
}
	
MGFog* MGFog::clone()const{
	return new MGFog(*this);
}

//assignment
MGFog& MGFog::operator=(const MGFog& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_density=gel2.m_density;
	m_start=gel2.m_start;
	m_end=gel2.m_end;
	m_color=gel2.m_color;
	return *this;
}
MGFog& MGFog::operator=(const MGGel& gel2){
	const MGFog* gel2_is_this=dynamic_cast<const MGFog*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGFog::operator<(const MGFog& gel2)const{
	if(m_flag==gel2.m_flag)
		return m_color<gel2.m_color;
	return m_flag<gel2.m_flag;
}
bool MGFog::operator<(const MGGel& gel2)const{
	const MGFog* gel2_is_this=dynamic_cast<const MGFog*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//render GLAttribute process.
void MGFog::exec()const{
	if(undefined()) return;
	if(disabled()){
		glDisable(GL_FOG);
		return;
	}
	glEnable(GL_FOG);
	glFogi(GL_FOG_MODE,GLfog_mode());
	glFogf(GL_FOG_START,m_start);
	glFogf(GL_FOG_END,m_end);
	glFogf(GL_FOG_DENSITY,m_density);
	const float* color=m_color.color();
	glFogfv(GL_COLOR,color);
}

// Output function.
std::ostream& MGFog::out(std::ostream& ostrm) const{
	ostrm<<std::endl<<"Fog="; MGGLAttrib::out(ostrm);
	if(enabled()){
		FOGMODE m=fog_mode();
		ostrm<<",fog mode=";
		switch(m){
			case LINEAR: ostrm<<"LINEAR"; break;
			case EXP: ostrm<<"EXP"; break;
			case EXP2: ostrm<<"EXP2"; break;
		}
		ostrm<<",m_density"<<m_density<<",m_start="<<m_start;
		ostrm<<",m_end="<<m_end<<std::endl;
		const float* color=m_color.color();
		ostrm<<"m_color=["<<color[0]<<","<<color[1]<<","<<color[2]<<","<<color[3]<<"]";
	}
	return ostrm;
}

//Write all member data
void MGFog::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	buf<<m_density;
	buf<<m_start;
	buf<<m_end;
	const float* color=m_color.color();
	for(size_t i=0; i<4; i++) buf<<color[i];
}
//Read all member data.
void MGFog::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	buf>>m_density;
	buf>>m_start;
	buf>>m_end;
	float* color=m_color.color();
	for(size_t i=0; i<4; i++) buf>>color[i];
}
