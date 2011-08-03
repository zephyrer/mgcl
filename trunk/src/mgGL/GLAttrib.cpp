/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <bitset>
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Color.h"
#include "mgGL/Lights.h"
#include "mgGL/TranspMode.h"
#include "mgGL/AlphaFunc.h"
#include "mgGL/BlendFunc.h"
#include "mgGL/ColorMask.h"
#include "mgGL/LineWidth.h"
#include "mgGL/DepthFunc.h"
#include "mgGL/DepthMask.h"
#include "mgGL/LightEnable.h"
#include "mgGL/LineStipple.h"
#include "mgGL/PolygonMode.h"
#include "mgGL/ShadeModel.h"
#include "mgGL/Texture.h"
#include "mgGL/Fog.h"
#include "mgGL/DirectionalLight.h"
#include "mgGL/SpotLight.h"
#include "mgGL/RenderAttr.h"
#include "mg/GelFactory.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Set or reset the bit of mask.
void set_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit){mask|=bit;}
void reset_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit){mask&=~bit;}

//////////////  MGGLAttrib  ////////////////

//Assignment
MGGLAttrib& MGGLAttrib::set_glattrib(const MGGLAttrib& gel2){
	if(this==&gel2)
		return *this;

	MGAttrib::operator=(gel2);
	m_flag=gel2.m_flag;
	return *this;
}

//Read all member data.
//Write all member data
void MGGLAttrib::WriteMembers(MGOfstream& buf)const{
	buf << m_flag;
}
void MGGLAttrib::ReadMembers(MGIfstream& buf){
	buf >> m_flag;
}

std::ostream& MGGLAttrib::out(
	std::ostream& ostrm
)const{
	ostrm<<this;
	if(undefined()) ostrm<<":UNDEFINED";
	else if(disabled()){
		ostrm<<":DISABLED";
	}
	return ostrm;
}

//Compare if this and at2 are the same leaf MGGLAttrib class.
bool MGGLAttrib::same_type(const MGGLAttrib& at2)const{
	return identify_type()==at2.identify_type();
}

//////////////  MGAlphaFunc  ////////////////
bool MGAlphaFunc::operator<(const MGAlphaFunc& gel2)const{
	return m_refval<gel2.m_refval;
}
bool MGAlphaFunc::operator<(const MGGel& gel2)const{
	const MGAlphaFunc* gel2_is_this=dynamic_cast<const MGAlphaFunc*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGAlphaFunc* MGAlphaFunc::clone()const{
	return new MGAlphaFunc(*this);
}

//assignment
MGAlphaFunc& MGAlphaFunc::operator=(const MGAlphaFunc& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_refval=gel2.m_refval;
	return *this;
}
MGAlphaFunc& MGAlphaFunc::operator=(const MGGel& gel2){
	const MGAlphaFunc* gel2_is_this=dynamic_cast<const MGAlphaFunc*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

std::ostream& MGAlphaFunc::out(std::ostream& ostrm) const{
	ostrm<<"AlphaFunc="; MGGLAttrib::out(ostrm);
	if(enabled()){
		FUNC m=get_func();
		ostrm<<",func=";
		switch(m){
			case NEVER: ostrm<<"NEVER"; break;
			case LESS: ostrm<<"LESS"; break;
			case EQUAL: ostrm<<"EQUAL"; break;
			case LEQUAL: ostrm<<"LEQUAL"; break;
			case GREATER: ostrm<<"GREATER"; break;
			case NOTEQUAL: ostrm<<"NOTEQUAL"; break;
			case GEQUAL: ostrm<<"GEQUAL"; break;
			case ALWAYS: ostrm<<"ALWAYS"; break;
		}
		ostrm<<",m_refval="<<m_refval;
	}
	return ostrm;
}

void MGAlphaFunc::exec()const{
	if(undefined()) return;
	if(disabled()){
		glDisable(GL_ALPHA_TEST);
	}else{
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GLfunc(),m_refval);
	}
}

// Serialization.
void MGAlphaFunc::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	buf<<m_refval;
}
void MGAlphaFunc::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	buf>>m_refval;
}

//////////////  MGBlendFunc  ////////////////

MGBlendFunc* MGBlendFunc::clone()const{
	return new MGBlendFunc(*this);
}

//assignment
MGBlendFunc& MGBlendFunc::operator=(const MGBlendFunc& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_dst_blend_func=gel2.m_dst_blend_func;
	return *this;
}
MGBlendFunc& MGBlendFunc::operator=(const MGGel& gel2){
	const MGBlendFunc* gel2_is_this=dynamic_cast<const MGBlendFunc*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGBlendFunc::operator<(const MGBlendFunc& gel2)const{
	if(m_flag==gel2.m_flag)
		return m_dst_blend_func<gel2.m_dst_blend_func;
	return m_flag<gel2.m_flag;
}
bool MGBlendFunc::operator<(const MGGel& gel2)const{
	const MGBlendFunc* gel2_is_this=dynamic_cast<const MGBlendFunc*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

// Output function.
std::ostream& MGBlendFunc::out(std::ostream& ostrm) const{
	ostrm<<"BlendFunc="; MGGLAttrib::out(ostrm);
	if(enabled()){
		SrcBlendFunc sf=get_sfunc(); DstBlendFunc df=get_dfunc();
		ostrm<<",source func=";
		switch(sf){
			case ZERO_SBLEND: ostrm<<"ZERO_SBLEND"; break;
			case ONE_SBLEND: ostrm<<"ONE_SBLEND"; break;
			case DST_COLOR_SBLEND: ostrm<<"DST_COLOR_SBLEND"; break;
			case ONE_MINUS_DST_COLOR_SBLEND: ostrm<<"ONE_MINUS_DST_COLOR_SBLEND"; break;
			case SRC_ALPHA_SATURATE_SBLEND: ostrm<<"SRC_ALPHA_SATURATE_SBLEND"; break;
			case SRC_ALPHA_SBLEND: ostrm<<"SRC_ALPHA_SBLEND"; break;
			case ONE_MINUS_SRC_ALPHA_SBLEND: ostrm<<"ONE_MINUS_SRC_ALPHA_SBLEND"; break;
			case DST_ALPHA_SBLEND: ostrm<<"DST_ALPHA_SBLEND"; break;
			case ONE_MINUS_DST_ALPHA_SBLEND: ostrm<<"ONE_MINUS_DST_ALPHA_SBLEND"; break;
		}
		ostrm<<",destination func=";
		switch(df){
			case ZERO_DBLEND: ostrm<<"ZERO_DBLEND"; break;
			case ONE_DBLEND: ostrm<<"ONE_DBLEND"; break;
			case SRC_COLOR_DBLEND: ostrm<<"SRC_COLOR_DBLEND"; break;
			case ONE_MINUS_SRC_COLOR_DBLEND: ostrm<<"ONE_MINUS_SRC_COLOR_DBLEND"; break;
			case SRC_ALPHA_DBLEND: ostrm<<"SRC_ALPHA_DBLEND"; break;
			case ONE_MINUS_SRC_ALPHA_SBLEND: ostrm<<"ONE_MINUS_SRC_ALPHA_SBLEND"; break;
			case DST_ALPHA_DBLEND: ostrm<<"DST_ALPHA_DBLEND"; break;
			case ONE_MINUS_DST_ALPHA_DBLEND: ostrm<<"ONE_MINUS_DST_ALPHA_DBLEND"; break;
		}
	}
	return ostrm;
}

void MGBlendFunc::exec()const{
	if(undefined()) return;
	if(disabled()){
		glDisable(GL_BLEND);
	}else{
		glEnable(GL_BLEND);
		glBlendFunc(GLsfunc(),GLdfunc());
	}
}

void MGBlendFunc::set_func(
	SrcBlendFunc sf,
	DstBlendFunc df
){
	m_flag=sf;
	m_dst_blend_func=df;
}

MGBlendFunc::MGBlendFunc(SrcBlendFunc sf, DstBlendFunc df)
:MGGLAttrib(sf),m_dst_blend_func(df){
}

void MGBlendFunc::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	buf<<static_cast<int>(m_dst_blend_func);
}
void MGBlendFunc::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	int data;
	buf>>data;
	m_dst_blend_func=static_cast<MGBlendFunc::DstBlendFunc>(data);
}

//////////////  MGColorMask  ////////////////
	
MGColorMask* MGColorMask::clone()const{
	return new MGColorMask(*this);
}

MGColorMask::MGColorMask(bool red, bool green, bool blue, bool alpha){
	set_mask(red, green, blue, alpha);
}

//assignment
MGColorMask& MGColorMask::operator=(const MGColorMask& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	return *this;
}
MGColorMask& MGColorMask::operator=(const MGGel& gel2){
	const MGColorMask* gel2_is_this=dynamic_cast<const MGColorMask*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGColorMask::operator<(const MGColorMask& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGColorMask::operator<(const MGGel& gel2)const{
	const MGColorMask* gel2_is_this=dynamic_cast<const MGColorMask*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

// Output function.
std::ostream& MGColorMask::out(std::ostream& ostrm) const{
	ostrm<<"ColorMask="; MGGLAttrib::out(ostrm);
	if(defined()){
		bool r,g,b,a;
		get_mask(r,g,b,a);
		ostrm<<",mask=["<<r<<","<<g<<","<<b<<","<<a;
		ostrm<<"]";
	}
	return ostrm;
}

void MGColorMask::exec()const{
	if(undefined()) return;
	bool r,g,b,a;
	get_mask(r,g,b,a);
	glColorMask(r,g,b,a);
}

void MGColorMask::get_mask(bool& red, bool& green, bool& blue, bool& alpha)const{
	std::bitset<32> mask(m_flag);
	red=mask[0];
	green=mask[1];
	blue=mask[2];
	alpha=mask[3];
}

void MGColorMask::set_mask(bool red, bool green, bool blue, bool alpha){
	std::bitset<32> mask;
	mask[0]=red;
	mask[1]=green;
	mask[2]=blue;
	mask[3]=alpha;
	m_flag=mask.to_ulong();
}

void MGColorMask::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
}
void MGColorMask::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
}

//////////////  MGDepthFunc  ////////////////
	
//assignment
MGDepthFunc& MGDepthFunc::operator=(const MGDepthFunc& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	return *this;
}
MGDepthFunc& MGDepthFunc::operator=(const MGGel& gel2){
	const MGDepthFunc* gel2_is_this=dynamic_cast<const MGDepthFunc*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGDepthFunc::operator<(const MGDepthFunc& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGDepthFunc::operator<(const MGGel& gel2)const{
	const MGDepthFunc* gel2_is_this=dynamic_cast<const MGDepthFunc*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGDepthFunc* MGDepthFunc::clone()const{
	return new MGDepthFunc(*this);
}

// Output function.
std::ostream& MGDepthFunc::out(std::ostream& ostrm) const{
	ostrm<<"DepthFunc="; MGGLAttrib::out(ostrm);
	if(enabled()){
		FUNC f=get_func();
		ostrm<<",func=";
		switch(f){
			case NEVER: ostrm<<"NEVER"; break;
			case LESS: ostrm<<"LESS"; break;
			case EQUAL: ostrm<<"EQUAL"; break;
			case LEQUAL: ostrm<<"LEQUAL"; break;
			case GREATER: ostrm<<"GREATER"; break;
			case NOTEQUAL: ostrm<<"NOTEQUAL"; break;
			case GEQUAL: ostrm<<"GEQUAL"; break;
			case ALWAYS: ostrm<<"ALWAYS"; break;
		}
	}
	return ostrm;
}

void MGDepthFunc::exec()const{
	if(undefined()) return;
	if(disabled()){
		glDisable(GL_DEPTH_TEST);
	}else{
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GLfunc());
	}
}

void MGDepthFunc::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
}
void MGDepthFunc::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
}

//////////////  MGDepthMask  ////////////////
	
MGDepthMask* MGDepthMask::clone()const{
	return new MGDepthMask(*this);
}

//assignment
MGDepthMask& MGDepthMask::operator=(const MGDepthMask& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	return *this;
}
MGDepthMask& MGDepthMask::operator=(const MGGel& gel2){
	const MGDepthMask* gel2_is_this=dynamic_cast<const MGDepthMask*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGDepthMask::operator<(const MGDepthMask& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGDepthMask::operator<(const MGGel& gel2)const{
	const MGDepthMask* gel2_is_this=dynamic_cast<const MGDepthMask*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

// Output function.
std::ostream& MGDepthMask::out(std::ostream& ostrm) const{
	ostrm<<"DepthMask="; MGGLAttrib::out(ostrm);
	if(enabled()) ostrm<<"ENABLED";
	return ostrm;
}

void MGDepthMask::exec()const{
	if(undefined()) return;
	glDepthMask(!disabled());
}

void MGDepthMask::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
}
void MGDepthMask::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
}

//////////////  MGLineWidth  ////////////////

//assignment
MGLineWidth& MGLineWidth::operator=(const MGLineWidth& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_line_width=gel2.m_line_width;
	return *this;
}
MGLineWidth& MGLineWidth::operator=(const MGGel& gel2){
	const MGLineWidth* gel2_is_this=dynamic_cast<const MGLineWidth*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLineWidth::operator<(const MGLineWidth& gel2)const{
	if(m_flag==gel2.m_flag)
		return m_line_width<gel2.m_line_width;
	return m_flag<gel2.m_flag;
}
bool MGLineWidth::operator<(const MGGel& gel2)const{
	const MGLineWidth* gel2_is_this=dynamic_cast<const MGLineWidth*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGLineWidth* MGLineWidth::clone()const{
	return new MGLineWidth(*this);
}

//get maximum line width
float MGLineWidth::get_maximum_width()const{
	float wd[2];
	glGetFloatv(GL_LINE_WIDTH_RANGE,wd);
	return wd[1];
}

// Output function.
std::ostream& MGLineWidth::out(std::ostream& ostrm) const{
	ostrm<<"LineWidth="; MGGLAttrib::out(ostrm);
	if(undefined()) return ostrm;
	ostrm<<"="<<m_line_width;
	return ostrm;
}

void MGLineWidth::set_width(float width){
	m_line_width=width;
	m_flag=ENABLED;
}

void MGLineWidth::exec()const{
	if(undefined()) return;
	glLineWidth(get_width());
}

void MGLineWidth::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	buf<<m_line_width;
}
void MGLineWidth::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	buf>>m_line_width;
}

//////////////  MGPolygonMode  ////////////////

//assignment
MGPolygonMode& MGPolygonMode::operator=(const MGPolygonMode& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	return *this;
}
MGPolygonMode& MGPolygonMode::operator=(const MGGel& gel2){
	const MGPolygonMode* gel2_is_this=dynamic_cast<const MGPolygonMode*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGPolygonMode::operator<(const MGPolygonMode& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGPolygonMode::operator<(const MGGel& gel2)const{
	const MGPolygonMode* gel2_is_this=dynamic_cast<const MGPolygonMode*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGPolygonMode* MGPolygonMode::clone()const{
	return new MGPolygonMode(*this);
}

// Output function.
std::ostream& MGPolygonMode::out(std::ostream& ostrm) const{
	ostrm<<"PolygonMode="; MGGLAttrib::out(ostrm);
	if(undefined()) return ostrm;

	MODE m=get_mode();
	switch(m){
		case POINT: ostrm<<"POINT"; break;
		case LINE: ostrm<<"LINE"; break;
		case FILL: ostrm<<"FILL"; break;
	}
	return ostrm;
}

void MGPolygonMode::exec()const{
	if(undefined()) return;
	glPolygonMode(GL_FRONT_AND_BACK,get_mode());
}

void MGPolygonMode::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
}
void MGPolygonMode::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
}

//////////////  MGShadeModel  ////////////////

// Output function.
std::ostream& MGShadeModel::out(std::ostream& ostrm) const{
	ostrm<<"ShadeModel=";
	MODEL m=get_model();
	switch(m){
		case SMOOTH: ostrm<<"SMOOTH"; break;
		case FLAT: ostrm<<"FLAT"; break;
	}
	return ostrm;
}

void MGShadeModel::exec()const{
	glShadeModel(GLmodel());
}

// Serialization fucntion.
MGOfstream& operator<< (MGOfstream& buf, const MGShadeModel& sm){
	buf<<static_cast<size_t>(sm.m_model); return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGShadeModel& sm){
	size_t sm2;
	buf>>sm2;
	sm.m_model=static_cast<MGShadeModel::MODEL>(sm2); return buf;
}

//////////////  MGTranspMode  ////////////////
	
MGTranspMode* MGTranspMode::clone()const{
	return new MGTranspMode(*this);
}

//assignment
MGTranspMode& MGTranspMode::operator=(const MGTranspMode& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	return *this;
}
MGTranspMode& MGTranspMode::operator=(const MGGel& gel2){
	const MGTranspMode* gel2_is_this=dynamic_cast<const MGTranspMode*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGTranspMode::operator<(const MGTranspMode& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGTranspMode::operator<(const MGGel& gel2)const{
	const MGTranspMode* gel2_is_this=dynamic_cast<const MGTranspMode*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

// Output function.
std::ostream& MGTranspMode::out(std::ostream& ostrm) const{
	ostrm<<"TranspMode="; MGGLAttrib::out(ostrm);
	if(undefined()) return ostrm;

	MODE m=get_transparency_mode();
	switch(m){
		case FAST: ostrm<<"FAST"; break;
		case NICE: ostrm<<"NICE"; break;
		case BLEND: ostrm<<"BLEND"; break;
		case SCREEN_DOOR: ostrm<<"SCREEN_DOOR"; break;
	}
	return ostrm;
}

MGTranspMode::MGTranspMode(MODE m):MGGLAttrib(m){;};

void MGTranspMode::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
}
void MGTranspMode::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
}

//////////////////MGLightEnable//////////////////
	
MGLightEnable* MGLightEnable::clone()const{
	return new MGLightEnable(*this);
}

//assignment
MGLightEnable& MGLightEnable::operator=(const MGLightEnable& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	return *this;
}
MGLightEnable& MGLightEnable::operator=(const MGGel& gel2){
	const MGLightEnable* gel2_is_this=dynamic_cast<const MGLightEnable*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLightEnable::operator<(const MGLightEnable& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGLightEnable::operator<(const MGGel& gel2)const{
	const MGLightEnable* gel2_is_this=dynamic_cast<const MGLightEnable*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

void MGLightEnable::exec()const{
	if(undefined()) return;
	if(disabled()){
		glDisable(GL_LIGHTING);
	}else{
		glEnable(GL_LIGHTING);
	}
}

// Output function.
std::ostream& MGLightEnable::out(std::ostream& ostrm) const{
	ostrm<<std::endl<<"LightEnable="; MGGLAttrib::out(ostrm);
	if(enabled()) ostrm<<"ENABLED";
	return ostrm;
}

// Serialization fucntion.
void MGLightEnable::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
}
void MGLightEnable::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
}

//////////////////MGLineStipple//////////////////

//assignment
MGLineStipple& MGLineStipple::operator=(const MGLineStipple& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_pattern=gel2.m_pattern;
	return *this;
}
MGLineStipple& MGLineStipple::operator=(const MGGel& gel2){
	const MGLineStipple* gel2_is_this=dynamic_cast<const MGLineStipple*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLineStipple::operator<(const MGLineStipple& gel2)const{
	if(m_flag==gel2.m_flag)
		return m_pattern<gel2.m_pattern;
	return m_flag<gel2.m_flag;
}
bool MGLineStipple::operator<(const MGGel& gel2)const{
	const MGLineStipple* gel2_is_this=dynamic_cast<const MGLineStipple*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGLineStipple::MGLineStipple(LineFont font):MGGLAttrib(2){
	switch(font){
//	case Solid: m_pattern=0xffff; break;
	case Dashed: m_pattern=0x3333; break;
	case Phantom: m_pattern=0x5757; break;
	case CenterLine: m_pattern=0x5f5f; break;
	case Dotted: m_pattern=0x1111; break;
	}
}

MGLineStipple* MGLineStipple::clone()const{
	return new MGLineStipple(*this);
}

void MGLineStipple::exec()const{
	if(undefined()) return;
	if(disabled()){
		glDisable(GL_LINE_STIPPLE);
	}else{
		glLineStipple(get_factor(),get_pattern());
		glEnable(GL_LINE_STIPPLE);
	}
}

//Get the font number
MGLineStipple::LineFont MGLineStipple::get_font_number()const{
	int factr=get_factor();
	if(factr!=2)
		return UndefinedFont;

	unsigned short font=get_pattern();
	switch(font){
	case 0x3333: return Dashed;
	case 0x5757: return Phantom;
	case 0x5f5f: return CenterLine;
	case 0x1111: return Dotted;
	}

	return UndefinedFont;
}

std::ostream& MGLineStipple::out(std::ostream& ostrm) const{
	ostrm<<"LineStipple="; MGGLAttrib::out(ostrm);
	if(undefined()) return ostrm;
	ostrm<<",pattern="<<m_pattern;
	return ostrm;
}

void MGLineStipple::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	buf<<m_pattern;
}
void MGLineStipple::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	buf>>m_pattern;
}

//Construct a null newed MGAttrib from the type id TID.
MGGLAttrib* MGNullGLAttrib(long TID){
	MGGelFactoryRegistry* reg = MGGelFactoryRegistry::get_instance();
	return static_cast<MGGLAttrib*>(reg->create_gel(TID));
}

AUTO_GEL_REGISTER(MGLights, MGLIGHTS_TID);
AUTO_GEL_REGISTER(MGFog, MGFOG_TID);

AUTO_GEL_REGISTER(MGLight, MGLIGHT_TID);
AUTO_GEL_REGISTER(MGDirectionalLight, MGDIRECTIONAL_LIGHT_TID);
AUTO_GEL_REGISTER(MGPointLight, MGPOINT_LIGHT_TID);
AUTO_GEL_REGISTER(MGSpotLight, MGSPOT_LIGHT_TID);

AUTO_GEL_REGISTER(MGAlphaFunc, MGALPHA_FUNC_TID);
AUTO_GEL_REGISTER(MGBlendFunc, MGBLEND_FUNC_TID);
AUTO_GEL_REGISTER(MGColor, MGCOLOR_TID);
AUTO_GEL_REGISTER(MGColorMask, MGCOLOR_MASK_TID);
AUTO_GEL_REGISTER(MGDepthFunc, MGDEPTH_FUNC_TID);
AUTO_GEL_REGISTER(MGDepthMask, MGDEPTH_MASK_TID);
AUTO_GEL_REGISTER(MGLightEnable, MGLIGHT_ENABLE_TID);
AUTO_GEL_REGISTER(MGLineStipple, MGLINE_STIPPLE_TID);
AUTO_GEL_REGISTER(MGLineWidth, MGLINE_WIDTH_TID);
AUTO_GEL_REGISTER(MGPolygonMode, MGPOLYGON_MODE_TID);
AUTO_GEL_REGISTER(MGRenderAttr, MGRENDER_ATTR_TID);
AUTO_GEL_REGISTER(MGTranspMode, MGTRANSP_MODE_TID);

