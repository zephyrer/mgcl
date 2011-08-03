/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Material.h"
#include "mgGL/RenderAttr.h"
#include "mgGL/Texture.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGRenderAttr Class.
//MGRenderAttr defines the attributes of rendering attributes.
//These attrubutes are not used for drawing(line drawing mode).

//copy constructor.
MGRenderAttr::MGRenderAttr(const MGRenderAttr& attr)
:MGGLAttrib(attr.m_flag),m_shade_model(attr.m_shade_model){
	MGRenderAttr* attr2=const_cast<MGRenderAttr*>(&attr);
	m_material=attr2->m_material; attr2->m_material=0;
	m_back_material=attr2->m_back_material; attr2->m_back_material=0;
	m_texture=attr2->m_texture; attr2->m_texture=0;
}

//Destructor.
MGRenderAttr::~MGRenderAttr(){
	delete m_material;
	delete m_back_material;
	delete m_texture;
}

//Assignment.
MGRenderAttr& MGRenderAttr::operator=(const MGRenderAttr& attr){
	if(this==&attr)
		return *this;

	MGGLAttrib::operator=(attr);
	MGRenderAttr* attr2=const_cast<MGRenderAttr*>(&attr);
	m_shade_model=attr.m_shade_model;
	delete m_material;
	m_material=attr2->m_material; attr2->m_material=0;
	delete m_back_material;
	m_back_material=attr2->m_back_material; attr2->m_back_material=0;
	delete m_texture;
	m_texture=attr2->m_texture; attr2->m_texture=0;
	return *this;
}
MGRenderAttr& MGRenderAttr::operator=(const MGGel& gel2){
	const MGRenderAttr* gel2_is_this=dynamic_cast<const MGRenderAttr*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGRenderAttr::operator<(const MGRenderAttr& gel2)const{
	return m_flag<gel2.m_flag;
}
bool MGRenderAttr::operator<(const MGGel& gel2)const{
	const MGRenderAttr* gel2_is_this=dynamic_cast<const MGRenderAttr*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//////////Member Function//////////

//Generate a newed clone object.
MGRenderAttr* MGRenderAttr::clone()const{
	return new MGRenderAttr(*this);
}

//Invoke appropriate OpenGL fucntion to this attribute.
void MGRenderAttr::exec()const{
	m_shade_model.exec();

	if(m_material){
		RENDERSIDE rs=render_side();
		if(rs==FRONT_AND_BACK)
			glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		else
			glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

		if(rs==BACK) glFrontFace(GL_CW);
		else glFrontFace(GL_CCW);

		if(rs==FRONT || rs==BACK) m_material->exec(FRONT);
		else if(m_back_material){
			m_material->exec(FRONT);
			m_back_material->exec(BACK);
		}else
			m_material->exec(FRONT_AND_BACK);
	}
	//*********TEXTURE PROCESS is necessary here.*********
}

//Set material.
void set_material_data(
	MGMaterial*& material,
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	if(!material){
		material=new MGMaterial(ambient,diffuse,specular,emission,shininess,transparency);
	}else{
		material->setAmbientColor(ambient);
		material->setDiffuseColor(diffuse);
		material->setSpecularColor(specular);
		material->setEmissiveColor(emission);
		material->setShininess(shininess);
		material->setTransparency(transparency);
	}
}

bool MGRenderAttr::material_defined()const{
	if(!m_material) return false;
	return true;
}

//Set the material. When rs=FRONT_AND_BACK and different material for the back side
//is used, set_back_material must be invoked after invoking set_material.
//Else the same material will be appllied for the both sides.
void MGRenderAttr::set_material(
	RENDERSIDE rs,
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	set_render_side(rs);
	set_material_data(m_material,
		ambient,diffuse,specular,emission,shininess,transparency);
	delete m_back_material;
	m_back_material=0;
}

//Set the back side material. Invoking set_back_material means two sided material
//and setting different material to the back side.
//Before use of set_back_material, set_material must be invoked first.
//set_back_material will set two sided material.
//set_material must be invoked before invoking set_back_material.
//Else the same material will be appllied for the both sides.
void MGRenderAttr::set_back_material(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	assert(m_material);//must be defined.
	set_render_side(FRONT_AND_BACK);
	set_material_data(m_back_material,
		ambient,diffuse,specular,emission,shininess,transparency);
}

//Return the material pointer of the back side when this has the two sided rendering
//(FRONT_AND_BACK) and the different materials are applied to each side.
//NULL will be returned when render_side()!=FRONT_AND_BACK. 
//Even when render_side()==FRONT_AND_BACK, if m_back_material==null(the same material
//is applied to both side), null will be returned.
const MGMaterial* MGRenderAttr::back_material()const{
	if(render_side()!=FRONT_AND_BACK) return 0;
	return m_back_material;
}

//Debug Function.
std::ostream& MGRenderAttr::out(std::ostream& ostrm)const{
	ostrm<<std::endl<<"RenderAttr="; MGGLAttrib::out(ostrm);
	if(render_side()==FRONT) ostrm<<"FRONT"; else if(render_side()==BACK) ostrm<<"BACK";
	else if(render_side()==FRONT_AND_BACK) ostrm<<"FRONT_AND_BACK";
	ostrm<<",ShadeModel=";m_shade_model.out(ostrm);
	ostrm<<",Material="<<m_material;if(m_material) ostrm<<"=";m_material->out(ostrm);
	if(back_material()) ostrm<<",Back Material=";m_back_material->out(ostrm);
	ostrm<<",Texture="<<m_texture;if(m_texture) ostrm<<"=";m_texture->out(ostrm);
	return ostrm;
}

// Serialization.
//Write all member data
void MGRenderAttr::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;

	buf<<m_shade_model;
	if(m_material){
		buf<<0xffffffffL;
		//The header that indicates an MGMaterial object is followed.
		buf<<(*m_material);
	}else{
		//When null pointer.
		buf<<0x00000000L;
	}
	if(m_back_material){
		buf<<0xffffffffL;
		//The header that indicates an MGMaterial object is followed.
		buf<<(*m_back_material);
	}else{
		//When null pointer.
		buf<<0x00000000L;
	}
	if(m_texture){
		buf<<0xffffffffL;
		//The header that indicates an MGTexture object is followed.
		buf<<(*m_texture);
	}else{
		//When null pointer.
		buf<<0x00000000L;
	}
}
void MGRenderAttr::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;

	buf>>m_shade_model;
	long addr;
	buf>>addr;
	if(addr){
		m_material=new MGMaterial();
		buf>>(*m_material);
	}else{
		//When null pointer.
		m_material=0;
	}

	buf>>addr;
	if(addr){
		m_back_material=new MGMaterial();
		buf>>(*m_back_material);
	}else{
		//When null pointer.
		m_back_material=0;
	}

	buf>>addr;
	if(addr){
		m_texture=new MGTexture();
		buf>>(*m_texture);
	}else{
		//When null pointer.
		m_texture=0;
	}
}

//Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void MGRenderAttr::set_render_attrib_mask(unsigned int& mask)const{
	set_Amask(mask,LIGHTING_BIT);
	if(m_texture) set_Amask(mask,TEXTURE_BIT);
}

//Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void MGRenderAttr::reset_render_attrib_mask(unsigned int& mask)const{
	reset_Amask(mask,LIGHTING_BIT);
	if(m_texture) reset_Amask(mask,TEXTURE_BIT);
}
