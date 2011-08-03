/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Material.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGMaterial defines OpenGL's Material attributes.
//See the private member data.

/////// Constructors and Destructors //////////
MGMaterial::MGMaterial():m_shininess(0.){
    m_ambientColor[0]=m_ambientColor[1]=m_ambientColor[2]=float(.2);
    m_diffuseColor[0]=m_diffuseColor[1]=m_diffuseColor[2]=float(.8);
    m_specularColor[0]=m_specularColor[1]=m_specularColor[2]=float(.0);
    m_emissiveColor[0]=m_emissiveColor[1]=m_emissiveColor[2]=float(.0);
	m_ambientColor[3]=1.;//Transparency is 1.
	m_diffuseColor[3]=1.;
	m_specularColor[3]=1.;
	m_emissiveColor[3]=1.;
}

MGMaterial::MGMaterial(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
):m_shininess(shininess){
	for(size_t i=0; i<3; i++){
		m_ambientColor[i]=ambient[i];
		m_diffuseColor[i]=diffuse[i];
		m_specularColor[i]=specular[i];
		m_emissiveColor[i]=emission[i];
	}
	m_ambientColor[3]=transparency;
	m_diffuseColor[3]=transparency;
	m_specularColor[3]=transparency;
	m_emissiveColor[3]=transparency;
}

///////////////Member function ///////////////

//render GLAttribute process.
void MGMaterial::exec(MGRenderAttr::RENDERSIDE rs)const{
	GLenum rs2=GLrender_side(rs);
	glMaterialfv(rs2,GL_AMBIENT,m_ambientColor);
	glMaterialfv(rs2,GL_DIFFUSE,m_diffuseColor);
	glMaterialfv(rs2,GL_SPECULAR,m_specularColor);
	glMaterialfv(rs2,GL_EMISSION,m_emissiveColor);
	glMaterialf (rs2,GL_SHININESS,m_shininess);
}

//Set the transparency.
void MGMaterial::setTransparency(float transparency){
	m_ambientColor[3]=transparency;
	m_diffuseColor[3]=transparency;
	m_specularColor[3]=transparency;
	m_emissiveColor[3]=transparency;
}

// Output function.
std::ostream& MGMaterial::out(std::ostream& ostrm) const{
	ostrm<<std::endl<<"Material="<<this;
	ostrm<<",ambientColor["<<m_ambientColor[0]<<","<<m_ambientColor[1]<<","<<m_ambientColor[2];
	ostrm<<"]"<<std::endl;
	ostrm<<",diffuseColor["<<m_diffuseColor[0]<<","<<m_diffuseColor[1]<<","<<m_diffuseColor[2];
	ostrm<<"]"<<std::endl;
	ostrm<<",specularColor["<<m_specularColor[0]<<","<<m_specularColor[1]<<","<<m_specularColor[2];
	ostrm<<"]"<<std::endl;
	ostrm<<",emissiveColor["<<m_emissiveColor[0]<<","<<m_emissiveColor[1]<<","<<m_emissiveColor[2];
	ostrm<<"]"<<std::endl;
	ostrm<<",shininess="<<m_shininess;
	return ostrm;
}
  
// Serialization fucntion.
MGOfstream& operator<< (MGOfstream& buf, const MGMaterial& mt){
	size_t i;
	for(i=0; i<4; i++) buf<<mt.m_ambientColor[i];
	for(i=0; i<4; i++) buf<<mt.m_diffuseColor[i];
	for(i=0; i<4; i++) buf<<mt.m_specularColor[i];
	for(i=0; i<4; i++) buf<<mt.m_emissiveColor[i];
	buf<<mt.m_shininess;
	return buf;
}
MGIfstream& operator>> (MGIfstream& buf, MGMaterial& mt){
	size_t i;
	for(i=0; i<4; i++) buf>>mt.m_ambientColor[i];
	for(i=0; i<4; i++) buf>>mt.m_diffuseColor[i];
	for(i=0; i<4; i++) buf>>mt.m_specularColor[i];
	for(i=0; i<4; i++) buf>>mt.m_emissiveColor[i];
	buf>>mt.m_shininess;
	return buf;
}
