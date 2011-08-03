/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGMaterial_HH_
#define _MGMaterial_HH_

#include "mgGL/GLAttrib.h"
#include "mgGL/RenderAttr.h"

class MGOfstream;
class MGIfstream;

/** @addtogroup GLAttrib
 *  @{
 */

///MGMaterial defines OpenGL's Material attributes.
///See the private member data.
class MGCLASS MGMaterial{
public:

/// Serialization fucntion.
MGDECL friend MGOfstream& operator<< (MGOfstream& buf, const MGMaterial& mt);
MGDECL friend MGIfstream& operator>> (MGIfstream& buf, MGMaterial& mt);

///////// Constructors and Destructors ////////////
MGMaterial();
MGMaterial(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
);

//////////////////Member function //////////////////

///render GLAttribute process.
void exec(MGRenderAttr::RENDERSIDE rs)const;

///////// Set and Get ///////////////////

////Ambient color
void setAmbientColor(const float ambientColor[3]){
	for(size_t i=0; i<3; i++) m_ambientColor[i]=ambientColor[i];
}
void getAmbientColor(float ambientColor[3])const{
	for(size_t i=0; i<3; i++) ambientColor[i]=m_ambientColor[i];
}
void setAmbientColor(float v0, float v1, float v2){
	m_ambientColor[0]=v0;
	m_ambientColor[1]=v1;
	m_ambientColor[2]=v2;
}

////Diffuse color
void setDiffuseColor(const float diffuseColor[3]){
	for(size_t i=0; i<3; i++) m_diffuseColor[i]=diffuseColor[i];
}
void getDiffuseColor(float diffuseColor[3])const{
	for(size_t i=0; i<3; i++) diffuseColor[i]=m_diffuseColor[i];
}
void setDiffuseColor(float v0, float v1, float v2){
	m_diffuseColor[0]=v0;
	m_diffuseColor[1]=v1;
	m_diffuseColor[2]=v2;
}

////Specular color
void setSpecularColor(const float specularColor[3]){
	for(size_t i=0; i<3; i++) m_specularColor[i]=specularColor[i];
}
void getSpecularColor(float specularColor[3])const{
	for(size_t i=0; i<3; i++) specularColor[i]=m_specularColor[i];
}
void setSpecularColor(float v0, float v1, float v2){
	m_specularColor[0]=v0;
	m_specularColor[1]=v1;
	m_specularColor[2]=v2;
}

////Emissive color
void setEmissiveColor(const float emissiveColor[3]){
	for(size_t i=0; i<3; i++) m_emissiveColor[i]=emissiveColor[i];
}
void getEmissiveColor(float emissiveColor[3])const{
	for(size_t i=0; i<3; i++) emissiveColor[i]=m_emissiveColor[i];
}
void setEmissiveColor(float v0, float v1, float v2){
	m_emissiveColor[0]=v0;
	m_emissiveColor[1]=v1;
	m_emissiveColor[2]=v2;
}

////Shiness
void setShininess(float shininess){m_shininess=shininess;};
float getShininess()const{return m_shininess;};

////Transparency
void setTransparency(float transparency);
float getTransparency()const{return m_ambientColor[3];};

/// Return This object's typeID
long identify_type() const{return MGMATERIAL_TID;};

/// Output function.
std::ostream& out(std::ostream& ostrm) const;

private:

    float m_ambientColor[4];
    float m_diffuseColor[4];
    float m_specularColor[4];
    float m_emissiveColor[4];
    float m_shininess;

};

/** @} */ // end of GLAttrib group
#endif // _MGMaterial_HH_
