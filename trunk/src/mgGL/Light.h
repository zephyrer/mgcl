/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGLIGHT_HH_
#define _MGLIGHT_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;
class MGMatrix;
class MGSpotLight;
class MGDirectionalLight;
class MGPointLight;

/** @addtogroup GLAttrib
 *  @{
 */

///an abstract base class for light sources(MGDirectionalLight,
///MGPointLight, or MGSpotLight).
class MGLight:public MGGLAttrib{

public:

MGDECL friend std::ostream& operator<< (std::ostream&, const MGLight&);

//////////Constructor///////////////

MGLight();
MGLight(
    float intensity,		///<applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	///<applied to GL_AMBIENT
    const float color[3]	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
);

///virtual ~MGLight();

public:

///Assignment

///Assignment
virtual MGLight& operator=(const MGGel& gel2);
virtual MGLight& operator=(const MGLight& gel2);

///comparison
virtual bool operator<(const MGLight& gel2)const;
virtual bool operator<(const MGGel& gel2)const;

///Generate a newed clone object.
virtual MGLight* clone()const;

///render process.
///Function's return value is the lightnumber of this light executed.
virtual size_t exec()const;

///draw GLAttribute process.
virtual void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{;};

///Obtain the light number of this.
GLenum get_light_num()const;

///render GLAttribute process.
virtual void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_draw_attrib_mask(unsigned int& mask)const{;};

///Set light number.
void set_light_number(size_t lnum){data()=lnum;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_draw_attrib_mask(unsigned int& mask)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LIGHTING_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LIGHTING_BIT);};

//////////// Set/Get  /////////////

void setIntensity(float intensity){m_intensity=intensity;};
float getIntensity()const{return m_intensity;};

void setAmbientIntensity(float ambientIntensity){
	m_ambientIntensity=ambientIntensity;
};
float getAmbientIntensity()const{return m_ambientIntensity;};

void setColor(const float color[3]){
	for(size_t i=0;i<3; i++) m_color[i]=color[i];
};
void getColor(float color[3]){
	for(size_t i=0;i<3; i++) color[i]=m_color[i];
};
void setColor(float v0, float v1, float v2){
	m_color[0]=v0;
	m_color[1]=v1;
	m_color[2]=v2;
}
void getColor(float& v0, float& v1, float& v2){
	v0=m_color[0];
	v1=m_color[1];
	v2=m_color[2];
}

virtual std::string whoami()const{return "Light";};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

protected:

    float m_intensity;	///<applied to GL_DIFFUSE and GL_SPECULAR
    float m_ambientIntensity;///<applied to GL_AMBIENT
    float m_color[3];	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR

///assignment
MGLight& set_light(const MGLight& gel2);

};

/** @} */ // end of GLAttrib group
#endif // _MGLIGHT_HH_
