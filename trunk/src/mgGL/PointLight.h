/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGPOINTLIGHT_HH_
#define _MGPOINTLIGHT_HH_

#include "mgGL/Light.h"

class MGOfstream;
class MGIfstream;
class MGPosition;

/** @addtogroup GLAttrib
 *  @{
 */

///MGPointLight is a point light source that radiates equally in all directions.
///The range of a MGPointLight's effect is localized to m_radius from
///the location(m_location).
class MGPointLight:public MGLight{

public:

///////// Constructors  ///////////////

MGPointLight();
MGPointLight(
    float intensity,		///<applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	///<applied to GL_AMBIENT
    const float color[3],	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGPosition& location,
	float radius,
	const float attenuation[3]///<[0]=GL_CONSTANT_ATTENUATION,
							///<[1]=GL_LINEAR_ATTENUATION
							///<[2]=GL_QUADRATIC_ATTENUATION
);

///virtual ~MGPointLight();

///Assignment
MGPointLight& operator=(const MGGel& gel2);
MGPointLight& operator=(const MGPointLight& gel2);

///comparison
bool operator<(const MGPointLight& gel2)const;
bool operator<(const MGGel& gel2)const;

///Generate a newed clone object.
virtual MGPointLight* clone()const;

//////////// Set/Get /////////////

void setLocation(const MGPosition& location);
void setLocation(const float location[3]){
	for(size_t i=0; i<3; i++) m_location[i]=location[i];
}
void setLocation(float x, float y, float z){
	m_location[0]=x;
	m_location[1]=y;
	m_location[2]=z;
}
void getLocation(MGPosition& location)const;
void getLocation(float location[3])const{
	for(size_t i=0; i<3; i++) location[i]=m_location[i];
}
void getLocation(float& x, float& y, float& z)const{
	x=m_location[0];
	y=m_location[1];
	z=m_location[2];
}

void setRadius(float radius){m_radius=radius;};
float getRadius()const{return m_radius;};

void setAttenuation(const float attenuation[3]){
	for(size_t i=0; i<3; i++) m_attenuation[i]=attenuation[i];
}
void setAttenuation(float const_att, float linear_att, float quadratic_att){
	m_attenuation[0]=const_att;
	m_attenuation[1]=linear_att;
	m_attenuation[2]=quadratic_att;
}
void getAttenuation(float attenuation[3])const{
	for(size_t i=0; i<3; i++) attenuation[i]=m_attenuation[i];
}
void getAttenuation(float& const_att, float& linear_att, float& quadratic_att)const{
	const_att=m_attenuation[0];
	linear_att=m_attenuation[1];
	quadratic_att=m_attenuation[2];
}

///exec Attribute process.
///Function's return value is the lightnumber of this light executed.
virtual size_t exec()const;

/// Return This object's typeID
long identify_type() const{return MGPOINT_LIGHT_TID;};

std::string whoami()const{return "PointLight";};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

///Transform the gel by the argument.

///translation
virtual void transform(const MGVector& v);

///scaling.
virtual void transform(double scale);

///matrix transformation.
virtual void transform(const MGMatrix& mat);

///general transformation.
virtual void transform(const MGTransf& tr);

private:

    float m_location[4];	///<GL_POSITION
    float m_radius;			///<radius this point light reaches.
	float m_attenuation[3];	///<[0]=GL_CONSTANT_ATTENUATION,
							///<[1]=GL_LINEAR_ATTENUATION
							///<[2]=GL_QUADRATIC_ATTENUATION

};

/** @} */ // end of GLAttrib group
#endif // _MGPOINTLIGHT_HH_
