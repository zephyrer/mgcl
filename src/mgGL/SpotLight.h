/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGSPOTLIGHT_HH_
#define _MGSPOTLIGHT_HH_

#include "mgGL/PointLight.h"

/** @addtogroup GLAttrib
 *  @{
 */

///A directional light source.
///Because MGSpotLight is subclassed from MGPointLight,
///MGSpotLight can be positioned in a world. 
///The light emanates as a cone; the axis of the cone(m_direction) specifies the direction
///of the spot light.
class MGCLASS MGSpotLight : public MGPointLight{

public:

/// Field enums
enum {
	DIRECTION=1,
	EXPONENT,
	CUT_OFF_ANGLE
};

//////////// Constructors ////////////////

MGSpotLight();
MGSpotLight(
    float intensity,		///<applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	///<applied to GL_AMBIENT
    const float color[3],	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGPosition& location,
	float radius,
	const float attenuation[3],///<[0]=GL_CONSTANT_ATTENUATION,
							///<[1]=GL_LINEAR_ATTENUATION
							///<[2]=GL_QUADRATIC_ATTENUATION
    const MGVector& direction,///<GL_SPOT_DIRECTION
    float exponent,		///<GL_SPOT_EXPONENT
    float cutOffAngle	///<GL_SPOT_CUTOFF
);

///~MGSpotLight();

///Assignment
MGSpotLight& operator=(const MGGel& gel2);
MGSpotLight& operator=(const MGSpotLight& gel2);

///comparison
bool operator<(const MGSpotLight& gel2)const;
bool operator<(const MGGel& gel2)const;

/////////////// Set/Get  /////////////

///Generate a newed clone object.
MGSpotLight* clone()const;

void setDirection(const MGVector& direction);
void setDirection(const float direction[3]){
	for(size_t i=0; i<3; i++) m_direction[i]=direction[i];
}
void setDirection(float v0, float v1, float v2){
	m_direction[0]=v0;
	m_direction[1]=v1;
	m_direction[2]=v2;
}
void getDirection(MGVector& direction);
void getDirection(float direction[3]){
	for(size_t i=0; i<3; i++) direction[i]=m_direction[i];
}
void getDirection(float& v0, float& v1, float& v2){
	v0=m_direction[0];
	v1=m_direction[1];
	v2=m_direction[2];
}

void setExponent(float exponent){m_exponent=exponent;};
float getExponent(){return m_exponent;};

void setCutOffAngle(float cutOffAngle){m_cutOffAngle=cutOffAngle;};
float getCutOffAngle(){return m_cutOffAngle;};

///exec Attribute process.
///Function's return value is the lightnumber of this light executed.
size_t exec()const;

/// Return This object's typeID
long identify_type() const{return MGSPOT_LIGHT_TID;};

std::string whoami()const{return "SpotLight";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

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

    float m_direction[3];	///<GL_SPOT_DIRECTION

    float m_exponent;		///<GL_SPOT_EXPONENT
    float m_cutOffAngle;	///<GL_SPOT_CUTOFF

};

/** @} */ // end of GLAttrib group
#endif // _MGSPOTLIGHT_HH_
