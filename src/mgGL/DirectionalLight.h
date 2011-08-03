/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGDIRECTIONAL_LIGHT_HH_
#define _MGDIRECTIONAL_LIGHT_HH_

#include "mgGL/Light.h"

class MGOfstream;
class MGIfstream;
class MGVector;

/** @addtogroup GLAttrib
 *  @{
 */

///MGDirectionalLight is a directional light source that approximates infinite light sources
///as the sun.
///Can improve rendering performance over local light sources,
///such as MGPointLight and MGSpotLight. Use MGDirectionalLight to set the direction
///of general lighting for a scene.
class MGDirectionalLight : public MGLight{

public:
	
///////////// Constructors & destructor /////////////

MGDirectionalLight();
MGDirectionalLight(
    float intensity,		///<applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	///<applied to GL_AMBIENT
    const float color[3],	///<applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGVector& direction
);

///~MGDirectionalLight();

///Assignment
MGDirectionalLight& operator=(const MGGel& gel2);
MGDirectionalLight& operator=(const MGDirectionalLight& gel2);

///comparison
bool operator<(const MGDirectionalLight& gel2)const;
bool operator<(const MGGel& gel2)const;

///Generate a newed clone object.
MGDirectionalLight* clone()const;

///render GLAttribute process.
///Function's return value is the lightnumber of this light executed.
size_t exec()const;
	
///////// Set/Get  /////////////

void setDirection(const MGVector& direction);
void setDirection(const float direction[3]){
	for(size_t i=0; i<3; i++) m_direction[i]=direction[i];
}
void setDirection(float v0, float v1, float v2){
	m_direction[0]=v0;
	m_direction[1]=v1;
	m_direction[2]=v2;
}

void getDirection(MGVector& direction)const;
void getDirection(float direction[3])const{
	for(size_t i=0; i<3; i++) direction[i]=m_direction[i];
}
void getDirection(float& v0, float& v1, float& v2)const{
	v0=m_direction[0];
	v1=m_direction[1];
	v2=m_direction[2];
}

/// Return This object's typeID
long identify_type() const{return MGDIRECTIONAL_LIGHT_TID;};

std::string whoami()const{return "DirectionalLight";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

///Transform the gel by the argument.

///translation
void transform(const MGVector& v);

///scaling.
void transform(double scale);

///matrix transformation.
void transform(const MGMatrix& mat);

///general transformation.
void transform(const MGTransf& tr);

private:

    float m_direction[3];	///<GL_SPOT_DIRECTION

};

/** @} */ // end of GLAttrib group
#endif // _MGDIRECTIONAL_LIGHT_HH_
