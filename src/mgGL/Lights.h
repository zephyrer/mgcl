/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGLIGHTS_HH_
#define _MGLIGHTS_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;
class MGLight;

/** @addtogroup GLAttrib
 *  @{
 */

///a container class for light sources(MGDirectionalLight, MGPointLight, or MGSpotLight).
class MGLights:public MGGLAttrib{

public:

	typedef MGPvector<MGLight>::iterator iterator ;
	typedef MGPvector<MGLight>::const_iterator const_iterator;
	typedef MGPvector<MGLight>::reverse_iterator reverse_iterator ;
	typedef MGPvector<MGLight>::const_reverse_iterator const_reverse_iterator;

// Constructors

MGLights():MGGLAttrib(1){;};

///Construct MGLights of one light.
MGLights(
	MGLight* light
):MGGLAttrib(1),m_lights(light){;};

///virtual ~MGLights();

public:

///Assignment
MGLights& operator=(const MGGel& gel2);
MGLights& operator=(const MGLights& gel2);

///comparison
bool operator<(const MGLights& gel2)const;
bool operator<(const MGGel& gel2)const;

const MGLight* back()const{return m_lights.back();};
MGLight* back(){return m_lights.back();};
iterator begin(){return m_lights.begin();};
const_iterator begin()const{return m_lights.begin();};
void clear(){m_lights.clear();};
bool empty()const{return m_lights.empty();};
iterator end(){return m_lights.end();};
const_iterator end()const{return m_lights.end();};
const MGLight* front()const{return m_lights.front();};
MGLight* front(){return m_lights.front();};
const MGLight* operator[](size_t i)const{return m_lights[i];};
MGLight* operator[](size_t i){return m_lights[i];};
void pop_back(){m_lights.pop_back();};

///Generate a newed clone object.
virtual MGLights* clone()const;

///render GLAttribute process.
virtual void exec()const;

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<<if true, color attribute will be neglected.
)const{;};

///add one light.
///Function's return value is the numbe of lights defined.
size_t push_back(MGLight* light);

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LIGHTING_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LIGHTING_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{;};

///Obtain the light number defined.
size_t size()const{return m_lights.size();};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{;};

/// Return This object's typeID
virtual long identify_type() const{return MGLIGHTS_TID;};

std::string whoami()const{return "Lights";};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

private:

	MGPvector<MGLight> m_lights;///<vector of MGLight.

};

/** @} */ // end of GLAttrib group
#endif // _MGLIGHTS_HH_
