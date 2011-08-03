/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGFog_HH_
#define _MGFog_HH_

#include "mgGL/GLAttrib.h"
#include "mgGL/Color.h"
#include <iosfwd>

class MGOfstream;
class MGIfstream;

/** @addtogroup GLAttrib
 *  @{
 */

///Defines OpenGL's fog data.
///We assume glHint(GL_FOG_HINT,xxx) is invoked outside MGFog.
///MGFog does not hold the GL_FOG_HINT attribute.
class MGCLASS MGFog : public MGGLAttrib{

public:

enum FOGMODE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	LINEAR=GL_LINEAR,
	EXP=GL_EXP,
	EXP2=GL_EXP2
};

MGFog();
MGFog(FOGMODE mode);
MGFog(FOGMODE mode,float density,float start,float end, const MGColor& color);

///Assignment
MGFog& operator=(const MGGel& gel2);
MGFog& operator=(const MGFog& gel2);

///comparison
bool operator<(const MGFog& gel2)const;
bool operator<(const MGGel& gel2)const;

///Generate a newed clone object.
MGFog* clone()const;

void setMode(FOGMODE mode){m_flag=mode;};
FOGMODE	fog_mode()const{return static_cast<FOGMODE>(m_flag);};
GLenum GLfog_mode()const{return static_cast<GLenum>(m_flag);};

void setDensity(float density){m_density=density;};
float getDensity()const{return m_density;};

void setStart(float start){m_start=start;};
float getStart()const{return m_start;};

void setEnd(float end){m_end=end;};
float getEnd()const{return m_end;};

void setColor(const MGColor& color){
	m_color=color;
}
void setColor(float v0, float v1, float v2, float v3){
	m_color=MGColor(v0,v1,v2,v3);
}
void getColor(float& v0, float& v1, float& v2, float& v3)const{
	const float* clrs=m_color.color();
	v0=clrs[0];v1=clrs[1];v2=clrs[2];v3=clrs[3];
}
const MGColor& getColor()const{return m_color;};

///render GLAttribute process.
void exec()const;

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{;};

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,FOG_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,FOG_BIT);};

/// Return This object's typeID
long identify_type() const{return MGFOG_TID;};

std::string whoami()const{return "Fog";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

private:

    float		m_density;		///<GL_FOG_DENSITY
 								///<Used when fog_mode()=EXP_FOG or EXP2_FOG.
	float		m_start, m_end;	///<GL_FOG_START, GL_FOG_END
								///<Used only when fog_mode()=LINEAR_FOG
    MGColor		m_color;		///<GL_FOG_COLOR

};

/** @} */ // end of GLAttrib group
#endif // _MGFog_HH_
