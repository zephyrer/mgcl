/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLightEnable_HH_
#define _MGLightEnable_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

/** @addtogroup GLAttrib
 *  @{
 */

///Define MGLightEnable Class.
///MGLightEnable defines the attributes of shading model.
class MGCLASS MGLightEnable:public MGGLAttrib{

public:

enum MODE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	ENABLED=1,
};

MGLightEnable(MODE m=UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};

///Assignment
MGLightEnable& operator=(const MGGel& gel2);
MGLightEnable& operator=(const MGLightEnable& gel2);

///comparison
bool operator<(const MGLightEnable& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGLightEnable* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_enabled(){data()=ENABLED;};

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
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LIGHTING_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LIGHTING_BIT);};

/// Return This object's typeID
long identify_type() const{return MGLIGHT_ENABLE_TID;};

std::string whoami()const{return "LightEnable";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
std::ostream& out(std::ostream&) const;

private:

};

/** @} */ // end of GLAttrib group
#endif
