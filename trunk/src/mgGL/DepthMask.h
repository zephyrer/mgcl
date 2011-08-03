/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDepthMask_HH_
#define _MGDepthMask_HH_
#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGDepthMask Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGDepthMask defines a depth mask, an attribute of shading model.
class MGCLASS MGDepthMask:public MGGLAttrib{

public:

enum MASK{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	ENABLED=1
};

MGDepthMask(MASK m=UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};

///Assignment
MGDepthMask& operator=(const MGGel& gel2);
MGDepthMask& operator=(const MGDepthMask& gel2);

///comparison
bool operator<(const MGDepthMask& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGDepthMask* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_mask_enabled(){data()=ENABLED;};
void set_mask_disabled(){data()=DISABLED;};
MASK get_mask_data()const{return static_cast<MASK>(data());};

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec();};

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,DEPTH_BUFFER_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{reset_Amask(mask,DEPTH_BUFFER_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,DEPTH_BUFFER_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,DEPTH_BUFFER_BIT);};

/// Return This object's typeID
long identify_type() const{return MGDEPTH_MASK_TID;};

std::string whoami()const{return "DepthMask";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

private:

};

/** @} */ // end of GLAttrib group
#endif
