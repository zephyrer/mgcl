/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGColorMask_HH_
#define _MGColorMask_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGColorMask Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGColorMask defines the color mask data.
class MGCLASS MGColorMask:public MGGLAttrib{

public:

enum MASK{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED
};

MGColorMask(MASK m=UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};
MGColorMask(bool red, bool green, bool blue, bool alpha);

///Assignment
MGColorMask& operator=(const MGGel& gel2);
MGColorMask& operator=(const MGColorMask& gel2);

///comparison
bool operator<(const MGColorMask& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function////////

///Generate a newed clone object.
MGColorMask* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void get_mask(bool& red, bool& green, bool& blue, bool& alpha)const;
void set_mask(bool red, bool green, bool blue, bool alpha);

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec();};

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,COLOR_BUFFER_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{reset_Amask(mask,COLOR_BUFFER_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,COLOR_BUFFER_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,COLOR_BUFFER_BIT);};

/// Return This object's typeID
long identify_type() const{return MGCOLOR_MASK_TID;};

std::string whoami()const{return "ColorMask";};

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
