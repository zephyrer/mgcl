/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTranspMode_HH_
#define _MGTranspMode_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGTranspMode Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGTranspMode defines the mode of transparency
///******Currently this function is not supported yet**********////
class MGCLASS MGTranspMode:public MGGLAttrib{

public:

enum MODE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	FAST=1,
	NICE,
	BLEND,
	SCREEN_DOOR
};

MGTranspMode(MODE m=UNDEFINED);

///Assignment
MGTranspMode& operator=(const MGGel& gel2);
MGTranspMode& operator=(const MGTranspMode& gel2);

///comparison
bool operator<(const MGTranspMode& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGTranspMode* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const{;};

void set_transparency_mode(MODE m=BLEND){m_flag=m;};

MODE get_transparency_mode()const{return static_cast<MODE>(data());};

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
void set_render_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{;};

/// Return This object's typeID
long identify_type() const{return MGTRANSP_MODE_TID;};

std::string whoami()const{return "TranspMode";};

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
