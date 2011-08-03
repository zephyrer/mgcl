/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGDepthFunc_HH_
#define _MGDepthFunc_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGDepthFunc Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGDepthFunc defines a comparison function for depth test.
///Depth test is used to judge a pixel drawn is visible or not,
///a kind of attributes of shading model.
class MGCLASS MGDepthFunc:public MGGLAttrib{

public:

enum FUNC{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	NEVER=GL_NEVER,
	LESS=GL_LESS,
	EQUAL=GL_EQUAL,
	LEQUAL=GL_LEQUAL,
	GREATER=GL_GREATER,
	NOTEQUAL=GL_NOTEQUAL,
	GEQUAL=GL_GEQUAL,
	ALWAYS=GL_ALWAYS
};

MGDepthFunc(FUNC m=UNDEFINED):MGGLAttrib(m){;};

///Assignment
MGDepthFunc& operator=(const MGGel& gel2);
MGDepthFunc& operator=(const MGDepthFunc& gel2);

///comparison
bool operator<(const MGDepthFunc& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function////////

///Generate a newed clone object.
MGDepthFunc* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_func(FUNC m=LESS){data()=static_cast<FUNC>(m);};
FUNC get_func()const{return static_cast<FUNC>(data());};
GLenum GLfunc()const{return static_cast<GLenum>(data());};

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///if true, color attribute will be neglected.
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
long identify_type() const{return MGDEPTH_FUNC_TID;};

std::string whoami()const{return "DepthFunc";};

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
