/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGAlphaFunc_HH_
#define _MGAlphaFunc_HH_

#include "mgGL/GLAttrib.h"

#include <iosfwd>

class MGOfstream;
class MGIfstream;

//
//Define MGAlphaFunc Class.

/** @defgroup GLAttrib OpenGL Attributes Class
 *  @ingroup DisplayHandling
 */

/** @addtogroup GLAttrib
 *  @{
 */

///MGAlphaFunc defines the attributes of alpha function.
class MGCLASS MGAlphaFunc:public MGGLAttrib{

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

MGAlphaFunc(FUNC f=UNDEFINED):MGGLAttrib(f){;};
MGAlphaFunc(FUNC f,float ref):MGGLAttrib(f),m_refval(ref){;};

///Assignment
MGAlphaFunc& operator=(const MGGel& gel2);
MGAlphaFunc& operator=(const MGAlphaFunc& gel2);

///comparison
bool operator<(const MGAlphaFunc& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function////////

///Generate a newed clone object.
MGAlphaFunc* clone()const;

///Execute this attribute using OpenGL function(s).
void exec()const;

void set_func(FUNC f=ALWAYS){m_flag=f;};

FUNC get_func()const{return static_cast<FUNC>(data());};
GLenum GLfunc()const{return static_cast<GLenum>(data());};

float ref_val()const{return m_refval;};
void set_ref_val(float ref){m_refval=ref;};

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
long identify_type() const{return MGALPHA_FUNC_TID;};

std::string whoami()const{return "AlphaFunc";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

private:

	float m_refval;	///reference value. clamped value between 0. and 1.

};

/** @} */ // end of GLAttrib group
#endif
