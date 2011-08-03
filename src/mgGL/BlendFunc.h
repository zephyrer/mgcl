/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGBlendFunc_HH_
#define _MGBlendFunc_HH_

#include "mgGL/GLAttrib.h"

//
//Define MGBlendFunc Class.

class MGOfstream;
class MGIfstream;

/** @addtogroup GLAttrib
 *  @{
 */

///MGBlendFunc defines the attributes of Blend function.
class MGCLASS MGBlendFunc:public MGGLAttrib{

public:

enum SrcBlendFunc{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	ZERO_SBLEND=GL_ZERO,
	ONE_SBLEND=GL_ONE,
	DST_COLOR_SBLEND=GL_DST_COLOR,
	ONE_MINUS_DST_COLOR_SBLEND=GL_ONE_MINUS_DST_COLOR,
	SRC_ALPHA_SATURATE_SBLEND=GL_SRC_ALPHA_SATURATE, 
	SRC_ALPHA_SBLEND=GL_SRC_ALPHA,
	ONE_MINUS_SRC_ALPHA_SBLEND=GL_ONE_MINUS_SRC_ALPHA,
	DST_ALPHA_SBLEND=GL_DST_ALPHA, 
	ONE_MINUS_DST_ALPHA_SBLEND=GL_ONE_MINUS_DST_ALPHA
};

enum DstBlendFunc{
	ZERO_DBLEND=GL_ZERO,
	ONE_DBLEND=GL_ONE,
	SRC_COLOR_DBLEND=GL_SRC_COLOR,
	ONE_MINUS_SRC_COLOR_DBLEND=GL_ONE_MINUS_SRC_COLOR,
	SRC_ALPHA_DBLEND=GL_SRC_ALPHA,
	ONE_MINUS_SRC_ALPHA_DBLEND=GL_ONE_MINUS_SRC_ALPHA,
	DST_ALPHA_DBLEND=GL_DST_ALPHA,
	ONE_MINUS_DST_ALPHA_DBLEND=GL_ONE_MINUS_DST_ALPHA
};

MGBlendFunc(SrcBlendFunc sf=UNDEFINED, DstBlendFunc df=ZERO_DBLEND);

///Assignment
MGBlendFunc& operator=(const MGGel& gel2);
MGBlendFunc& operator=(const MGBlendFunc& gel2);

///comparison
bool operator<(const MGBlendFunc& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function////////

///Generate a newed clone object.
MGBlendFunc* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_func(
	SrcBlendFunc sf=ONE_SBLEND,
	DstBlendFunc df=ZERO_DBLEND
);

SrcBlendFunc get_sfunc()const{return static_cast<SrcBlendFunc>(data());};
DstBlendFunc get_dfunc()const{return m_dst_blend_func;};
GLenum GLsfunc()const{return static_cast<GLenum>(data());};
GLenum GLdfunc()const{return static_cast<GLenum>(m_dst_blend_func);};

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
long identify_type() const{return MGBLEND_FUNC_TID;};

std::string whoami()const{return "BlendFunc";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

private:

	DstBlendFunc m_dst_blend_func;///<destination blending func.

};

/** @} */ // end of GLAttrib group
#endif
