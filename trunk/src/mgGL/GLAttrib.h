/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifdef MGCL_NO_MFC
#include "mgGL/GLAttribNOMFC.h"
#else

#ifndef _MGGLAttrib_HH_
#define _MGGLAttrib_HH_

#include <afxwin.h>
#include <gl\gl.h>
#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/Attrib.h"

//
//Define MGGLAttrib Class.

class MGOfstream;
class MGIfstream;
class MGAlphaFunc;
class MGDepthFunc;
class MGColor;
class MGBlendFunc;
class MGDepthMask;
class MGColorMask;
class MGFog;
class MGLightEnable;
class MGTranspMode;
class MGLight;
class MGRenderAttr;
class MGPolygonMode;
class MGLineStipple;
class MGLights;
class MGLineWidth;

/** @addtogroup GLAttrib
 *  @{
 */

///MGGLAttrib is an abstract class which defines the enum of undefined, disabled.
///Subclass of MGGLAttrib can use m_flag as a part of its own class attributes.
///In this case, -3 and -2 must be avoided. -3 and -2 do not appear in OpenGL attributes.
class MGCLASS MGGLAttrib:public MGAttrib{

public:

enum FLAG{
	UNDEFINED=-3,
	DISABLED=-2
};

enum ATTRIB_MASK{
	CURRENT_BIT= GL_CURRENT_BIT,
	POINT_BIT= GL_POINT_BIT,
	LINE_BIT= GL_LINE_BIT,
	POLYGON_BIT= GL_POLYGON_BIT,
	POLYGON_STIPPLE_BIT= GL_POLYGON_STIPPLE_BIT,
	PIXEL_MODE_BIT= GL_PIXEL_MODE_BIT,
	LIGHTING_BIT= GL_LIGHTING_BIT,
	FOG_BIT= GL_FOG_BIT,
	DEPTH_BUFFER_BIT= GL_DEPTH_BUFFER_BIT,
	ACCUM_BUFFER_BIT= GL_ACCUM_BUFFER_BIT,
	STENCIL_BUFFER_BIT= GL_STENCIL_BUFFER_BIT,
	VIEWPORT_BIT= GL_VIEWPORT_BIT,
	TRANSFORM_BIT= GL_TRANSFORM_BIT,
	ENABLE_BIT= GL_ENABLE_BIT,
	COLOR_BUFFER_BIT= GL_COLOR_BUFFER_BIT,
	HINT_BIT= GL_HINT_BIT,
	EVAL_BIT= GL_EVAL_BIT,
	LIST_BIT= GL_LIST_BIT,
	TEXTURE_BIT= GL_TEXTURE_BIT,
	SCISSOR_BIT= GL_SCISSOR_BIT
};

MGGLAttrib(int flag=UNDEFINED):m_flag(flag){;};

////////////Assignment operator////////////
virtual MGGLAttrib& operator=(const MGGLAttrib& gel2){set_glattrib(gel2);return *this;};

////////////Member Function////////////

///Generate a newed clone object.
virtual MGGLAttrib* clone()const=0;

bool undefined()const{return m_flag==UNDEFINED;};
bool defined()const{return m_flag!=UNDEFINED;};
bool disabled()const{return m_flag==DISABLED;};
bool enabled()const{return m_flag!=UNDEFINED && m_flag!=DISABLED;};

void set_undefined(){m_flag=UNDEFINED;};
void set_disabled(){m_flag=DISABLED;};

///retrieve the data.
int data()const{return m_flag;};
int& data(){return m_flag;};

///draw GLAttribute process.
virtual void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const=0;

///render GLAttribute process.
virtual void render()const=0;

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_draw_attrib_mask(unsigned int& mask)const=0;

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_draw_attrib_mask(unsigned int& mask)const=0;

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_render_attrib_mask(unsigned int& mask)const=0;

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_render_attrib_mask(unsigned int& mask)const=0;

/// Return This object's typeID
virtual long identify_type() const{return MGGLATTRIBUTE_TID;};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const;

///Compare if this and at2 are the same leaf MGGLAttrib class.
bool same_type(const MGGLAttrib& at2)const;

protected:

int m_flag;	///< =-3:undefined, will be inheritted.
			///< =-2:disabled.
			///< =other:each subclass's data(enabled).

///Assignment
MGGLAttrib& set_glattrib(const MGGLAttrib& gel2);

};

///Set or reset the bit of mask.
void set_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit);
void reset_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit);

///Construct a null newed MGAttrib from the type id TID.
MGGLAttrib* MGNullGLAttrib(long TID);

/** @} */ // end of GLAttrib group
#endif //#ifndef _MGGLAttrib_HH_
#endif //#ifdef MGCL_NO_MFC

