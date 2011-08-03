/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPolygonMode_HH_
#define _MGPolygonMode_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGPolygonMode Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGPolygonMode defines polygon mode.
class MGCLASS MGPolygonMode:public MGGLAttrib{

public:

enum MODE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	POINT=GL_POINT,
	LINE=GL_LINE,
	FILL=GL_FILL
};

MGPolygonMode(MODE m=UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};

///Assignment
MGPolygonMode& operator=(const MGGel& gel2);
MGPolygonMode& operator=(const MGPolygonMode& gel2);

///comparison
bool operator<(const MGPolygonMode& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGPolygonMode* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_mode(MODE m=FILL){data()=m;};
MODE get_mode()const{return static_cast<MODE>(data());};
GLenum GLmode()const{return static_cast<GLenum>(data());};

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec();};

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,POLYGON_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{reset_Amask(mask,POLYGON_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,POLYGON_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,POLYGON_BIT);};

/// Return This object's typeID
long identify_type() const{return MGPOLYGON_MODE_TID;};

std::string whoami()const{return "PolygonMode";};

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
