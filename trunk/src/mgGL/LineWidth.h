/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifdef MGCL_NO_MFC
#include "mgGL/LineWidthNOMFC.h"
#else

#ifndef _MGLineWidth_HH_
#define _MGLineWidth_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//
//Define MGLineWidth Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGLineWidth defines line width.
class MGCLASS MGLineWidth:public MGGLAttrib{

public:

enum MODE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	ENABLED=1
};

MGLineWidth(MODE m=UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};
MGLineWidth(float width):MGGLAttrib(ENABLED),m_line_width(width){;};

///Assignment
MGLineWidth& operator=(const MGGel& gel2);
MGLineWidth& operator=(const MGLineWidth& gel2);

///comparison
bool operator<(const MGLineWidth& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGLineWidth* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_width(float width);
float get_width()const{return m_line_width;};

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{exec();};

///get maximum line width
float get_maximum_width()const;

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{set_Amask(mask,LINE_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LINE_BIT);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{set_Amask(mask,LINE_BIT);};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{reset_Amask(mask,LINE_BIT);};

/// Return This object's typeID
long identify_type() const{return MGLINE_WIDTH_TID;};

std::string whoami()const{return "LineWidth";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

private:

	float m_line_width;	///<line width.

};

/** @} */ // end of GLAttrib group
#endif //#ifndef _MGLineWidth_HH_
#endif //#ifdef MGCL_NO_MFC
