/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGColor_HH_
#define _MGColor_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;
class MGColors;

//
//Define MGColor Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGColor defines the GL color.
class MGCLASS MGColor:public MGGLAttrib{

public:

enum MODE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	ENABLED=1,
};

/// Color enumuration.
/// This is completely compatible to GDI+ color enumuration.

///The color id from 1 to 8(black to white) is the id of IGES.
enum ColorID{
	dummy                = 0,
    Black                = 1,
    Red                  = 2,
    Green                = 3,
    Blue                 = 4,
    Yellow               = 5,
    Magenta              = 6,
    Cyan                 = 7,
    White                = 8,
	endID
};

MGColor(MODE m=UNDEFINED):MGGLAttrib(static_cast<int>(m)){;};

/// Construct a color from float values.
MGColor(float red, float green, float blue, float alpha=1.){;};

/// Construct a color from ARGB value. In argb, each 8bits of (a,r,g,b) is the value of
///the range 0 to 255 .
MGColor(unsigned int argb){;};

///Assignment
MGColor& operator=(const MGGel& gel2){return *this;};
MGColor& operator=(const MGColor& gel2){return *this;};

///Scaling the color values by the factor scale.
///All of the elements except trancparency element will be multiplied
///by the scale and be clamped between 0. and 1.
MGColor& operator*=(float scale){return *this;};

///Add a color values to RGB data..
///All of the elements except trancparency element will be added by value
///and be clamped between 0. and 1.
MGColor& operator+=(float value){return *this;};

///comparison
bool operator<(const MGColor& gel2)const{return true;};
bool operator<(const MGGel& gel2)const{return true;};

////////Member Function////////

/// helper function
void argb_to_float(unsigned int argb, float out[4]) const{;};

///Generate a newed clone object.
MGColor* clone()const{return new MGColor();};

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const{;};

float* color(){return m_color;};
const float* color()const{return m_color;};

///get the color id of GDI+ color enumeration.
///If not found, 0 will be returned.
///maxID is the max id of ColorID enumeration.
int get_ColorID(int maxID)const{return 0;};
///Get the color instance reference by the color id.
static const MGColor& get_instance(ColorID id){return *(MGColor*)0;};

///Get unsigned integer valur of a color id.
///Function's return value is (A, R, G, B) data of GDI+.
static size_t get_ARGBinstance(MGColor::ColorID id){return 0;};

void set_color(const float color[4]){;};
void set_color(float red, float green, float blue, float alpha=1.){;};
void get_color(float color[4])const{;};
void get_color(float& red, float& green, float& blue, float& alpha)const{;};

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{;};

///render GLAttribute process.
void render()const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const{;};

/// Return This object's typeID
long identify_type() const{return MGCOLOR_TID;};

std::string whoami()const{return "Color";};

///Read all member data.
void ReadMembers(MGIfstream& buf){;};
///Write all member data
void WriteMembers(MGOfstream& buf)const{;};

/// Output function.
std::ostream& out(std::ostream& ostrm) const{
	ostrm<<"Color="; MGGLAttrib::out(ostrm);return ostrm;
};

///Output to IGES stream file(Color=PD314).
///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

private:

	float m_color[4];///<color data.
	static MGColors m_colors;///<color instance;

	friend class MGColors;

};

///@cond

///MGColors defines standard MGColor array that can be accessed through
///ColorID.
class MGColors{
public:
	const MGColor& color(size_t i){return *(m_colors[i]);};
	~MGColors();
private:
	MGColor** m_colors;
	MGColors();
	friend class MGColor;
};

namespace Gdiplus{
	typedef long GdiplusStartupInput;
	int GdiplusStartup( unsigned long*, GdiplusStartupInput*, unsigned* );
	int GdiplusShutdown( unsigned long);
};

///@endcond

/** @} */ // end of GLAttrib group

#endif
