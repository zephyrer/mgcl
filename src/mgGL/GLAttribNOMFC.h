/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGLAttrib_HH_
#define _MGGLAttrib_HH_

#include "mg/MGCL.h"
#include "mg/Attrib.h"

//
//Define MGGLAttrib Class.

class MGOfstream;
class MGIfstream;

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

enum ATTRIB_MASK{};

MGGLAttrib(int flag=UNDEFINED):m_flag(flag){;};

////////////Assignment operator////////////
virtual MGGLAttrib& operator=(const MGGLAttrib& gel2){return *this;};

////////////Member Function////////////

///Generate a newed clone object.
virtual MGGLAttrib* clone()const{return new MGGLAttrib;};

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
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
	)const{;};

///render GLAttribute process.
virtual void render()const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_draw_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_draw_attrib_mask(unsigned int& mask)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
virtual void set_render_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
virtual void reset_render_attrib_mask(unsigned int& mask)const{;};

/// Return This object's typeID
virtual long identify_type() const{return MGGLATTRIBUTE_TID;};

///Read all member data.
virtual void ReadMembers(MGIfstream& buf){;};
///Write all member data
virtual void WriteMembers(MGOfstream& buf)const{;};

/// Output virtual function.
virtual std::ostream& out(std::ostream& ostrm) const{
	ostrm<<this;
	return ostrm;
};

///Compare if this and at2 are the same leaf MGGLAttrib class.
bool same_type(const MGGLAttrib& at2)const{return false;};

virtual std::string whoami()const{return "MGGLAttrib";};

protected:

int m_flag;	///< =-3:undefined, will be inheritted.
			///< =-2:disabled.
			///< =other:each subclass's data(enabled).

};

///Set or reset the bit of mask.
void set_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit);
void reset_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit);

///Construct a null newed MGAttrib from the type id TID.
MGGLAttrib* MGNullGLAttrib(long TID);

/** @} */ // end of GLAttrib group
#endif
