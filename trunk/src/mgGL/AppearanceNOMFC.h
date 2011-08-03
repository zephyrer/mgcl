/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGAppearance_HH_
#define _MGAppearance_HH_

#include "mg/MGCL.h"
#include "mg/Gel.h"
#include "mg/Attrib.h"
#include "mg/Group.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"

//
//Define MGAppearance Class.
class MGOfstream;
class MGIfstream;
class MGGLAttrib;
class MGFog;
class MGLight;
class MGColor;

/** @addtogroup DisplayHandling
 *  @{
 */

///MGAppearance is a class to contain MGGLAttrib objects.
///MGAppearance acts just like as std::auto_ptr.
///That is, MGAppearance holds newed object pointers of MGGLAttrib,
///and when copy constructor or assignment operator is invoked,
///the pointer ownership is transfered to the new MGAppearance object.
///A list of newed MGGLAttrib object pointer will be stored in parent MGGroup.
///No two same leaf type MGGLAttrib objects are included in this list.
class MGCLASS MGAppearance: public MGGroup{

public:

typedef MGGroup::iterator iterator;
typedef MGGroup::const_iterator const_iterator;

MGAppearance():m_no_display(false){;};

MGAppearance& operator=(const MGGel& gel2){return *this;};
MGAppearance& operator=(const MGAppearance& gel2){return *this;};

///comparison
bool operator<(const MGAppearance& gel2)const{return true;};
bool operator<(const MGGel& gel2)const{return true;};

////////Member Function////////

std::ostream& out(
	std::ostream& ostrm
)const{
	ostrm<<"MGAppearance="<<this<<",no_display="<<0<<",";
	MGGroup::out(ostrm);
	return ostrm;
};

MGAttrib* back(){return static_cast<MGAttrib*>(MGGroup::back());};
const MGAttrib* back()const{return static_cast<const MGAttrib*>(MGGroup::back());};
iterator begin(){return MGGroup::begin();};
const_iterator begin()const{return MGGroup::begin();};
void clear(){MGGroup::clear();};
bool empty()const{return MGGroup::empty();};
iterator end(){return MGGroup::end();};
const_iterator end()const{return MGGroup::end();};
MGAttrib* front(){return static_cast<MGAttrib*>(MGGroup::front());};
const MGAttrib* front()const{return static_cast<const MGAttrib*>(MGGroup::front());};
void pop_back(){MGGroup::pop_back();};
void pop_front(){MGGroup::pop_front();};
size_t size()const{return MGGroup::size();};

///Add a light. light must be a newed object, and the ownership will be
///transfered to this object.
///Function's return value is the number of lights after added.
size_t add_light(MGLight* light){return 0;};

///Test if this MGAppearance can be removed or not.
bool can_be_removed()const{return true;};

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
MGAppearance* clone()const{return new MGAppearance;};

///draw GLAttributes process.
void drawAttrib(
	bool no_color=false	///if true, color attribute will be neglected.
)const{;};

///Erase the specified attribute.
///Function's return value is the iterator after the erased data.
iterator erase(iterator i){return MGGroup::erase(i);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
size_t get_draw_attrib_mask()const{return 0;};
size_t get_render_attrib_mask()const{return 0;};

/// Return This object's typeID
long identify_type() const{return MGAPPEARANCE_TID;};

///Test this is no_display MGAppearance.
bool no_display()const{return m_no_display;};

///Release the attribute of specified type.
///Function's return value is the MGGLAttrib* that is released.
MGGLAttrib* release_attrib(long tid){return 0;};

///render GLAttributes process.
void render()const{;};

///Search the same type MGGLAttrib leaf class object in this list.
///Function's return value is the iterator found.
///If not found, end() will be returned.
iterator search(const MGGLAttrib* atr){return end();};
iterator search_by_id(MGGEL_TID tid){return end();};
const_iterator search(const MGGLAttrib* atr)const{return end();};
const_iterator search_by_id(MGGEL_TID tid)const{return end();};

///Set the attribute in this list. attr must be a newed object, and the
///ownership will be transfered to this MGAppearance.
void set_attrib(MGGLAttrib* attr){;};
void set_attrib(MGGLAttribs& attrs){;};

///Set the attribute in this list. attr must be a newed object, and the
///ownership will be transfered to this MGAppearance.
///When the appearance held an attribute, the old one will be returned
///as the function's return value. Users must delete it.
MGGLAttrib* set_attrib_with_old(MGGLAttrib* attr){;};

void set_color(const MGColor& color){;};
void set_color(const float color[4]){;};
void set_color(float red, float green, float blue, float alpha=1.){;};

void set_display(){	m_no_display=false;};
void set_no_display(){	m_no_display=true;}
virtual std::string whoami()const{return "Appearance";};

protected:

	bool m_no_display;///True if not to display, false if to display.

///メンバデータを読み出す関数
void ReadMembers(MGIfstream& buf){
	MGGroup::ReadMembers(buf);
	int no_disp;
	buf>>no_disp;
	m_no_display=false; if(no_disp) m_no_display=true;
};

///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const{
	MGGroup::WriteMembers(buf);
	int no_disp=1;if(!no_display()) no_disp=0;
	buf<<no_disp;
};

};

/** @} */ // end of DisplayHandling group
#endif
