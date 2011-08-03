/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifdef MGCL_NO_MFC
#include "mgGL/AppearanceNOMFC.h"
#else

#ifndef _MGAppearance_HH_
#define _MGAppearance_HH_

#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/Group.h"
#include "mgGL/RenderAttr.h"
#include "mgGL/AlphaFunc.h"
#include "mgGL/BlendFunc.h"
#include "mgGL/DepthFunc.h"
#include "mgGL/PolygonMode.h"

//
//Define MGAppearance Class.
class MGOfstream;
class MGIfstream;
class MGGLAttrib;
class MGFog;
class MGLight;

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

MGAppearance& operator=(const MGGel& gel2);
MGAppearance& operator=(const MGAppearance& gel2);

///comparison
bool operator<(const MGAppearance& gel2)const;
bool operator<(const MGGel& gel2)const;

////////Member Function////////

std::ostream& out(
	std::ostream& ostrm
)const;

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
size_t add_light(MGLight* light);

///Test if this MGAppearance can be removed or not.
bool can_be_removed()const;

///Generate copied gel of this gel.
///Returned is a newed object. User must delete the object.
MGAppearance* clone()const;

///draw GLAttributes process.
void drawAttrib(
	bool no_color=false	///if true, color attribute will be neglected.
)const;

///Erase the specified attribute.
///Function's return value is the iterator after the erased data.
iterator erase(iterator i){return MGGroup::erase(i);};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
size_t get_draw_attrib_mask()const;
size_t get_render_attrib_mask()const;

/// Return This object's typeID
long identify_type() const{return MGAPPEARANCE_TID;};

///Test this is no_display MGAppearance.
bool no_display()const{return m_no_display;};

///Release the attribute of specified type.
///Function's return value is the MGGLAttrib* that is released.
MGGLAttrib* release_attrib(long tid);

///render GLAttributes process.
void render()const;

///Search the same type MGGLAttrib leaf class object in this list.
///Function's return value is the iterator found.
///If not found, end() will be returned.
iterator search(const MGGLAttrib* atr);
iterator search_by_id(MGGEL_TID tid);
const_iterator search(const MGGLAttrib* atr)const;
const_iterator search_by_id(MGGEL_TID tid)const;

///Set the attribute in this list. attr must be a newed object, and the
///ownership will be transfered to this MGAppearance.
void set_attrib(MGGLAttrib* attr);
void set_attrib(MGGLAttribs& attrs);

///Set the attribute in this list. attr must be a newed object, and the
///ownership will be transfered to this MGAppearance.
///When the appearance held an attribute, the old one will be returned
///as the function's return value. Users must delete it.
MGGLAttrib* set_attrib_with_old(MGGLAttrib* attr);

void set_color(const MGColor& color);
void set_color(const float color[4]);
void set_color(float red, float green, float blue, float alpha=1.);

///Set the fog data. fog must be a newed object and the ownership will be transfered
///to this Object.
void set_fog(MGFog* fog);

void set_display(){	m_no_display=false;};
void set_no_display(){	m_no_display=true;}

void set_light_disabled();
void set_light_enabled();

///Set the material. When rs=FRONT_AND_BACK and different material for the back side
///is used, set_back_material must be invoked after invoking set_material.
///Else the same material will be appllied for the both sides.
void set_material(
	MGRenderAttr::RENDERSIDE rs,
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess=0.,
	float transparency=0.
);

///Set the back side material. Invoking set_back_material means two sided material
///and setting different material to the back side.
///Before use of set_back_material, set_material must be invoked first.
///set_back_material will set two sided material.
void set_back_material(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess=0.,
	float transparency=0.
);

void set_shade_model(MGShadeModel::MODEL model);

void set_texture(
	MGTexture* texture	///<texture must be newed one.
		///<The ownership will be transfered to this MGAppearance.
		///<Texture will be enabled.
);

void set_transparency_disabled();
void set_transparency_mode(MGTranspMode::MODE tr=MGTranspMode::BLEND);
void set_alpha_func_disabled();
void setAlphaFunc(MGAlphaFunc::FUNC func,float ref);
void set_blending_func_disabled();

///By invoking setBlendFunc, Blending function will be enabled.
void setBlendFunc(
	MGBlendFunc::SrcBlendFunc srcBlendFunc=MGBlendFunc::ONE_SBLEND,
	MGBlendFunc::DstBlendFunc dstBlendFunc=MGBlendFunc::ZERO_DBLEND
);

void set_depth_func_disabled();
void setDepthFunc(MGDepthFunc::FUNC depthFunc=MGDepthFunc::LESS);///Depth test will be enabled.
void set_depth_mask_enabled();
void set_depth_mask_disabled();
void setColorMask(bool R, bool G, bool B, bool A);
void setPolyMode(MGPolygonMode::MODE polyMode);
void setLineStipple(unsigned factor,unsigned short pattern);
void setLineWidth(float width);

virtual std::string whoami()const{return "Appearance";};

protected:

	bool m_no_display;///True if not to display, false if to display.

///メンバデータを読み出す関数
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
void WriteMembers(MGOfstream& buf) const;

private:
iterator insert(iterator i, MGAttrib* atr){return MGGroup::insert(i,atr);};
void push_back(MGAttrib* atr){MGGroup::push_back(atr);};
void push_front(MGAttrib* atr){MGGroup::push_front(atr);};

};

/** @} */ // end of DisplayHandling group
#endif //#ifndef _MGAppearance_HH_
#endif //#ifdef MGCL_NO_MFC
