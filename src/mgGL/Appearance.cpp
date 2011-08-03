/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Color.h"
#include "mgGL/Appearance.h"
#include "mgGL/Fog.h"
#include "mgGL/Light.h"
#include "mgGL/LightEnable.h"
#include "mgGL/DepthMask.h"
#include "mgGL/ColorMask.h"
#include "mgGL/LineStipple.h"
#include "mgGL/LineWidth.h"
#include "mgGL/Lights.h"

//
//Define MGAppearance Class.
//MGAppearance is a class to contain MGGLAttrib objects.
//MGAppearance acts just like as std::auto_ptr.
//That is, MGAppearance holds newed object pointers of MGGLAttrib,
//and when copy constructor or assignment operator is invoked,
//the pointer ownership is transfered to the new MGAppearance object.

MGAppearance& MGAppearance::operator=(const MGAppearance& gel2){
	if(this==&gel2)
		return *this;

	MGGroup::operator=(gel2);
	m_no_display=gel2.m_no_display;
	return *this;
}
MGAppearance& MGAppearance::operator=(const MGGel& gel2){
	const MGAppearance* gel2_is_this=dynamic_cast<const MGAppearance*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGAppearance::operator<(const MGAppearance& gel2)const{
	size_t n1=size(), n2=gel2.size();
	if(n1==n2){
		return size_t(this)<size_t(&gel2);
	}else
		return n1<n2;
}
bool MGAppearance::operator<(const MGGel& gel2)const{
	const MGAppearance* gel2_is_this=dynamic_cast<const MGAppearance*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//Add a light. light must be a newed object, and the ownership will be
//transfered to this object.
//Function's return value is the number of lights after added.
size_t MGAppearance::add_light(MGLight* light){
	iterator i=search_by_id(MGLIGHTS_TID);
	MGLights* lights;
	if(i==end()){
		lights=new MGLights;
		push_back(lights);
	}else
		lights=static_cast<MGLights*>(*i);
	lights->push_back(light);
	return lights->size();
}

//Test if this MGAppearance can be removed or not.
bool MGAppearance::can_be_removed()const{
	return (size()==0 && !no_display());
}

//Generate copied gel of this gel.
//Returned is a newed object. User must delete the object.
MGAppearance* MGAppearance::clone()const{
	MGAppearance* gel=new MGAppearance;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		gel->m_gels.push_back((**i).clone());
	}
	gel->m_no_display=m_no_display;
	return gel;
}

std::ostream& MGAppearance::out(
	std::ostream& ostrm
)const{
	ostrm<<"MGAppearance="<<this<<",no_display="<<m_no_display;
	ostrm<<", number of attribs="<<size()<<std::endl;
	const_iterator i=begin(), ie=end();	
	for(size_t j=0; i!=ie; i++, j++){
		ostrm<<"  attr["<<j<<"]::"<<(**i)<<std::endl;
	}
	return ostrm;
}

//////////Member Function//////////

//draw GLAttributes process.
void MGAppearance::drawAttrib(
	bool no_color	//if true, color attribute will be neglected.
)const{
	if(no_display())
		return;

	const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		if(no_color){
			const MGColor* colr=dynamic_cast<const MGColor*>(*i);
			if(colr)
				continue;
		}
		const MGGLAttrib* gla=dynamic_cast<const MGGLAttrib*>(*i);
		if(gla)
			gla->drawAttrib();
	}
}

//Release the specified attribute.
//Function's return value is the MGGLAttrib* that is released.
MGGLAttrib* MGAppearance::release_attrib(long tid){
	iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		if((*i)->identify_type()==tid){
			MGGLAttrib* gla=static_cast<MGGLAttrib*>(*i);
			MGGroup::release(i);
			return gla;
		}
	}
	return 0;
}

//render GLAttributes process.
void MGAppearance::render()const{
	if(no_display())
		return;

	const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=dynamic_cast<const MGGLAttrib*>(*i);
		if(gla)
			gla->render();
	}
}

//Set the attribute in this list. attr must be a newed object, and the
//ownership will be transfered to this MGAppearance.
void MGAppearance::set_attrib(MGGLAttrib* attr){
	MGGLAttrib* olda=set_attrib_with_old(attr);
	if(olda)
		delete olda;
}
void MGAppearance::set_attrib(MGGLAttribs& attrs){
	size_t n=attrs.size();
	for(size_t i=0; i<n; i++){
		MGGLAttrib* attr=attrs.release(i);
		MGGLAttrib* olda=set_attrib_with_old(attr);
		if(olda)
			delete olda;
	}
}

//Set the attribute in this list. attr must be a newed object, and the
//ownership will be transfered to this MGAppearance.
//When the appearance held an attribute, the old one will be returned
//as the function's return value. Users must delete it.
MGGLAttrib* MGAppearance::set_attrib_with_old(MGGLAttrib* attr){
	if(!attr)
		return 0;

	iterator i=search(attr);
	if(i!=end()){//If found.
		MGGLAttrib* gla=static_cast<MGGLAttrib*>(*i);
		iterator j=MGGroup::release(i); insert(j,attr);
		return gla;
	}else{//If not found.
		push_back(attr);
		return 0;
	}
}

//Functional object for find_if.
class MGAppearanceSearch{
public:
	MGAppearanceSearch(const MGGLAttrib* atr):m_attr(atr){;};
	bool operator()(const MGGel* atr2){
		return m_attr->same_type(*(static_cast<const MGGLAttrib*>(atr2)));
	};
	const MGGLAttrib* m_attr;
};
//Search the same MGGLAttrib leaf class object in this list.
//If not found, end() will be returned.
MGAppearance::iterator MGAppearance::search(const MGGLAttrib* atr){
	return std::find_if(begin(), end(), MGAppearanceSearch(atr));
}
MGAppearance::const_iterator MGAppearance::search(const MGGLAttrib* atr)const{
	return std::find_if(begin(), end(), MGAppearanceSearch(atr));
}
MGAppearance::iterator MGAppearance::search_by_id(MGGEL_TID tid){
	MGAppearance::iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(*i);
		if(gla->identify_type()==tid) return i;
	}
	return ie;
}
MGAppearance::const_iterator MGAppearance::search_by_id(MGGEL_TID tid)const{
	MGAppearance::const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(*i);
		if(gla->identify_type()==tid) return i;
	}
	return ie;
}
//Turn on the appropriate mask bit for this attribute. See glPushAttrib().
size_t MGAppearance::get_draw_attrib_mask()const{
	size_t mask=0;
	MGAppearance::const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(*i);
		gla->set_draw_attrib_mask(mask);
	}
	return mask;
}

//Turn on the appropriate mask bit for this attribute. See glPushAttrib().
size_t MGAppearance::get_render_attrib_mask()const{
	size_t mask=0;
	MGAppearance::const_iterator i=begin(), ie=end();
	for(;i!=ie; i++){
		const MGGLAttrib* gla=static_cast<const MGGLAttrib*>(*i);
		gla->set_render_attrib_mask(mask);
	}
	return mask;
}

//メンバデータを書き込む関数
void MGAppearance::WriteMembers(MGOfstream& buf)const{
	MGGroup::WriteMembers(buf);
	int no_disp=1;if(!no_display()) no_disp=0;
	buf<<no_disp;
}

//メンバデータを読み出す関数
void MGAppearance::ReadMembers(MGIfstream& buf){
	MGGroup::ReadMembers(buf);
	int no_disp;
	buf>>no_disp;
	m_no_display=false; if(no_disp) m_no_display=true;
}

void MGAppearance::set_color(const MGColor& color){
	MGColor* colr=new MGColor(color);
	set_attrib(colr);
}
void MGAppearance::set_color(const float color[4]){
	MGColor* colr=new MGColor(color[0],color[1],color[2],color[3]);
	set_attrib(colr);
}
void MGAppearance::set_color(float red, float green, float blue, float alpha){
	MGColor* colr=new MGColor(red,green,blue,alpha);
	set_attrib(colr);
}

void MGAppearance::set_light_disabled(){
	MGLightEnable* ld=new MGLightEnable(MGLightEnable::DISABLED);
	set_attrib(ld);
}
void MGAppearance::set_light_enabled(){
	MGLightEnable* le=new MGLightEnable(MGLightEnable::ENABLED);
	set_attrib(le);
}

//Set the material. When rs=FRONT_AND_BACK and different material for the back side
//is used, set_back_material must be invoked after invoking set_material.
//Else the same material will be appllied for the both sides.
void MGAppearance::set_material(
	MGRenderAttr::RENDERSIDE rs,
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	MGRenderAttr* ra;
	iterator i=search_by_id(MGRENDER_ATTR_TID);
	if(i==end()){
		ra=new MGRenderAttr();
		push_back(ra);
	}else ra=static_cast<MGRenderAttr*>(*i);
	ra->set_material(rs,ambient,diffuse,specular,emission,shininess,transparency);
}

//Set the back side material. Invoking set_back_material means two sided material
//and setting different material to the back side.
//Before use of set_back_material, set_material must be invoked first.
//set_back_material will set two sided material.
void MGAppearance::set_back_material(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess,
	float transparency
){
	MGRenderAttr* ra;
	iterator i=search_by_id(MGRENDER_ATTR_TID);
	if(i==end()){
		ra=new MGRenderAttr();
		push_back(ra);
	}else ra=static_cast<MGRenderAttr*>(*i);
	ra->set_back_material(
		ambient,diffuse,specular,emission,shininess,transparency);
}

//Set the fog data. fog must be a newed object and the ownership will be transfered
//to this Object.
void MGAppearance::set_fog(MGFog* fog){
	set_attrib(fog);
}

void MGAppearance::set_shade_model(MGShadeModel::MODEL model){
	MGRenderAttr* ra;
	iterator i=search_by_id(MGRENDER_ATTR_TID);
	if(i==end()){
		ra=new MGRenderAttr();
		push_back(ra);
	}else ra=static_cast<MGRenderAttr*>(*i);
	set_light_enabled();
	ra->set_shade_model(model);
}

void MGAppearance::set_texture(
	MGTexture* texture	//texture must be newed one.
		//The ownership will be transfered to this MGAppearance.
)		//Texture will be enabled.
{
	MGRenderAttr* ra;
	iterator i=search_by_id(MGRENDER_ATTR_TID);
	if(i==end()){
		ra=new MGRenderAttr();
		push_back(ra);
	}else ra=static_cast<MGRenderAttr*>(*i);
	//***************UNDONE*****************//
}

void MGAppearance::set_transparency_disabled(){
	MGTranspMode* td=new MGTranspMode(MGTranspMode::DISABLED);
	set_attrib(td);
}
void MGAppearance::set_transparency_mode(MGTranspMode::MODE tr){
	MGTranspMode* tm=new MGTranspMode(tr);
	set_attrib(tm);
}

void MGAppearance::set_alpha_func_disabled(){
	MGAlphaFunc* af=new MGAlphaFunc(MGAlphaFunc::DISABLED);
	set_attrib(af);
}
void MGAppearance::setAlphaFunc(MGAlphaFunc::FUNC func,float ref){
	MGAlphaFunc* af=new MGAlphaFunc(func,ref);
	set_attrib(af);
}

void MGAppearance::set_blending_func_disabled(){
	MGBlendFunc* bf=new MGBlendFunc(MGBlendFunc::DISABLED);
	set_attrib(bf);
}

//By invoking setBlendFunc, Blending function will be enabled.
void MGAppearance::setBlendFunc(
	MGBlendFunc::SrcBlendFunc srcBlendFunc,
	MGBlendFunc::DstBlendFunc dstBlendFunc
){
	MGBlendFunc* bf=new MGBlendFunc(srcBlendFunc,dstBlendFunc);
	set_attrib(bf);
}

void MGAppearance::set_depth_func_disabled(){
	MGDepthFunc* df=new MGDepthFunc(MGDepthFunc::DISABLED);
	set_attrib(df);
}

void MGAppearance::setDepthFunc(MGDepthFunc::FUNC depthFunc)//Depth test will be enabled.
{
	MGDepthFunc* df=new MGDepthFunc(depthFunc);
	set_attrib(df);
}

void MGAppearance::set_depth_mask_enabled(){
	MGDepthMask* dm=new MGDepthMask(MGDepthMask::ENABLED);
	set_attrib(dm);
}

void MGAppearance::set_depth_mask_disabled(){
	MGDepthMask* dm=new MGDepthMask(MGDepthMask::DISABLED);
	set_attrib(dm);
}

void MGAppearance::setColorMask(bool R, bool G, bool B, bool A){
	MGColorMask* colrm=new MGColorMask(R,G,B,A);
	set_attrib(colrm);
}

void MGAppearance::setPolyMode(MGPolygonMode::MODE polyMode){
	MGPolygonMode* pm=new MGPolygonMode(polyMode);
	set_attrib(pm);
}

void MGAppearance::setLineStipple(unsigned factor,unsigned short pattern){
	MGLineStipple* ls=new MGLineStipple(factor,pattern);
	set_attrib(ls);
}

void MGAppearance::setLineWidth(float width){
	MGLineWidth* lw=new MGLineWidth(width);
	set_attrib(lw);
}
