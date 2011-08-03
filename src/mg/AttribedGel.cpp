/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/AttribedGel.h"
#include "mgGL/Appearance.h"
#include "mgGL/GLAttrib.h"

//
//Define MGAttribedGel Class.
//MGAttribedGel is an abstract class which provides function interfaces of
//MGGroup and MGObject that have MGAppearance. 
//*******MGAttribedGel does not have any member data.************
//

MGAttribedGel& MGAttribedGel::operator=(const MGAttribedGel& gel2){
	MGGel::operator=(gel2);
	return *this;
}

//Process of draw or render attributes.
void MGAttribedGel::draw_attribute(
	bool no_color	//if true, color attribute will be neglected.
)const{
	const MGAppearance* app=appearance();
	if(app)
		app->drawAttrib(no_color);
}
void MGAttribedGel::render_attribute()const{
	const MGAppearance* app=appearance();
	if(app)
		app->render();
}

//Obtain attribute mask for glPushAttrib().
size_t MGAttribedGel::get_draw_attrib_mask()const{
	size_t mask=0;
	const MGAppearance* app=appearance();
	if(app)
		mask=app->get_draw_attrib_mask();
	return mask;
}
size_t MGAttribedGel::get_render_attrib_mask()const{
	size_t mask=0;
	const MGAppearance* app=appearance();
	if(app)
		mask=app->get_render_attrib_mask();
	return mask;
}

//Test if this group is attributed  as no display.
//true if attributed as no display.
bool MGAttribedGel::no_display()const{
	const MGAppearance* app=appearance();
	if(app)
		return app->no_display();
	return false;
}

//Removed the attribute of specified type.
void MGAttribedGel::remove_GLattrib(long tid){
	MGAppearance* app=appearance();
	if(app){
		delete app->release_attrib(tid);
		if(app->can_be_removed())
			remove_appearance();
	}
}

//Set the attribute in this list. attr must be a newed object, and the
//ownership will be transfered to this MGAppearance.
void MGAttribedGel::set_GLattrib(MGGLAttrib* attr){
	if(!attr)
		return;
	MGAppearance* app=ensure_appearance();
	app->set_attrib(attr);
}

//Set this group as display or no display group.
void MGAttribedGel::set_display(){
	MGAppearance* app=appearance();
	if(!app)
		return;
	app->set_display();
}
void MGAttribedGel::set_no_display(){
	MGAppearance* app=ensure_appearance();
	app->set_no_display();
}
