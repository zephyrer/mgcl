/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGAttribedGel_HH_
#define _MGAttribedGel_HH_

#include "mg/MGCL.h"
#include "mg/Gel.h"

//
//Define MGAttribedGel Class.

class MGAppearance;
class MGObject;
class MGGroup;
class MGGLAttrib;

/** @addtogroup GelRelated
 *  @{
 */

///MGAttribedGel is an abstract class which provides function interfaces of
///MGGroup and MGObject that have MGAppearance. 
///*******MGAttribedGel does not have any member data.************
class MGCLASS MGAttribedGel:public MGGel{

public:

///Assignment.
///When the leaf objects of this and gel2 are not equal, this assignment
///does nothing.
virtual MGAttribedGel& operator=(const MGAttribedGel& gel2);

///Get the MGAppearance pointer in this group. If not defined, null will be
///returned.
virtual MGAppearance* appearance()=0;
virtual const MGAppearance* appearance()const=0;

///Process of draw or render attributes.
virtual void draw_attribute(
	bool no_color=false	///if true, color attribute will be neglected.
)const;
virtual void render_attribute()const;

///Make sure that this MGAttribedGel has appearance,
/// and get the MGAppearance pointer.
virtual MGAppearance* ensure_appearance()=0;

///Obtain attribute mask for glPushAttrib().
virtual size_t get_draw_attrib_mask()const;
virtual size_t get_render_attrib_mask()const;

///Test if this group is attributed  as no display.
///true if attributed as no display.
virtual bool no_display()const;

///Remove the MGAppearance of this MGAttribedGel.
virtual void remove_appearance()=0;

///Removed the attribute of specified type.
void remove_GLattrib(long tid);

///Set the attribute in this list. attr must be a newed object, and the
///ownership will be transfered to this MGAppearance.
virtual void set_GLattrib(MGGLAttrib* attr);

///Set this group as display or no display group.
virtual void set_display();
virtual void set_no_display();
bool visible()const{return !no_display();};

private:

friend class MGGroup;
friend class MGObject;

};

/** @} */ // end of GelRelated group
#endif
