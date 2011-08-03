/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSYSGL_HH_
#define _MGSYSGL_HH_

#include <iosfwd>

class MGGel;
class MGOpenGLView;

/** @addtogroup DisplayHandling
 *  @{
 */

///mgSysGL is a class to provide a facility to draw temporal pictures.
///MGOpenGLView holds a list of mgSysGL and draws the pictures by invoking
///display list drawer.
///As long as the codes are unique, function codes can be
///any numbers. Usually the id is a command id.
class mgSysGL{

public:

///////////// mgSysGL /////////////
friend std::ostream& operator<< (std::ostream& outp, const mgSysGL& sysgl);

///////////////Constructor/////////////

mgSysGL():m_fucntion_id(0),m_gel(0){;};
mgSysGL(size_t fucntion_code,const MGGel* object_id)
:m_fucntion_id(fucntion_code),m_gel(const_cast<MGGel*>(object_id)){;};

///Copy constructor, replacing gel_old to gel_new.
mgSysGL(mgSysGL& glold,const MGGel* gel_old, const MGGel* gel_new);

virtual ~mgSysGL(){;};

////////////////Operator overload//////////////////

/////////////

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
virtual mgSysGL* clone()const;

///Return display list name of this.
virtual size_t dlist_name()const{return size_t(this);};

///Draw this Sysgl.
///This draw is used to draw the pictures for Undo(, Redo) operations.
///When this draw is invoked, all of the functions of mgGDL can be
///used in that context.
virtual void draw(MGOpenGLView& glv)const{;};

///Get the i-th ellement's object_id or fucntion_code.
size_t function_code()const{return m_fucntion_id;};

///Test if this mgSysGL includes gel(return true) or not.
///The default includes tests if the input gel is m_gel of this
///member data.
virtual bool includes(const MGGel* gel)const;

///Make system display list in glv.
///This must be a newed object and the ownership will be transfered to
///glv(glv.m_sysgllist).
void make_display_list(
	MGOpenGLView& glv
);

/// Output virtual function.
///Output to stream file:メンバデータを標準出力に出力する。
virtual std::ostream& out(std::ostream& ostrm) const;

///replace gel_old to gel_new.
///If gel_old is not included in this, do nothing.
virtual void replace(const MGGel* gel_old, const MGGel* gel_new);

void set_function_code(size_t fc){m_fucntion_id=fc;};

protected:
MGGel* object_id()const{return m_gel;};
void set_object_id(MGGel* oi){m_gel=oi;};

private:

	size_t m_fucntion_id;///<fucntion code.
	MGGel* m_gel;		///<Object id. When more than 1 objects are concerned, 
						///<mgSysGL class will be inheritted and the subclass
						///<will retain the objcets.
};

/** @} */ // end of DisplayHandling group
#endif
