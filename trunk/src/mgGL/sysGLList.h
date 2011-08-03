/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSYSGLList_HH_
#define _MGSYSGLList_HH_

#include <vector>
#include <utility>
#include <bitset>
#include "mg/MGCL.h"
#include "mg/Pvector.h"
#include "mgGL/sysGL.h"
#include "mgGL/SysZebraGL.h"

class MGGel;
class MGCurve;

/** @addtogroup DisplayHandling
 *  @{
 */

#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<mgSysGL*>;
#pragma warning( pop )
#else
#include <list>
#endif

///mgSysGLList is a class to constrol system generated display list.
///System generated display list means the display lists not generated
///from the document, m_group of the document.
///mgSysGLList returns a display list name when an mgSysGL(fucntion code, object id)
///is added(push_back, or push_front). Using the display list name, invoke
///	glNewList(name, GL_COMPILE), ..., glEndList().
///Then mgSysGLList will manage to delete the display list,
/// giving the fucntion code or the object id.
class MGCLASS mgSysGLList{

public:

#if defined(MGCL_DLL)
	typedef MGListProxy<mgSysGL*> container_type;
	typedef MGListProxy<mgSysZebraGL*> zebra_container_type;
#else
	typedef std::list<mgSysGL*> container_type;
	typedef std::list<mgSysZebraGL*> zebra_container_type;
#endif

///In m_sysgls mgSysGL pointers are stored.
container_type m_sysgls;			///<list of mgSysGL.
zebra_container_type m_sysgls_zebra;///<list of mgSysZebraGL.

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;
typedef zebra_container_type::iterator zebra_iterator;
typedef zebra_container_type::const_iterator const_zebra_iterator;

///////////////Constructor/////////////

mgSysGLList();

///Copy constructor.
mgSysGLList(const mgSysGLList& list2);

////////////Destructor////////////////
~mgSysGLList();

////////////////Operator overload//////////////////

///Assignment operator.
mgSysGLList& operator=(const mgSysGLList& list2);

/// オペレーション

///const_iterator begin()const{return m_sysgls.begin();};
///iterator begin(){return m_sysgls.begin();};

///Clear the list. Delete all the display lists.
///make_RC_current() of MGOpenGLView is necessary to use.
void clear();

///Delete all the display lists that have the fucntion_code fc.
///make_RC_current() of MGOpenGLView is necessary to use.
///function's retrurn value is true if any one is deleted.
bool delete_lists_by_function_code(size_t fc);

///Delete all the display lists that have the object_id oi.
///make_RC_current() of MGOpenGLView is necessary to use.
///function's retrurn value is true if any one is deleted.
bool delete_lists_by_object_id(
	const MGGel* oi,
	MGPvector<mgSysGL>& functions	///</mgSysGL pointer will be appended.
);

///Delete the display list that have the fucntion_code fc and the object id gel.
///make_RC_current() of MGOpenGLView is necessary to use.
///function's retrurn value is true if any one is deleted.
bool delete_lists_by_function_object_code(size_t fc, const MGGel* gel);

///Draw all the objects by calling glCallList in this list.
///make_RC_current() of MGOpenGLView is necessary to use.
void draw_list()const;

///Test if this list includes the fucntion code fc's SysGL or not.
bool includes(size_t fc)const;

///const_iterator end()const{return m_sysgls.end();};
///iterator end(){return m_sysgls.end();};

///Invoke pre_transform_process of m_sysgls_zebra;
void pre_transform_process();

///Push back (function_code, object id) to the end or the beginning of the list.
///Function's return value is the display list name to use for glNewList.
///size_t push_back(size_t fc, MGGel* oi);
///size_t push_front(size_t fc, MGGel* oi);
size_t push_back(size_t fc, const MGGel* oi);
size_t push_front(size_t fc, const MGGel* oi);

///sysgl must be a newed object, and the ownership will be 
///transfered to this.
size_t push_back(mgSysGL* sysgl);

///Get the size of this list.
size_t size()const{return m_sysgls.size();};

private:

};

/** @} */ // end of DisplayHandling group
#endif
