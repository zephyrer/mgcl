/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGelPosition_HH_
#define _MGGelPosition_HH_

#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/Group.h"

//
//Define MGGelPosition Class.

/** @addtogroup GelRelated
 *  @{
 */

///MGGelPosition is a class which expresses which group a gel belongs to.
class MGCLASS MGGelPosition{

public:

///String stream function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGGelPosition&);

////////////Constructor////////////

///Void constructor.
MGGelPosition():m_gel(0), m_group(0){;};

///Constructor.
explicit MGGelPosition(MGGel* gel, MGGroup* group)
:m_gel(gel), m_group(group){;};

///Copy constructor.
MGGelPosition(const MGGelPosition& obj2){
	m_gel=obj2.m_gel;m_group=obj2.m_group;
}

///Destructor
virtual ~MGGelPosition(){;};

///Assignment.
virtual MGGelPosition& operator= (const MGGelPosition& GelPosition2);

///Equal operator
bool operator== (const MGGelPosition& gelp2) const;
bool operator!= (const MGGelPosition& gelp2) const{return !(*this==gelp2);}

bool operator<(const MGGelPosition& gp2)const;
bool operator>(const MGGelPosition& gp2)const{return gp2<(*this);};
bool operator<=(const MGGelPosition& gp2)const{return !((*this)>gp2);};
bool operator>=(const MGGelPosition& gp2)const{return !(gp2>(*this));};

////////////Member Function////////////

///Return the MGGel
const MGGel* gel()const{return m_gel;};
MGGel* gel(){return m_gel;};

///Return if the bottom gel is a MGGroup or MGObject.
bool gel_is_object()const;

///Clear the content.
void clear(){ m_gel = 0; m_group=0;}

///Generate a newed clone object.
virtual MGGelPosition* clone()const;

///perform add operation of this gel position.
///(push back the gel of m_gel to the group of m_group)
void do_add();

///perform remove operation of this gel position.
///(Release the gel of m_gel from the group of m_group, but does not delete the gel).
void do_remove();

///Get the group pointer.
const MGGroup* group()const{return m_group;};
MGGroup* group(){return m_group;};

///Test if this is null.
bool is_null() const{ return m_gel==0;}

///Get the object pointer of this.
const MGObject* object()const;
MGObject* object();

///Set the gel data.
void set_gel(MGGel* gel){m_gel=gel;};

///Set the group data.
void set_group(MGGroup* group){m_group=group;};

///Set this as null.
void set_null(){m_gel=0; m_group=0;};

///Test if this is symmetric to gel2.
///Symmetric means:
///both gels are MGObject and they have the same manifold dimension.
bool symmetric(const MGGelPosition& gp2)const;

protected:
	MGGel*	m_gel;	///<A Gel pointer.
	MGGroup* m_group;///<The group pointer which (is supposed to) includes the gel m_gel.
};

/** @} */ // end of GelRelated group
#endif
