/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGAttrib_HH_
#define _MGAttrib_HH_

#include "mg/Pvector.h"
#include "mg/Gel.h"
#include "mg/types.h"

//
//Define MGAttrib Class.

class MGIfstream;
class MGOfstream;
class MGAttrib;
class MGGLAttrib;
class MGContext;

/** @addtogroup GelRelated
 *  @{
 */

///MGAttrib is an abstract class which represents a system attribute element.
///Currently main subclasses of MGAttrib are MGGLAttrib for OpenGL attributes.
class MGCLASS MGAttrib: public MGGel{

public:

///Virtual Destructor
virtual ~MGAttrib();

///Assignment.
///When the leaf objects of this and gel2 are not equal, this assignment does nothing.
virtual MGAttrib& operator=(const MGAttrib& gel2){MGGel::operator=(gel2); return *this;};

////////Member Function////////

/// Output virtual function.
virtual std::ostream& out(std::ostream&) const=0;

///Return MGAttrib pointer if this MGGel is an MGAttrib, else return null.
MGAttrib* attrib(){return this;};
const MGAttrib* attrib()const{return this;};

///Obtain display list name. 0(null) means this gel need not to be displayed.
///dlist_name() of MGAttrib is null.
size_t dlist_name()const{return 0;};

/// Return This object's typeID
///Sub class of MGAttrib must have the id as 0x02nnnnnnL.
virtual long identify_type() const{return MGATTRIB_TID;};

///Test if this gel includes an object.
const MGObject* includes_object()const{return 0;};

///Test if this gel includes an object.
MGObject* includes_object(){return 0;};

protected:

///Read all member data.
virtual void ReadMembers(MGIfstream& buf);

///Write all member data
virtual void WriteMembers(MGOfstream& buf)const;

private:

friend class MGIfstream;
friend class MGOfstream;

};

typedef MGPvector<MGGLAttrib> MGGLAttribs;///<attributes.

///Construct a null newed MGAttrib from the type id TID.
MGDECL MGAttrib* MGNullAttrib(long TID);

/** @} */ // end of GelRelated group
#endif
