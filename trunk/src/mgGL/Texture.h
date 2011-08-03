/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTexture_HH_
#define _MGTexture_HH_

#include "mg/MGCL.h"
#include <iosfwd>
#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//Define MGTexture Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGTexture defines the attributes of texture.
class MGCLASS MGTexture{

public:
  
/// Serialization fucntion.
MGDECL friend MGOfstream& operator<< (MGOfstream& buf, const MGTexture& txtr);
MGDECL friend MGIfstream& operator>> (MGIfstream& buf, MGTexture& txtr);

MGTexture(){;};

/*
///copy constructor.
MGTexture(const MGTexture& txtr);

///Destructor.
~MGTexture();

///Assignment.
MGTexture& operator=(const MGTexture& attr);
*/

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

/// Return This object's typeID
long identify_type() const{return MGTEXTURE_TID;};

/// Output function.
std::ostream& out(std::ostream&) const;

};

/** @} */ // end of GLAttrib group
#endif
