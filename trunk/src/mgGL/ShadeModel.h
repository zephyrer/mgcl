/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGShadeModel_HH_
#define _MGShadeModel_HH_

#include "mgGL/GLAttrib.h"

class MGOfstream;
class MGIfstream;

//Define MGShadeModel Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGShadeModel defines shading model.
class MGCLASS MGShadeModel{

public:

enum MODEL{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	SMOOTH=GL_SMOOTH,
	FLAT=GL_FLAT
};

/// Serialization fucntion.
MGDECL friend MGOfstream& operator<< (MGOfstream& buf, const MGShadeModel& sm);
MGDECL friend MGIfstream& operator>> (MGIfstream& buf, MGShadeModel& sm);

MGShadeModel(MODEL m=UNDEFINED):m_model(m){;};

////////////Member Function////////////

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

void set_model(MODEL m=SMOOTH){m_model=m;};
MODEL get_model()const{return m_model;};
GLenum GLmodel()const{return static_cast<GLenum>(m_model);};

/// Return This object's typeID
long identify_type() const{return MGSHADE_MODEL_TID;};

/// Output function.
std::ostream& out(std::ostream&) const;

private:
	MODEL m_model;///Shade model
};

/** @} */ // end of GLAttrib group
#endif
