/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTextureImage2_HH_
#define _MGTextureImage2_HH_

#include <iosfwd>
#include <memory>
#include "mg/MGCL.h"
#include "mgGL/TextureImage.h"

/** @addtogroup GLAttrib
 *  @{
 */

///Define MGTextureImage2 Class.
class MGCLASS MGTextureImage2:public MGTextureImage{

public:

MGTextureImage2(
):m_repeatS(true),m_repeatT(true){;};

MGTextureImage2(
	MGImage& image,
	double total_width,	///<width of the image in the world coordinate.
	double total_height,///<height of the image in the world coordinate.
	bool RepeatS=true,
	bool RepeatT=true
):MGTextureImage(image,0.,0.,total_width, total_height),
	m_repeatS(m_repeatS),m_repeatT(RepeatT){;};

MGTextureImage2(
	Gdiplus::Bitmap& gdi_image,
	double total_width,	///<width of the image in the world coordinate.
	bool Repeat=true
);

///render process.
void make_texture_object()const;

bool repeatS()const{return m_repeatS;};
bool repeatT()const{return m_repeatT;};

private:
	bool m_repeatS, m_repeatT;	///<repeat flag for the texture coordinate (s,t) each.
				///<If true, GL_REPEAT will be applied, if false, GL_CLAMP be applied.

};

/** @} */ // end of GLAttrib group
#endif
