/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mgGL/TextureImage2.h"

//
//Implements MGTextureImage2 Class.

MGTextureImage2::MGTextureImage2(
	Gdiplus::Bitmap& gdi_image,
	double total_width,	//width of the image in the world coordinate.
	bool Repeat
):m_repeatS(Repeat),m_repeatT(Repeat){
	m_image=MGImage(gdi_image);

	int nwidth = m_image.width(), nheight = m_image.height();
	double pixel_length=total_width/double(nwidth);
	double height=pixel_length*double(nheight);

	int nwidth2,nheight2;
	MGImageCompute_2spower(nwidth,nheight,nwidth2,nheight2);
	if(Repeat){
		m_image.resize(nwidth2,nheight2);
	}else{
		m_image.resize_and_add_zero_border(nwidth2,nheight2);
	}
	set_address(0.,0.,total_width,height);
}

//render process.
//Function's return value is the lightnumber of this light executed.
void MGTextureImage2::make_texture_object()const{
	int tex_repeatS=GL_REPEAT;
	if(!m_repeatS){
		tex_repeatS=GL_CLAMP;
	}
	int tex_repeatT=GL_REPEAT;
	if(!m_repeatT){
		tex_repeatT=GL_CLAMP;
	}

	bind_texture_object();
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,tex_repeatS);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,tex_repeatT);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,nwidth(),nheight(),0,GL_RGBA,GL_UNSIGNED_BYTE,image());
}
