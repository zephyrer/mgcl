/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mgGL/TextureImages.h"

//
//Define MGTextureImage Class.
//MGTextureImage defines the image data for OpenGL's glTexImage2D.
//That is, the image size will be rounded to 2's power.

MGTextureImage::MGTextureImage(
	Gdiplus::Bitmap& bitmap,	//texture image data.
	double left, double bottom,	//texture coordinates of left and bottom point of the image.
	double right, double top,	//texture coordinates of right and top point of the image.
	bool scalable				//indicates if the image can be scalabe or not.
			//if scalable=false, (right, top) will be changed according to the expanded image area.
):m_image(bitmap),m_box(MGPosition(left,bottom), MGPosition(right,top)){
	round(scalable);
}

MGTextureImage::MGTextureImage(
	MGImage& bitmap,		//image data.
	double left, double bottom,//texture coordinates of left and bottom point of the image.
	double right, double top,//texture coordinates of right and top point of the image.
	bool scalable				//indicates if the image can be scalabe or not.
			//if scalable=false, (right, top) will be changed according to the expanded image area.
):m_image(bitmap),m_box(MGPosition(left,bottom), MGPosition(right,top)){
	round(scalable);
}

//copy constructor.
MGTextureImage::MGTextureImage(
	const MGTextureImage& txtr
	):m_box(txtr.m_box){
	MGTextureImage* timage=const_cast<MGTextureImage*>(&txtr);
	m_image=timage->m_image;
}

//Assignment.
MGTextureImage& MGTextureImage::operator=(const MGTextureImage& txtr){
	m_box=txtr.m_box;
	MGTextureImage* timage=const_cast<MGTextureImage*>(&txtr);
	m_image=timage->m_image;
	return *this;
}

//render process.
void MGTextureImage::bind_texture_object()const{
	size_t tname=dlist_name();
	glDeleteTextures(1,&tname);
	glBindTexture(GL_TEXTURE_2D,tname);
}

void MGTextureImage::apply_texture_object()const{
	size_t tname=dlist_name();
	glBindTexture(GL_TEXTURE_2D,tname);
}

void MGTextureImage::round(bool scalable){
	if(!m_image.image())
		return;

	int w1=m_image.width(), h1=m_image.height();
	int w2,h2;
	MGImageCompute_2spower(w1,h1,w2,h2);
	if(w1==w2 && h1==h2)
		return;

	if(scalable){
		m_image.resize(w2,h2);
		return;
	}

	MGImage::pixel onepixel;onepixel.uint=0;
	m_image.resize_with_fill_color(w2,h2,onepixel);
	if(w1!=w2)
		m_box[0].set_high_point(left()+width()*double(w2)/double(w1));
	if(h1!=h2)
		m_box[1].set_high_point(bottom()+height()*double(h2)/double(h1));
}

//Set texture coordinates of this.
void MGTextureImage::set_address(
	double left, double bottom,
	double right, double top
){
	m_box=MGBox(MGPosition(left,bottom), MGPosition(right,top));
}
