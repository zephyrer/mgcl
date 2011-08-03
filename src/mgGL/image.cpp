/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mgGL/Image.h"

//
//Implements MGTextureImages Class.
//MGTextureImages defines the attributes of TextureImages.

MGImage::MGImage(
	Gdiplus::Bitmap& bitmap
):m_width(bitmap.GetWidth()),m_height(bitmap.GetHeight()),
m_image(new pixel[m_width*m_height]){
	extract(bitmap,0,0,m_width,m_height);
}

//Extract a part of bitmap.
MGImage::MGImage(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	//left bottom address of bitmap.
	int width, int height
):m_width(width),m_height(height),
m_image(new pixel[width*height]){
	extract(bitmap,x,y,width,height);
}

MGImage::pixel& MGImage::operator()(int i, int j){
	pixel* pixelP=image();
	return pixelP[j*m_width+i];
}
const MGImage::pixel& MGImage::operator()(int i, int j)const{
	const pixel* pixelP=image();
	return pixelP[j*m_width+i];
}

//Extract a part of bitmap into this, from(x,y) to (x+width, y+height).
void MGImage::extract(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	//left bottom address of bitmap.
	int width, int height
){
	assert(x+width<=int(bitmap.GetWidth()));
	assert(y+height<=int(bitmap.GetHeight()));

	int totalhm1=bitmap.GetHeight()-1;//total height -1.
	// Change the pixel data format
	int hm1=m_height-1;
	pixel* pixelP=m_image.get();
	if(!pixelP)
		return;

	for(int j=0; j<height; j++){
		int jpy=j+y;
		int jrow_s=j*m_width;
		int ybitmap=totalhm1-(j+y);
			//This is because MFC's row address is reverse order to MGImage's.
		for(int i=0; i<width; i++){
			Gdiplus::Color gdipixel;
			bitmap.GetPixel(i+x, ybitmap, &gdipixel);
			pixel& pixelij=pixelP[jrow_s+i];
			pixelij.uchar[0]=gdipixel.GetR();
			pixelij.uchar[1]=gdipixel.GetG();
			pixelij.uchar[2]=gdipixel.GetB();
			pixelij.uchar[3]=gdipixel.GetA();
		}
	}
}

//resize the image size to (width,height).
void MGImage::resize(
	int width,
	int height
){
	std::auto_ptr<pixel> glimage2ptr(new pixel[width*height]);
	if(!glimage2ptr.get())
		return;

	int error=gluScaleImage(GL_RGBA,
		m_width,m_height,GL_UNSIGNED_BYTE,m_image.get(),
		width,height,GL_UNSIGNED_BYTE,glimage2ptr.get());
	if(!error){
		m_image=glimage2ptr;
		m_width=width;
		m_height=height;
	}
}

//Resize the image size to (width,height) filling the color to
//the extra part for the size(width,height).
//resize_with_fill_color() does not perform scaling to the image.
//width and height can be less than the original length. In this case,
//image trimming will be done.
void MGImage::resize_with_fill_color(
	int width,
	int height,
	const pixel& pdata
){
	std::auto_ptr<pixel> glimage2ptr(new pixel[width*height]);
	pixel* pixelP2=glimage2ptr.get();
	if(!pixelP2)
		return;

	int w2=m_width;
	if(width<w2) w2=width;
	int h2=m_height;
	if(height<h2) w2=height;

	int i,j;
	pixel* pixelP1=m_image.get();
	for(j=0; j<h2; j++){
		int jw1=j*m_width, jw2=j*width;
		for(i=0; i<w2; i++)
			pixelP2[jw2+i].uint=pixelP1[jw1+i].uint;
		for(;i<width; i++)
			pixelP2[jw2+i].uint=pdata.uint;

	}
	for(; j<height; j++){
		int jw2=j*width;
		for(i=0; i<width; i++)
			pixelP2[jw2+i].uint=pdata.uint;

	}
	m_image=glimage2ptr;
	m_width=width;
	m_height=height;
}

//Add border color(0,0,0,0)=(alfa=0) of pixel size2 for each perimeter.
void MGImage::resize_and_add_zero_border(int nwidth2,int nheight2){
	resize(nwidth2-4,nheight2-4);
	int i, j;
	std::auto_ptr<pixel> image3ptr(new pixel[nwidth2*nheight2]);
	if(!image3ptr.get())
		return;

	int width2by2=nwidth2*2;
	int width2Mborder=nwidth2-4;

	pixel* image2=m_image.get();
	pixel* image3=image3ptr.get();
	for(i=0; i<width2by2; i++) image3[i].uint=0;
	for(j=0; j<nheight2-4; j++){
		int jrow_s3=(j+2)*nwidth2;
		int jrow_s2=j*width2Mborder;
		image3[jrow_s3].uint=0;image3[jrow_s3+1].uint=0;
		int jrow_s3p2=jrow_s3+2;
		for(i=0; i<width2Mborder; i++){
			image3[jrow_s3p2+i].uint=image2[jrow_s2+i].uint;
		}
		image3[jrow_s3p2+width2Mborder].uint=0;
		image3[jrow_s3p2+width2Mborder+1].uint=0;
	}
	size_t erow_s3=(nheight2-2)*nwidth2;
	for(i=0; i<width2by2; i++) image3[erow_s3+i].uint=0;
	m_image=image3ptr;
	m_width=nwidth2;
	m_height=nheight2;
}

//Compute 2's power of width and height.
void MGImageCompute_2spower(
	int width, int height,
	int& width2, int& height2//The smallest 2's power of width and height will be output.
){
	size_t nshift=5;
	width2=64;
	while(width2<width && nshift<32){
		width2=width2<<1;nshift++;
	}

	nshift=5;
	height2=64;
	while(height2<height && nshift<32){
		height2=height2<<1 ;nshift++;
	}
}
