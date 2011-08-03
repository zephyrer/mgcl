/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mgGL/TextureImages.h"

//
//Implements MGTextureImages Class.
//MGTextureImages defines the attributes of TextureImages.

MGTextureImages::MGTextureImages(
	Gdiplus::Bitmap& bitmap,
	int texture_size,	//texture size to use.
	double pixelw, double pixelh,//One pixel size, width and height.
	bool repeatS,bool repeatT//repaet to s or t coordinate.
):m_repeatS(repeatS), m_repeatT(repeatT),m_texture_size(texture_size){
	build_textures(bitmap,pixelw,pixelh);
}

MGTextureImages::MGTextureImages(
	Gdiplus::Bitmap& bitmap,
	int texture_size,	//texture size to use.
	double pixel_length,//One pixel size in world coordinate for width and height,
						//which are the same.
	bool repeat			//repeat to s or t coordinate, which are the same.
):m_repeatS(repeat), m_repeatT(repeat),m_texture_size(texture_size){
	build_textures(bitmap,pixel_length,pixel_length);
}

//Build textures after m_repeatS, m_repeatT, and m_texture_size are set up.
void MGTextureImages::build_textures(
	Gdiplus::Bitmap& bitmap,
	double pixelw, double pixelh//One pixel size, width and height.
){
	m_total_pixel_w=bitmap.GetWidth();
	m_total_pixel_h=bitmap.GetHeight();
	m_total_width=double(m_total_pixel_w)*pixelw;
	m_total_height=double(m_total_pixel_h)*pixelh;

	int nwspan=m_n_columns=m_total_pixel_w/m_texture_size;
	int nhspan=m_n_rows=m_total_pixel_h/m_texture_size;
	m_remainder_w=m_total_pixel_w%m_texture_size;
	m_remainder_h=m_total_pixel_h%m_texture_size;
	if(m_remainder_w) m_n_columns++;
	if(m_remainder_h) m_n_rows++;
	m_remainder_wd=m_remainder_w*pixelw;
	m_remainder_hd=m_remainder_h*pixelh;

	m_1width=pixelw*double(m_texture_size);
	m_1height=pixelh*double(m_texture_size);
	double wremainder=pixelw*double(m_remainder_w);
	double hremainder=pixelh*double(m_remainder_h);

	int ix, iy=0;//pixel address toward width and height.
	double x,xend, y=0.,yend;
	for(int j=0; j<nhspan; j++){
		ix=0;
		x=0.;
		yend=y+m_1height;
		for(int i=0; i<nwspan; i++){
			xend=x+m_1width;
			MGImage image(bitmap,ix,iy,m_texture_size,m_texture_size);
			push_back(MGTextureImage(image,x,y,xend,yend,false));
			ix+=m_texture_size;
			x=xend;
		}
		if(m_remainder_w){
			MGImage image(bitmap,ix,iy,m_remainder_w,m_texture_size);
			push_back(MGTextureImage(image,x,y,x+wremainder,yend,false));
		}
		iy+=m_texture_size;
		y=yend;
	}
	if(m_remainder_h){
		ix=0;
		x=0.;
		yend=y+hremainder;
		for(int i=0; i<nwspan; i++){
			xend=x+m_1width;
			MGImage image(bitmap,ix,iy,m_texture_size,m_remainder_h);
			push_back(MGTextureImage(image,x,y,xend,yend,false));
			ix+=m_texture_size;
			x=xend;
		}
		if(m_remainder_w){
			MGImage image(bitmap,ix,iy,m_remainder_w,m_remainder_h);
			push_back(MGTextureImage(image,x,y,x+wremainder,yend,false));
		}
	}
}

//Refer to (i,j)th MGTextureImage.
const MGTextureImage& MGTextureImages::operator()(int i, int j)const{
	assert(0<=i&&i<m_n_columns && 0<=j&&j<m_n_rows);
	return m_images[m_n_columns*j+i];
}

//Refer to (i,j)th MGTextureImage.
MGTextureImage& MGTextureImages::operator()(int i, int j){
	assert(0<=i&&i<m_n_columns && 0<=j&&j<m_n_rows);
	return m_images[m_n_columns*j+i];
}

//compute the center of the image in the world (s,t) coordinage.
MGPosition MGTextureImages::center()const{
	return MGPosition(m_total_width*.5,m_total_height*.5);
}

//Locate at which index of the mesh the value belongs to.
//Function's return value is the id.
int MGTextureImage_get_ijwh2(
  double span,	//input the span length of the mesh.
  double value,	//the value to get the index.
  bool maximum	//true if the value is on the maximum side.
 ){
	double error=MGTolerance::wc_zero();
	int i=int(value/span);
	double ivalue=double(i)*span;
	if(value>=0){
		if(maximum){
			if(value<ivalue+error) i--;
		}else{
			if(value>ivalue+span-error) i++;
		}
	}else{
		if(maximum){
			i--;
			if(value<ivalue-span+error) i--;
		}else{
			if(value<=ivalue-error) i--;
		}
	}
	return i;
}

void MGTextureImage_get_ijwh(
  double width, double height,
  double s0,double s1,
  double t0,double t1,
  int&i, int&j,
  int& nwidth, int&nheight
 ){
	i=MGTextureImage_get_ijwh2(width,s0,false);
	j=MGTextureImage_get_ijwh2(height,t0,false);
	int i1=MGTextureImage_get_ijwh2(width,s1,true);
	int j1=MGTextureImage_get_ijwh2(height,t1,true);
	nwidth=i1-i+1;
	nheight=j1-j+1;
}

//Given stbox, find which MGTextureImage this stbox belongs to.
//Function's return value is:
//0: stbox is outside all the MGTextureImage's.
//1: stbox is inside from (*this)(i,j) to (*this)(i+nwidth-1, j+nheight-1).
//2: stbox is inside from (*this)(i,j) to (*this)(i+nwidth-1, j+nheight-1) when
//	the output base is subtracted from stbox.
//3: a part of stbox is inside total_box and a part is outside total_box.
int MGTextureImages::compute_image_range(
	const MGBox& stbox,	//Box of the target triangle.
	int& iTotal, int& jTotal,//id of the rectangle of this MGTextureImages.
	int& nwTotal, int& nhTotal,//number of repitition of this MGTextureImages.
	int& i, int&j,		//id in the 1st MGTextureImages will be output.
	int& nwidth, int& nheight,//number of repitition of MGTextureImage in the 1st MGTextureImage.
		//i,j, nwidth, nheight are valid only when return=1 or 2.
	MGVector& base		//Base(left, bottom) address of the (i,j) MGTextureImage.
)const{
	double error=MGTolerance::wc_zero();

	MGPosition minST=stbox.low();
	MGPosition maxST=stbox.high();
	double t0=minST[1], t1=maxST[1];
	if(!repeatT()){
		if(t0>=m_total_height-error || t1<=error)
			return 0;
	}
	double s0=minST[0], s1=maxST[0];
	if(!repeatS()){
		if(s0>=m_total_width-error || s1<=error)
			return 0;
	}

	base=MGVector(0.,0.);
	if(s0>-error && s1<=m_total_width+error && t0>-error && t1<=m_total_height+error){
		//if stbox is inside total_box.
		MGTextureImage_get_ijwh(m_1width,m_1height,s0,s1,t0,t1,i,j,nwidth,nheight);
		return 1;
	}

	MGTextureImage_get_ijwh(m_total_width,m_total_height,s0,s1,t0,t1,iTotal,jTotal,nwTotal,nhTotal);
	double& sbase=base(0);
	if(repeatS()){
		sbase=m_total_width*double(iTotal);
		s0-=sbase; s1-=sbase;
	}else{
		iTotal=0; nwTotal=1;
	}

	double& tbase=base(1);
	if(repeatT()){
		tbase=m_total_height*double(jTotal);
		t0-=tbase; t1-=tbase;
	}else{
		jTotal=0; nhTotal=1;
	}

	if(s0>-error && s1<m_total_width+error && t0>-error && t1<=m_total_height+error){
		//if stbox is inside total_box.
		MGTextureImage_get_ijwh(m_1width,m_1height,s0,s1,t0,t1,i,j,nwidth,nheight);
		return 2;
	}

	return 3;
}

double MGTextureImages::pixel_width()const{
	return total_width()/double(m_total_pixel_w);
}
double MGTextureImages::pixel_height()const{
	return total_height()/double(m_total_pixel_h);
}

double MGTextureImages::total_width()const{
	return m_total_width;
}
double MGTextureImages::total_height()const{
	return m_total_height;
}

//delete texture objects made by make_texture_object();
void MGTextureImages::delete_texture_object(){
	for(int j=0; j<m_n_rows; j++){
	for(int i=0; i<m_n_columns; i++){
		const MGTextureImage& timage=(*this)(i,j);
		size_t tname=timage.dlist_name();
		glDeleteTextures(1,&tname);//Erase texture object.
	}
	}
}

//render process.
//Function's return value is the lightnumber of this light executed.
void MGTextureImages::make_texture_object()const{
	for(int j=0; j<m_n_rows; j++){
	for(int i=0; i<m_n_columns; i++){
		const MGTextureImage& timage=(*this)(i,j);
		timage.bind_texture_object();
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,timage.nwidth(),timage.nheight(),0,GL_RGBA,
			GL_UNSIGNED_BYTE,timage.image());
	}
	}
}
